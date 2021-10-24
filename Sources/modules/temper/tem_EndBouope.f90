subroutine tem_EndBouope(a)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!    for the moment the routine only calculates forces on bodies
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use Mod_iofile
   use Mod_int2str
   use Mod_Temperature
   use MPI
   implicit none
   class(TemperatureProblem) :: a   
   class(FiniteElement), pointer :: e => NULL() 
   real(rp)   , allocatable :: heatfR(:)
   real(rp)   , allocatable :: eltem(:),grate(:)
   character(150) :: fil_force  !Output data file
   real(rp)    :: weightfactor  !Parallel factor
   real(rp)    :: acden,acsph,actco,acrea,acsou,dsurf,dsurf0
   integer(ip) :: idime,iboun,ndime,igaub,jdime,nbody,nboun,nelem,ibody,ielem
   integer(ip), pointer :: lbody => NULL()
   integer(ip) :: icount,npoinLocal,inodb,ierr
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)   
   call a%Mesh%GetNbody(nbody)
   call a%Mesh%GetNpoinLocal(npoinLocal)     
   
   if(a%kfl_openforce) then
      if (a%MPIrank == a%MPIroot) then 
         do ibody=1,nbody
            fil_force = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//adjustl(trim(int2str(ibody)))//'.frc'
            call iofile(zero,a%lun_force(ibody),fil_force,adjustl(trim(a%exmod))//'FORCES & MOMENTS')
            if (ndime==2) write (a%lun_force(ibody),12)
            if (ndime==3) write (a%lun_force(ibody),13)
         end do
         a%kfl_openforce = .false.
      endif      
   end if   
  

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','tem_EndBouope')      
   
   !Allocation
   call a%Memor%alloc(e%mnode,eltem,'eltem','tem_EndBouope')
   call a%Memor%alloc(ndime,grate,'grate','tem_EndBouope')
   !working variables (reuced variables)
   call a%Memor%alloc(nbody,heatfR,'heatfR','tem_EndBouope')
  
   !Inializations
   !Force and moments
   a%heatf=0.0_rp
   !Reduced values
   heatfR=0.0_rp
   !Workin variables
   ibody=0
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)
      
      call a%Mesh%GetLbody(lbody,iboun)
       
      ibody=lbody
      
      !to delete repeated entities 
      icount = 0
      do inodb = 1,e%pnodb
         if (e%lnodb(inodb) <= npoinLocal) then
            icount = icount +1
         endif
      enddo
      weightfactor = real(icount)/real(e%pnodb)
      
      if(ibody>0)then
         !Physical Parameters
         call a%Mesh%GetBoundaryIelem(iboun,ielem)
         call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
         
         call e%elmdel
         
         !Gathers in boundary elements 
         call e%gather(1,eltem,a%tempe(:,1))
        
         
         dsurf0 = 0.0_rp           

         !Gauss-Point Loop
         do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            dsurf0 = dsurf0 + dsurf
            
            !Derivatives at the boundary
            call e%elmderb         
            call e%gradientb(1,eltem,grate)
            call e%gradientb(1,eltem,grate)
            
            !----------------------------------------------------------------
            !Heat flux
            !  k*n*grad(T)              
            heatfR(ibody) = heatfR(ibody) + weightfactor*dot_product(grate,e%baloc(:,ndime))*dsurf*actco
             
         end do 
      
      end if
      
   end do boundaries 

   !Collect Values
   call MPI_REDUCE(heatfR,a%heatf,nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
   
   ! Output forces
   
   if (a%MPIrank == a%MPIroot) then  
      do ibody = 1,nbody
         write(a%lun_force(ibody),20) a%ctime, a%heatf(ibody)
      end do
   end if

   if (a%kfl_flush == 1) then
      do ibody = 1,nbody
         call flush(a%lun_force(ibody))
      end do
   endif
   
   !Dellocation
   call a%Memor%dealloc(e%mnode,eltem,'eltem','tem_EndBouope')
   call a%Memor%dealloc(ndime,grate,'grate','tem_EndBouope')
   !working variables (reduced variables)
   call a%Memor%dealloc(nbody,heatfR,'heatfR','tem_EndBouope')
   
   
   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','tem_EndBouope')
   
   !Formats

   12   format( &
          '$        Heat flux over bodies          ',/,          &
          '$  Time',15x,'GraT')
   13   format( &
          '$          Heat flux over bodies            ',/,      &
          '$ Time  GraT')

   20   format(1x,e12.6,2x,1(2x,e15.8))
   30   format(1x,e12.6,2x,1(2x,e15.8))

   50   format(a1,i1)
   51   format(a1,i2)
   52   format(a1,i3)
   53   format(5x,'WARNING: CANNOT WRITE MORE THAN 999 bodies')   
   
end subroutine













subroutine lmn_EndBouope(a)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!    for the moment the routine only calculates forces on bodies
!
!-----------------------------------------------------------------------
   use MPI
   use Mod_iofile
   use Mod_int2str
   use typre 
   use Mod_Element
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a   
   class(FiniteElement), pointer :: e => NULL() 
   real(rp)   , allocatable :: tract(:),forceR(:,:),momenR(:,:),heatfR(:)
   real(rp)   , allocatable :: grvel(:,:),elvel(:,:,:),elpre(:,:),gpcod(:),bopre(:),eltem(:),grate(:)
   character(150) :: fil_force  !Output data file
   real(rp)    :: weightfactor  !Parallel factor
   real(rp)    :: prgau(1),divu,acden,acvis,actco,dsurf,dsurf0,vnor
   integer(ip) :: idime,iboun,ndime,igaub,jdime,nbody,nboun,nelem,ibody
   integer(ip) :: icount,npoinLocal,inodb,ierr,ibopo,ipoin
   integer(ip), pointer :: lbody => NULL()
   real(rp), pointer    :: exnor(:,:) => NULL()
   logical  :: cycleflag
   
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
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_EndBouope')      
   
   !Allocation
   call a%Memor%alloc(ndime,tract,'tract','lmn_EndBouope')
   call a%Memor%alloc(ndime,e%mnode,1,elvel,'elvel','lmn_EndBouope')
   call a%Memor%alloc(e%mnode,eltem,'eltem','lmn_EndBouope')
   call a%Memor%alloc(      e%mnode,1,elpre,'elpre','lmn_EndBouope')  
   call a%Memor%alloc(ndime,ndime,grvel,'grvel','lmn_EndBouope')
   call a%Memor%alloc(ndime,gpcod,'gpcod','lmn_EndBouope')
   call a%Memor%alloc(      e%mnodb,bopre,'bopre','lmn_EndBouope')
   call a%Memor%alloc(ndime,grate,'grate','lmn_EndBouope')
   !working variables (reuced variables)
   call a%Memor%alloc(ndime,nbody,forceR,'forceR','lmn_EndBouope')
   call a%Memor%alloc(3,nbody,momenR,'momenR','lmn_EndBouope')
   call a%Memor%alloc(nbody,heatfR,'heatfR','lmn_EndBouope')
  
   !Inializations
   !Force and moments
   a%force=0.0_rp
   a%momen=0.0_rp
   a%heatf=0.0_rp
   !Reduced values
   forceR=0.0_rp
   momenR=0.0_rp
   heatfR=0.0_rp
   !Workin variables
   ibody=0

   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      cycleflag = .false.
      do inodb = 1,e%pnodb
         ipoin = e%lnodb(inodb)
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if (ibopo == 0) then
            cycleflag =  .true.
         else
            call vecnor(exnor(:,1),e%ndime,vnor,2)
            if (vnor == 0.0_rp) cycleflag =  .true. 
         end if
      end do

      if (cycleflag) cycle

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
         call a%GetPhysicalParameters(acden=acden,acvis=acvis,actco=actco) 
         
         call e%elmdel
         
         !Gathers in boundary elements 
         call e%gatherb(1_ip   ,bopre,a%press(1:e%pnode,1))       
         call e%gather(ndime,elvel,a%veloc(:,:,1))      
         call e%gather(1,eltem,a%tempe(:,1))
         call e%divergence(elvel,divu)
        
         
         dsurf0 = 0.0_rp           

         !Gauss-Point Loop
         do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Traction at gauss point
            tract=0.0_rp

            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            dsurf0 = dsurf0 + dsurf
            
            call e%interpb(ndime,e%bocod,gpcod)
            call e%interpb(1,bopre,prgau)
            
            !Derivatives at the boundary
            call e%elmderb         
            call e%gradientb(ndime,elvel,grvel)
            call e%gradientb(1,eltem,grate)
            call e%gradientb(1,eltem,grate)
      
            !----------------------------------------------------------------
            !Traction value
            !  t=n*S  S=-p*I+2mu*gradS(u)                  
            do idime=1,ndime
                tract(idime)= tract(idime) + e%baloc(idime,ndime) * (0.66666666_rp*acvis*divu + prgau(1))
               
               do jdime=1,ndime
                  tract(idime) = tract(idime) - acvis*e%baloc(jdime,ndime)*(grvel(idime,jdime) + grvel(jdime,idime))
               end do
            end do
            
            !---------------------------------------------------------------
            !Forces
            do idime=1,ndime
               forceR(idime,ibody)=forceR(idime,ibody) + tract(idime)*dsurf*weightfactor
            end do
            
            !Moments
            if(ndime==2)then
               momenR(3,ibody) = momenR(3,ibody) &
                  + tract(2)*(gpcod(1)-a%origm(1)) &
                  - tract(1)*(gpcod(2)-a%origm(2)) 
            elseif(ndime==3)then
               momenR(1,ibody) = momenR(1,ibody) &
                  + tract(3)*(gpcod(2)-a%origm(2)) &
                  - tract(2)*(gpcod(3)-a%origm(3))          
               momenR(2,ibody) = momenR(1,ibody) &
                  + tract(1)*(gpcod(3)-a%origm(3)) &
                  - tract(3)*(gpcod(1)-a%origm(1))  
               momenR(3,ibody) = momenR(3,ibody) &
                  + tract(2)*(gpcod(1)-a%origm(1)) &
                  - tract(1)*(gpcod(2)-a%origm(2))
            end if
            
            !----------------------------------------------------------------
            !Heat flux
            !  n*gradT              
            heatfR(ibody) = heatfR(ibody) + weightfactor*dot_product(grate,e%baloc(:,ndime))*dsurf*actco
             
         end do 
      
      end if
      
   end do boundaries 

   if (nbody > 0) then
      !Collect Values
      call MPI_REDUCE(forceR,a%force,ndime*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr) 
      call MPI_REDUCE(momenR,a%momen,3*nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr) 
      call MPI_REDUCE(heatfR,a%heatf,nbody,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
   end if
   
   ! Output forces
   
   if (a%MPIrank == a%MPIroot) then  
      
      if(ndime==2) then
         do ibody = 1,nbody
            write(a%lun_force(ibody),20) a%ctime, &
               -a%force(1,ibody)*a%adimf(1),      &
               -a%force(2,ibody)*a%adimf(2),      &
               -a%momen(3,ibody)*a%adimm(3),      &
               a%heatf(ibody)*a%adimh,            &
               a%pther(1)/a%itpre
         end do
      else if(ndime==3) then
         do ibody = 1,nbody
            write(a%lun_force(ibody),30) a%ctime, &
               -a%force(1,ibody)*a%adimf(1),      &
               -a%force(2,ibody)*a%adimf(2),      &
               -a%force(3,ibody)*a%adimf(3),      &
               -a%momen(1,ibody)*a%adimm(1),      &
               -a%momen(2,ibody)*a%adimm(2),      &
               -a%momen(3,ibody)*a%adimm(3),      &
               a%heatf(ibody)*a%adimh,            &
               a%pther(1)/a%itpre
         end do
      end if   
   
   end if

   if (a%kfl_flush == 1) then
      do ibody = 1,nbody
         call flush(a%lun_force(ibody))
      end do
   endif
   
   !Dellocation
   call a%Memor%dealloc(ndime,tract,'tract','lmn_EndBouope')   
   call a%Memor%dealloc(ndime,e%mnode,1,elvel,'elvel','lmn_EndBouope')
   call a%Memor%dealloc(e%mnode,eltem,'eltem','lmn_EndBouope')
   call a%Memor%dealloc(      e%mnode,1,elpre,'elpre','lmn_EndBouope')  
   call a%Memor%dealloc(ndime,ndime,grvel,'grvel','lmn_EndBouope') 
   call a%Memor%dealloc(ndime,gpcod,'gpcod','lmn_EndBouope')  
   call a%Memor%dealloc(      e%mnodb,bopre,'bopre','lmn_EndBouope')
   call a%Memor%dealloc(ndime,grate,'grate','lmn_EndBouope')
   !working variables (reduced variables)
   call a%Memor%dealloc(ndime,nbody,forceR,'forceR','lmn_EndBouope')
   call a%Memor%dealloc(3,nbody,momenR,'momenR','lmn_EndBouope')   
   call a%Memor%dealloc(nbody,heatfR,'heatfR','lmn_EndBouope')
   
   
   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','lmn_froces')
   
   !Formats

   12   format( &
          '$        Forces, Moments and Heat flux over bodies          ',/,          &
          '$  Time',15x,'Fx',15x,'Fy',15x,'Mz',15x,'GraT',15x,'Pther')
   13   format( &
          '$          Forces, Moments and Heat flux over bodies            ',/,      &
          '$ Time  Fx  Fy  Fz  Mx  My  Mz  GraT  Pther')

   20   format(1x,e12.6,2x,5(2x,e15.8))
   30   format(1x,e12.6,2x,8(2x,e15.8))

   50   format(a1,i1)
   51   format(a1,i2)
   52   format(a1,i3)
   53   format(5x,'WARNING: CANNOT WRITE MORE THAN 999 bodies')   
   
end subroutine












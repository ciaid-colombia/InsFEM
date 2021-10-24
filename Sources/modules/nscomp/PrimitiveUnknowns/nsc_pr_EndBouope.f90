subroutine nsc_pr_EndBouope(a)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!    The routine calculates forces on bodies
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_Mesh
   use Mod_Element
   use Mod_Memor
   use Mod_iofile
   use Mod_int2str
   use Mod_NSCompressiblePrimitive
   use Mod_NSCompressiblePrimitiveNBC
   use MPI
   
   implicit none
   class(NSCompressiblePrimitiveProblem) :: a   
   class(FiniteElement) , pointer:: e => NULL()
 

   real(rp), allocatable :: tract(:),forceR(:,:),momenR(:,:)
   real(rp), allocatable :: gpcod(:)
   real(rp), allocatable :: bopre(:,:),gpbp(:)
   real(rp), allocatable :: bovel(:,:,:),gpbv(:,:)
   real(rp), allocatable :: botem(:,:),gpbt(:)
   real(rp), allocatable :: elpre(:), elvel(:,:), eltem(:)
   real(rp), allocatable :: grpre(:), grvel(:,:), grtem(:)
   real(rp), allocatable :: qflx(:)
   real(rp), allocatable :: incre(:)
   real(rp), allocatable :: borhs(:,:),TimeIntegratedFluxes(:)

   character(150) :: fil_force
   real(rp)    :: weightfactor
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: gpbd,gpbe,divv
   real(rp)    :: tempn0(a%ndofn)
   integer(ip) :: idime,iboun,ndime,igaub,jdime,nbody,nboun,nelem,ibody
   integer(ip), pointer :: lbody => NULL()
   integer(ip) :: icount,npoinLocal,npoin,inodb,ipoin,ierr,idofn
   
   logical            :: isALE
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)   
   call a%Mesh%GetNbody(nbody)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)     
   call a%Mesh%GetALE(isALE)
   
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
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_EndBouope')      
   
   !Allocation

   call a%Memor%alloc(ndime+1,tract,'tract','nsc_pr_EndBouope')   
   call a%Memor%alloc(ndime,gpcod,'gpcod','nsc_pr_EndBouope')
   call a%Memor%alloc(e%pnodb,2,bopre,'bopre','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%pnodb,2,bovel,'bovel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%pnodb,2,botem,'botem','nsc_pr_EndBouope')
   call a%Memor%alloc(e%mnode,elpre,'elpre','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%mnode,eltem,'eltem','nsc_pr_EndBouope')
   call a%Memor%alloc(2,gpbp,'gpbp','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,2,gpbv,'gpbv','nsc_pr_EndBouope')
   call a%Memor%alloc(2,gpbt,'gpbt','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,grpre,'grpre','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,grtem,'grtem','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,qflx,'qflx','nsc_pr_EndBouope')
   call a%Memor%alloc(a%ndofn,incre,'incre','nsc_pr_EndBouope')
   !nscnbc(boundary problem)
   call a%Memor%alloc(a%ndofn,npoin,borhs,'borhs','nsc_pr_EndBouope')
   !working variables (reuced variables)
   call a%Memor%alloc(ndime+1,nbody,forceR,'forceR','nsc_pr_EndBouope')
   call a%Memor%alloc(3,nbody,momenR,'momenR','nsc_pr_EndBouope')
  
   !Inializations
   a%btraction=0.0_rp
   !Force and moments fields
   if(a%npp_stepi(11) /= 0) then
      a%forfl=0.0_rp
      a%momfl=0.0_rp
      a%hfxfl=0.0_rp
   end if
   if(a%kfl_outfm==1) then   
   !Force and moments
      a%force=0.0_rp
      a%formo=0.0_rp
   end if
   !Reduced values
   forceR=0.0_rp
   momenR=0.0_rp
   !Working variables
   ibody=0
   !nscnbc
   borhs=0.0_rp
   
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
      
      if((a%kfl_fixbo(iboun)==10) .or. (a%kfl_fixbo(iboun)==12) .or. ((ibody>0).and.(a%kfl_outfm==1)) .or. (a%kfl_restrictions==1))then
         !Physical Parameters
         call a%GetPhysicalParameters(acvis,actco,accph,accvh)         
         
         call e%elmdel
         
         !Gathers in boundary elements 
         call e%gatherb(1_ip,bopre(:,1),a%press(:,1))
         call e%gatherb(1_ip,bopre(:,2),a%press(:,3))
         call e%gatherb(e%ndime,bovel(:,:,1),a%veloc(:,:,1))
         call e%gatherb(e%ndime,bovel(:,:,2),a%veloc(:,:,3))
         call e%gatherb(1_ip,botem(:,1),a%tempe(:,1))
         call e%gatherb(1_ip,botem(:,2),a%tempe(:,3))
         call e%gather(1_ip,elpre,a%press(:,1))       
         call e%gather(ndime,elvel,a%veloc(:,:,1))      
         call e%gather(1_ip,eltem,a%tempe(:,1))       
        
         dsurf0 = 0.0_rp           

         !Gauss-Point Loop
         do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Traction at gauss point
            tract=0.0_rp
            divv=0.0_rp

            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            dsurf0 = dsurf0 + dsurf
            
            call e%interpb(ndime,e%bocod,gpcod)
            call e%interpb(1_ip,bopre(:,1),gpbp(1))
            call e%interpb(1_ip,bopre(:,2),gpbp(2))
            call e%interpb(e%ndime,bovel(:,:,1),gpbv(:,1))
            call e%interpb(e%ndime,bovel(:,:,2),gpbv(:,2))
            call e%interpb(1_ip,botem(:,1),gpbt(1))
            call e%interpb(1_ip,botem(:,2),gpbt(2))
            
            !Derivatives at the boundary
            call e%elmderb         
            call e%gradientb(1_ip,elpre,grpre)
            call e%gradientb(ndime,elvel,grvel)
            call e%gradientb(1_ip,eltem,grtem)

            !Velocity divergence
            do idime =1,e%ndime
              divv = divv + grvel(idime,idime)
            end do
  
            !Traction value (also for momentum conservation restrictions)
            !  t=n·S  S=-p*I+tau
            !  tau=(2mu)*gradS(vel)-(2mu/3)*div(vel)*I

            do idime=1,ndime
               tract(idime)= tract(idime) - e%baloc(idime,ndime)*(gpbp(1)+a%relpre)&
                -(2/3)*acvis*divv*e%baloc(idime,ndime)
               do jdime=1,ndime
                  tract(idime) = tract(idime) &
                     + acvis*e%baloc(jdime,ndime)*(grvel(idime,jdime)+grvel(jdime,idime))
               end do
            end do

            !Now we fill the traction vector
            do inodb = 1,e%pnodb
                ipoin= e%lnodB(inodb)
                a%btraction(1:ndime,ipoin) = a%btraction(1:ndime,ipoin) + tract(1:ndime)*e%shapb(inodb,e%igaub)*dsurf*weightfactor
            end do

            !Heat flux (also for energy conservation restriction) 
            !   Q=q·n 
            !   q=-actco*grad(tempe)
         
            qflx = -actco*grtem

            tract(ndime+1) = dot_product(qflx,e%baloc(:,e%ndime))

            !Forces and moments to global fields
            if(a%npp_stepi(11) /= 0) then
               do inodb=1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  do idime=1,ndime
                   a%forfl(idime,ipoin)=a%forfl(idime,ipoin)+dsurf*tract(idime)*e%shapb(inodb,e%igaub)*weightfactor
                  end do
                  if(ndime==2)then
                     a%momfl(1,ipoin) = a%momfl(1,ipoin) &
                        + tract(2)*(gpcod(1)-a%origm(1))*dsurf*e%shapb(inodb,e%igaub)*weightfactor &
                        - tract(1)*(gpcod(2)-a%origm(2))*dsurf*e%shapb(inodb,e%igaub)*weightfactor 
                  elseif(ndime==3)then
                     a%momfl(1,ipoin) = a%momfl(1,ipoin) &
                        + tract(3)*(gpcod(2)-a%origm(2))*dsurf*e%shapb(inodb,e%igaub)*weightfactor &
                        - tract(2)*(gpcod(3)-a%origm(3))*dsurf*e%shapb(inodb,e%igaub)*weightfactor          
                     a%momfl(2,ipoin) = a%momfl(2,ipoin) &
                        + tract(1)*(gpcod(3)-a%origm(3))*dsurf*e%shapb(inodb,e%igaub)*weightfactor &
                        - tract(3)*(gpcod(1)-a%origm(1))*dsurf*e%shapb(inodb,e%igaub)*weightfactor  
                     a%momfl(3,ipoin) = a%momfl(3,ipoin) &
                        + tract(2)*(gpcod(1)-a%origm(1))*dsurf*e%shapb(inodb,e%igaub)*weightfactor &
                        - tract(1)*(gpcod(2)-a%origm(2))*dsurf*e%shapb(inodb,e%igaub)*weightfactor
                  end if
                  a%hfxfl(ipoin)=a%hfxfl(ipoin)+dsurf*tract(ndime+1)*e%shapb(inodb,e%igaub)*weightfactor
               end do
            end if

            !---------------------------------------------------------------

            tract = tract*dsurf*weightfactor

            !Forces
            if(ibody>0)then

               do idime=1,ndime+1
                  a%force(idime,ibody)=a%force(idime,ibody) + tract(idime)
               end do
               
               !Moments
               if(ndime==2)then
                  a%formo(3,ibody) = a%formo(3,ibody) &
                     + tract(2)*(gpcod(1)-a%origm(1)) &
                     - tract(1)*(gpcod(2)-a%origm(2)) 
               elseif(ndime==3)then
                  a%formo(1,ibody) = a%formo(1,ibody) &
                     + tract(3)*(gpcod(2)-a%origm(2)) &
                     - tract(2)*(gpcod(3)-a%origm(3))          
                  a%formo(2,ibody) = a%formo(2,ibody) &
                     + tract(1)*(gpcod(3)-a%origm(3)) &
                     - tract(3)*(gpcod(1)-a%origm(1))  
                  a%formo(3,ibody) = a%formo(3,ibody) &
                     + tract(2)*(gpcod(1)-a%origm(1)) &
                     - tract(1)*(gpcod(2)-a%origm(2))
               end if

            end if

            
            if(a%kfl_fixbo(iboun)==10) then
               incre = 0.0_rp
               call nsc_pr_nscbc_ex(e,accph,accvh,gpbp(1),gpbv,gpbt(1),grpre,grvel,grtem,a%relpre,a%reltem,a%bvnat(iboun)%a,incre)

               do inodb=1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  borhs(1:a%ndofn,ipoin) = borhs(1:a%ndofn,ipoin) + incre(1:a%ndofn)*e%shapb(inodb,e%igaub)*weightfactor*dsurf
               end do

            elseif(a%kfl_fixbo(iboun)==12) then
               incre = 0.0_rp
               call nsc_pr_nscbc_press_ex(e,a%dtinv,accph,accvh,gpbp(1),gpbv,gpbt(1),gpbv(:,2),grpre,grvel,a%relpre,a%reltem,incre)

               do inodb=1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  borhs(1,ipoin) = borhs(1,ipoin) + incre(1)*e%shapb(inodb,e%igaub)*weightfactor*dsurf
               end do

            end if

         end do 
      
      end if
      
   end do boundaries 

  if( a%kfl_outfm==1) then   
      !Collect Values
      call  MPI_REDUCE(a%force, forceR, (ndime+1)*nbody, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr ) 
      a%force = forceR
      call  MPI_REDUCE(a%formo, momenR, 3*nbody, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr ) 
      a%formo = momenR
  end if   
   
   ! Output forces
   ! Applied to body => Negative sign.
   if(a%kfl_outfm==1) then   
      if (a%MPIrank == a%MPIroot) then  
         
         if(ndime==2) then
            do ibody = 1,nbody
               write(a%lun_force(ibody),20) a%ctime,      &
                  -a%force(1,ibody)/a%adimf(1),    &
                  -a%force(2,ibody)/a%adimf(2),    &
                  -a%formo(3,ibody)/a%adimm(3),    &
                  -a%force(3,ibody)/a%adimf(4)
            end do
         else if(ndime==3) then
            do ibody = 1,nbody
               write(a%lun_force(ibody),30) a%ctime,      &
                  -a%force(1,ibody)/a%adimf(1),    &
                  -a%force(2,ibody)/a%adimf(2),    &
                  -a%force(3,ibody)/a%adimf(3),    &
                  -a%formo(1,ibody)/a%adimm(1),    &
                  -a%formo(2,ibody)/a%adimm(2),    &
                  -a%formo(3,ibody)/a%adimm(3),    &
                  -a%force(4,ibody)/a%adimf(4)
            end do
         end if   
      
      end if
   end if
   
   if (a%kfl_flush == 1) then
      do ibody = 1,nbody
         call flush(a%lun_force(ibody))
      end do
   endif
   
   if( a%kfl_nscnbc==1) then   
      ! Loop over points
      do ipoin = 1,npoin
         tempn0(1:a%ndofn)=a%bvess(1:a%ndofn,ipoin,2)
         a%bvess(1:a%ndofn,ipoin,2)=(4.0_rp*a%bvess(1:a%ndofn,ipoin,1)-tempn0(1:a%ndofn))/3.0_rp
         a%bvess(1:a%ndofn,ipoin,1)=a%bvess(1:a%ndofn,ipoin,2) - a%bvmass(ipoin)*borhs(1:a%ndofn,ipoin)/a%dtinv
      enddo
      !Communicate between subdomains
      call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%bvess(:,:,1))
   end if


   !Dellocation

   call a%Memor%dealloc(ndime+1,tract,'tract','nsc_pr_EndBouope')   
   call a%Memor%dealloc(ndime,gpcod,'gpcod','nsc_pr_EndBouope')  
   call a%Memor%dealloc(e%pnodb,2,bopre,'bopre','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%pnodb,2,bovel,'bovel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%pnodb,2,botem,'botem','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%mnode,elpre,'elpre','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%mnode,eltem,'eltem','nsc_pr_EndBouope')
   call a%Memor%dealloc(2,gpbp,'gpbp','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,2,gpbv,'gpbv','nsc_pr_EndBouope')
   call a%Memor%dealloc(2,gpbt,'gpbt','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,grpre,'grpre','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,grtem,'grtem','nsc_pr_EndBouope')
   call a%Memor%dealloc(a%ndofn,incre,'incre','nsc_pr_EndBouope')   
   call a%Memor%dealloc(e%ndime,qflx,'qflx','nsc_pr_EndBouope')
   !nscnbc(boundary problem)
   call a%Memor%dealloc(a%ndofn,npoin,borhs,'borhs','nsc_pr_EndBouope')
   !working variables (reduced variables)
   call a%Memor%dealloc(ndime+1,nbody,forceR,'forceR','nsc_pr_EndBouope')
   call a%Memor%dealloc(3,nbody,momenR,'momenR','nsc_pr_EndBouope')   
   
   
   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_froces')
   
   !Formats

   12   format( &
          '$        Forces and Moments over bodies          ',/,          &
          '$    Time',15x,'Fx',15x,'Fy',15x,'Mz',15x'Q')
   13   format( &
          '$          Forces and Moments over bodies            ',/,      &
          '$ Time  Fx  Fy  Fz  Mx  My  Mz  Q')

   20   format(1x,e12.6,2x,4(2x,e15.8))
   30   format(1x,e12.6,2x,7(2x,e15.8))

   50   format(a1,i1)
   51   format(a1,i2)
   52   format(a1,i3)
   53   format(5x,'WARNING: CANNOT WRITE MORE THAN 999 bodies')   
   
end subroutine


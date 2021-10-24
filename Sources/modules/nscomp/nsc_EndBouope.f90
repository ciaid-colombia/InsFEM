subroutine nsc_EndBouope(a)
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
   use Mod_NSCompressible
   use Mod_nsc_nscbc_ex
   use MPI
   
   implicit none
   class(NSCompressibleProblem) :: a   
   class(FiniteElement), pointer :: e => NULL() 
 

   real(rp), allocatable :: tract(:),forceR(:,:),momenR(:,:)
   real(rp), allocatable :: gpcod(:)
   real(rp), allocatable :: boden(:), bomom(:,:), boene(:)
   real(rp), allocatable :: elden(:), elmom(:,:), elene(:)
   real(rp), allocatable :: gpmom(:)
   real(rp), allocatable :: grden(:), grmom(:,:), grene(:)
   real(rp), allocatable :: qflx(:), mgrm(:)
   real(rp), allocatable :: incre(:)
   real(rp), allocatable :: bomom0(:,:), gpmom0(:)
   real(rp), allocatable :: boden0(:),boene0(:)
   real(rp), allocatable :: borhs(:,:)

   character(150) :: fil_force
   real(rp)    :: weightfactor
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: aux,invcvh,invgpd
   real(rp)    :: gpden(1),gpene(1),gppre,divm,mgrd
   real(rp)    :: gpden0(1),gpene0(1)
   real(rp)    :: tempn0(a%ndofn)
   integer(ip) :: idime,iboun,ndime,igaub,jdime,nbody,nboun,nelem,ibody
   integer(ip), pointer :: lbody => NULL()
   integer(ip) :: icount,npoinLocal,inodb,ipoin,idofn,ierr,npoin
   
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
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_EndBouope')      
   
   !Allocation

   call a%Memor%alloc(ndime+1,tract,'tract','nsc_EndBouope')   
   call a%Memor%alloc(ndime,gpcod,'gpcod','nsc_EndBouope')
   call a%Memor%alloc(e%mnodb,boden,'boden','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,e%mnodb,bomom,'bomom','nsc_EndBouope')
   call a%Memor%alloc(e%mnodb,boene,'boene','nsc_EndBouope')
   call a%Memor%alloc(e%mnode,elden,'elden','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,e%mnode,elmom,'elmom','nsc_EndBouope')
   call a%Memor%alloc(e%mnode,elene,'elene','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,gpmom,'gpmom','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,grden,'grden','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,e%ndime,grmom,'grmom','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,grene,'grene','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,qflx,'qflx','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,mgrm,'mgrm','nsc_EndBouope')
   call a%Memor%alloc(a%ndofn,incre,'incre','nsc_bouope_im')
   call a%Memor%alloc(e%ndime,e%mnodb,bomom0,'bomom0','nsc_EndBouope')
   call a%Memor%alloc(e%ndime,gpmom0,'gpmom0','nsc_EndBouope')
   call a%Memor%alloc(e%mnodb,boden0,'boden0','nsc_EndBouope')
   call a%Memor%alloc(e%mnodb,boene0,'boene0','nsc_EndBouope')
   !nscnbc(boundary problem)
   call a%Memor%alloc(a%ndofn,npoin,borhs,'borhs','nsc_EndBouope') 
   !working variables (reuced variables)
   call a%Memor%alloc(ndime+1,nbody,forceR,'forceR','nsc_EndBouope')
   call a%Memor%alloc(3,nbody,momenR,'momenR','nsc_EndBouope')
  
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
      
      if((a%kfl_fixbo(iboun)==10) .or. (a%kfl_fixbo(iboun)==12) .or. ((ibody>0).and.(a%kfl_outfm==1)))then
         !Physical Parameters
         call a%GetPhysicalParameters(acvis,actco,accph,accvh)         
         
         call e%elmdel
         
         !Gathers in boundary elements 
         call e%gatherb(1_ip,boden,a%densf(:,1))       
         call e%gatherb(ndime,bomom,a%momen(:,:,1))       
         call e%gatherb(1_ip,boene,a%energ(:,1))       
         call e%gather(1_ip,elden,a%densf(:,1))       
         call e%gather(ndime,elmom,a%momen(:,:,1))      
         call e%gather(1_ip,elene,a%energ(:,1))      
         call e%gatherb(1_ip,boden0,a%densf(:,3))       
         call e%gatherb(ndime,bomom0,a%momen(:,:,3))  
         call e%gatherb(1_ip,boene0,a%energ(:,3))       
        
         dsurf0 = 0.0_rp           

         !Gauss-Point Loop
         do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Traction at gauss point
            tract=0.0_rp
            divm=0.0_rp

            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            dsurf0 = dsurf0 + dsurf
            
            call e%interpb(ndime,e%bocod,gpcod)
            call e%interpb(1_ip,boden,gpden)
            call e%interpb(ndime,bomom,gpmom)
            call e%interpb(1_ip,boene,gpene)
            
            !Derivatives at the boundary
            call e%elmderb         
            call e%gradientb(1_ip,elden,grden)
            call e%gradientb(ndime,elmom,grmom)
            call e%gradientb(1_ip,elene,grene)

            !Momentum divergence
            do idime =1,e%ndime
              divm = divm + grmom(idime,idime)
            end do

            !----------------------------------------------------------------
   
            !Auxiliary variables
            invcvh = 1.0_rp / accvh 
            invgpd = 1.0_rp / gpden(1)
            aux = dot_product(gpmom,gpmom)
            
            !----------------------------------------------------------------
   
            !Traction value
            !  t=n·S  S=-p*I+tau
            if (a%lawde == 1) then
               !  Ideal state law pressure 
               !  p=(gamma-1)*(s- m·m/2rho)

               gppre = (accph*invcvh - 1.0_rp)*(gpene(1) - (aux*invgpd*0.5_rp))
    
            else if (a%lawde /= 1) then
               call runend('Nsc_EndBouope: Non-ideal state law not ready')
            endif

            !  tau=(2mu/rho)*gradS(m)-(2mu/3rho)*div(m)*I
            !      -(2mu/rho^2)*gradS(rho)+(2mu/3rho^2)*m·grad(rho)*I

            mgrd = dot_product(gpmom,grden)

            do idime=1,ndime
                tract(idime)= tract(idime) - e%baloc(idime,ndime)*gppre&
                -(2/3)*acvis*invgpd*divm*e%baloc(idime,ndime)&
                +(2/3)*acvis*invgpd*invgpd*mgrd*e%baloc(idime,ndime)
               do jdime=1,ndime
                  tract(idime) = tract(idime) &
                     + acvis*invgpd*e%baloc(jdime,ndime)*(grmom(idime,jdime)+grmom(jdime,idime))&
                     - acvis*invgpd*invgpd*e%baloc(jdime,ndime)*(gpmom(idime)*grden(jdime)&
                           + gpmom(jdime)*grden(idime))
               end do
            end do

            !Now we fill the traction vector
            do inodb = 1,e%pnodb
                ipoin= e%lnodB(inodb)
                a%btraction(1:ndime,ipoin) = a%btraction(1:ndime,ipoin) + tract(1:ndime)*e%shapb(inodb,e%igaub)*dsurf*weightfactor
            end do

            !Heat flux: 
            !   Q=q·n 
            if (a%lawde == 1) then
               !   Ideal state law heat flux 
               !   q=-actco*(1/rhoCv)*[grad(s)-(1/rho)m·grad(m)^T
               !             +(m·m/rho^2)*grad(rho)-(s/rho)*grad(rho)]
         
               mgrm = matmul(gpmom,grmom)
               qflx = -actco*invgpd*invcvh*(grene - invgpd*mgrm&
                    + (aux*invgpd*invgpd - gpene(1)*invgpd)*grden)

            else if (a%lawde /= 1) then
               call runend('Nsc_EndBouope: Non-ideal state law not ready')
            endif

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
               call nsc_nscbc_ex(e,accph,accvh,gpden(1),gpmom,gpene(1),grden,grmom,grene,a%bvnat(iboun)%a,incre)

               do inodb=1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  borhs(1:a%ndofn,ipoin) = borhs(1:a%ndofn,ipoin) + incre(1:a%ndofn)*e%shapb(inodb,e%igaub)*weightfactor*dsurf
               end do

            elseif(a%kfl_fixbo(iboun)==12) then
               call e%interpb(1_ip,boden0,gpden0)
               call e%interpb(ndime,bomom0,gpmom0)
               call e%interpb(1_ip,boene0,gpene0)
               incre = 0.0_rp
               call nsc_nscbc_density_ex(e,a%dtinv,accph,accvh,gpden(1),gpmom,gpene(1),gpden0(1),gpmom0,gpene0(1),grden,grmom,grene,incre)

               do inodb=1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  borhs(1,ipoin) = borhs(1,ipoin) + incre(1)*e%shapb(inodb,e%igaub)*weightfactor*dsurf
               end do

            end if

         end do 
      
      end if
      
   end do boundaries 

  if(nbody>0 .and. a%kfl_outfm==1) then   
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

   call a%Memor%dealloc(ndime+1,tract,'tract','nsc_EndBouope')   
   call a%Memor%dealloc(ndime,gpcod,'gpcod','nsc_EndBouope')  
   call a%Memor%dealloc(e%mnodb,boden,'boden','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,e%mnodb,bomom,'bomom','nsc_EndBouope')
   call a%Memor%dealloc(e%mnodb,boene,'boene','nsc_EndBouope')
   call a%Memor%dealloc(e%mnode,elden,'elden','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,e%mnode,elmom,'elmom','nsc_EndBouope')
   call a%Memor%dealloc(e%mnode,elene,'elene','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,gpmom,'gpmom','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,grden,'grden','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,e%ndime,grmom,'grmom','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,grene,'grene','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,qflx,'qflx','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,mgrm,'mgrm','nsc_EndBouope')
   call a%Memor%dealloc(a%ndofn,incre,'incre','nsc_EndBouope')   
   call a%Memor%dealloc(e%ndime,e%mnodb,bomom0,'bomom0','nsc_EndBouope')
   call a%Memor%dealloc(e%ndime,gpmom0,'gpmom0','nsc_EndBouope')
   call a%Memor%dealloc(e%mnodb,boden0,'boden0','nsc_EndBouope')
   call a%Memor%dealloc(e%mnodb,boene0,'boene0','nsc_EndBouope')
   !nscnbc(boundary problem)
   call a%Memor%dealloc(a%ndofn,npoin,borhs,'borhs','nsc_EndBouope') 
   !working variables (reduced variables)
   call a%Memor%dealloc(ndime+1,nbody,forceR,'forceR','nsc_EndBouope')
   call a%Memor%dealloc(3,nbody,momenR,'momenR','nsc_EndBouope')   
   
   
   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsc_EndBouope')
   
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


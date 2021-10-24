subroutine ale_updbcs(a)
   !-----------------------------------------------------------------------
   !> This routine updates the free surface boundary condition of the mesh 
   !! at the beginning of the time step by integrating the Navier-Stokes 
   !! velocity. 
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Alemov
   implicit none
   class(AlemovProblem),target :: a
   type(TimeIntegratorDt1) :: Integrator

   real(rp) :: LHSDtinv,valRHS(1)
   integer(ip) :: ipoin,npoin,idime,nsteps,i,ndime
   
   real(rp) :: auxbvess(10),w_r
   
   real(rp) :: center(3), w(3), waux(3)
   real(rp) :: coord(3),vector(3),alphas(3),rotX(3,3), rotY(3,3), rotZ(3,3),coordnew(3)
   real(rp), pointer :: coordpointer(:) => NULL()
   
   integer(ip) :: funno, funty(2)
   
   
   
   call a%Mesh%GetNpoin(npoin)
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)

   if (a%istep < nsteps) then
      nsteps = 2
   end if

   if (a%istep > 1) then
       !Update bvess with NSVeloc only if it is allocated (FMALE and Free Surface)
       if (associated(a%NSveloc)) then
           do ipoin=1,npoin
               do idime=1,a%ndofbc
                   if(a%kfl_fixno(idime,ipoin)<=0 .or. a%kfl_fixno(idime,ipoin)==3) then
                       auxbvess(1:nsteps-1) = a%displacement(idime,ipoin,1:nsteps-1)
                       call Integrator%GetRHS(1,auxbvess,valRHS)
                       a%bvess(idime,ipoin,1)=(1/LHSDtinv)*a%NSveloc(idime,ipoin)+valRHS(1)*a%dtinv/LHSDtinv
                   end if
               end do
           end do  
       endif

       !call a%FilePostpr%postpr(a%bvess(:,:,1),'ToEnforceALEDisp',a%istep,a%ctime,a%Mesh)

   endif

   !For solid interaction
   !bvess is updated in sldale_begite

   !Rotation boundary conditions
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   do ipoin = 1,npoin
      if(a%kfl_fixno(1,ipoin)==1) then
         funno = a%kfl_funno(ipoin)
         if(funno /= 0) then
            funty = a%kfl_funty(funno,:)
            if (funty(1) == -1) then
               
               center = a%funpa(funno)%a(1:3)
               w = a%funpa(funno)%a(4:6)
               coord = 0;
               call a%Mesh%GetPointCoord(ipoin,coordpointer)
               coord(1:ndime) = coordpointer
               
               vector = coord-center
               
               !waux = 0.0_rp
               !waux(3) = w(3)*cos(40.0*a%ctime) 
               !alphas = waux*a%ctime
               alphas = w*a%ctime
               
               rotX = transpose(reshape((/ 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, cos(alphas(1)), -sin(alphas(1)), 0.0_rp, sin(alphas(1)), cos(alphas(1)) /), shape(rotX)))
               rotY = transpose(reshape((/ cos(alphas(2)), 0.0_rp, sin(alphas(2)), 0.0_rp, 1.0_rp, 0.0_rp, -sin(alphas(2)), 0.0_rp, cos(alphas(2)) /), shape(rotY)))
               rotZ = transpose(reshape((/ cos(alphas(3)), -sin(alphas(3)), 0.0_rp, sin(alphas(3)), cos(alphas(3)), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp /), shape(rotZ)))
               
               vector = matmul(rotX,matmul(rotY,matmul(rotZ,vector)))
               coordNew = center+vector
               
               a%bvess(:,ipoin,1) = CoordNew(1:ndime) - coord(1:ndime)
            endif
         end if
      end if
   end do
  
end subroutine ale_updbcs

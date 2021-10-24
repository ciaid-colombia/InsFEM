subroutine nsi_begste(a)
   ! NAME 
   !    nsi_begste
   ! DESCRIPTION
   !    This routine prepares for a new time step of the incompressible NS
   !    equations.
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   integer(ip) :: icomp,ndime,npoin
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy
   
   interface
      subroutine nsi_coupling_outerr(a)
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
   end interface
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)   
   
   !!Check if the coupling between modules is done
   call nsi_coupling_outerr(a)
   
   if (a%kfl_timei==1) then
      call Integrator%Init(a%kfl_tsche_1st_current)
      call Integrator%GetNumberOfTimeSteps(nsteps)
      call Integrator%GetAccuracy(Accuracy)
      
      if (Accuracy == 1) then
         a%veloc(:,:,2) = a%veloc(:,:,3)
      !Higher order schemes
      elseif (Accuracy <= 3 .and. a%ncomp >= Accuracy + 2) then
         call GetTimeProjection(Accuracy,ndime,npoin,a%ncomp,a%veloc,a%veloc(:,:,2)) 
       !Default
      else
         a%veloc(:,:,2) = a%veloc(:,:,3)
      endif   
      
      !Set the first one
      a%veloc(:,:,1) = a%veloc(:,:,2)

   end if
   
   !Check if flow is confined (only if non-constant boundary conditions)
   if (a%kfl_conbc /= 1) call a%Ifconf
   
   ! Update frame of reference
   !call nsi_updfor
   
   !We put the initial residual projection to zero
   if (a%kfl_StabilizeFreeSurface == 1) then
      !a%StabilizeFreeSurface_VelocityGradient = 0.0_rp
      !a%StabilizeFreeSurface_PressureGradientAndForce = 0.0_rp
   endif
   
end subroutine nsi_begste

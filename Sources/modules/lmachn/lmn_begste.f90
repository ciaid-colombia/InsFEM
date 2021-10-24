subroutine lmn_begste(a)
   ! DESCRIPTION
   !    This routine prepares for a new time step of the Low Mach number   
   !    equations.
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator     
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy
   integer(ip) :: ndime,npoin

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)   
   
   if (a%kfl_timei==1) then
      call Integrator%Init(a%kfl_tsche_1st_current)
      call Integrator%GetNumberOfTimeSteps(nsteps)
      call Integrator%GetAccuracy(Accuracy)
      
      if (Accuracy == 1) then
         a%veloc(:,:,2) = a%veloc(:,:,3)
         a%tempe(:,2) = a%tempe(:,3)
         a%press(:,2) = a%press(:,3)
         a%pther(2) = a%pther(3)
      !Higher order schemes
      !elseif (Accuracy <= 3 .and. a%ncomp >= Accuracy + 2) then
         !call GetTimeProjection(Accuracy,ndime,npoin,a%ncomp,a%veloc,a%veloc(:,:,2)) 
         !call GetTimeProjection(Accuracy,1,npoin,a%ncomp,a%tempe,a%tempe(:,2)) 
         !call GetTimeProjection(Accuracy,1,npoin,a%ncomp,a%press,a%press(:,2)) 
         !call GetTimeProjection(Accuracy,1,1,a%ncomp,a%pther,a%pther(:,2)) 
      !Default
      else
         a%veloc(:,:,2) = a%veloc(:,:,3)
         a%tempe(:,2) = a%tempe(:,3)
         a%press(:,2) = a%press(:,3)
         a%pther(2) = a%pther(3)
      endif   
      
      !Set the first one
      a%veloc(:,:,1) = a%veloc(:,:,2)
      a%tempe(:,1) = a%tempe(:,2)
      a%press(:,1) = a%press(:,2)
      a%pther(1) = a%pther(2)

   end if
   
   !Check if flow is confined (only if non-constant boundary conditions)
   if (a%kfl_conbc /= 1) call a%Ifconf
   
end subroutine lmn_begste


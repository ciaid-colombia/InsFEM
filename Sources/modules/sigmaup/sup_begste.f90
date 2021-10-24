subroutine sup_begste(a)
   ! NAME 
   !    sup_begste
   ! DESCRIPTION
   !    This routine prepares for a new time step of the incompressible NS
   !    equations.
   !-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   use Mod_Mesh   
   use Mod_TimeIntegrator   
   implicit none  
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: icomp,ndime,ntens,npoin
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy   
   
   interface
      subroutine sup_coupling_outerr(a)
         use Mod_ThreeField
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine
   end interface
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   !Check if the coupling between modules is done
   call sup_coupling_outerr(a)
   
   
   ntens=(ndime-1)*(ndime-1)+2
   
   if(a%kfl_timei==1) then
      call Integrator%Init(a%kfl_tsche_1st_current)
      call Integrator%GetNumberOfTimeSteps(nsteps)
      call Integrator%GetAccuracy(Accuracy)
      
      if (Accuracy == 1) then
         a%veloc(:,:,2) = a%veloc(:,:,3)
         a%sigma(:,:,2) = a%sigma(:,:,3)
      !Higher order schemes
      elseif (Accuracy <= 3 .and. a%ncomp >= Accuracy + 2) then
         call GetTimeProjection(Accuracy,ndime,npoin,a%ncomp,a%veloc,a%veloc(:,:,2))         
         call GetTimeProjection(Accuracy,ntens,npoin,a%ncomp,a%sigma,a%sigma(:,:,2))         
       !Default
      else
         a%veloc(:,:,2) = a%veloc(:,:,3)
         a%sigma(:,:,2) = a%sigma(:,:,3)
      endif   
      
      !Set the first one
      a%veloc(:,:,1) = a%veloc(:,:,2)
      a%sigma(:,:,1) = a%sigma(:,:,2)

   end if   
   
   
   
   

   
   !Check if flow is confined (only if non-constant boundary conditions)
   if (a%kfl_conbc /= 1) call a%Ifconf     
   
   
end subroutine sup_begste

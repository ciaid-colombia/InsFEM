subroutine supf_begite(a)
!-----------------------------------------------------------------------
! NAME 
!    sup_begite
! DESCRIPTION
!    This routine starts an internal iteration for the incompressible NS
!    equations. 
!-----------------------------------------------------------------------
   use typre
   use Mod_SUPFractionalStep  
   implicit none
   class(SUPFractionalStepProblem) :: a
   
   !Initializations
   a%kfl_goiteS = 1
   a%iteraS     = 0
   
   a%kfl_goiteY = 1
   a%iteraY     = 0   
   
   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   a%sigma(:,:,1) = a%sigma(:,:,2)
   a%veloc(:,:,1) = a%veloc(:,:,2)
   a%press(:,1) = a%press(:,2)

end subroutine supf_begite

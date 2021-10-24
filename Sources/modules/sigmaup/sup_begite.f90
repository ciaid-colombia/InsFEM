subroutine sup_begite(a)
!-----------------------------------------------------------------------
! NAME 
!    sup_begite
! DESCRIPTION
!    This routine starts an internal iteration for the incompressible NS
!    equations. 
!-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   
   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   a%sigma(:,:,1) = a%sigma(:,:,2)
   a%veloc(:,:,1) = a%veloc(:,:,2)
   a%press(:,1) = a%press(:,2)

end subroutine sup_begite

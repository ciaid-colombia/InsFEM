subroutine nsi_begite(a)
!-----------------------------------------------------------------------
! NAME 
!    nsi_begite
! DESCRIPTION
!    This routine starts an internal iteration for the incompressible NS
!    equations. 
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   use Mod_nsi_BouwalStats
   implicit none
   class(NavierStokesProblem) :: a
   
   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   a%veloc(:,:,1) = a%veloc(:,:,2)
   a%press(:,1) = a%press(:,2)
   
   !If wall law, update statistics
   call nsi_BouwalStats(a,0)

end subroutine nsi_begite

subroutine lmn_begite(a)
!-----------------------------------------------------------------------
! NAME 
!    lmn_begite
! DESCRIPTION
!    This routine starts an internal iteration for the incompressible NS
!    equations. 
!-----------------------------------------------------------------------
   use typre
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
  
   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   !velocity, pressure and temperature
   a%veloc(:,:,1) = a%veloc(:,:,2)
   a%press(:,1) = a%press(:,2)
   a%tempe(:,1) = a%tempe(:,2)
   a%pther(1) = a%pther(2)
   a%rilmn = 1.0_rp
   a%rhsnorm = 1.0_rp
   a%kfl_acite = 1
   a%kfl_linop = 1

end subroutine lmn_begite

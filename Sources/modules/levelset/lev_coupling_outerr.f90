subroutine lev_coupling_outerr(a)
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem) :: a

!    !Check Advection from outerr Physical problem has not been set
!    if(a%kfl_advec==1.and.(a%kfl_ExternalAdvection==0)) &
!       call runend('EXTERNAL CONVECTION IN LEVELSET IS IMPOSSIBLE IF EXTERNAL VELOCITY IS NOT SET')
! 

end subroutine
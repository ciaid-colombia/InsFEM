subroutine lev_CrankNicolsonEndste(a)
   use typre
   use Mod_LevelSet
   use def_parame
   implicit none
   class(LevelSetProblem) :: a
  
   a%level(:,1) = 2.0_rp*a%level(:,1)-a%level(:,3)
   !write(*,*) 'tem_CrankNicolsonEndste: executing'
end subroutine

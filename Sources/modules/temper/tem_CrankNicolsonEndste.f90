subroutine tem_CrankNicolsonEndste(a)
   use typre
   use Mod_Temperature
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
  
   a%tempe(:,1) = 2.0_rp*a%tempe(:,1)-a%tempe(:,3)
   !write(*,*) 'tem_CrankNicolsonEndste: executing'
end subroutine

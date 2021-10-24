subroutine tem_updbcs(a)
   use typre
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   
   call runend('tem_exacso: exact solution not ready for mod_temperature')   
end subroutine tem_updbcs



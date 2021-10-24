subroutine tem_turnon(a)
   use Mod_Temperature
   use Mod_tem_DustTransport
   implicit none
   class(TemperatureProblem) :: a
   
   if (a%DT_kfl_DustTransport == 1) then
      call ComputeSettlingVelocity(a)
   endif
end subroutine
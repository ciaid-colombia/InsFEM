subroutine tem_coupling_outerr(a)
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a

   !Check Advection from outerr Physical problem has not been set
   if(a%kfl_advec==1.and.(a%kfl_ExternalAdvection==0)) &
      call runend('EXTERNAL CONVECTION IN TEMPER IS IMPOSSIBLE IF EXTERNAL VELOCITY IS NOT SET')


end subroutine
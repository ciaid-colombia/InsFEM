subroutine opt_coupling_outerr(a)
   use Mod_Optics
   implicit none
   class(OpticsProblem) :: a

   !Check Advection from outerr Physical problem has not been set
   if(a%kfl_velocity==1.and.(a%kfl_Externalvelocity==0)) &
      call runend('EXTERNAL VELOCITY IN OPTICS IS IMPOSSIBLE IF EXTERNAL VELOCITY IS NOT SET')
   
   if(a%kfl_pressure==1.and.(a%kfl_ExternalPressure==0)) &
      call runend('EXTERNAL PRESSURE IN OPTICS IS IMPOSSIBLE IF EXTERNAL VELOCITY IS NOT SET')
   
   if(a%kfl_temperature==1.and.(a%kfl_ExternalTemperature==0)) &
      call runend('EXTERNALTEMPERATURE IN OPTICS IS IMPOSSIBLE IF EXTERNAL TEMPERATURE IS NOT SET')
   

end subroutine
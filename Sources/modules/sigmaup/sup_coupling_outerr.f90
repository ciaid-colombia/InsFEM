subroutine sup_coupling_outerr(a)
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
    
    !Check Temperature from outerr Physical problem has not been set
   if(a%kfl_cotem==1 .and.a%kfl_ExternalTemperature==0) &
      call runend('EXTERNAL TEMPERATURE IN NSTINC IS IMPOSSIBLE IF EXTERNAL TEMPERATURE IS NOT SET')
end subroutine

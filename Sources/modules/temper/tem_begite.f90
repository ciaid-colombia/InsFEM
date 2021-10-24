subroutine tem_begite(a)
   use typre
   use Mod_Temperature
   use def_parame
   use Mod_SmoothedFieldGradient
   use Mod_Memor
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: ndime

   if (a%kfl_CouplingThreeField==1) then
      call a%Mesh%GetNdime(ndime)
      call ComputeSmoothedFieldGradient(a%Mesh,a%Memor,ndime,a%veloc,a%SmoothedVelocityGradient) 
   endif

end subroutine

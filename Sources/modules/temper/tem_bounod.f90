subroutine tem_bounod(a,itask)
   use typre
   use Mod_Temperature
   implicit none

   class(TemperatureProblem) :: a
   integer(ip) :: itask

   !Nothing to be done
   !Originally the subroutine passed Neumann and Robin Conditions on nodes to 
   !conditions on boundaries, but the correct thing to do is define this conditions
   !as conditions on boundaries directly
   
end subroutine

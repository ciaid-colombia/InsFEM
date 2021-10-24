subroutine opt_SetExmod(a)
   use typre
   use Mod_NavierStokes
   implicit none
   
   class(NavierStokesProblem) :: a
   
   a%exmod = 'opt'
   a%namod = 'OPTICS'
   
end subroutine


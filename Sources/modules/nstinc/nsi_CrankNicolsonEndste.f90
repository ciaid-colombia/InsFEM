subroutine nsi_CrankNicolsonEndste(a)
   !This routine ends a time step of the incoma%pressible NS equations.
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   a%veloc(:,:,1) = 2.0_rp*a%veloc(:,:,1)-a%veloc(:,:,3)
   
end subroutine

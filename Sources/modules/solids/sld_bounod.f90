subroutine sld_bounod(a,itask)
   !This subroutine passes boundary conditions from boundary to nodes
   use typre
   use Mod_Solids
   implicit none

   class(SolidsProblem) :: a
   integer(ip) :: itask
   
   if (itask == 1) then
      a%kfl_fixrs(a%gipoin) = a%kfl_bours(a%giboun)
      
   endif

end subroutine

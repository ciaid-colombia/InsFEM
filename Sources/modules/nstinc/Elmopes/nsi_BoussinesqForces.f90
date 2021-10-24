module Mod_nsi_BoussinesqForces
contains
   subroutine nsi_BoussinesqForces(e,acden,gravi,boube,boutr,gptem,elext)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      real(rp) :: acden, grnor, gravi(e%ndime), boube,boutr,gptem,elext(e%ndime)
   
        elext(1:e%ndime) = elext(1:e%ndime) - acden*boube*(gptem-boutr)*gravi(1:e%ndime)
   end subroutine
end module

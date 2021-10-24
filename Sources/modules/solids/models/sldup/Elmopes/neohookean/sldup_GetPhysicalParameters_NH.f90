subroutine sldup_GetPhysicalParameters_NH(a,ndime,sz,densi,lam,G,ielem)
   use typre
   use Mod_UPSolids_NH
   implicit none
   class(UPSolidsProblem_NH), target :: a
   integer(ip),intent(in) :: ndime,sz
   real(rp),   intent(out):: densi,G,lam
   integer(ip),optional   :: ielem
   
   densi = a%densi
   lam   = a%lambda
   G     = a%mu
end subroutine

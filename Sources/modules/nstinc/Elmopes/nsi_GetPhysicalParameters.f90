subroutine nsi_GetPhysicalParameters(a,imat,acden,acvis)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   real(rp)              :: acden, acvis
   integer(ip)           :: imat 
   
   acden = a%MatProp(imat)%densi
   acvis = a%MatProp(imat)%visco
   
end subroutine
   
subroutine GetPhysicalParametersPointers(a,imat,acden,acvis)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem), target :: a
   real(rp), pointer     :: acden, acvis
   integer(ip)           :: imat 
   
   acden => a%MatProp(imat)%densi
   acvis => a%MatProp(imat)%visco
   
end subroutine

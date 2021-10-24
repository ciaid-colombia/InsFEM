subroutine tem_GetPhysicalParameters(a,ielem,acden,acsph,actco,acrea,acsou)
   use typre
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: ielem, auxmaterial = 1
   real(rp) :: acden, acsph,actco,acrea,acsou
   
   if (a%NumberOfMaterials > 1) auxmaterial = a%ElementMaterials(ielem)
   call a%Materials(auxmaterial)%GetDensity(acden)
   call a%Materials(auxmaterial)%GetSpecificHeat(acsph)
   call a%Materials(auxmaterial)%GetThermalConductivity(actco)
   
   acrea = a%react
   acsou = a%sourc(1)
   
   
end subroutine

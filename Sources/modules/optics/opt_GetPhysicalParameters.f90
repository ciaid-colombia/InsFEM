subroutine opt_GetPhysicalParameters(a,acden,acvis,acsph,actco)
   use typre
   use Mod_Optics
   implicit none
   class(OpticsProblem) :: a
   real(rp) :: acden, acsph,actco,acvis
   
   acden = a%densi
   acsph = a%sphea
   actco = a%tcond
   acvis = a%visco
   
   
end subroutine
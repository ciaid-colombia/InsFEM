subroutine nsc_GetPhysicalParameters(a,acvis,actco,accph,accvh)
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   real(rp) :: acvis, actco, accph, accvh
   
   acvis = a%visco
   actco = a%tcond
   accph = a%cphea
   accvh = a%cvhea
end subroutine
   

subroutine plcd_Getste(a,dtinv)
   use typre
   use Mod_PLCD
   use Mod_plcd_Stages
   implicit none
   class(PLCDProblem), target :: a
   real(rp) :: dtinv

   
   !If arclength, this should be implemented here

   dtinv = 1.0_rp/a%css%TimeStep

end subroutine
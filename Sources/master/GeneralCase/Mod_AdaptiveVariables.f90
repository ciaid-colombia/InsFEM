module Mod_AdaptiveVariables
   use typre
   use Mod_DCHashCharSize
   use Mod_AdaptiveInterface
   use Mod_IsotropicAdaptiveRefiner
   use Mod_Timer
   implicit none

   type AdaptiveVariables
      !Adaptive Mesh Refinement
      integer(ip) :: kfl_AdaptiveRefinement = 0           !Default is no adaptive refinement
      integer(ip) :: kfl_AdaptiveMulticomm = 0            !Default is no multicomm
         
      character(5) :: RefinementLeader = 'TEMPE'          !Default is temperature leads

      class(AdaptiveRefinerInterface), pointer :: Refiner => NULL()
      integer(ip), allocatable                 :: RefinerMarkel(:)
      real(rp)    :: RefinerRebalancingRatio = 1.5
      integer(ip) :: kfl_21Balancing = 0
      integer(ip) :: RefinerMaxLevel = 0
      integer(ip) :: ProgressiveMaxLevel = 0
      integer(ip) :: NumberOfInitialUniformRefinementSteps
      integer(ip) :: AdaptiveDelay = -1   !Number of delay steps before starting adaptive

      logical     :: UpdateIO = .false.   !Update IO
      integer(ip) :: UpdateFreq = -1      !Update IO frequency

      type(Timer) :: cpu_adaptive
      integer(ip) :: kfl_perio = 0      !To deactivate and reactivate periodic boundary conditions when initial refinement
   end type

end module

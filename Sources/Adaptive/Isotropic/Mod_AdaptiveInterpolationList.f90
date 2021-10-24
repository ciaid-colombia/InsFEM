module Mod_AdaptiveInterpolationList
   use typre
   implicit none

   type :: InterpolationList
      integer(ip) :: n
      integer(ip) :: childElementList(10)
      integer(ip) :: childObjectiveNodeList(10)
      integer(ip) :: nInterpolatingNodes(10)
      integer(ip) :: InterpolatingNodesList(10,10)
      real(rp)    :: InterpolatingCoefficientsList(10,10)
   end type

end module
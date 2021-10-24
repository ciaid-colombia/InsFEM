module Mod_Debugging
   use typre
   use Mod_timer
   implicit none
   
   type, extends(Timer) :: AuxTimer
   
      integer(ip) :: kfl_isSet = 0

contains      
   
      procedure :: Tic =>auxTic
      procedure :: IsSet
   
   end type
   
   type(AuxTimer) :: AuxTimers(200)
   integer(ip) :: deb_PostprocessMatrix = 0

contains
   subroutine auxTic(a)
      implicit none
      class(AuxTimer) :: a
      
      a%kfl_isSet = 1
      call a%Timer%Tic
   end subroutine
   
   subroutine IsSet(a,kfl_isSet)
      implicit none
      class(AuxTimer) :: a
      integer(ip) :: kfl_isSet
      
      kfl_isSet = a%kfl_isSet
   end subroutine
end module

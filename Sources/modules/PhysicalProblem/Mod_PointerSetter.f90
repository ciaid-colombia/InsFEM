module Mod_PointerSetter
   use typre
   implicit none
   
   type, abstract :: PointerSetter
   
      integer(ip) :: kfl_isSet = -2
   
contains
      procedure :: Initialize
      procedure :: Set
      procedure(SpecificSet), deferred :: SpecificSet
      procedure :: Finalize
   end type

   abstract interface
      subroutine SpecificSet(d)
         import PointerSetter
         implicit none
         class(PointerSetter) :: d
      end subroutine
   end interface
   
contains
   subroutine Initialize(a)
      class(PointerSetter) :: a
      if (a%kfl_isSet /= -2) call runend('Pevious PointerSetter was not finalized')
      a%kfl_isSet = -1
   end subroutine
   
   subroutine Set(a)
      class(PointerSetter) :: a
      if (a%kfl_isSet == -2) then
         call runend('PointerSetter was not initialized')
      elseif (a%kfl_isSet == -1) then
         a%kfl_isSet = 1
         call a%SpecificSet
      endif
   end subroutine
   
   subroutine Finalize(a)
      class(PointerSetter) :: a
      a%kfl_isSet = -2
   end subroutine
   
end module

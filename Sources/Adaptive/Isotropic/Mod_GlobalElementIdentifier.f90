module Mod_GlobalElementIdentifier
   use typre
   implicit none
   private
   public GlobalElementIdentifier, GEI_SetupOriginal, GEI_SetupFromParent, GEI_GetPositionInElementForLevel, GEI_ElementMaxLevels,ListGEI
   
   !This is not going to be a class so that we can make arrays of it
   !No possibility to extend etc.

   !This structure allows to identify refined elements Globally
   integer, parameter :: GEI_ElementBitsPerLevel = 3
   integer, parameter :: GEI_ElementMaxLevels = 21  !   21*3 = 63, position can be stored in an integer(8)
   type :: GlobalElementIdentifier
      integer(ip) :: GlobalOriginalElement = -1
      integer(1)  :: Level = 0
      integer(8)  :: PositionInElement = 0
   end type
   type :: ListGEI
      type(GlobalElementIdentifier), allocatable :: a(:)
   end type

contains

   subroutine GEI_SetupOriginal(e,GlobalOriginalElement)
      implicit none
      type(GlobalElementIdentifier) :: e
      integer(ip) :: GlobalOriginalElement
      
      e%GlobalOriginalElement = GlobalOriginalElement
      e%Level = 0
   end subroutine
   
   subroutine GEI_SetupFromParent(e,eParent,LocalPositionInElement)
      implicit none
      type(GlobalElementIdentifier) :: e, eParent
      integer(ip) :: LocalPositionInElement
      
      integer(8) :: auxLocalPositionInElement
      integer(1) :: ibit, baselen
      
      e%GlobalOriginalElement = eParent%GlobalOriginalElement
      e%Level = eParent%Level+1
      if (e%Level > GEI_ElementMaxLevels) call runend('Maximum number of element levels reached in GlobalElementIdentifier')
      
      baselen = eParent%Level*GEI_ElementBitsPerLevel
      !Copy the position in element from parent
      e%PositionInElement = IBITS(eParent%PositionInElement,0,baselen)
      
      auxLocalPositionInElement = LocalPositionInElement-1
      !Add the bits for the new position in element
      call MVBITS(auxLocalPositionInElement, 0, GEI_ElementBitsPerLevel, e%PositionInElement, baselen) 
   end subroutine  
   
   subroutine GEI_GetPositionInElementForLevel(e,ilevel,ichild)
      implicit none
      type(GlobalElementIdentifier) :: e
      integer(ip) :: ilevel,ichild
      
      ichild = IBITS(e%PositionInElement,GEI_ElementBitsPerLevel*(ilevel-1),GEI_ElementBitsPerLevel)+1
   end subroutine
   
   subroutine GEI_Compare(e1,e2,AreEqual)
      implicit none
      type(GlobalElementIdentifier) e1,e2
      logical :: AreEqual
      
      AreEqual = .false.
      
      if (e1%GlobalOriginalElement /= e2%GlobalOriginalElement) return
      if (e1%Level /= e2%Level) return
      if (e1%PositionInElement /= e2%PositionInElement) return
      
      AreEqual = .true.
   end subroutine
      
end module
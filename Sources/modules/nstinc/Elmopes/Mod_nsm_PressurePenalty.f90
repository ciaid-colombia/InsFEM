module Mod_nsm_PressurePenalty
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use typre
   implicit none
   private
   public SetPointersPressurePenalty 
   
   !For elmopes
   type, extends(PointerSetter) :: SPPressurePenalty
contains
      procedure :: SpecificSet => SpecificSetPP
   end type
   type(SPPressurePenalty) :: SetPointersPressurePenalty

contains
   subroutine SpecificSetPP(d)
      implicit none
      class(SPPressurePenalty) :: d
      
      !PenaltyTerm for the pressure
      if (a%kfl_penal == 1) call ConcatenateProcedures(ProcHook_InGaussElmats,PenaltyPressure)
   end subroutine
   
   !PenaltyTerm
   subroutine PenaltyPressure
      implicit none
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elmpq(1,inode,1,inode) = elmpq(1,inode,1,inode) + a%penal*dvol
      enddo
   end subroutine
   
end module

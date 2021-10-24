module Mod_plcd_HangingNodes
   use typre
   use Mod_plcd_BaseElmope
   implicit none
   private
   public SetPointersHangingNodes
   
   type, extends(PointerSetter) :: SPHangingNodes
contains
      procedure :: SpecificSet => SpecificSetHangingNodes
   end type
   type(SPHangingNodes) :: SetPointersHangingNodes
   
   
   integer(ip) :: kfl_HangingNodes

   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetHangingNodes(d)
      implicit none
      class(SPHangingNodes) :: d
   
      !-------------------------------------------------------
      !Hanging Nodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) then
         call ConcatenateProcedures(ProcHook%PreDirichlet,ModifyMatricesHanging)

      endif

   end subroutine   
   
   
   !For actual computations
   !------------------------------------------------------
   !Hanging nodes
   subroutine ModifyMatricesHanging
      implicit none
      
      call a%Mesh%PrepareHangingMatrices(e,elmat,elrhs)
   end subroutine



   
end module

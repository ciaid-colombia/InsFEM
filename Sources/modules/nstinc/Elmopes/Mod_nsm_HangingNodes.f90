module Mod_nsm_HangingNodes
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersHangingNodes,ModifyMatricesHanging
   
   type, extends(PointerSetter) :: SPHangingNodes
contains
      procedure :: SpecificSet => SpecificSetHangingNodes
   end type
   type(SPHangingNodes) :: SetPointersHangingNodes

   class(FiniteElement) , pointer :: whe => NULL()
   real(rp), allocatable :: whelmat(:,:,:,:), whelrhs(:,:)
   integer(ip) :: kfl_HangingNodes
   integer(ip) :: HangingCount
   
contains
   
   subroutine SpecificSetHangingNodes(d)
      implicit none
      class(SPHangingNodes) :: d
         
      !Hanging Nodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) then
         call PrependProcedure(ProcHook_PreDirichlet,ModifyMatricesHanging)
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

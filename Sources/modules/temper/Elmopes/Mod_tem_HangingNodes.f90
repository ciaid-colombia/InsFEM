module Mod_tem_HangingNodes
   use typre
   use mod_tem_BaseElmope
   implicit none
   private
   public SetPointersHangingNodes

   integer(ip), allocatable :: kfl_IsSet
   
   integer(ip) :: kfl_HangingNodes


   
contains
   
   !Set Pointers
   subroutine SetPointersHangingNodes(itask)
      implicit none
      integer(ip) :: itask
      
      integer(ip) :: kfl_nonlinear
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            call a%Mesh%GetHanging(kfl_HangingNodes)
            if (kfl_HangingNodes == 1) then
               call ConcatenateProcedures(ProcHook%PreDirichlet,ModifyMatricesHanging)

            endif 
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
    
   !------------------------------------------------------
   !Hanging nodes
   subroutine ModifyMatricesHanging
      implicit none
      
      call a%Mesh%PrepareHangingMatrices(e,elmat,elrhs)
   end subroutine


end module
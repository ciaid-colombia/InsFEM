module Mod_nsc_pr_NonLinearDerivatives
   use typre
   use mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersNonLinearDerivatives

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !Set Pointers
   subroutine SetPointersNonLinearDerivatives(itask)
      implicit none
      integer(ip) :: itask
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            call ConcatenateProcedures(ProcHook_nsc_pr_InGauss,InGaussNonLinearDerivatives)

         endif 
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
   !Computation Subroutine
   subroutine InGaussNonLinearDerivatives
      implicit none
      
      call e%elmder
      call e%elmhes
   end subroutine

end module

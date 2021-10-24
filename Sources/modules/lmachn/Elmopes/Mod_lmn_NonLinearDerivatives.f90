module Mod_lmn_NonLinearDerivatives
   use typre
   use mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersNonLinearDerivatives

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !Set Pointers
   subroutine SetPointersNonLinearDerivatives(itask)
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
         
            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) then
               call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinearDerivatives)
            endif
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

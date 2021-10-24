 module Mod_supm_InGaussNonLinear
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersInGaussNonLinearDerivatives
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersInGaussNonLinearDerivatives(itask)
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
            call ConcatenateProcedures(ProcHook_InGauss,InGaussNonLinearDerivatives)
         endif  
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES     
   subroutine InGaussNonLinearDerivatives
      implicit none
      call e%elmder 
      call e%elmhes
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
   end subroutine
   
end module

module Mod_supm_ComputeTestf
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersComputeTestf

   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeTestf(itask)
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
            
            !Allocations
            call ConcatenateProcedures(ProcHook_Initializations, AllocTestf)
            call ConcatenateProcedures(ProcHook_Finalizations, DeallocTestf)
         
            ProcHook_ComputeTestf => ComputeTestf
            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) then
               call ConcatenateProcedures(ProcHook_ComputeTestf,ComputeTestfNonlinear)
            endif
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute Testf values
   subroutine AllocTestf
      implicit none
      
      call a%Memor%alloc(e%mnode,testf,'testf','supm_EnditeElmope')
   end subroutine
   
   subroutine DeallocTestf
      implicit none
      
      call a%Memor%dealloc(e%mnode,testf,'testf','supm_EnditeElmope')
   end subroutine
   
   subroutine ComputeTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms : -tau L*v
      call nsm_ComputeTestf(e,acden,timom,AGradV,testf)
   end subroutine
      
   subroutine ComputeTestfNonLinear
      implicit none
      
      call nsm_ComputeTestfNonLinear(e,acvis,timom,a%kfl_stabm,testf)
   end subroutine
   
end module 

module Mod_nsm_ComputeTestf
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersComputeTestf

   type, extends(PointerSetter) :: SPComputeTestf
contains
      procedure :: SpecificSet => SpecificSetComputeTestf
   end type
   type(SPComputeTestf) :: SetPointersComputeTestf
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetComputeTestf(d)
      implicit none
      class(SPComputeTestf) :: d
      integer(ip) :: kfl_nonlinear
 
      !Allocations
      call ConcatenateProcedures(ProcHook_Initializations, AllocTestf)
      call ConcatenateProcedures(ProcHook_Finalizations, DeallocTestf)
      
      ProcHook_ComputeTestf => ComputeTestf
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook_ComputeTestf,ComputeTestfNonlinear)
      endif

   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute Testf values
   subroutine AllocTestf
      implicit none
      
      call a%Memor%alloc(e%mnode,testf,'testf','nsm_EnditeElmope')
   end subroutine
   
   subroutine DeallocTestf
      implicit none
      
      call a%Memor%dealloc(e%mnode,testf,'testf','nsm_EnditeElmope')
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

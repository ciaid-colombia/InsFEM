module Mod_plcd_ComputeSmoothDisplacementGradients
   use typre
   use Mod_plcd_BaseElmope
   implicit none
   private
   public SetPointersComputeSmoothDisplacementGradients
   
   type, extends(PointerSetter) :: SPComputeSmoothDisplacementGradients
contains
      procedure :: SpecificSet => SpecificSetComputeSmoothDisplacementGradients
   end type
   
   type(SPComputeSmoothDisplacementGradients) :: SetPointersComputeSmoothDisplacementGradients
   
   
   real(rp), allocatable, target :: elSmoothGradient(:,:,:),auxGradDisp(:,:)
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetComputeSmoothDisplacementGradients(d)
      implicit none
      class(SPComputeSmoothDisplacementGradients) :: d
      
      !Element size required by EMDS
      if (a%UseSmoothedDisplacementGradient) then
         call ConcatenateProcedures(ProcHook%Initializations,Alloc)
         call ConcatenateProcedures(ProcHook%Finalizations,DeAlloc)
         call ConcatenateProcedures(ProcHook%PreGauss,Gathers)
         
         call ConcatenateProcedures(ProcHook%InGauss,ComputeSmoothDisplacementGradients)
      endif
   end subroutine   
   
   
   !For actual computations
   !------------------------------------------------------
   subroutine Alloc
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elSmoothGradient,'elSmoothGradient','ComputeDisplacementGradients')
      call a%Memor%Alloc(e%ndime,e%ndime,auxGradDisp,'auxGradDisp','plcd_EnditeElmope')
   
   
   end subroutine
   
   subroutine Gathers
      call e%gather(e%ndime*e%ndime,elSmoothGradient,a%SmoothedDisplacementGradient)
   end subroutine
   
   subroutine Dealloc  
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elSmoothGradient,'elSmoothGradient','ComputeDisplacementGradients')
      call a%Memor%dealloc(e%ndime,e%ndime,auxGradDisp,'auxGradDisp','plcd_EnditeElmope')
   end subroutine
   
   
   subroutine ComputeSmoothDisplacementGradients
      implicit none
      
      !Real smooth computation
      call e%interpg(e%ndime*e%ndime,elSmoothGradient,auxGradDisp)
      gradDispHistory => auxGradDisp
   end subroutine
end module

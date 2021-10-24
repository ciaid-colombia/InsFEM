module Mod_nsm_InterpolateResidualProjection
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_ComputeResidualProjection
   implicit none
   private
   public SetPointersInterpolateResidualProjection, gprep
   
   type, extends(PointerSetter) :: SPInterpolateResidualProjection
contains
      procedure :: SpecificSet => SpecificSetInterpolateResidualProjection
   end type
   type(SPInterpolateResidualProjection) :: SetPointersInterpolateResidualProjection
   
   real(rp), allocatable :: elrep(:,:),gprep(:) !The residual projection
   
contains
   
   subroutine SpecificSetInterpolateResidualProjection(d)
      implicit none
      class(SPInterpolateResidualProjection) :: d

      if (a%kfl_repro /= 0) then
         call ConcatenateProcedures(ProcHook_Initializations,AllocRep)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocRep)
         call ConcatenateProcedures(ProcHook_Gathers,GatherRep)
         call ConcatenateProcedures(ProcHook_Interpolates,InterpolateRep)
      endif
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocRep
      implicit none
      
      call a%Memor%alloc(a%ResidualSize,e%mnode,elrep,'elrep','nsm_EnditeElmope')
      call a%Memor%alloc(a%ResidualSize,gprep,'gprep','nsm_EnditeElmope')
   end subroutine
   
   subroutine DeallocRep
      implicit none
      
      call a%Memor%dealloc(a%ResidualSize,e%mnode,elrep,'elrep','nsm_EnditeElmope')
      call a%Memor%dealloc(a%ResidualSize,gprep,'gprep','nsm_EnditeElmope')
   end subroutine
   
   !Residual Interpolation Subroutines
   subroutine GatherRep
      implicit none
      
      call e%gather(a%ResidualSize,elrep,a%repro)
   end subroutine

   subroutine InterpolateRep
      implicit none
      
      !Interpolate
      call e%interpg(a%ResidualSize,elrep,gprep)
   end subroutine
   
end module

module Mod_nsm_SubgridSpaceResidual
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_ComputeGpResidual
   use Mod_nsm_InterpolateResidualProjection
   implicit none
   private
   public SetPointersComputeSubgridSpaceResidual, gpSubscaleSpaceResidual
   
   type, extends(PointerSetter) :: SPComputeSubigridSpaceResidual
contains
      procedure :: SpecificSet => SpecificSetComputeSubgridSpaceResidual
   end type
   type(SPComputeSubigridSpaceResidual) :: SetPointersComputeSubgridSpaceResidual

   real(rp), allocatable :: gpSubscaleSpaceResidual(:)

contains
   
   subroutine SpecificSetComputeSubgridSpaceResidual(d)
      implicit none
      class(SPComputeSubigridSpaceResidual) :: d
         
      !Allocations
      call ConcatenateProcedures(ProcHook_Initializations,AllocSGRes)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocSGres)
            
      !Compute the Residual at the Gauss point 
      call SetPointersComputeGpResidual%Set
            
      !If residual projection, gather, interpolate and substract the residual projection
      if (a%kfl_repro == 0) then
         !We need the residual in the subscale space
         !The subscale space is obtained with and Identity matrix projector (ASGS) (=gpres)
         call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceResidualASGS)
      elseif (a%kfl_repro /= 0) then
         !We need to gather and interpolate the residual projection
         call SetPointersInterpolateResidualProjection%Set
         !Now we can obtain the projection of the residual onto the subscale space (gpres-gprep)
         call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceResidualOSS)
      endif 
   end subroutine   

   subroutine AllocSGRes
      implicit none
      
      call a%Memor%alloc(a%ResidualSize,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','nsm_EnditeElmope')
   end subroutine
   
   subroutine DeallocSGRes
      implicit none
      
      call a%Memor%dealloc(a%ResidualSize,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','nsm_EnditeElmope')
   end subroutine

   subroutine SubgridSpaceResidualOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres- gprep
   end subroutine
   
   subroutine SubgridSpaceResidualASGS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres
   end subroutine
end module

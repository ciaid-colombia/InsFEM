module Mod_sldsup_SubgridSpaceResidual
   use typre
   use Mod_sld_BaseElmope
   use Mod_sldsup_ComputeGpResidual
   use Mod_sldsup_InterpolateResidualProjection

   implicit none
   private
   public SetPointersComputeSubgridSpaceResidual, gpSubscaleSpaceResidual

   real(rp), allocatable :: gpSubscaleSpaceResidual(:) 
   integer(ip), allocatable :: kfl_IsSet

contains
 
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubgridSpaceResidual(itask)
      implicit none
      integer(ip) :: itask

      !External Procedures
      procedure() :: NULLSUB

      select case (itask) 
 
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
 
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
 
         !Allocations
         call ConcatenateProcedures(ProcHook%Initializations,AllocSGRes)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocSGRes)
 
         !Compute the Residual at the Gauss point 
            call SetPointersComputeGpResidual(1)
 
            !If residual projection, gather, interpolate and substract the residual projection
            if (sup%kfl_repro == 0) then
               !We need the residual in the subscale space
               !The subscale space is obtained with and Identity matrix projector (ASGS) (=gpres)
               call ConcatenateProcedures(ProcHook%InGaussElmats,SubgridSpaceResidualSS)
            elseif (sup%kfl_repro /= 0) then
               !We need to gather and interpolate the residual projection
               !ProcPointer%TimeIntegrationToElext => NULLSUB
               call SetPointersInterpolateResidualProjection(1)
               !Now we can obtain the projection of the residual onto the subscale space (gpres-gprep)
               call ConcatenateProcedures(ProcHook%InGaussElmats,SubgridSpaceResidualOSS)
            endif
         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select 
   end subroutine

   subroutine AllocSGRes
      implicit none
      
      call a%Memor%alloc(a%ndofn,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','sldsup_EndElmope')
   end subroutine
   
   subroutine DeallocSGRes
      implicit none
      
      call a%Memor%dealloc(a%ndofn,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','sldsup_EndElmope')
   end subroutine

   subroutine SubgridSpaceResidualOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres - gprep
   end subroutine
   
   subroutine SubgridSpaceResidualSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres
   end subroutine
end module

module Mod_nsc_pr_InterpolateResidualProjection
   
   use typre
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_ComputeResidualProjection
   implicit none
   private
   public SetPointersInterpolateResidualProjection, gprep
   
   real(rp), allocatable :: elrep(:,:),gprep(:) !The residual projection

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersInterpolateResidualProjection(itask)
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
         
         if (a%kfl_repro /= 0) then
               call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocRep)
               call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocRep)
               call ConcatenateProcedures(ProcHook_nsc_pr_Gathers,GatherRep)
               call ConcatenateProcedures(ProcHook_nsc_pr_Interpolates,InterpolateRep)
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocRep
      implicit none
      
      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','nsc_pr_Elmope')
      call a%Memor%alloc(a%ndofn,gprep,'gprep','nsc_pr_Elmope')
   end subroutine
   
   subroutine DeallocRep
      implicit none
      
      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','nsc_pr_Elmope')
      call a%Memor%dealloc(a%ndofn,gprep,'gprep','nsc_pr_Elmope')
   end subroutine
   
   !Residual Interpolation Subroutines
   subroutine GatherRep
      implicit none
      
      call e%gather(a%ndofn,elrep,a%repro)
   end subroutine

   subroutine InterpolateRep
      implicit none
      
      !Interpolate
      call e%interpg(a%ndofn,elrep,gprep)
   end subroutine
   
end module



module Mod_nsc_pr_SubgridSpaceResidual
   use typre
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_ComputeGpResidual
   use Mod_nsc_pr_InterpolateResidualProjection
   implicit none
   private
   public SetPointersComputeSubgridSpaceResidual, gpdSGSpaceResidual,gpmSGSpaceResidual,gpeSGSpaceResidual

   real(rp), allocatable :: gpdSGSpaceResidual(:),gpmSGSpaceResidual(:),gpeSGSpaceResidual(:)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubgridSpaceResidual(itask)
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
         
         !Allocations
         call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocSGRes)
         call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocSGres)
               
         !Compute the Residual at the Gauss point 
            call SetPointersComputeGpResidual(1)
                  
            !If residual projection, gather, interpolate and substract the residual projection
            if (a%kfl_repro == 0) then
               !We need the residual in the subscale space
               !The subscale space is obtained with and Identity matrix projector (ASGS) (=gpres)
               call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,SubgridSpaceResidualASGS)
            elseif (a%kfl_repro /= 0) then
               !We need to gather and interpolate the residual projection
               call SetPointersInterpolateResidualProjection(1)
               !Now we can obtain the projection of the residual onto the subscale space (gpres-gprep)
               call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,SubgridSpaceResidualOSS)
            endif 
         endif 
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select   
   end subroutine   

   subroutine AllocSGRes
      implicit none
      
      call a%Memor%alloc(1,gpdSGSpaceResidual,'gpdSGSpaceResidual','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,gpmSGSpaceResidual,'gpmSGSpaceResidual','nsc_pr_EndElmope')
      call a%Memor%alloc(1,gpeSGSpaceResidual,'gpeSGSpaceResidual','nsc_pr_EndElmope')
   end subroutine
   
   subroutine DeallocSGRes
      implicit none
      
      call a%Memor%dealloc(1,gpdSGSpaceResidual,'gpdSGSpaceResidual','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,gpmSGSpaceResidual,'gpmSGSpaceResidual','nsc_pr_EndElmope')
      call a%Memor%dealloc(1,gpeSGSpaceResidual,'gpeSGSpaceResidual','nsc_pr_EndElmope')
   end subroutine

   subroutine SubgridSpaceResidualOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpdSGSpaceResidual(1) = gpred(1) - gprep(1)
      gpmSGSpaceResidual(1:e%ndime) = gprem(1:e%ndime) - gprep(2:e%ndime+1)
      gpeSGSpaceResidual(1) = gpree(1) - gprep(e%ndime+2)
   end subroutine
   
   subroutine SubgridSpaceResidualASGS
      implicit none
      
      !Substract the residual projection to gpres
      gpdSGSpaceResidual(1) = gpred(1)
      gpmSGSpaceResidual(1:e%ndime) = gprem(1:e%ndime)
      gpeSGSpaceResidual(1) = gpree(1)
   end subroutine
end module



   
  

module Mod_supm_InterpolateResidualProjection
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersInterpolateResidualProjection

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
               call ConcatenateProcedures(ProcHook_Initializations,AllocRep)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocRep)
               call ConcatenateProcedures(ProcHook_Gathers,GatherRep)
               call ConcatenateProcedures(ProcHook_Interpolates,InterpolateRep)
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

      auxndim=a%ndofn

      call a%Memor%alloc(a%ResidualSize,e%mnode,elrep,'elrep','supm_EnditeElmope')
      call a%Memor%alloc(a%ResidualSize,gprep,'gprep','supm_EnditeElmope')
   end subroutine
   
   subroutine DeallocRep
      implicit none
      
      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','supm_EnditeElmope')
      call a%Memor%dealloc(a%ndofn,gprep,'gprep','supm_EnditeElmope')
      
      
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

!**********************************************************************************
! SUBGRID SPACE RESIDUAL 
!**********************************************************************************

module Mod_supm_SubgridSpaceResidual
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_ComputeGpResidual
   use Mod_supm_InterpolateResidualProjection
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
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
         !Allocations
         call ConcatenateProcedures(ProcHook_Initializations,AllocSGRes)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocSGres)
               
            !Compute the Residual at the Gauss point 
            call SetPointersComputeGpResidual(1)
                  
            !If residual projection, gather, interpolate and substract the residual projection
            if (a%kfl_repro == 0 .or. a%kfl_repro == 3) then
               !We need the residual in the subscale space
               !The subscale space is obtained with and Identity matrix projector (ASGS) (=gpres)
               call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceResidualASGS)
            elseif (a%kfl_repro == 1 .or. a%kfl_repro==2) then
               !We need to gather and interpolate the residual projection
               call SetPointersInterpolateResidualProjection(1)
               !Now we can obtain the projection of the residual onto the subscale space (gpres-gprep)
               call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceResidualOSS)
            endif 
         endif 
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select   
   end subroutine   

   subroutine AllocSGRes
      implicit none
      
      call a%Memor%alloc(a%ResidualSize,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','supm_EnditeElmope')
   end subroutine
   
   subroutine DeallocSGRes
      implicit none
      
      call a%Memor%dealloc(a%ResidualSize,gpSubscaleSpaceResidual,'gpSubscaleSpaceResidual','supm_EnditeElmope')
   end subroutine

   subroutine SubgridSpaceResidualOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres - gprep
   end subroutine
   
   subroutine SubgridSpaceResidualASGS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres
   end subroutine
end module



   
  

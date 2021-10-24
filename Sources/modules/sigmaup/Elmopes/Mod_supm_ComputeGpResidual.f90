  module Mod_supm_ComputeGpResidual
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateGradients  
   implicit none
   private
   public SetPointersComputeGpResidual, gpres

   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SetPointersComputeGpResidual(itask) 
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
         
            !We need the gradients in the gauss point
            call SetPointersInterpolateGradients(1)
            call ConcatenateProcedures(ProcHook_Initializations,AllocGpRes)
           
            !if dynamic subscales are not used, SplitOSS in viscoelastic case is not exactly a common "split",
            !-only for tau1 terms- and the total residual is employed in both cases (OSS, SplitOSS).
            call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpRes)  

            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeRHSNonLinearRep)

             call ConcatenateProcedures(ProcHook_Finalizations,DeallocGpRes)
         end if

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------
   !Residual computation
   subroutine AllocGpRes
      implicit none
      
      call a%Memor%alloc(a%ResidualSize,gpres,'gpres','supm_EnditeElmope')
   end subroutine
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(a%ResidualSize,gpres,'gpres','supm_EnditeElmope')
   end subroutine
   
   subroutine ComputeGpRes
      implicit none    
      call  ProcPointer%FEResPro 
   end subroutine
   
   subroutine ComputeRHSNonLinearRep
      implicit none   
      call supm_elmrfe_oto_nonlinear(e,acvis,auxtens,beta,elvel,gpres)
   end subroutine
   
end module



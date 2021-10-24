module Mod_tem_ComputeGpSubscaleSpaceResidual
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_TempeGradient
   use Mod_tem_InterpolateResidualProjection
   use Mod_tem_ComputeGpResidual
   implicit none
   private
   public SetPointersComputeGpSubscaleSpaceResidual,gpSubscaleSpaceResidual

   
   integer(ip), allocatable :: kfl_IsSet

   real(rp) :: gpSubscaleSpaceResidual
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeGpSubscaleSpaceResidual(itask)
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
            
            call SetPointersComputeGpResidual(1)
            
            !Now we need to compute the residual in the subscales space
            if (a%kfl_repro == 0) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpSubscaleSpaceResidualASGS)
            else 
               call SetPointersInterpolateResidualProjection(1)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpSubscaleSpaceResidualOSS)
            endif
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
  !----------------------------------------------------------------------
   subroutine ComputeGpSubscaleSpaceResidualASGS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres
   end subroutine
   
    subroutine ComputeGpSubscaleSpaceResidualOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceResidual = gpres- gprep(1)
   end subroutine
   
   
   
end module 

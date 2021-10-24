module Mod_tem_ComputeGpResidual
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_TempeGradient
   use Mod_tem_InterpolateResidualProjection
   implicit none
   private
   public SetPointersComputeGpResidual,gpres

   
   integer(ip), allocatable :: kfl_IsSet

   real(rp) :: gpres,gpSubscaleSpaceResidual
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeGPResidual(itask)
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
            
            !We need to interpolate the temperature gradients
            call SetPointersComputeTempeGradient(1)
            
            !we need to compute the gauss point residual, and substract the residual projection
            !if subscales are orthogonal
            if (a%kfl_repro == 0) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpResidualASGS)
            else
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpResidualOSS)
            endif
            !contribution of nonlinear terms to the residual  
            call a%Mesh%isNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpNonLinearRes)
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
  !----------------------------------------------------------------------
   !FOR COMPUTING THE GAUSS POINT RESIDUAL
   subroutine ComputeGpResidualASGS
      implicit none
      
      !Compute
      call tem_elmrfe_osa(e,LHSDtinv,acden,acsph,acrcp,acsou,gpadv, &
         gptem(1),grtem,elext,gpres)
   end subroutine
   
   subroutine ComputeGpResidualOSS
      implicit none
      
      !Compute
      call tem_elmrfe_oto(e,LHSDtinv,acden,acsph,acrcp,acsou,gpadv, &
         gptem(1),grtem,elext,gpres)
   end subroutine
      
   subroutine ComputeGpNonLinearRes
      implicit none
   
      call tem_elmrfe_oto_nonlinear(e,acvis,eltem,gpres)
   end subroutine
end module 

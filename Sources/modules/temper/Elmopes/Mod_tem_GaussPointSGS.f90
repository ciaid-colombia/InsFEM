module Mod_tem_GaussPointSGS
   use typre
   use mod_tem_BaseElmope
   use Mod_tem_ComputeTaus
   use Mod_tem_ComputeGPSubscaleSpaceResidual
   use Mod_tem_InterpolateResidualProjection
   implicit none
   private
   public SetPointersComputeGpSGS,gptempesgs

   integer(ip), allocatable :: kfl_IsSet
   
   real(rp) :: gptempesgs(2)
   
contains
   
   !Set Pointers
   subroutine SetPointersComputeGpSGS(itask)
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
         
            
            !Sets all the pointers so that the residual at the gauss point is computed
            call SetPointersComputeGpSubscaleSpaceResidual(1)
            !Static subscales
            if (a%kfl_tacsg == 0) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeTempeSgsQSS)
                  
            !Transient subgrid scales
            elseif (a%kfl_tacsg == 1) then
            
               call SetPointersComputeTaus(1)
               !Computes the temperature subscale value
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeTempeSgsDSS)
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
   !-------------------------------------------------------------------
   !TRACKING OF THE SUBSCALE
   subroutine ComputeTempeSgsDSS
      implicit none
      
      !Uses a Backward Euler scheme
      a%tesgs(ielem)%a(1,e%igaus) = timom*(-gpSubscaleSpaceResidual + acden*a%tesgs(ielem)%a(2,e%igaus)*ReferenceDtinv)
      GpTempeSgs(1) = a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
      
   !Static subscales   
   subroutine ComputeTempeSgsQSS
      implicit none
      
      GpTempeSgs(1) = -timom*gpSubscaleSpaceResidual
      if (a%kfl_trasg /= 0) then
         a%tesgs(ielem)%a(1,e%igaus) = GpTempeSgs(1)
      endif
   end subroutine

end module
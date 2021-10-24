module Mod_lmn_ComputeNonLinearSubscales
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_SubgridSpaceResidual
   use Mod_lmn_ComputeSubscales
   implicit none
   private
   public SetPointersComputeNonLinearSubscales

   !Non-linear subgrid scales
   integer(ip) :: iiterNolSGS,iiterNolSGS_print
   real(rp)    :: SGSMaxResidual(2)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeNonLinearSubscales(itask)
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
            
            call ConcatenateProcedures(ProcHook%Initializations,NonLinearSubscalesInitializations)
            
            !Set Iteration Counter To Zero
            call ConcatenateProcedures(ProcHook%InGauss,NonLinearSubscalesIiterToZero)
            !We need to compute the subscales at the GaussPoints
            call SetPointersComputeSubscales(1)
            call ConcatenateProcedures(ProcHook%InGaussElmats,NonLinearSubscaleConvergence)   
            call ConcatenateProcedures(ProcHook%Finalizations,NonLinearSubscalesFinalizations)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute 
   !Non-linear subscales
   subroutine NonLinearSubscalesInitializations
      implicit none
      
      iiterNolSGS_print = 0.0_rp
      SGSMaxResidual = 0.0_rp
   end subroutine

   subroutine NonLinearSubscalesFinalizations
      implicit none
     
      if (a%MPIrank == a%MPIroot) then 
         write(a%lun_nolin,301) a%istep, a%itera, iiterNolSGS_print, SGSMaxResidual(1), SGSMaxResidual(2)
         if (a%kfl_flush == 1) call flush(a%lun_nolin)
      endif

301 format(2x,i9,2x,i9,2x,i9,4x,2(4x,e12.6))

   end subroutine
   
   subroutine NonLinearSubscalesIiterToZero
      implicit none
      
      iiterNolSGS = 0
      
      !We copy the subscales at the previous iteration
      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+2) = a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
   
   subroutine NonLinearSubscaleConvergence
     implicit none

      real(rp) :: SgsNorm_mom,SgsNorm_ene,tolsgs
      real(rp) :: SGSResidual(e%ndime+1)
      real(rp) :: SGSResidualNorm_mom,SGSResidualNorm_ene,SGSResidualNorm(2)
     
      SgsNorm_mom = 0.0_rp
      SgsNorm_ene = 0.0_rp
      tolsgs = a%tosgs*a%rilmn

      !Convergence check
      SGSResidual(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) - GpSGS(1:e%ndime)
      SGSResidual(e%ndime+1) = a%tesgs(ielem)%a(1,e%igaus) - GpSGS(e%ndime+2)

      call vecnor(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,SgsNorm_mom,2)
      call vecnor(a%tesgs(ielem)%a(1,e%igaus),e%ndime,SgsNorm_ene,2)
      call vecnor(SGSResidual(1:e%ndime),e%ndime,SgsResidualNorm_mom,2)
      call vecnor(SGSResidual(e%ndime+1),1,SgsResidualNorm_ene,2)
      
      if (SgsNorm_mom == 0.0_rp) then
         SgsResidualNorm_mom = 0.0_rp
      else
         SgsResidualNorm_mom = SgsResidualNorm_mom/SgsNorm_mom
      end if
      SgsResidualNorm(1) = SgsresidualNorm_mom
      if (SgsNorm_ene == 0.0_rp) then
         SgsResidualNorm_ene = 0.0_rp 
      else
         SgsResidualNorm_ene = SgsResidualNorm_ene/SgsNorm_ene
      end if
      SgsResidualNorm(2) = SgsResidualNorm_ene
      
      !Tolerance criteria
      if (SgsResidualNorm(1) < tolsgs .and. SgsResidualNorm(2) < tolsgs) then
         !Exit the non-linear subscales loop
         kfl_GoIteInGauss = 0
      endif
      
      !Maximum number of iterations criteria
      iiterNolSGS = iiterNolSGS + 1
      if (iiterNolSGS >= a%mtrit) then
         kfl_GoIteInGauss = 0
      endif
      
      !If not done, continue
      if (kfl_goiteInGauss /= 0) then 
         kfl_goiteInGauss = 2
      else
         !Maximum residual for the subscales
         SGSMaxResidual(1) = max(SGSMaxResidual(1),SgsResidualNorm(1))
         SGSMaxResidual(2) = max(SGSMaxResidual(2),SgsResidualNorm(2))
         iiterNolSGS_print = max(iiterNolSGS_print,iiterNolSGS)
      endif
      
      !Keep the subscales value for the next convergence check
      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+2) = a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
end module

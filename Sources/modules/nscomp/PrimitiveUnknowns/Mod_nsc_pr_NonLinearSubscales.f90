module Mod_nsc_pr_ComputeNonLinearSubscales
   use typre
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_SubgridSpaceResidual
   use Mod_nsc_pr_ComputeSubscales
   implicit none
   private
   public SetPointersComputeNonLinearSubscales

   !Non-linear subgrid scales
   integer(ip) :: iiterNolSGS
   real(rp)    :: GpSGS_i(5)
   real(rp)    :: SGSMaxResidual
   
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
            
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,NonLinearSubscalesInitializations)
            
            !Set Iteration Counter To Zero
            call ConcatenateProcedures(ProcHook_nsc_pr_InGauss,NonLinearSubscalesIiterToZero)
            !We need to compute the subscales at the GaussPoints
            call SetPointersComputeSubscales(1)
            !Compute residual and decide if we need to continue iterating
            call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,NonLinearSubscalesInGaussElmats)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,NonLinearSubscalesFinalizations)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute 

   subroutine NonLinearSubscalesInitializations
      implicit none
      
      SGSMaxResidual = 0.0_rp
   end subroutine
   
   subroutine NonLinearSubscalesFinalizations
      implicit none
      
      if (a%MPIrank == a%MPIroot) write(a%lun_nolin,*) 'Step: ', a%istep, ' Iteration: ', a%itera, ' Non Linear Subscales Maximum Residual: ', SGSMaxResidual
      !if (a%kfl_flush == 1) call flush(a%lun_nolin)  !memory error
   end subroutine
   
   subroutine NonLinearSubscalesIiterToZero
      implicit none
      
      iiterNolSGS = 0
      
      !We copy the subscales at the previous iteration
     GpSGS_i(1:e%ndime)   = a%mosgs(ielem)%a(1:e%ndime,1,e%igaus)
     GpSGS_i(4)   = a%cosgs(ielem)%a(1,e%igaus)
     GpSGS_i(5)   = a%ensgs(ielem)%a(1,e%igaus)


   end subroutine
   
   subroutine NonLinearSubscalesInGaussElmats
      implicit none
      
      a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime) 
      a%cosgs(ielem)%a(1,e%igaus) = GpSGS(4) 
      a%ensgs(ielem)%a(1,e%igaus) = GpSGS(5) 

      !Relaxation parameter
      if (a%subrelax /= 1.0_rp) then
         a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = a%subrelax*GpSGS(1:e%ndime) + (1-a%subrelax)*GpSGS_i(1:e%ndime)
         a%cosgs(ielem)%a(1,e%igaus) = a%subrelax*GpSGS(4) + (1-a%subrelax)*GpSGS_i(4)
         a%ensgs(ielem)%a(1,e%igaus) = a%subrelax*GpSGS(5) + (1-a%subrelax)*GpSGS_i(5)
      endif
   end subroutine
end module 

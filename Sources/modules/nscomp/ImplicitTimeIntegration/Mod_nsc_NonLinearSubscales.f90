module Mod_nsc_ComputeNonLinearSubscales
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_SubgridSpaceResidual
   use Mod_nsc_ComputeSubscales
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
            
            call ConcatenateProcedures(ProcHook_nsc_Initializations,NonLinearSubscalesInitializations)
            
            !Set Iteration Counter To Zero
            call ConcatenateProcedures(ProcHook_nsc_InGauss,NonLinearSubscalesIiterToZero)
            !We need to compute the subscales at the GaussPoints
            call SetPointersComputeSubscales(1)
            !Compute residual and decide if we need to continue iterating
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,NonLinearSubscalesInGaussElmats)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,NonLinearSubscalesFinalizations)
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
      if (a%kfl_flush == 1) call flush(a%lun_nolin)
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
      
!      real(rp) :: SgsNorm
!      real(rp) :: SGSResidual(e%ndime)
!      real(rp) :: SGSResidualNorm
!      
!      real(rp) :: Matrix(e%ndime,e%ndime),InvMatrix(e%ndime,e%ndime), RHS(e%ndime),deter
!      integer(ip) :: idime,jdime
!      real(rp) :: ResidualTerm(e%ndime), Tau_1, TauDeriv(e%ndime),Tau_1Deriv(e%ndime)
      
      
!      if (a%kfl_nolsgNewtonRaphson == 1 .and. gpvno > 0.0_rp) then
         !If Newton-Raphson
         !\delta \tilde u +tau_t*\delta \tilde u · \nabla u = \tau_t (R(\tilde u) - \ro \tilde u *dtinv)) - \tilde u
!         ResidualTerm = (gpmSGSpaceResidual(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv) FOR DSS
!         ResidualTerm = gpmSGSpaceResidual(1:e%ndime) 
!         Tau_1Deriv = gpden*a%staco(2)/chale(1)*gpadv/gpvno
!         TauDeriv = -timom(2)*timom(2)*Tau_1Deriv
!         Tau_1 = 1/timom(2)
!         
!         !System Matrix
!         Matrix = 0.0_rp
!         
!         !Identity term
!         do idime = 1,e%ndime
!            Matrix(idime,idime) = 1.0_rp
!         enddo
!         
!         !Tau_t times convective term
!         do idime = 1,e%ndime
!            do jdime = 1,e%ndime
!               Matrix(idime,jdime) = Matrix(idime,jdime) + gpden*timom(2)*grvel(idime,jdime)
!            enddo
!         enddo
!         
!         !Tau_t times residual term
!         do idime = 1,e%ndime
!            do jdime = 1,e%ndime
!               Matrix(idime,jdime) = Matrix(idime,jdime) -TauDeriv(jdime)*ResidualTerm(idime)
!            enddo
!         enddo
!         
!         !RHS
!         RHS = timom*ResidualTerm - GpSGS_i(1:e%ndime)
!         
!         !Now we invert and solve
!         call invmtx(Matrix,InvMatrix,deter,e%ndime)
!         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = matmul(InvMatrix,RHS) + GpSGS_i(1:e%ndime)
!
!      endif
!      
      a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime) 
      a%cosgs(ielem)%a(1,e%igaus) = GpSGS(4) 
      a%ensgs(ielem)%a(1,e%igaus) = GpSGS(5) 

      !Relaxation parameter
      if (a%subrelax /= 1.0_rp) then
!         a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = a%subrelax*a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) + (1-a%subrelax)*GpSGS_i(1:e%ndime)
!         a%cosgs(ielem)%a(1,e%igaus) = a%subrelax*a%cosgs(ielem)%a(1,e%igaus) + (1-a%subrelax)*GpSGS_i(4)
!         a%ensgs(ielem)%a(1,e%igaus) = a%subrelax*a%ensgs(ielem)%a(1,e%igaus) + (1-a%subrelax)*GpSGS_i(5)
         a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = a%subrelax*GpSGS(1:e%ndime) + (1-a%subrelax)*GpSGS_i(1:e%ndime)
         a%cosgs(ielem)%a(1,e%igaus) = a%subrelax*GpSGS(4) + (1-a%subrelax)*GpSGS_i(4)
         a%ensgs(ielem)%a(1,e%igaus) = a%subrelax*GpSGS(5) + (1-a%subrelax)*GpSGS_i(5)
      endif
!      
!      !Convergence check
!      SGSResidual(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) - GpSGS_i(1:e%ndime) 
!      call vecnor(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,SgsNorm,2)
!      call vecnor(SGSResidual(1:e%ndime),e%ndime,SgsResidualNorm,2)
!      
!      if (SgsNorm == 0.0_rp) then
!         SgsResidualNorm = 0.0_rp
!      else
!         SgsResidualNorm = SgsResidualNorm/SgsNorm
!      end if
!      
!      !Tolerance criteria
!      if (SgsResidualNorm < a%tosgs) then
!         !Exit the non-linear subscales loop
!         kfl_goiteInGauss = 0
!      endif
!      
!      !Maximum number of iterations criteria
!      iiterNolSGS = iiterNolSGS + 1
!      if (iiterNolSGS >= a%mtrit) then
!         kfl_goiteInGauss = 0
!      endif
!      
!      !If not done, continue
!      if (kfl_goiteInGauss /= 0) then 
!         kfl_goiteInGauss = 2
!      else
!         !Maximum residual for the subscales
!         SGSMaxResidual = max(SGSMaxResidual,SgsResidualNorm)
!      endif
!      
!      !Keep the subscales value for the next convergence check
!      GpSGS_i(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
   end subroutine
end module 

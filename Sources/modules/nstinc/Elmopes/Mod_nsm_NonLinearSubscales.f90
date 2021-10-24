module Mod_nsm_ComputeNonLinearSubscales
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_SubgridSpaceResidual
   use Mod_nsm_ComputeSubscales
   implicit none
   private
   public SetPointersNonLinearSubscales
   
   type, extends(PointerSetter) :: SPNonLinearSubscales
contains
      procedure :: SpecificSet => SpecificSetNonLinearSubscales
   end type
   type(SPNonLinearSubscales) :: SetPointersNonLinearSubscales

   integer(ip) :: iiterNolSGS
   real(rp)    :: GpSGS_i(3)
   real(rp)    :: SGSMaxResidual

contains
   
   subroutine SpecificSetNonLinearSubscales(d)
      implicit none
      class(SPNonLinearSubscales) :: d
            
      call ConcatenateProcedures(ProcHook_Initializations,NonLinearSubscalesInitializations)
      !Set Iteration Counter To Zero
      call ConcatenateProcedures(ProcHook_InGauss,NonLinearSubscalesIiterToZero)
      !We need to compute the subscales at the GaussPoints
      call SetPointersComputeSubscales%Set
      !Compute residual and decide if we need to continue iterating
      call ConcatenateProcedures(ProcHook_InGaussElmats,NonLinearSubscalesInGaussElmats)
      call ConcatenateProcedures(ProcHook_Finalizations,NonLinearSubscalesFinalizations)
   end subroutine   
   
   !-------------------------------------------------------------------
   subroutine NonLinearSGSNonLinearSubscales
      implicit none
      
      gpadv = gpadv + a%vesgs(ielem)%a(:,1,e%igaus)
   end subroutine
  
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
      GpSGS_i(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
   end subroutine
   
   subroutine NonLinearSubscalesInGaussElmats
      implicit none
      real(rp)    :: SgsNorm
      real(rp)    :: SGSResidual(e%ndime)
      real(rp)    :: SGSResidualNorm
      real(rp)    :: Matrix(e%ndime,e%ndime),InvMatrix(e%ndime,e%ndime), RHS(e%ndime),deter
      integer(ip) :: idime,jdime
      real(rp)    :: ResidualTerm(e%ndime), Tau_1, TauDeriv(e%ndime),Tau_1Deriv(e%ndime)
      
      if (a%kfl_nolsgNewtonRaphson == 1 .and. gpvno > 0.0_rp) then
         !If Newton-Raphson
         !\delta \tilde u +tau_t*\delta \tilde u · \nabla u = \tau_t (R(\tilde u) - \ro \tilde u *dtinv)) - \tilde u

         ResidualTerm = (-gpSubscaleSpaceResidual(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
         Tau_1Deriv = acden*a%staco(2)/chale(1)*gpadv/gpvno
         TauDeriv = -timom*timom*Tau_1Deriv
         Tau_1 = 1/timom
         
         !System Matrix
         Matrix = 0.0_rp
         !Implementation 1
         !Identity term
         do idime = 1,e%ndime
            Matrix(idime,idime) = 1.0_rp
         enddo
         
         !Tau_t times convective term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) + acden*timom*grvel(idime,jdime)
            enddo
         enddo
         
         !Tau_t times residual term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) -TauDeriv(jdime)*ResidualTerm(idime)
            enddo
         enddo
         
         !RHS
         RHS = timom*ResidualTerm - GpSGS_i(1:e%ndime)
         
         !Now we invert and solve
         call invmtx(Matrix,InvMatrix,deter,e%ndime)
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = matmul(InvMatrix,RHS) + GpSGS_i(1:e%ndime)
      endif
      
      !Relaxation parameter
      if (a%relsg /= 1.0_rp) then
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = a%relsg*a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) + (1-a%relsg)*GpSGS_i(1:e%ndime)
      endif
      
      !Convergence check
      SGSResidual(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) - GpSGS_i(1:e%ndime) 
      call vecnor(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,SgsNorm,2)
      call vecnor(SGSResidual(1:e%ndime),e%ndime,SgsResidualNorm,2)
      
      if (SgsNorm == 0.0_rp) then
         SgsResidualNorm = 0.0_rp
      else
         SgsResidualNorm = SgsResidualNorm/SgsNorm
      end if
      
      !Tolerance criteria
      if (SgsResidualNorm < a%tosgs) then
         !Exit the non-linear subscales loop
         kfl_goiteInGauss = 0
      endif
      
      !Maximum number of iterations criteria
      iiterNolSGS = iiterNolSGS + 1
      if (iiterNolSGS >= a%mtrit) then
         kfl_goiteInGauss = 0
      endif
      
      !If not done, continue
      if (kfl_goiteInGauss /= 0) then 
         kfl_goiteInGauss = 2
      else
         !Maximum residual for the subscales
         SGSMaxResidual = max(SGSMaxResidual,SgsResidualNorm)
      endif
      
      !Keep the subscales value for the next convergence check
      GpSGS_i(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
   end subroutine

end module 

module Mod_supm_LogarithmicProblem
   use typre
   use Mod_supm_BaseElmope
   use Mod_LogOperations
   implicit none
   private
   public SetPointersLogarithmicComponents
!   procedure()                :: NULL()
   integer(ip), allocatable   :: kfl_IsSet   
  
contains

   subroutine SetPointersLogarithmicComponents(itask)
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
            if (a%LogFormulation/=0) then
            
               !Allocations
               call ConcatenateProcedures(ProcHook_Initializations,AllocComp)
               if (a%kfl_repro==4) call ConcatenateProcedures(ProcHook_Initializations,AllocCompSplitConstitutive)
               call ConcatenateProcedures(ProcHook_Initializations,AllocCompGiesekus)
               
               !Compute Exponential Matrices
               if (a%kfl_exacs/=0) then
                  call ConcatenateProcedures(ProcPointer%ComputeAdvectionVelocity_sup,SetGpExactSol)
                  call ConcatenateProcedures (ProcPointer%ComputeAGradV_sup,ComputeExponentialMatricesExactSolution)
               else
                  call ConcatenateProcedures(ProcPointer%ComputeAGradV_sup,ComputeExponentialMatrices)
               end if
               if (a%Giesekus_model==1) call ConcatenateProcedures(ProcPointer%ComputeAGradV_sup,ComputeExponentialProductsGiesekus)
               
               !Compute RHS linearization terms
               call ConcatenateProcedures(ProcPointer%TwoFieldTerms,ComputeRHSLinearizations)
               if ((a%kfl_repro==1 .or. a%kfl_repro==2) .and. a%kfl_tacsg==1) then
                  call ConcatenateProcedures(ProcPointer%TwoFieldTerms,ComputeRHSLinearizationsDynamicSubscales)
               end if   
               if(a%kfl_repro==4) call ConcatenateProcedures(ProcPointer%TwoFieldTerms,ComputeRHSLinearizationsSplit)
               if(a%Giesekus_model==1) call ConcatenateProcedures(ProcPointer%TwoFieldTerms,ComputeRHSLinearizationsGiesekus) 
               
               !Compute Deallocations
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocComp)
               if (a%kfl_repro==4) call ConcatenateProcedures(ProcHook_Finalizations,DeallocCompSplitConstitutive)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocCompGiesekus)
               
               
            end if   
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   !-----------------------------------------------------------------------------------------------
   !Computation Subroutines
   
   !Allocations ------------------------------------------------------------------------------------ 
   subroutine AllocComp
      implicit none
      call a%Memor%alloc(auxtens,ExpGpPsi,'ExpGpPsi','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,ExpGpPsi_Matrix,'ExpGpPsi_Matrix','supm_elmope')
      call a%Memor%alloc(auxtens,e%ndime,GrExpPsi,'GrExpPsi','supm_elmope')
      call a%Memor%alloc(e%ndime,DivExpPsi,'DivExpPsi','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%ndime,GrExpMatrix,'GrExpMatrix','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,ConvExpMatrix,'ConvExpMatrix','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,GrvelExp,'GrvelExp','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%ndime,GrPsiMatrix,'GrPsiMatrix','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,ExpGradV,'ExpGradV','supm_elmope') 
      
      call a%Memor%alloc(e%ndime,e%ndime,RHSconst,'RHSconst','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,RHSconstEstab,'RHSconst','supm_elmope')
      call a%Memor%alloc(e%ndime,RHSmom,'RHSmom','supm_elmope')
      call a%Memor%alloc(e%ndime,RHS1momGal,'RHS1momGal','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,RHS2momGal,'RHS2momGal','supm_elmope')
      
   end subroutine
   
   subroutine AllocCompGiesekus
      implicit none
      call a%Memor%alloc(e%ndime,e%ndime,ExpGiesekus_Matrix,'ExpGiesekus_Matrix','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,RHSGiesekus,'RHSGiesekus','supm_elmope')

   end subroutine
   
   subroutine AllocCompSplitConstitutive
      implicit none
      call a%Memor%alloc(e%ndime,e%ndime,RHSconst_conv,'RHSconst_conv','supm_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,RHSconst_deform,'RHSconst_deform','supm_elmope')
   end subroutine   
   
   !Compute Exponential Matrices ----------------------------------------------------------------
   subroutine SetGpExactSol
      implicit none
      real(rp)    :: gpcod(e%ndime), exgpsig(auxtens), exgpsigr(auxtens,e%ndime)
      real(rp)    :: exvel(e%ndime), exveg(e%ndime,e%ndime),  expre, exprg(e%ndime)

      call e%interpg(e%ndime,e%elcod,gpcod)
      call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)
      call exacso%sup_GetVelocity(e%ndime,exvel,exveg)
      call exacso%sup_GetPsi(e%ndime,exgpsig,exgpsigr) 

      gpsig(:,1)=exgpsig(:)
      gpvel(:,1)=exvel(:)
      
      grsig=exgpsigr
      grvel=exveg
         
      if (a%kfl_repro /=0) then
         call exacso%sup_GetPressure(e%ndime,expre,exprg)
         grpre(1,:)=exprg(:)
      end if 

   end subroutine
   
   
  subroutine ComputeExponentialMatricesExactSolution  
      implicit none
       real(rp)  :: gpcod(e%ndime),mnode, gradpsi(auxtens,e%ndime)
       real(rp)  :: psi(auxtens)
      !Initializations
      ExpGpPsi=0.0_rp
      ExpGpPsi_Matrix=0.0_rp
      GrExpPsi=0.0_rp
      DivExpPsi=0.0_rp
      GrExpMatrix=0.0_rp
      GrPsiMatrix=0.0_rp
      
      call e%interpg(e%ndime,e%elcod,gpcod) 
      call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)
      call exacso%sup_GetExponential(e%ndime,ExpGpPsi)
      call sup_SigmaMatrix(e%ndime,auxtens,ExpGpPsi,ExpGpPsi_Matrix)
      call exacso%sup_GetExponentialGradient(e%ndime,GrExpPsi)
      call exacso%sup_GetExponentialDivergence(e%ndime,DivExpPsi)
      call exacso%sup_GetPsi(e%ndime,psi,gradpsi)
      call PassGradStressToTensor3(e,auxtens,gradpsi,GrPsiMatrix)  
      
      call ComputeExponentialProducts
      
  end subroutine    
   
   
   subroutine ComputeExponentialMatrices 
      implicit none
      real(rp)  :: elesig(auxtens,e%mnode)
      real(rp)  :: psi(auxtens), GradStress(auxtens,e%ndime)
      
      !Initializations
      ExpGpPsi=0.0_rp
      ExpGpPsi_Matrix=0.0_rp
      GrExpPsi=0.0_rp
      DivExpPsi=0.0_rp
      GrPsiMatrix=0.0_rp

      psi(:)=gpsig(:,1)
      elesig(1:auxtens,:)=elsig(1:auxtens,:,1)
      
      call sup_ComputeExponential(e%ndime,auxtens,psi,ExpGpPsi,ExpGpPsi_Matrix)
      call ComputeGradientExponential(e,auxtens,elesig,GrExpPsi)
      !call ComputeGradientExponentialGaussPoints(e,auxtens,grsig,ExpGpPsi_Matrix,GrExpPsi)
      call sup_GetDivS(e%ndime,auxtens,GrExpPsi,DivExpPsi)
      call PassGradStressToTensor3(e,auxtens,grsig,GrPsiMatrix)
      call ComputeExponentialProducts
       
   end subroutine
   
   
   subroutine ComputeExponentialProducts 
      implicit none
      GrExpMatrix=0.0_rp
      ExpGradV=0.0_rp
      ExpGiesekus_Matrix=0.0_rp
      call PassGradStressToTensor3(e,auxtens,GrExpPsi,GrExpMatrix)
      call ComputeExponentialGradV(e,auxL,ExpGpPsi_Matrix,ExpGradV)
      call ComputeConvectiveExponentialTerm(e,gpvel,GrExpMatrix,ConvExpMatrix)
      
      GrvelExp=matmul(grvel,ExpGpPsi_Matrix)
   end subroutine  
   
   subroutine ComputeExponentialProductsGiesekus
      implicit none
      real(rp)  :: ExpPosMatrix2(e%ndime,e%ndime)
      ExpGiesekus_Matrix=0.0_rp
      ExpPosMatrix2=0.0_rp
      ExpPosMatrix2=matmul(ExpGpPsi_Matrix,ExpGpPsi_Matrix)
      ExpGiesekus_Matrix=ExpPosMatrix2-ExpGpPsi_Matrix
   end subroutine   
   
   !Compute linearizations -----------------------------------------------------------------------
   subroutine ComputeRHSLinearizations
      implicit none
      call ComputeMomentumGalerkinRHS(e,auxtens,auxL,acden,elext,grvel,gpvel,gpsig,ExpGpPsi_Matrix,a%kfl_linearConvectiveTerm,RHS1momGal,RHS2momGal)
      call ComputeMomentumRHS(e,auxtens,auxL,acden,gpsig,GrPsiMatrix,DivExpPsi,ExpGpPsi_Matrix,RHSmom)
      call ComputeConstitutiveRHS(e,lambda,lambda0,auxtens,LHSdtinv,gpsig,GrPsiMatrix,gpvel,grvel,ExpGpPsi_Matrix,GrExpMatrix,a%kfl_linearConstitutiveTerms,RHSconst)
      RHSconstEstab=RHSconst 
   end subroutine
   
   subroutine ComputeRHSLinearizationsDynamicSubscales
      implicit none
       call ComputeConstitutiveRHS(e,lambda,lambda0,auxtens,0.0_rp,gpsig,GrPsiMatrix,gpvel,grvel,ExpGpPsi_Matrix,GrExpMatrix,a%kfl_linearConstitutiveTerms,RHSconstEstab)
   end subroutine   
   
   subroutine ComputeRHSLinearizationsSplit
      implicit none
      !compute RHSconst_conv and RHSconst_def for Split-OSS dynamic model
      call ComputeConstitutiveRHSconvective(e,lambda,lambda0,auxtens,gpsig,GrPsiMatrix,gpvel,ExpGpPsi_Matrix,GrExpMatrix,ConvExpMatrix,a%kfl_linearConstitutiveTerms,RHSconst_conv) 
!          call ComputeConstitutiveRHSconvective(e,lambda,lambda0,auxtens,gpsig,GrPsiMatrix,gpvel,ExpGpPsi_Matrix,GrExpMatrix,a%kfl_linearConstitutiveTerms,RHSconst_conv)
      call ComputeConstitutiveRHSdeformation(e,lambda,lambda0,auxtens,gpsig,grvel,ExpGpPsi_Matrix,a%kfl_linearConstitutiveTerms,RHSconst_deform) 
      
   end subroutine   
   
   subroutine ComputeRHSLinearizationsGiesekus
      implicit none
       call ComputeTerm_rhsGiesekus(e,auxtens,auxG,gpsig,ExpGiesekus_Matrix,ExpGpPsi_Matrix,RHSGiesekus)
   end subroutine    

!---------------------------------------------------------------------------------------------------
! Subroutines that computes exponential elements.

   subroutine ComputeGradientExponential(e,auxtens,elsig,GrExpPsi)
      !Compute  Grad(exp(S))
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: auxtens,i
      real(rp) :: elsig(auxtens,e%mnode),sig(auxtens)
      real(rp) :: Expelsig(auxtens,e%mnode)
      real(rp) :: GrExpPsi(auxtens,e%ndime), Expo(auxtens)
      real(rp) :: ExpMatrix(e%ndime,e%ndime)
      
      do i=1,e%mnode
         sig=0.0_rp
         sig(:)=elsig(:,i)
         call sup_ComputeExponential(e%ndime,auxtens,sig,Expo,ExpMatrix)
         Expelsig(:,i)=Expo(:)
      end do
   
      call e%gradient(auxtens,Expelsig,GrExpPsi)
   end subroutine 
   
   
   subroutine ComputeGradientExponentialGaussPoints(e,auxtens,grsig,ExpoSigMatrix,GrExpPsi)
      !Compute  Grad(exp(S))
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: auxtens,i,j,k,l
      real(rp) :: ExpoSigMatrix(e%ndime,e%ndime), grsig(auxtens,e%ndime)
      real(rp) :: GrExpPsi(auxtens,e%ndime), GrPsiMatrix(e%ndime,e%ndime,e%ndime)
      real(rp) :: GrExpPsi_aux(e%ndime,e%ndime,e%ndime)
      
      GrExpPsi_aux=0.0_rp
      call PassGradStressToTensor3(e,auxtens,grsig,GrPsiMatrix)
     
      do concurrent (l=1:e%ndime, k=1:e%ndime, j=1:e%ndime)
         GrExpPsi_aux(:,k,l)=ExpoSigMatrix(:,j)*GrPsiMatrix(j,k,l) + GrExpPsi_aux(:,k,l)
      end do   
  
      call PassTensor3ToTensor2Left_sym(e,auxtens,GrExpPsi_aux,GrExpPsi)
    
  end subroutine 
  
  subroutine ComputeExponentialGradV(e,auxL,ExpPosMatrix, ExpGradV)
      !Compute Exp(S)*Grad(V)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: i,j,inode
      real(rp)  :: ExpGradV(e%ndime,e%mnode),ExpPosMatrix(e%ndime,e%ndime)
      real(rp)  :: auxL
      
      ExpGradV=0.0_rp

      do concurrent (i=1:e%ndime, j=1:e%ndime)
         ExpGradV(i,:)=ExpPosMatrix(j,i)*e%cartd(j,:) + ExpGradV(i,:)
      end do

  end subroutine
  
  subroutine ComputeConvectiveExponentialTerm(e,gpvel,GrExposMatrix, ConvectiveTerm)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: i, j, k
      real(rp) :: gpvel(e%ndime), GrExposMatrix(e%ndime, e%ndime, e%ndime)
      real(rp) :: ConvectiveTerm(e%ndime,e%ndime)
      
      ConvectiveTerm=0.0_rp
      do k=1,e%ndime
         ConvectiveTerm(:,1:e%ndime)=gpvel(k)*GrExposMatrix(:,1:e%ndime,k) + ConvectiveTerm(:,1:e%ndime)
      end do
   
   end subroutine
   
   !--------------------------------------------------------------------------------------------
   !Deallocation
   
    subroutine DeallocComp
      implicit none
      call a%Memor%dealloc(auxtens,ExpGpPsi,'ExpGpPsi','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,ExpGpPsi_Matrix,'ExpGpPsi_Matrix','supm_elmope')
      call a%Memor%dealloc(auxtens,e%ndime,GrExpPsi,'GrExpPsi','supm_elmope')
      call a%Memor%dealloc(e%ndime,DivExpPsi,'DivExpPsi','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%ndime,GrExpMatrix,'GrExpMatrix','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,ConvExpMatrix,'ConvExpMatrix','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,GrvelExp,'GrvelExp','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%ndime,GrPsiMatrix,'GrPsiMatrix','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,ExpGradV,'ExpGradV','supm_elmope')
      
      call a%Memor%dealloc(e%ndime,e%ndime,RHSconst,'RHSconst','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,RHSconstEstab,'RHSconst','supm_elmope')
      call a%Memor%dealloc(e%ndime,RHSmom,'RHSmom','supm_elmope')
      call a%Memor%dealloc(e%ndime,RHS1momGal,'RHS1momGal','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,RHS2momGal,'RHS2momGal','supm_elmope')
      
   end subroutine
   
   subroutine DeallocCompGiesekus
      implicit none
      call a%Memor%dealloc(e%ndime,e%ndime,ExpGiesekus_Matrix,'ExpGiesekus_Matrix','supm_elmope') 
      call a%Memor%dealloc(e%ndime,e%ndime,RHSGiesekus,'RHSGiesekus','supm_elmope')
   end subroutine   
   
   subroutine DeallocCompSplitConstitutive
      implicit none
       call a%Memor%dealloc(e%ndime,e%ndime,RHSconst_conv,'RHSconst_conv','supm_elmope')
       call a%Memor%dealloc(e%ndime,e%ndime,RHSconst_deform,'RHSconst_deform','supm_elmope')
   
   end subroutine

end module

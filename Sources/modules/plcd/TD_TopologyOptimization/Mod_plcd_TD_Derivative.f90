module Mod_plcd_TD_Derivative
   use typre
   use Mod_Element
   use Mod_Mesh
   use Mod_Postpr
   use Mod_int2str
   use Mod_plcd_TD_PolarizationTensor
   use Mod_plcd_SIMPMaterial
   use Mod_plcd_Material
   use Mod_Memor
   use Mod_plcd_TDData
   use Mod_MpiObject
   use Mod_CutMesh
   use MPI
   use Mod_plcd_ConstitutiveTensor
   use Mod_plcd_Isotrop2DPlainStress
   use Mod_plcd_IsotropConstitutiveTensor
   use Mod_debugging
   implicit none
   private
   public :: TD_Derivative, PostprocessTopologicalDerivative

   type, extends(MpiObject) :: TD_Derivative

      real(rp), allocatable :: PW_Derivative(:,:), ElementVolumeFraction(:),NodalChi(:), LocalKappa(:),PW_Compliance(:)
      real(rp) :: MeanDerivative = 0.0_rp, MinDerivative = 0.0_rp, MaxDerivative = 0.0_rp
      real(rp) :: TempMaxDerivative = -1e12_rp, TempMinDerivative = 1e12_rp, TempMeanDerivative = 0.0_rp, TempVolume = 0.0_rp
      class(FiniteElement), pointer :: eClosedRule

      real(rp) :: HistoryLambda(3) = 0.0_rp, HistoryConstraint(3) = 0.0_rp
      real(rp) :: HistoryCompliance(3) = 1e32_rp, PostprocessCompliance = 0.0_rp
      integer(ip) :: itercompliance = 0

      class(FemMesh), pointer :: Mesh => NULL()
      type(MemoryMan), pointer :: MemorP => NULL()
      type(TDDataType), pointer :: TDData => NULL()

      class(PolarizationTensor), pointer :: PT => NULL()
      integer :: IterCounter = 0

      real(rp) :: GaussTD(1,50), GaussCompliance(1,50)
      real(rp) :: kappa = 1.0_rp
      logical :: kfl_ComputePW_Compliance = .false.
      real(rp) :: CurrentVolumeFraction
      
      integer(ip) :: kfl_LargeStrains = 0 !Linear Elastic is 0. Updated Formulation Nonlinear Elastic 1. Total Formulation Nonlinear Elastic 2.

contains
      procedure :: SetMesh
      procedure :: SetMemor
      procedure :: SetTDData
      procedure :: SetPTType
      procedure :: SetLargeStrains
      procedure :: SetComputeCompliance
      procedure :: Initialize
      procedure :: Adaptive
      procedure :: Finalize
      procedure :: AddGaussPointData
      procedure :: PostGaussData
      procedure :: FinalizeDataComputations
      procedure :: GetPointerToTopologicalDerivativeArray
      procedure :: GetPointerToComplianceArray
      procedure :: GetCompliance
      procedure :: Reset
      procedure :: UpdateTopologicalDerivative
      procedure :: GetGaussPointChi
      procedure :: PostprocessChi
      procedure :: GetNodalChi
      procedure :: WriteLogData
   end type

contains

   subroutine SetMesh(a,Mesh)
      class(TD_Derivative) :: a
      class(FemMesh), target :: Mesh
      a%Mesh => Mesh
   end subroutine

   subroutine SetMemor(a,Memor)
      class(TD_Derivative) :: a
      type(MemoryMan), target :: Memor
      a%MemorP => Memor
   end subroutine

   subroutine SetTDData(a,TDData)
      implicit none
      class(TD_Derivative) :: a
      type(TDDataType), target :: TDData

      a%TDData => TDData
   end subroutine
   
   subroutine SetLargeStrains(a,kfl_LargeStrains)
      implicit none 
      class(TD_Derivative) :: a
      integer(ip) :: kfl_LargeStrains
      
      a%kfl_LargeStrains = kfl_LargeStrains
   end subroutine
   
   subroutine WriteLogData(a,lun_outpu)
      implicit none
      class(TD_Derivative) :: a
      integer(ip) :: lun_outpu
      write(lun_outpu,100,advance='no') ' LocalKappa: ', a%LocalKappa(1), ' VolumeFraction: ',a%CurrentVolumeFraction
      100 format (A11, E12.5, A11, E12.5)
   end subroutine


   subroutine SetPTType(a,CT)
      class(TD_Derivative) :: a
      class(ConstitutiveTensor) :: CT

      select type(CT)
         type is (Isotrop2dPlainStress)
            if (a%kfl_LargeStrains == 0) then
               allocate(PlainStressPolarizationTensor::a%PT)
            else
               allocate(PolarizationTensor2DLargeStrains::a%PT)
            endif
         type is (IsotropCT)
            if (a%kfl_LargeStrains == 0) then
               allocate(PolarizationTensor3D::a%PT)
            else
               allocate(PolarizationTensor3DLargeStrains::a%PT)
            endif
         !type default
         !   call runend('PolarizationTensor not implemented')
      end select
   end subroutine
   
   subroutine SetComputeCompliance(a,kfl_ComputePW_Compliance)
      class(TD_Derivative) :: a
      logical :: kfl_ComputePW_Compliance
      
      a%kfl_ComputePW_Compliance = kfl_ComputePW_Compliance
   end subroutine

   subroutine Initialize(a)
      class(TD_Derivative) :: a
      integer(ip) :: npoin, nelem
      real(rp) :: direction(3) = 1.0_rp, E, nu

      integer(ip) :: imaterial
      real(rp) :: gamma

      nu = a%TDData%nu
      gamma = a%TDData%chimax/a%TDData%chimin
      call a%PT%ComputePTensor(nu,gamma)

      call a%Mesh%Getnpoin(npoin)
      call a%Mesh%GetNelem(nelem)

      a%IterCounter = 0
      call a%MemorP%alloc(npoin,4,a%PW_Derivative,'PW_Derivative','PW')
      if (a%kfl_ComputePW_Compliance) call a%MemorP%alloc(npoin,a%PW_Compliance,'PW_Derivative','PW')
      call a%MemorP%alloc(npoin,a%LocalKappa,'LocalKappa','PW')
      a%LocalKappa = 0.01_rp
      call a%MemorP%alloc(npoin,a%NodalChi,'NodalChi','PW')

      a%PW_Derivative = 1.0_rp
      a%NodalChi = a%TDData%chimax

      call a%MemorP%alloc(nelem,a%ElementVolumeFraction,'ElementVolumeFraction','PW')
      call a%Mesh%ElementAlloc(a%eClosedRule,a%MemorP,'ForceClosedRule','SIMP_TopOpt')
      
      a%HistoryCompliance(1) = 0.0_rp
      a%HistoryCompliance(2:3) = 1e32_rp
   end subroutine

    subroutine Adaptive(a,b,itask)
      use Mod_AdaptiveInterface
      use Mod_PhysicalProblem
      use Mod_phpRefineArrays
      class(TD_Derivative) :: a
      character(6) :: itask
      class(PhysicalProblem) :: b


      integer(ip) :: nelem,npoin,npoinLocal,ierr
      real(rp) :: norm

      call a%MemorP%dealloc(size(a%ElementVolumeFraction),a%ElementVolumeFraction,'ElementVolumeFraction','PW')

      call a%Mesh%GetNelem(nelem)
      call a%MemorP%alloc(nelem,a%ElementVolumeFraction,'ElementVolumeFraction','PW')

      call a%MemorP%dealloc(size(a%NodalChi),a%NodalChi,'NodalChi','PW')
      call a%Mesh%GetNpoin(npoin)
      call a%MemorP%alloc(npoin,a%NodalChi,'NodalChi','PW')
      call php_RefineArrays(b,itask,a%PW_Derivative,'PW_Derivative')
      if (a%kfl_ComputePW_Compliance) call php_RefineArrays(b,itask,a%PW_Compliance,'PW_Derivative')
      call php_RefineArrays(b,itask,a%LocalKappa,'LocalKappa')

      !ReNormalize the topological derivative for the new mesh
      call a%Mesh%GetNpoinLocal(npoinlocal)
      call vecnorMPI(a%PW_Derivative(:,2),npoinLocal,norm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
      call MPI_BCAST(norm,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      a%PW_Derivative(:,2) = a%PW_Derivative(:,2)/norm

      !Allocate again so that npoinlocal is correct
      call a%Mesh%ElementDeAlloc(a%eClosedRule,a%MemorP,'ForceClosedRule','SIMP_TopOpt')
      call a%Mesh%ElementAlloc(a%eClosedRule,a%MemorP,'ForceClosedRule','SIMP_TopOpt')
      

   end subroutine

   subroutine Finalize(a)
      class(TD_Derivative) :: a
      integer(ip) :: npoin,nelem

      call a%Mesh%Getnpoin(npoin)
      call a%Mesh%GetNelem(nelem)
      call a%MemorP%dealloc(size(a%PW_Derivative,1),size(a%PW_Derivative,2),a%PW_Derivative,'PW_Derivative','PW')
      if (a%kfl_ComputePW_Compliance) call a%MemorP%dealloc(size(a%PW_Compliance),a%PW_Compliance,'PW_Derivative','PW')
      call a%MemorP%dealloc(size(a%LocalKappa,1),a%LocalKappa,'LocalKappa','PW')
      call a%MemorP%dealloc(nelem,a%ElementVolumeFraction,'ElementVolumeFraction','PW')
      call a%MemorP%dealloc(size(a%NodalChi),a%NodalChi,'NodalChi','PW')

      call a%Mesh%ElementDeAlloc(a%eClosedRule,a%MemorP,'ForceClosedRule','SIMP_TopOpt')
      deallocate(a%PT)

      a%IterCounter = 0
   end subroutine

   subroutine AddGaussPointData(a,ielem,e,strain,stress,chi)
      class(TD_Derivative) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: strain(:), stress(:), chi

      real(rp) :: PStress(size(stress)),TDerivative,weightfactor
      real(rp) :: relativeChi,GaussPointCompliance
      
      relativeChi = (chi-a%TdData%Chimin)/(a%TdData%Chimax-a%TdData%Chimin)

      call a%PT%ApplyPTensor(relativeChi,stress,PStress)
      TDerivative = dot_product(PStress,strain(1:size(stress)))
      
      !call a%PT%ApplyPTensor(relativeChi,strain(1:size(stress)),PStress)
      !TDerivative = dot_product(PStress,stress)
  
      !TDerivative = dot_product(stress,strain(1:size(stress)))


      !HistoryCompliance
      call e%GetWeightFactor(weightfactor)
      GaussPointCompliance = dot_product(stress,strain(1:size(stress)))
      a%HistoryCompliance(1) = a%HistoryCompliance(1) + e%detjm*e%weigp(e%igaus)*GaussPointCompliance*weightfactor
      
      
      
      a%GaussTD(1,e%igaus) = TDerivative
      if (a%kfl_ComputePW_Compliance) a%GaussCompliance(1,e%igaus) = GaussPointCompliance
   end subroutine

   subroutine PostGaussData(a,ielem,e)
      class(TD_Derivative) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e

      real(rp) :: dvol,elTD(e%pnode),TD,interp_TD(1),elCompliance(e%pnode),GaussPointCompliance,interp_Compliance(1)
      integer(ip) :: igaus
      real(rp) :: NodeTD(1,e%pnode),NodeCompliance(1,e%pnode)
      real(rp) :: weightfactor,Compliance

      call e%GaussToNodes(a%GaussTD,NodeTD)
      call e%GetWeightFactor(weightfactor)
      elTD = 0.0_rp
      
      if (a%kfl_ComputePW_Compliance) then
         call e%GaussToNodes(a%GaussCompliance,NodeCompliance)
         elCompliance = 0.0_rp
      endif

      call a%Mesh%ElementLoad(ielem,a%eClosedRule)
      do igaus = 1,a%eClosedRule%pgaus
         a%eClosedRule%igaus = igaus

         call a%eClosedRule%elmder
         dvol = a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%detjm

         call a%eClosedRule%interpg(1_ip,NodeTD,interp_TD)
         TD = interp_TD(1)

         elTD(1:a%eClosedRule%pnode) = elTD(1:a%eClosedRule%pnode) + a%eClosedRule%detjm*a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%shape(1:a%eClosedRule%pnode,a%eClosedRule%igaus)*TD
         
         
         if (a%kfl_ComputePW_Compliance) then
            call a%eClosedRule%interpg(1_ip,NodeCompliance,interp_Compliance)
            Compliance = interp_Compliance(1)

            elCompliance(1:a%eClosedRule%pnode) = elCompliance(1:a%eClosedRule%pnode) + a%eClosedRule%detjm*a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%shape(1:a%eClosedRule%pnode,a%eClosedRule%igaus)*Compliance
         endif
         


      enddo
      call a%Mesh%AssemblyToArray(a%eClosedRule,1_ip,elTD,a%PW_Derivative(:,1))
      if (a%kfl_ComputePW_Compliance) call a%Mesh%AssemblyToArray(a%eClosedRule,1_ip,elCompliance,a%PW_Compliance(:))
      
   end subroutine

   subroutine FinalizeDataComputations(a)
      class(TD_Derivative) :: a

      real(rp) :: sumderivative,HistoryCompliance
      integer(ip) :: ierr, gnpoin, npoinLocal

      call a%Mesh%Smooth(1,a%PW_Derivative(:,1))
      !where(a%PW_Derivative(:,1) <= 1e-8_rp) a%PW_Derivative(:,1) = 1e-8
      
      if (a%kfl_ComputePW_Compliance) call a%Mesh%Smooth(1,a%PW_Compliance)
      
      
      call  MPI_AllREDUCE( a%HistoryCompliance(1), HistoryCompliance, 1_ip, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
      a%HistoryCompliance(1) = HistoryCompliance
      a%PostprocessCompliance = HistoryCompliance
      
   end subroutine
   
   subroutine GetPointerToTopologicalDerivativeArray(a,TopologicalDerivative)
      class(TD_Derivative), target :: a
      real(rp), pointer :: TopologicalDerivative(:)
      
      TopologicalDerivative => a%PW_Derivative(:,1)
   end subroutine
   
   subroutine GetPointerToComplianceArray(a,Compliance)
      class(TD_Derivative), target :: a
      real(rp), pointer :: Compliance(:)
      
      Compliance => a%PW_Compliance
   end subroutine
   
   subroutine GetCompliance(a,HistoryCompliance)
      class(TD_Derivative) :: a
      real(rp) :: HistoryCompliance 
      
      HistoryCompliance = a%PostprocessCompliance
   end subroutine
   
   

   subroutine PostprocessTopologicalDerivative(a,FilePostpr,itera,istep,ctime)
      class(TD_Derivative) :: a
      class(PostprFile) :: FilePostpr
      integer(ip) :: itera, istep
      real(rp) :: ctime

      call FilePostpr%postpr(a%PW_Derivative(:,1),'PW_Derivative'//adjustl(trim(int2str(itera))),istep,ctime,a%Mesh)
   end subroutine

   subroutine Reset(a)
      class(TD_Derivative) :: a

      a%PW_Derivative(:,1) = 0.0_rp
      if (a%kfl_ComputePW_Compliance) a%PW_Compliance = 0.0_rp
   end subroutine

   subroutine UpdateTopologicalDerivative(a,FilePostpr,itera,istep,ctime)
      use Mod_Postpr
      class(TD_Derivative) :: a
      class(PostprFile) :: FilePostpr
      integer(ip) :: itera,istep
      real(rp) :: ctime

      real(rp) :: VolumeFraction

      integer(ip) :: npoinLocal,ipoin
      type(CutMesh) :: CMesh

      real(rp) :: weigp(50),fractionminus, fractionplus, totalfraction,norm
      integer(ip) :: eStatus,ndime,ielem,nelem,ierr,ngauss_minus,ngauss_plus,ngauss_total,npoin
      real(rp) :: theta,vdot, constraint, lambda, slope
      real(rp) :: gvolume, gtvolume, tvolume, volume
      integer(ip) :: iiter
      real(rp) :: sumderivative
      integer(ip) :: gnpoin
      real(rp) :: diff,ratio
      integer(ip) :: iDetectOscillation
      real(rp)    :: DetectOscillation,auxtransfer(2), auxtransfer2(2),fractio

      a%IterCounter = a%IterCounter + 1
      call auxtimers(1)%Tic


      !If Progressive Volume, we update it here
      if (a%TDData%kfl_ProgressiveVolume == 1) then
         if (a%IterCounter <= a%TDData%ProgressiveVolume_IniDelay) then
            a%TDData%VolumeFraction = a%TDData%ProgressiveVolume_inivol
         elseif (a%IterCounter < a%TDData%ProgressiveVolume_IniDelay +a%TDData%ProgressiveVolume_nsteps) then
            fractio = real(a%IterCounter-a%TDData%ProgressiveVolume_IniDelay)/real(a%TDData%ProgressiveVolume_nsteps)
            a%TDData%VolumeFraction = a%TDData%ProgressiveVolume_inivol + fractio*(a%TDData%ProgressiveVolume_endvol - a%TDData%ProgressiveVolume_inivol)
         else
            a%TDData%VolumeFraction = a%TDData%ProgressiveVolume_endvol
         endif 
      endif
      

      !Normalize the topological derivative
      call a%Mesh%GetNpoinLocal(npoinlocal)
      call vecnorMPI(a%PW_Derivative(:,1),npoinLocal,norm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
      call MPI_BCAST(norm,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      a%PW_Derivative(:,1) = a%PW_Derivative(:,1)/norm
      
      
      a%HistoryCompliance(1) = 0.0_rp
      
      !We advance with a fixed point method
      if (a%IterCounter > 1 ) then
! 
!              !Oscillation capturing fixed Point method
!              a%kappa = 1
!              call a%Mesh%GetNpoin(npoin)   
!              !Additional local relaxation
!              if (a%IterCounter > 3) then
!                 do concurrent (ipoin = 1:npoin)
!                    DetectOscillation = (a%PW_Derivative(ipoin,1)-a%Pw_Derivative(ipoin,2) )/(a%PW_Derivative(ipoin,2) - a%PW_Derivative(ipoin,3) )
!                    iDetectOscillation = -nint(DetectOscillation/abs(DetectOscillation))
!                    if (iDetectOscillation == 1) then
!                       a%LocalKappa(ipoin) = min(a%kappa,0.5_rp*a%LocalKappa(ipoin) )
!                       !a%LocalKappa(ipoin) = 0.5_rp*a%LocalKappa(ipoin)
!                    else
!                       a%LocalKappa(ipoin) = min(a%kappa,1.5_rp*a%LocalKappa(ipoin) )
!                       !a%LocalKappa(ipoin) = 1.5_rp*a%LocalKappa(ipoin) )
!                    endif
!                 enddo
!              else
!                  a%LocalKappa = a%kappa
!              endif
!              !write(*,*) 'Mean kappa: ', sum(a%LocalKappa)/npoin
!              
!              do concurrent (ipoin = 1:npoin)
!                a%PW_Derivative(ipoin,1) = (1-a%LocalKappa(ipoin))*a%PW_Derivative(ipoin,2)+a%LocalKappa(ipoin)*a%PW_Derivative(ipoin,1)
!              enddo

    !Oscillation capturing fixed Point method
             a%kappa = 1
             call a%Mesh%GetNpoin(npoin)   
             !Additional local relaxation
             if (a%IterCounter > 3) then
                do concurrent (ipoin = 1:npoin)
                   DetectOscillation = (a%PW_Derivative(ipoin,1)-a%Pw_Derivative(ipoin,2) )/(a%PW_Derivative(ipoin,2) - a%PW_Derivative(ipoin,3) )
                   iDetectOscillation = -nint(DetectOscillation/abs(DetectOscillation))
                   if (iDetectOscillation == 1) then
                      a%LocalKappa(ipoin) = min(a%kappa,a%TDData%KappaIncreaseCoefficient(2)*a%LocalKappa(ipoin) )
                      !a%LocalKappa(ipoin) = 0.5_rp*a%LocalKappa(ipoin)
                   else
                      a%LocalKappa(ipoin) = min(a%kappa,a%TDData%KappaIncreaseCoefficient(1)*a%LocalKappa(ipoin) )
                      !a%LocalKappa(ipoin) = 1.5_rp*a%LocalKappa(ipoin) )
                   endif
                enddo
             else
                 a%LocalKappa = a%kappa
             endif

             
             call a%Mesh%GetNpoinLocal(npoinLocal)
             auxtransfer(1) = sum(a%LocalKappa(1:npoinLocal)**a%TDData%KappaPowerCoefficient)
             auxtransfer(2) = npoinLocal
             call MPI_REDUCE(auxtransfer, auxtransfer2, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
             auxtransfer2(1) = auxtransfer2(1)/auxtransfer2(2)
             call MPI_BCAST(auxtransfer2(1), 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
             a%LocalKappa = auxtransfer2(1)**(1.0_rp/a%TDData%KappaPowerCoefficient)
             a%LocalKappa = max(a%LocalKappa(1),a%TDData%KappaMin)
             
             
             
             do concurrent (ipoin = 1:npoin)
               a%PW_Derivative(ipoin,1) = (1-a%LocalKappa(ipoin))*a%PW_Derivative(ipoin,2)+a%LocalKappa(ipoin)*a%PW_Derivative(ipoin,1)
             enddo


             
!              
!             !Now we do the relaxation algorithm (SLERP) 
!             call MPIdot(a%PW_Derivative(:,1),a%PW_Derivative(:,2),npoinLocal,vdot,a%MPIcomm,a%MPIroot)
!             call MPI_BCAST(vdot,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
!             vdot = max(-1.0_rp,min(1.0_rp,vdot))
!             theta = acos(vdot)
!             a%kappa =((vdot+1)/2.0_rp*0.8)**5*3
!             a%kappa = min(1.0_rp,a%kappa)
!             if (a%MPirank == a%MPIroot) then
!                write(*,*) 'a%kappa', a%kappa
!             endif
!             a%LocalKappa = a%kappa
!              
!             !Apply relaxation
!             if (abs(sin(theta)) > 1e-12_rp) then
!                !We advance towards the topological derivative
!                call a%Mesh%GetNpoin(npoin)   
!                do concurrent (ipoin = 1:npoin)
!                   a%PW_Derivative(ipoin,1) = 1/sin(theta)*(sin((1-a%LocalKappa(ipoin))*theta)*a%PW_Derivative(ipoin,2)+sin(a%LocalKappa(ipoin)*theta)*a%PW_Derivative(ipoin,1))
!                enddo
!             endif
          
      end if
      !We update the value of the previous iteration so that we can use it in the following iteration
      a%PW_Derivative(:,3) = a%PW_Derivative(:,2)
      a%PW_Derivative(:,2) = a%PW_Derivative(:,1)

      !Just make sure that the initial guess are not two equal values of lambda
      if (abs(a%HistoryLambda(1)-a%HistoryLambda(2)) < 1e-3_rp*abs(a%HistoryLambda(1))) then
         a%HistoryLambda(2) = 0.99*a%HistoryLambda(2)
      end if

      !We initially compute the volume fraction here
      a%PW_Derivative(:,1) = a%PW_Derivative(:,1) - a%HistoryLambda(1)
      call ComputeVolumeFraction


      !Here we do an update of the volume fraction
      lambdado : do iiter = 1,a%TDData%NVolumeIterations
         Constraint = VolumeFraction - a%TDData%VolumeFraction
         if (abs(Constraint)/a%TDData%VolumeFraction <= a%TDData%VolumeTolerance) then
            exit
         endif

         a%HistoryConstraint(3) = a%HistoryConstraint(2)
         a%HistoryLambda(3) = a%HistoryLambda(2)
         a%HistoryConstraint(2) = Constraint
         a%HistoryLambda(2)         = a%HistoryLambda(1)

         if (a%IterCounter == 1 .and. iiter ==1) then
            !Approximate mean derivative for starting the iterations
            call a%Mesh%GetNpoinLocal(npoinLocal)
            sumderivative = sum(a%PW_Derivative(1:npoinLocal,1))
            call  MPI_AllREDUCE( sumderivative, a%MeanDerivative, 1_ip, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
            gnpoin = 0_ip
            call  MPI_AllREDUCE( npoinLocal, gnpoin, 1_ip, MPI_INTEGER4, MPI_SUM,a%MPIcomm, ierr )
            a%MeanDerivative = a%MeanDerivative/gnpoin

            lambda = a%MeanDerivative
         elseif (a%IterCounter == 1 .and. iiter ==2) then
            lambda = a%HistoryLambda(1)*0.9
         else
            !Secant method
            if (abs(a%HistoryLambda(2)-a%HistoryLambda(3)) < 1e-3_rp*abs(a%HistoryLambda(2))) exit lambdado
            slope = (a%HistoryConstraint(2)- a%HistoryConstraint(3) ) /( a%HistoryLambda(2) - a%HistoryLambda(3) )
            lambda = a%HistoryLambda(2) - a%HistoryConstraint(2)/slope
            lambda = max(a%HistoryLambda(2)*0.5_rp,min(a%HistoryLambda(2)*2.0_rp,lambda))
         endif

         a%HistoryLambda(1) = lambda
         a%PW_Derivative(:,1) = a%PW_Derivative(:,2) - lambda

         call ComputeVolumeFraction
      enddo lambdado

      call auxtimers(1)%Toc
      
      !if (a%MPIrank == a%MPIroot) write(*,*) 'Iters: ', iiter, ' Volume Fraction: ', VolumeFraction
contains

      subroutine ComputeVolumeFraction
         real(rp) :: weightfactor

         !Compute the fraction volume for each point
         call CMesh%allocCutMesh(a%MemorP,a%Mesh)
         call a%Mesh%GetNpoin(npoin)
         call a%Mesh%GetNdime(ndime)
         call CMesh%ComputeIntersectionPoints(a%PW_Derivative(:,1))
         call CMesh%ComputeElementType(a%PW_Derivative(:,1))
         call CMesh%ComputePoinType(a%PW_Derivative(:,1))
         call CMesh%ComputeSubelements(ndime)


         call a%Mesh%GetNelem(nelem)
         call a%Mesh%GetNdime(ndime)
         Volume = 0.0_rp
         TVolume = 0.0_rp
         do ielem = 1,nelem
            call a%Mesh%ElementLoad(ielem,a%eClosedRule)
            call a%eClosedRule%elmdcg
            call CMesh%GetElementType(ielem,eStatus)
            if (eStatus >= 1) a%ElementVolumeFraction(ielem) = 1.0_rp
            if (eStatus <= -1) a%ElementVolumeFraction(ielem) = 0.0_rp
            if (eStatus == 0) then
               call CMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
               call CMesh%GetWeigpCut(ielem,ndime,weigp)
               fractionminus = sum(weigp(1:ngauss_minus))
               fractionplus  = sum(weigp(ngauss_minus+1:ngauss_minus+ngauss_plus))
               totalfraction = fractionminus + fractionplus
               fractionminus = fractionminus/totalfraction
               fractionplus  = fractionplus/totalfraction
               a%ElementVolumeFraction(ielem) = fractionplus
            end if
            call a%eClosedRule%GetWeightFactor(weightfactor)
            if (weightfactor == 0.0_rp) then
               call a%eClosedRule%GetWeightFactor(weightfactor)
            end if
            TVolume = TVolume + a%eClosedRule%detjm*weightfactor
            Volume = Volume + a%ElementVolumeFraction(ielem)*a%eClosedRule%detjm*weightfactor
         end do
         !Dealloc the cut element
         call CMesh%deallocCutElement
         call CMesh%deallocCutMesh
         call  MPI_AllREDUCE( Volume, gVolume, 1_ip, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
         call  MPI_AllREDUCE( TVolume, gTVolume, 1_ip, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
         VolumeFraction = gVolume/gTVolume
         a%CurrentVolumeFraction = VolumeFraction
         !write(*,*) 'VolumeFraction: ', VolumeFraction
      end subroutine

      
   end subroutine

   subroutine GetGaussPointChi(a,ielem,e,Chi)
      class(TD_Derivative) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: Chi

      chi = a%ElementVolumeFraction(ielem)*a%TDData%chimax + (1.0_rp-a%ElementVolumeFraction(ielem))*a%TDData%chimin
   end subroutine

   subroutine PostprocessChi(a,FilePostpr,itera,istep,ctime)
      class(TD_Derivative) :: a
      class(PostprFile) :: FilePostpr
      integer(ip) :: itera,istep
      real(rp) :: ctime

      call FilePostpr%postpr(a%PW_Derivative(:,1),'TopologicalDerivative'//adjustl(trim(int2str(itera))),istep,ctime,a%Mesh)
   end subroutine

   subroutine GetNodalChi(a,NodalChi)
      class(TD_Derivative),target :: a
      real(rp), pointer :: NodalChi(:)
      where(a%PW_Derivative(:,1) > 0.0_rp) a%NodalChi = a%TDData%chimax
      where(a%PW_Derivative(:,1) <= 0.0_rp) a%NodalChi = a%TDData%chimin

      NodalChi => a%NodalChi
   end subroutine


end module


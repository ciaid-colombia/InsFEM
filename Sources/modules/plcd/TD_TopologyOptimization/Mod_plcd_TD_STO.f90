module Mod_plcd_TD_StochasticTopologyOptimization
   use typre
   use Mod_PLCD
   use MPI
   use Mod_plcd_TD_StochasticTopologyOptimizationData
   use Mod_plcd_Stages
   use Mod_sparsegridFHS
   implicit none
   
contains

   subroutine plcd_TDSTO_turnon(a)
      use Mod_int2str
      implicit none
      class(PLCDProblem) :: a
      
      integer(ip) :: WorldMPIrank,WorldMPIsize,ierr,iint,sqrtNIntPoints,jint,ndime,n,ifun,adddimension
      
      !Initialize a transversal communicator
      !Color is by rank of the domain decomposition communicator
      call MPI_Comm_split(MPI_COMM_WORLD,a%MPIrank,a%MPIrank,a%STOData%TransversalMPIcomm,ierr)
      !Now I need to infere which problem I am solving
      !This info is is stored in kfl_multicomm
      !So I recompute it here
      call MPI_COMM_RANK( MPI_COMM_WORLD, WorldMPIrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, WorldMPIsize, ierr )
      a%STOData%NTransversalParallelProcesses = a%kfl_Multicomm
      
      
      
      call a%Mesh%GetNdime(ndime)
      a%STOData%StochasticDimensions = 0
      a%STOData%NForcePoints = 0
      a%STOData%FunctionIntegrationPointIniDimension = 0
      a%STOData%FunctionForcePointIniDimension = 0
      call a%Mesh%GetNdime(ndime)
      do ifun = 1,size(a%kfl_funty,1)
         if (a%kfl_funty(ifun,1) == -1) then
            if (a%STOData%kfl_SimultaneousLoads == 1) then
               a%STOData%FunctionForcePointIniDimension(ifun) = 1
               a%STOData%FunctionIntegrationPointIniDimension(ifun) = a%STOData%StochasticDimensions+1
               adddimension = count(a%funpa(ifun)%a(1:1)/=0.0_rp) + (ndime-1)*count(a%funpa(ifun)%a(2:2)/=0.0_rp)
               a%STOData%StochasticDimensions = a%STOData%StochasticDimensions + adddimension
            else
               a%STOData%FunctionIntegrationPointIniDimension(ifun) = 1
               adddimension = count(a%funpa(ifun)%a(1:1)/=0.0_rp) + (ndime-1)*count(a%funpa(ifun)%a(2:2)/=0.0_rp)
               a%STOData%StochasticDimensions = max(a%STOData%StochasticDimensions,adddimension)
               a%STOData%FunctionForcePointIniDimension(ifun) = a%STOData%NForcePoints +1
               a%STOData%NForcePoints = a%STOData%NForcePoints + 1
            endif
         elseif (a%kfl_funty(ifun,1) == -3) then
            !For 3d, just a rotation around z
            if (a%STOData%kfl_SimultaneousLoads == 1) then
               call runend('not ready for simulataneous loads')
            else
               a%STOData%FunctionIntegrationPointIniDimension(ifun) = 1
               adddimension = 1
               a%STOData%StochasticDimensions = a%STOData%StochasticDimensions + adddimension
               a%STOData%FunctionForcePointIniDimension(ifun) = a%STOData%NForcePoints +1
               a%STOData%NForcePoints = a%STOData%NForcePoints + 1
            endif   
         elseif (a%kfl_funty(ifun,1) == -2) then
            if (a%STOData%kfl_SimultaneousLoads == 1) then
               a%STOData%FunctionForcePointIniDimension(ifun) = 1
               a%STOData%FunctionIntegrationPointIniDimension(ifun) = a%STOData%StochasticDimensions+1
               adddimension = count(a%funpa(ifun)%a(1:1)/=0.0_rp) + (ndime-1)*count(a%funpa(ifun)%a(2:2)/=0.0_rp)
               a%STOData%StochasticDimensions = a%STOData%StochasticDimensions + adddimension*a%STOData%NFunctionForces(ifun)
            else
               a%STOData%FunctionIntegrationPointIniDimension(ifun) = 1
               adddimension = count(a%funpa(ifun)%a(1:1)/=0.0_rp) + (ndime-1)*count(a%funpa(ifun)%a(2:2)/=0.0_rp)
               a%STOData%StochasticDimensions = max(a%STOData%StochasticDimensions,adddimension)
               a%STOData%FunctionForcePointIniDimension(ifun) = a%STOData%NForcePoints +1
               a%STOData%NForcePoints = a%STOData%NForcePoints + a%STOData%NFunctionForces(ifun)
            endif   
         endif
      enddo
      if (a%STOData%kfl_SimultaneousLoads == 1) a%STOData%NForcePoints = 1
      
      call nwspgr_size ( gqn_order, a%STOData%StochasticDimensions, a%STOData%level, n )
      !if (n > a%STOData%NTransversalParallelProcesses) call runend('Wrong number of integration points multicomm')
      call a%Memor%alloc(a%StoData%StochasticDimensions,n,a%STOData%IntegrationPointPositions,'IntegrationPointPositions','sto')
      call a%Memor%alloc(n,a%STOData%IntegrationPointWeights,'IntegrationPointPositions','sto')
      
      
      call nwspgr ( gqn, gqn_order, a%StoData%StochasticDimensions, a%StoData%Level, n, a%STOData%NStochasticIntegrationPoints, a%STODAta%IntegrationPointPositions,a%STOData%IntegrationPointWeights )
      if (a%STOData%NStochasticIntegrationPoints /= a%STOData%NTransversalParallelProcesses/a%STOData%NForcePoints)  &
         call runend('Wrong number of Transversal Parallel Processes, should be '//int2str(a%STOData%NStochasticIntegrationPoints*a%STOData%NForcePoints))
         
         
         
      !Which is my integration Point
      !a%STOData%TransversalPoint = mod(WorldMPIrank,a%STOData%NTransversalParallelProcesses)+1    !old, worse performance
      a%STOData%TransversalPoint = a%MulticommColor+1
      
      
      a%STOData%IntegrationPoint = mod(a%STOData%TransversalPoint-1,a%STOData%NStochasticIntegrationPoints)+1
      a%STOData%ForcePoint = (a%STOData%TransversalPoint-1)/a%STOData%NStochasticIntegrationPoints+1
     
     
      
      !write(*,*) 'MPIrank: ', a%MPIrank, 'WorldMPIrank: ', WorldMPIrank, 'IntegrationPoint: ', a%STODATA%IntegrationPoint, 'ForcePoint: ', a%STOData%ForcePoint
   end subroutine
   
   
   subroutine plcd_TDSTO_Turnof(a)
      implicit none
      class(PLCDProblem) :: a
   
      call a%Memor%dealloc(size(a%STOData%IntegrationPointPositions,1),size(a%STOData%IntegrationPointPositions,2),a%STOData%IntegrationPointPositions,'IntegrationPointPositions','sto')
      call a%Memor%dealloc(size(a%STOData%IntegrationPointWeights),a%STOData%IntegrationPointWeights,'IntegrationPointPositions','sto')
      
      
   
   
   end subroutine
   
   subroutine plcd_TDSTO_Initial(a)
      use Mod_int2str
      implicit none
      class(PLCDProblem) :: a
      
      if (a%STOData%Alpha > 0) call a%TopologicalDerivative%SetComputeCompliance(.true.)
   end subroutine
   
   subroutine plcd_TDSTO_modifybcs(a,ipoin,Extvector,ifun)
      implicit none
      class(PLCDProblem), target :: a
      integer(ip) :: ipoin
      real(rp) :: ExtVector(:)
      
      integer(ip) :: npoin
      integer(ip) :: ifun, ifuntype
      real(rp),  parameter :: pi  = 4 * atan (1.0_8)
      real(rp) :: Vector(3) = 0.0_rp,Ortovec(3,2),vnorort,rotmat(3,3), CrossProdMatrix(3,3)
      real(rp) :: costeta,sinteta,rot
      
      real(rp) :: AngleStandardDeviation, newangle,vnor,angle
      integer(ip) :: ndime,i,j,ivec
      real(rp) :: newvnor, VecnorStandardDeviation
      
      integer(ip) :: iforce,inidimension,nforces,counter,IntegrationPointDimension,IntegrationPointForce
      real(rp) :: ForceVectors(3,3),ModifiedVector(3),IntegrationPointPositions(3)
      
      
      type(Stage), pointer :: s => NULL()
         
      s => a%CurrentReadingStage
      
      
      
      !ifun = a%Listener%param(2+a%ndofbc+1)
      if(ifun>0) then
         ifuntype = a%kfl_funty(ifun,1)
         if (ifuntype == -1) then
            if (a%STOData%FunctionForcePointIniDimension(ifun) == a%STOData%ForcePoint) then
            
               Vector(1:size(s%bvess,1)) = ExtVector
               VecnorStandardDeviation = a%funpa(ifun)%a(1)
               AngleStandardDeviation = a%funpa(ifun)%a(2)
               call a%Mesh%GetNdime(ndime)
               if (ndime == 2) then
                  call vecnor(Vector,ndime,vnor,2)
                  angle = atan2(Vector(2),Vector(1))
                                       
                  !We start with a first approximation where we 
                  !integrate between minus 2 times the standard dev
                  !and 2 times the standard dev
                  counter = 0
                  newangle = angle
                  newvnor  = vnor
                  if (VecnorStandardDeviation /= 0.0_rp) then
                     newvnor = vnor + VecnorStandardDeviation*a%STOData%IntegrationPointPositions(a%STOData%FunctionIntegrationPointIniDimension(ifun)+counter,a%STOData%IntegrationPoint)
                     counter = counter +1
                  endif
                  if (AngleStandardDeviation /= 0.0_rp) then
                     newangle = angle + AngleStandardDeviation*a%STOData%IntegrationPointPositions(a%STOData%FunctionIntegrationPointIniDimension(ifun)+counter,a%STOData%IntegrationPoint)
                     counter = counter +1
                  endif
                  
                  if (newvnor < 0.0_rp) call runend('negative vnor')
                  !New vector
                  ExtVector(1) = newvnor*cos(newangle)
                  ExtVector(2) = newvnor*sin(newangle)
               elseif (ndime == 3) then
                  counter = 0
                  IntegrationPointDimension = a%STOData%FunctionIntegrationPointIniDimension(ifun)
                  IntegrationPointForce = a%STODAta%FunctionForcePointIniDimension(ifun)
            
            
                  ForceVectors(:,1) = ExtVector(:)
                  
                  ExtVector = 0.0_rp
                  do iforce = 1,1
                     IntegrationPointPositions = 0.0_rp
                     if ((IntegrationPointForce == a%STOData%ForcePoint) .or. (a%STOData%kfl_SimultaneousLoads == 1)) then
                        VecnorStandardDeviation = a%funpa(ifun)%a(1)
                        if (VecnorStandardDeviation /= 0.0_rp) then
                           IntegrationPointPositions(1) = a%STOData%IntegrationPointPositions(IntegrationPointDimension,a%STOData%IntegrationPoint)
                           IntegrationPointDimension = IntegrationPointDimension+1
                        endif   
                        
                        AngleStandardDeviation = a%funpa(ifun)%a(2)
                        if (AngleStandardDeviation /= 0.0_rp) then
                           IntegrationPointPositions(2:3) = a%STOData%IntegrationPointPositions(IntegrationPointDimension:IntegrationPointDimension+1,a%STOData%IntegrationPoint)
                           IntegrationPointDimension = IntegrationPointDimension+2
                        endif
                     
                        call Stochastic3dVectorModification(ForceVectors(:,iforce),AngleStandardDeviation,VecnorStandardDeviation,IntegrationPointPositions,ModifiedVector)
                        ExtVector = ExtVector + ModifiedVector
                     endif
                     IntegrationPointForce = IntegrationPointForce +1
                  enddo

               endif
            endif
         
         elseif (ifuntype == -2) then
            nforces = a%STOData%NFunctionForces(ifun)
            counter = 0
            IntegrationPointDimension = a%STOData%FunctionIntegrationPointIniDimension(ifun)
            IntegrationPointForce = a%STODAta%FunctionForcePointIniDimension(ifun)
            
            
            ForceVectors(:,1:nforces) = a%STOData%FunctionForces(:,1:nforces,ifun)
            
            ExtVector = 0.0_rp
            do iforce = 1,nforces
               IntegrationPointPositions = 0.0_rp
               if ((IntegrationPointForce == a%STOData%ForcePoint) .or. (a%STOData%kfl_SimultaneousLoads == 1)) then
                  VecnorStandardDeviation = a%funpa(ifun)%a(1)
                  if (VecnorStandardDeviation /= 0.0_rp) then
                     IntegrationPointPositions(1) = a%STOData%IntegrationPointPositions(IntegrationPointDimension,a%STOData%IntegrationPoint)
                     IntegrationPointDimension = IntegrationPointDimension+1
                  endif   
                  
                  AngleStandardDeviation = a%funpa(ifun)%a(2)
                  if (AngleStandardDeviation /= 0.0_rp) then
                     IntegrationPointPositions(2:3) = a%STOData%IntegrationPointPositions(IntegrationPointDimension:IntegrationPointDimension+1,a%STOData%IntegrationPoint)
                     IntegrationPointDimension = IntegrationPointDimension+2
                  endif
               
                  call Stochastic3dVectorModification(ForceVectors(:,iforce),AngleStandardDeviation,VecnorStandardDeviation,IntegrationPointPositions,ModifiedVector)
                  ExtVector = ExtVector + ModifiedVector
               endif
               IntegrationPointForce = IntegrationPointForce +1
            enddo
         elseif(ifuntype == -3) then
            Vector(1:size(s%bvess,1)) = ExtVector
            VecnorStandardDeviation = a%funpa(ifun)%a(1)
            AngleStandardDeviation = a%funpa(ifun)%a(2)
            call a%Mesh%GetNdime(ndime)
            
            if (ndime == 3) then
               counter = 0
               IntegrationPointDimension = a%STOData%FunctionIntegrationPointIniDimension(ifun)
               IntegrationPointForce = a%STODAta%FunctionForcePointIniDimension(ifun)
         
         
               ForceVectors(:,1) = ExtVector(:)
               
               ExtVector = 0.0_rp
               do iforce = 1,1
                  IntegrationPointPositions = 0.0_rp
                  if ((IntegrationPointForce == a%STOData%ForcePoint)) then
                     AngleStandardDeviation = a%funpa(ifun)%a(2)
                     if (AngleStandardDeviation /= 0.0_rp) then
                        IntegrationPointPositions(1) = a%STOData%IntegrationPointPositions(IntegrationPointDimension,a%STOData%IntegrationPoint)
                     endif
                     call Stochastic3dVectorRotationZ(ForceVectors(:,iforce),AngleStandardDeviation,IntegrationPointPositions,ModifiedVector)
                     ExtVector = ExtVector + ModifiedVector
                  endif
                  IntegrationPointForce = IntegrationPointForce +1
               enddo
            else
               call runend('only for 3d in here')
            
            endif
            
      
         
         
         
         
         
         
         
         endif
      endif

   end subroutine   
   
   
   subroutine Stochastic3dVectorModification(OutVector,AngleStandardDeviation,VecnorStandardDeviation,StochasticNodalPoint,ModifiedVector)
      implicit none
      real(rp) :: OutVector(3), AngleStandardDeviation, VecnorStandardDeviation,StochasticNodalPoint(3),ModifiedVector(3)

      
      real(rp),  parameter :: pi  = 4 * atan (1.0_8)
      real(rp) :: Ortovec(3,2),vnorort,rotmat(3,3), CrossProdMatrix(3,3)
      real(rp) :: costeta,sinteta,rot
      
      real(rp) :: newangle,vnor,angle
      integer(ip), parameter :: ndime = 3
      integer(ip) ::i,j,ivec
      real(rp) :: newvnor,Vector(3),auxVector(3,2)
      integer(ip) :: itest,itestorder(2,2),ivecaux
      
      
      
      Vector = OutVector
      
      !1rst, we build 2 vectors orthogonal to the original one (and to each other)
               
      call vecnor(Vector,ndime,vnor,2)
      Vector = Vector/vnor
      OrtoVec(:,:) = 0.0_rp
      if (Vector(1) > 1e-6_rp) then
         OrtoVec(2,1) = 1.0_rp
      else
         OrtoVec(1,1) = 1.0_rp
      endif
      
      !Ortogonalize ortovec 1
      OrtoVec(1:ndime,1) = OrtoVec(1:ndime,1) - dot_product(Vector(1:ndime),OrtoVec(1:ndime,1))*Vector(1:ndime)
      call vecnor(Ortovec(1,1),ndime,vnorOrt,2)
      Ortovec(1:ndime,1) = Ortovec(1:ndime,1)/vnorOrt
      
      !Find ortovec2 by using a cross product
      call vecpro(Vector,Ortovec(1,1),Ortovec(1,2),ndime)
      
      
      !Since there are 2 possible rotation orders, we are going to do both of them, and then take the mean value
      auxVector(:,1) = Vector
      auxVector(:,2) = Vector
      
      itestorder(1,1) = 1
      itestorder(2,1) = 2
      itestorder(1,2) = 2
      itestorder(2,2) = 1
      
      do itest = 1,2
         do ivecaux = 1,2
            ivec = itestorder(ivecaux,itest)
            !Rotation1
            rot = AngleStandardDeviation*StochasticNodalPoint(1+ivec)
            costeta = cos(rot)
            sinteta = sin(rot)
            
            
            rotmat = 0.0_rp
            forall(i=1:3)
               rotmat(i,i) = costeta
            end forall
            CrossProdMatrix = 0.0_rp
            CrossProdMatrix(1,2) = -Ortovec(3,ivec)
            CrossProdMatrix(1,3) = Ortovec(2,ivec)
            CrossProdMatrix(2,1) = Ortovec(3,ivec)
            CrossProdMatrix(2,3) = -Ortovec(1,ivec)
            CrossProdMatrix(3,1) = -Ortovec(2,ivec)
            CrossProdMatrix(3,2) = Ortovec(1,ivec)
            rotmat = rotmat + sinteta*CrossProdMatrix
            forall(i=1:3,j=1:3)
            rotmat(i,j) = rotmat(i,j) +  Ortovec(i,ivec)*Ortovec(j,ivec)*(1-costeta)
            end forall
            
            !New vector
            auxVector(:,itest) = matmul(rotmat,auxVector(:,itest))
         enddo
      enddo
      
      !We combine both
      Vector = auxVector(:,1)+auxVector(:,2)
      Vector = Vector/sqrt(dot_product(Vector,Vector))
      
      !Modify the vector modulus
      Vector = Vector*vnor*(1+VecnorStandardDeviation*StochasticNodalPoint(1))
               
               
      ModifiedVector = Vector
               
               
   
   
   
   
   
   
   end subroutine
   
   subroutine Stochastic3dVectorRotationZ(OutVector,AngleStandardDeviation,StochasticNodalPoint,ModifiedVector)
      implicit none
      real(rp) :: OutVector(3), AngleStandardDeviation, VecnorStandardDeviation,StochasticNodalPoint(3),ModifiedVector(3)

      
      real(rp),  parameter :: pi  = 4 * atan (1.0_8)
      real(rp) :: Ortovec(3,2),vnorort,rotmat(3,3), CrossProdMatrix(3,3)
      real(rp) :: costeta,sinteta,rot
      
      real(rp) :: newangle,vnor,angle
      integer(ip), parameter :: ndime = 3
      integer(ip) ::i,j,ivec
      real(rp) :: newvnor,Vector(3),auxVector(3,2)
      integer(ip) :: itest,itestorder(2,2),ivecaux
      
      Vector = OutVector
      
      rot = AngleStandardDeviation*StochasticNodalPoint(1)
      costeta = cos(rot)
      sinteta = sin(rot)
      
      rotmat = 0.0_rp
      rotmat(1,1) = costeta
      rotmat(1,2) = -sinteta
      rotmat(2,1) = sinteta
      rotmat(2,2) = costeta
      rotmat(3,3) = 1.0_rp
      
      ModifiedVector = matmul(rotmat,Vector)
   end subroutine


   
   subroutine plcd_TDSTO_AddStochasticTopologicalDerivative(a)
      implicit none
      class(PLCDProblem), target :: a
      
      real(rp), pointer :: TopologicalDerivative(:),ComplianceArray(:)
      
      real(rp), allocatable :: TopoMean(:,:), TopoVariance(:,:), TopoStandardDeviation(:),auxTopoComplianceMean(:,:),auxTopoComplianceVariance(:,:)
      
      integer(ip) :: ierr,idime,ndime
      real(rp) :: Compliance, auxCompliance,TotalWeight,aux
      integer(ip) :: MPIrank,ipoin,npoin
      
      type(Stage), pointer :: s => NULL()
      real(rp) :: auxComplianceVariance,ComplianceMean, ComplianceVariance,ComplianceStandardDeviation
      integer(ip) :: npoinLocal
     
      !if (a%istep < 50) return
      
      !call a%TopologicalDerivative%PostprocessChi(a%FilePostpr,0,a%istep,a%ctime)
      
      call a%TopologicalDerivative%GetPointerToTopologicalDerivativeArray(TopologicalDerivative)
      
      call a%Memor%Alloc(size(TopologicalDerivative),2,TopoMean,'TopoMean','STO')
      call a%Memor%Alloc(size(TopologicalDerivative),2,TopoVariance,'Topovariance','STO')
      call a%Memor%Alloc(size(TopologicalDerivative),TopoStandardDeviation,'Topovariance','STO')
      
       call a%Memor%Alloc(size(TopologicalDerivative),2,auxTopoComplianceMean,'TopoMean','STO')
       call a%Memor%Alloc(size(TopologicalDerivative),2,auxTopoComplianceVariance,'Topovariance','STO')
      
      
      call a%Mesh%GetNdime(ndime)
      
      !TopologicalDerivativeMean
      TotalWeight= a%STOData%IntegrationPointWeights(a%STOData%IntegrationPoint)
      TopoMean(:,1) = TopologicalDerivative*TotalWeight
      
      call  MPI_AllREDUCE( TopoMean(:,1), TopoMean(:,2), size(TopologicalDerivative), MPI_REAL8, MPI_SUM,a%STOData%TransversalMPIcomm, ierr ) 
      
      
      !heyyyy
      !TopoMean(:,2) = TopologicalDerivative
      !CALL MPI_BCAST(TopoMean(:,2),size(TopoMean,1),MPI_REAL8,1,a%STOData%TransversalMPIcomm,ierr)
      
      
      !Total Compliance Mean  (ComplianceMean (E[J])
      call a%TopologicalDerivative%GetCompliance(Compliance)
      auxCompliance = Compliance*TotalWeight
      call MPI_AllREDUCE( auxCompliance, ComplianceMean, 1, MPI_REAL8, MPI_SUM, a%STOData%TransversalMPIcomm, ierr )
      if (a%MPIrank == 0) then
         call MPI_COMM_RANK(a%STOData%TransversalMPIcomm, MPIrank, ierr )
         if (MPIrank == 0) then
            write(a%lun_outph,*) '                       Total Compliance Mean: ', ComplianceMean
         endif
      endif
      !Total Compliance Variance (Var(J))
      auxComplianceVariance  = (Compliance - ComplianceMean)**2*TotalWeight
      call MPI_ALLREDUCE( auxComplianceVariance, ComplianceVariance, 1, MPI_REAL8, MPI_SUM,a%STOData%TransversalMPIcomm, ierr )
      ComplianceStandardDeviation = sqrt(ComplianceVariance)
      if (a%MPIrank == 0) then
         call MPI_COMM_RANK(a%STOData%TransversalMPIcomm, MPIrank, ierr )
         if (MPIrank == 0) then
            write(a%lun_outph,*) '                       Total Compliance Standard Deviation: ', ComplianceStandardDeviation
         endif
      endif
         
      !endif
      
      
      
      if (a%STOData%Alpha == 0 .or. a%istep < a%STOData%VarianceDelaySteps  ) then
         TopologicalDerivative = TopoMean(:,2)
      elseif (a%STOData%Alpha /= 0) then
         !TopologicalDerivative of the  Variance
         call a%Mesh%GetNpoin(npoin)
         TopoVariance(:,1) = +(2*Compliance*TopologicalDerivative(:)*TotalWeight-2*ComplianceMean*TopoMean(:,2)*TotalWeight)
         call  MPI_AllREDUCE( TopoVariance(:,1), TopoVariance(:,2), size(TopologicalDerivative), MPI_REAL8, MPI_SUM,a%STOData%TransversalMPIcomm, ierr )
         
         if (ComplianceStandardDeviation /= 0.0_rp) then
            TopoStandardDeviation(:) = (0.5_rp/ComplianceStandardDeviation)*TopoVariance(:,2)
         else
            TopoStandardDeviation(:) = 0.0_rp
         endif
         
         TopologicalDerivative = TopoMean(:,2) + a%STOData%Alpha*TopoStandardDeviation
      endif
      
      
      call a%Memor%DeAlloc(size(TopologicalDerivative),2,TopoMean,'TopoMean','STO')
      call a%Memor%DeAlloc(size(TopologicalDerivative),2,TopoVariance,'Topovariance','STO')
      
      call a%Memor%DeAlloc(size(TopologicalDerivative),2,auxTopoComplianceMean,'TopoMean','STO')
      call a%Memor%DeAlloc(size(TopologicalDerivative),2,auxTopoComplianceVariance,'Topovariance','STO')
      
      call a%Memor%dealloc(size(TopologicalDerivative),TopoStandardDeviation,'Topovariance','STO')
      
      
      
      !Postprocess of the total mean compliance and total compliance variance
      !For this we use the nodal values directly and use vmass as a measure of the volume associated to each node
      call a%Mesh%GetNpoinLocal(npoinLocal)
      
      
      
      
      
      
      
      
      
         
      
      
   end subroutine
end module
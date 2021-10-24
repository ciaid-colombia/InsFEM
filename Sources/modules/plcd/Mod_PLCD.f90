module Mod_PLCD
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_BroadCastBuffer
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_plcd_Material
   use Mod_plcd_Stages
   use Mod_plcd_SIMPData
   use Mod_plcd_TDData
   use Mod_plcd_TD_Derivative
   use Mod_plcd_TD_StochasticTopologyOptimizationData
   implicit none
   private
   public PLCDProblem, PLCDProblem_Const,zeplcd,MatArray,MatElementArray

   real(rp)                :: zeplcd = 1e-12

   type :: MatArray
      class(PLCDMaterial), pointer :: p => NULL()
   end type

   type :: MatElementArray
      class(ElementMaterialData), pointer :: p => NULL()
   end type

   type, extends(PhysicalProblem) ::  PLCDProblem

      integer(ip) :: NumberOfMaterials
      type(MatArray), allocatable :: Materials(:)  !This array is of dimension the number of materials

      type(MatElementArray), pointer :: ElementMaterialsData(:) => NULL()

      !Auxiliary, read only
      character(5),allocatable :: AuxRead_MaterialTypeList(:)
      !integer(ip), allocatable :: AuxRead_ElementMaterialList(:)

      real(rp), allocatable :: Displacement(:,:,:)
      real(rp), allocatable :: Velocity(:,:,:)
      real(rp), allocatable :: Acceleration(:,:,:)

      real(rp), allocatable :: InternalForcesVector(:,:), ExternalForcesVector(:,:), ResidualForcesVector(:,:)
      
      real(rp) :: ForcesResidualNorm = 1.0_rp

      integer(ip) :: NumberOfStages = 1, CurrentStage = 1

      type(Stage), allocatable :: Stages(:)
      type(Stage), pointer :: CurrentReadingStage => NULL()
      type(Stage), pointer :: cs => NULL()
      type(Substage), pointer :: css => NULL()

      type(r2p), allocatable :: Stress(:)
      real(rp), allocatable :: MeanStress(:)
      real(rp), allocatable :: GMeanStress(:)
      real(rp), allocatable :: MeanConst(:,:)
      logical :: ComputeMeanflag = .false.

      !Postprocess
      integer(ip) :: lun_trap
      character(150) :: fil_trap

      !Strain
      type(r2p), allocatable :: Strain(:)
      logical :: PostprocessStrain = .false.

      !PostprocessMatData at each ieration
      integer(ip) :: kfl_PostprocessMatDataAtEachIteration = 0
      
      !Postprocess Displacements at each iteration
      integer(ip) :: kfl_PostprocessDisplacementAtEachIteration = 0

      !Is the element size required by some element material datas
      logical :: IsElementSizeRequiredByEMDs = .false.
      
      !Element rotator flag
      integer(ip) :: kfl_ElementRotators = 0

      !SmoothedDisplacementGradient
      logical :: UseSmoothedDisplacementGradient = .false.
      real(rp), allocatable :: SmoothedDisplacementGradient(:,:,:)

      !UpFormulation
      logical :: UseUPFormulation = .false., UPStoreSubscales = .false.
      real(rp), allocatable :: pressure(:,:)
      real(rp) :: staco(1)
      real(rp), allocatable :: UPResidualProjection(:,:)
      type(r2p), allocatable :: UPSubscales(:)
      integer(ip) :: ErrorEstimatorTypeOfSubscales = 0  !0: Orthogonal, 1: ASGS
      integer(ip) :: kfl_confi = 0, nodpr = 0

      !Topoly Optimization
      integer(ip) :: kfl_TopologyOptimization = 0 !Default is do not do it

      type(SIMPDataType) :: SIMPData
      type(TDDataType)   :: TDData
      class(TD_Derivative), pointer :: TopologicalDerivative
      !Stochastic analysis of topology optimization
      integer(ip) :: kfl_StochasticTopologyOptimization = 0
      type(TD_StochasticData) :: STOData
      
      !Adaptive mesh refinement, layers
      integer(ip) :: NLayersRefinement = 0, GeneralRefinementLevels = 0, InterfaceRefinementLevels = 0
      
      !Large Strain Deformations
      integer(ip) :: kfl_LargeStrains = 0 !0: Linear Elastic, 1: Updated Formulation NonLinear Elastic, 2: Total Formulation Nonlinear Elastic
      
      !Time Integration Scheme: Newmark
      integer(ip) :: kfl_TransientProblem = 0   !0: Stationary problem, 1:Transient problem ( Newmark scheme)
      real(rp) :: Beta = 0.0_rp                 !Beta parameter for Newmark scheme (usually 0.5)
      real(rp) :: Gamma = 0.0_rp                !Gamma parameter for Newmark scheme (usually 0.5)
      
      !Line Search 
      integer(ip) :: kfl_LineSearch = 0                     !0: Disabled, 1:Enabled
      integer(ip) :: kfl_InternalLineSearchLoop = 0         !0: Out of the loop, 1: Inside the loop
      real(rp) :: LineSearchRadiusOfConvergence = 0.5_rp    !Radius of convergence
      integer(ip) :: LineSearchMaxIterations = 10           !Maximum Line Search iterations
      
      !Gravity Force
      integer(ip) :: kfl_GravityForce = 0 !0:Disabled, 1:Enabled
      real(rp) :: gravity(3) = 0.0_rp     ! Gravity vector 
      
      !Rotating Frame of Reference
      integer(ip) :: kfl_RotatingFrame = 0   !0:Disabled, 1:Newmark Implicit Time scheme, 2:Explicit Time integrator scheme
      real(rp) :: angvelocitynorm2 = 0.0_rp  !Angular Velocity Norm ^2
      real(rp) :: angvelocity(3) = 0.0_rp    !Angular Velocity Vector
      
      !Fluid-Structure Interaction
      integer(ip) :: kfl_FSI = 0                               !0:Disabled 1:Enabled
      real(rp), pointer     :: fluidtraction(:,:) => NULL()    !Tractions Given by the fluid
      real(rp), allocatable :: btraction(:,:)                  !Traction in the boundary to check convergence
                                     

contains

      !specific deferred procedures from PhysicalProblemr
      procedure :: SetExmod          => plcd_SetExmod
      procedure :: SetNdofn          => plcd_SetNdofn
      procedure :: SetNdofbc         => plcd_SetNdofbc

      procedure :: SpecificReaphy    => plcd_reaphy   !PHYSICAL_PROBLEM
      procedure :: SpecificReanut    => plcd_reanut   !NUMERICAL_TREATMENT
      procedure :: SpecificReaous    => plcd_reaous   !OUTPUT_&_POST-PROCESS
      procedure :: SpecificReaMPI    => plcd_reampi   !Distribuye a todos los treads

      procedure :: SpecificReabcs =>  plcd_Reabcs
      procedure :: SpecificReadOnNodes =>  plcd_NULLSUB
      procedure :: SpecificReadSource =>  plcd_NULLSUB
      procedure :: SpecificReadOnBoundaries =>  plcd_NULLSUB
      procedure :: SpecificReadOnElements =>  plcd_ReadOnElements

      procedure :: SpecificBounod    => plcd_bounod
      procedure :: SpecificSolite    => plcd_solite
      procedure :: SpecificTurnon    => plcd_turnon
      procedure :: SpecificTurnof    => plcd_turnof
      procedure :: SpecificBegite    => plcd_begite
      procedure :: SpecificBegste    => plcd_begste
      procedure :: SpecificUpdbcs    => plcd_updbcs
      procedure :: SpecificIniunk    => plcd_iniunk
      procedure :: SpecificEndste    => plcd_endste
      procedure :: SpecificEndite    => plcd_endite
      procedure :: SpecificCrankNicolsonEndste => plcd_NULLSUB
      procedure :: SpecificExaerr    => plcd_NULLSUB
      procedure :: SpecificReadOnFunctions => plcd_ReadOnFunctions
      procedure :: SpecificScatterOnFunctions => plcd_ScatterOnFunctions


      procedure :: SpecificOuterr    => plcd_outerr
      procedure :: Memall            => plcd_memall
      procedure :: InnerResiduals    => plcd_InnerResiduals
      procedure :: Getste            => plcd_Getste
      procedure :: Cvgunk            => plcd_cvgunk
      procedure :: Output            => plcd_output
      procedure :: plcd_outtpo
      procedure :: Elmope            => plcd_Elmope
      procedure :: EnditeElmope      => plcd_EnditeElmope
      procedure :: EndsteElmope      => plcd_EndsteElmope
      procedure :: Restart           => plcd_NULLSUBitaskintentin
      procedure :: SpecificRestart   => plcd_NULLSUBitaskintentin

      !Boundary operations
      procedure :: Bouope            => plcd_bouope
      procedure :: EndBouope         => plcd_EndBouope

      !Point Operations
      procedure :: PoinOpe           => plcd_poinope


      !OverWrittenProcedures
      procedure :: ReadOnNodes       => plcd_ReadOnNodes
      procedure :: ReadOnBoundaries  => plcd_ReadOnBoundaries
      procedure :: ReadSource        => plcd_ReadSource

      procedure :: OTAReadOnNodesDo => plcd_OTAReadOnNodesDo
      procedure :: OTAReadOnBoundariesDo => plcd_OTAReadOnBoundariesDo
      procedure :: OTAReadSourceDo => plcd_OTAReadSourceDo

      !Adaptive Mesh Refinement
      procedure :: PreRefine => plcd_PreRefine
      procedure :: SpecificRefine => plcd_Refine
      procedure :: GetRefinementCriteria => plcd_GetRefinementCriteria

      !Mean Computes
      procedure :: SetComputeMean
      procedure :: GetMeanStress
      procedure :: GetMeanConst

      !Stages
      procedure :: DeallocStages
      procedure :: SetStages
      
      !Getters
      procedure :: GetDisplacementArray
      procedure :: GetVelocityArray
      procedure :: GetTractionArray
      procedure :: SetExternalTractionArray

   end type

   interface

      subroutine plcd_SetExmod(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_SetNdofn(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_SetNdofbc(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_reaphy(a,itask)
         use typre
         import PLCDProblem
         implicit none
         integer(ip) :: itask
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_reanut(a,itask)
         use typre
         import PLCDProblem
         implicit none
         integer(ip) :: itask
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_reaous(a,itask)
         use typre
         import PLCDProblem
         implicit none
         integer(ip) :: itask
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_reaMPI(a)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_Reabcs(a,itask,kflag)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine
      
      subroutine plcd_ReadOnFunctions(a,ifunc)
         import
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: ifunc
      end subroutine
      
      subroutine plcd_ScatterOnFunctions(a)
         import
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_updbcs(a)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_ReadOnNodes(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_ReadSource(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_ReadOnBoundaries(a,knodb)
         import
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: knodb(:)
      end subroutine

      subroutine plcd_ReadOnElements(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_bounod(a,itask)
         use typre
         import PLCDProblem
         implicit none
         integer(ip) :: itask
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_outerr(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_memall(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_begste(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_begite(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_endite(a,itask)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: itask
      end subroutine

      subroutine plcd_cvgunk(a,itask)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine

      subroutine plcd_outtpo(a)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_output(a,itask)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: itask
      end subroutine

      subroutine plcd_endste(a,itask)
         use typre
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: itask
      end subroutine

      subroutine plcd_iniunk(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_Elmope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_bouope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_Poinope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_turnon(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_turnof(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_EnditeElmope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         character(6) :: task
      end subroutine

      subroutine plcd_EndsteElmope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
         character(6) :: task
      end subroutine

      subroutine plcd_EndBouope(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_OTAReadOnNodesDo(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_OTAReadOnBoundariesDo(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine plcd_OTAReadSourceDo(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a
      end subroutine


      subroutine plcd_Getste(a,dtinv)
         import
         implicit none
         class(PLCDProblem) :: a
         real(rp) :: dtinv

      end subroutine

      subroutine plcd_InnerResiduals(a,rinsi)
         use typre
         import
         implicit none
         class(PLCDProblem) :: a
         real(rp) :: rinsi
      end subroutine

      subroutine plcd_solite(a,itask)
         import
         implicit none
         class(PLCDProblem) :: a
         integer(ip) :: itask
      end subroutine

      subroutine plcd_Refine(a,itask)
         import
         implicit none
         class(PLCDProblem) :: a
         character(6) :: itask
      end subroutine


      subroutine plcd_GetRefinementCriteria(a,markel)
         import
         class(PLCDProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine plcd_Prerefine(a)
         import
         implicit none
         class(PLCDProblem) :: a
      end subroutine

      subroutine DeallocStages(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a

      end subroutine

      subroutine SetStages(a)
         import PLCDProblem
         implicit none
         class(PLCDProblem) :: a

      end subroutine

  end interface


  interface PLCDProblem_Const
      procedure constructor
  end interface PLCDProblem_Const

contains

   function constructor()
         class(PLCDProblem), pointer :: constructor

         allocate(constructor)

   end function constructor

   subroutine plcd_NULLSUB(a)
      implicit none
      class(PLCDProblem) :: a

   end subroutine

   subroutine plcd_NULLSUB_knode(a,knode)
      implicit none
      class(PLCDProblem) :: a
      integer(ip) :: knode(:)

   end subroutine

   subroutine plcd_NULLSUBitask(a,itask)
      implicit none
      class(PLCDProblem) :: a
      integer(ip) :: itask

   end subroutine

   subroutine plcd_NULLSUBitaskc(a,itask)
      implicit none
      class(PLCDProblem) :: a
      character(6) :: itask

   end subroutine

   subroutine plcd_NULLSUBitaskintentin(a,itask)
      implicit none
      class(PLCDProblem) :: a
      integer(ip), intent(in) :: itask

   end subroutine

   subroutine plcd_NULLSUBitaskkflag(a,itask,kflag)
      implicit none
      class(PLCDProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag

   end subroutine

   subroutine plcd_NULLSUBdtinv(a,dtinv)
      implicit none
      class(PLCDProblem) :: a
      real(rp) :: dtinv

   end subroutine

   subroutine SetComputeMean(a,Meanflag)
      implicit none
      class(PLCDProblem) :: a
      logical :: Meanflag

      a%ComputeMeanflag = Meanflag

   end subroutine

   subroutine GetMeanStress(a,StressMean)
      implicit none
      class(PLCDProblem), target :: a
      real(rp), pointer :: StressMean(:)

      StressMean => a%GMeanStress

   end subroutine

   subroutine GetMeanConst(a,RVEConst)
      implicit none
      class(PLCDProblem), target :: a
      real(rp), pointer :: RVEConst(:,:)

      RVEConst => a%MeanConst

   end subroutine
   
   subroutine GetDisplacementArray(a,displacement)
      implicit none
      class(PLCDProblem), target :: a
      real(rp), pointer :: displacement(:,:)
      
      displacement => a%Displacement(:,:,1)
   end subroutine

   subroutine GetVelocityArray(a,velocity)
      implicit none
      class(PLCDProblem), target :: a
      real(rp), pointer :: velocity(:,:)
      
      velocity => a%Velocity(:,:,1)
   end subroutine

   subroutine GetTractionArray(a,traction)
      implicit none
      class(PLCDProblem), target :: a
      real(rp), pointer :: traction(:,:)

      traction => a%btraction

   end subroutine
   
   subroutine SetExternalTractionArray(a,traction)
      implicit none
      class(PLCDProblem) :: a
      real(rp), target :: traction(:,:)

      a%fluidtraction => traction

   end subroutine

end module Mod_PLCD

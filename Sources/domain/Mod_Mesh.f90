module Mod_mesh
   use typre
   use Mod_MpiObject
   use Mod_Memor
   use Mod_Listen
   use Mod_ParallelLibraryInterface
   use Mod_ParallelOrderingInterface
   use Mod_ParallelCommunicatorInterface
   use Mod_ParallelSystemInterface
   use Mod_Element
   use Mod_Timer
   use Mod_Octree
   use Mod_AdaptiveInterface
   use Mod_ElementDataStructures
   implicit none

   private
   public :: FemMesh  , Global2Local, Local2Global
   
   type :: MeshTimer
      type(Timer) :: &
         Total,&
         reaLnods,&
         Partitionate,&
         LocalDimensionsAndIRenumbering,&
         GhostPoints,&
         BuildLocalOrdering,&
         reaCoord,&
         LnodsToLocal,&
         LocalGraph,&
         BuildArrayCommunicator,&
         GhostCoord,&
         PeriodicBoundaryConditions,&
         Dombou,&
         ExtnorLpoty,&
         ComputeVmass,&
         Blbopo,&
         Refine,&
            TestGhostCommunicate,&
            TestAllToAll,&
            TestReduce,&
            TestBcast,&
            TestBandwidth
          
   end type

   
   type, extends(MpiObject) :: FemMesh 
      
      !Timers
      type(MeshTimer)       :: Timer
      
      
      !Flags
      integer(ip) :: &
         kfl_mnelt,&               !Nelty is given by user/automatically
         kfl_exnpo,&               !Post-process normals 
         kfl_TestCommunications, & !Set to 1 if we want to test communications
         kfl_MeshInputType = 0, &  !0: reads from disk, 1: Manufactured 
         kfl_ManuMeshType = 0      !0: Cartesian domains, 1: Polar domain, 2: Cylindrical domain 
      
     
      ! Logical units
      integer(ip) ::&
         lun_outpu_dom,&       !Output for the Mesh Information   
         lun_solve_dom
   
      !Global variables (to be used by MPIroot only)
      integer(ip) :: &
         gnelem,&                      !# of elements
         gnpoin,&                      !# of points
         gnpoinghost,&                 !# of ghost points
         MaxNpoinGhostPerProc, &       
         gnskew,&                      !# of skew systems
         MaxNelemPerProc,&
         MaxNpoinPerProc
         
      !Temporary arrays for reading data   
      integer(ip), allocatable :: &
         gpnods(:),glnods(:)    !Element connectivities
        
      integer(ip), allocatable :: & 
         gPointGlobNumber(:),&        !Global array of renumbering, to be used by MPIroot
         igPointGlobNumber(:),&       !Inverse Global array of renumbering
         gPointProcNumber(:),&        !Global array of processor to which a point belongs, to be used by MPIroot 
         lPointGlobNumber(:),&        !Local array of renumbering
         ilPointGlobNumber(:),&       !Inverse Local array of renumbering
         lPointProcNumber(:)          !Local array of processor to which a point belongs, to be used by MPIroot 
         
      !Temporary arrays for reading partitioned (RP)
      !Also RP_LocalOrdering which is written lower because of INTEL
      integer(ip) :: RP_npoinLocal, RP_npoinGhost,RP_npoin
      integer(ip), allocatable :: RP_NpoinOthersNeed(:), RP_NpoinINeed(:)
      integer(ip) :: RP_poin0
      
      !Temporary arrays for testing bandwidth
      integer(ip) :: BandwidthTestSize
      real(rp)    :: Bandwidth
     
      !Local variables
      integer(ip) :: &
         nelem,&                   ! # of elements
         elem0,&                   ! First element for current proc
         npoin,&                   ! # of nodal points
         poin0,&                   ! First point for current proc (global numbering) -1
         npoinGhost,&              ! Number of ghost points
         npoinLocal,&              ! Number of local points
         numGeoBlocks,&            ! # of blocks for output
         ndime,&                   ! # of space dimensions
         nelty,&                   ! # of element types
         mnode,&                   ! Max # of nodes per element
         mnodb,&                   ! Max # of nodes per boundary
         mface,&                   ! Max # of faces per element
         mgaus,&                   ! Max # of Gauss points per element
         mgaub,&                   ! Max # of Gauss points per boundary
         mlapl,&                   ! Max # of llapl(ielty), ielty=1,nelty
         ntens,&                   ! # of components of symmetric tensors
         nrule,&                   ! # of integration rule
         nbopo,&                   ! # of boundary points
         nboun,&                   ! # of boundary elements
         ndimb,&                   ! # of space dimensions-1
         nskew,&                   ! # of skew systems
         mitsm,&                   ! Max # of smoothing iterations
         nfic,&                    ! Number of fictitious nodes
         nbody                     ! # of bodies

      logical :: kfl_ReadBlock= .false.
      
      logical :: kfl_UseElementDataStructures = .false.
      type(ElementDataStructure), allocatable :: ElementDataStructure(:),ElementDataStructureClosedRule(:)  
      type(BoundaryDataStructure), allocatable :: BoundaryDataStructure(:), BoundaryDataStructureClosedRule(:)

      real(rp) :: &
         tolsm,&                   ! Smoothing tolerance
         vodom,&                   ! Measure of the domain
         scale(3)               ! Geometric scale factors   
      
      !Sizes for Element types   
      integer(ip),allocatable :: &
         nnode(:),&                ! # of nodes per element     
         ngaus(:),&                ! # of Gauss points per element  
         nnodb(:),&                ! # of nodes per boundary element 
         nface(:),&                ! # of faces per element 
         ngaub(:)                  ! # of Gauss points per boundary element     

      ! Lists
      integer(ip), allocatable :: &
         ltype(:),&                ! List of element types
         lquad(:),&                ! List of quadrature
         ltopo(:),&                ! List of topology
         llapl(:),&                ! List of Laplacian
         lrule(:),&                ! List of integration rules
         ltypb(:),&                ! List of boundary types
         lpoty(:)                  ! List of point types
         
      !Shape functions and derivatives  
      real(rp), allocatable  :: &
            shape(:,:,:),&
            deriv(:,:,:,:),&
            heslo(:,:,:,:),&
            hnatu(:),&                    ! Natural element length
            weigp(:,:),&
            posgp(:,:,:)                  !position of the gauss point coordinates
      
      !>Shape functions and derivatives for a closed integration rule (projections, vmass)
      real(rp), allocatable :: &
            shape_clo(:,:,:), &
            deriv_clo(:,:,:,:), &
            heslo_clo(:,:,:,:), &
            weigp_clo(:,:)
            
      ! Shape functions and derivatives center of gravity  
      real(rp), allocatable  :: &       
            shacg(:,:),&
            dercg(:,:,:),&
            weicg(:)
      ! Boundary shape functions and derivatives  
      real(rp), allocatable  :: &
            shapb(:,:,:),&
            derib(:,:,:,:),&
            weigb(:,:)
      real(rp), allocatable  :: &
            shaga(:,:,:),&
            shagb(:,:,:)   
      
      !Temporary arrays for reading boundary conditions
      integer(ip), allocatable :: &
         pbopo(:), lbopo(:)             ! pbopo(npoin), lbopo() list of iboun for ipoin (used for reading boundary conditions)      
      
      !Arrays with mesh info 
      integer(ip), allocatable :: & 
         cfael(:,:,:),&                 ! List of nodes in element boundaries
         lnods(:),pnods(:),&            ! Interior element connectivity   
         lboel(:),pboel(:),&            ! List of boundary elements
         pelpo(:),lelpo(:),&            ! List of elements to which a node belongs
         ia(:),ja(:),&                  ! Local graph
         lbody(:)                       ! List of Bodies on boundaries         

      real(rp), allocatable  :: &
         coord(:,:),&                  ! Coordinates
         exnor(:,:,:),&                ! Exterior normal  
         geoblock(:),&                 ! List of Blocks for output
         weight(:)                     ! # processors per boundary
         
      !For Fixing External normal
      logical, allocatable :: IsExnor(:)
      type(r1p), allocatable :: ExternalNormal(:)
         
      real(rp) :: grange(2,3), range(2,3)
      
      integer(ip),allocatable :: &
         fictio(:),                &   ! Nodal numeration for periodic bc's
         fic_nodes(:,:)                ! List of fictitious nodes
      logical, allocatable :: &
         mslave(:)                     ! Master-Slave nodes
         
      real(rp), allocatable  :: &
         skcos(:,:,:)                ! Cosine matrices of skew systems   
      
      !For projections
      real(rp),allocatable  ::   &
         vmass(:)                      ! Lumped Mass matrix   
      !Linear System
      class(ParallelSystemInterface), pointer   :: L2ProjectorSystem => NULL()
      real(rp), allocatable                     :: L2ProjectorUnkno(:)
      integer(ip)                               :: kfl_ProjectorReady = 0
      integer(ip)                               :: kfl_ProjectionFull = 0
      
      
      !External ParallelLibrary
      class(ParallelLibraryInterface),          pointer :: ParallelLibrary => NULL()
      class(ParallelOrderingInterface),         pointer :: LocalOrdering => NULL()
      class(ParallelCommunicatorInterface),     pointer :: ArrayCommunicator => NULL()
      !Do not move RP_LocalOrdering from here, otherwise intel crashes
      class(ParallelOrderingInterface),         pointer :: RP_LocalOrdering => NULL(), Map_LocalOrdering => NULL()
      class(ParallelOrderingInterface),         pointer :: RP_InitialOrdering => NULL(), Map_InitialOrdering => NULL()
      
      !For Periodic Boundary Conditions
      integer(ip) :: kfl_perio       !Mesh has periodic boundary conditions
      integer(ip) :: kfl_ForcePeriodicMesh = 0
      real(rp)    :: ForcePeriodicMeshTol = 1e-2
      logical     :: PerDime(3)      !do we have periodicity in x, y, z
      integer(ip) :: nslave = 0
      integer(ip), allocatable :: MasterSlave(:)
      integer(ip), allocatable :: SlaveList(:)
      
      integer(ip), allocatable :: gMasterSlave(:) !Temporary, for reading boundary conditions
      
      !For ALE Mesh
      integer(ip) :: kfl_alemov = 0
      real(rp), pointer :: &
         displ(:,:,:) => NULL(),&
         veloc(:,:,:) => NULL()
      logical :: AreFoldedALE = .false.
      
      
      !For Manufactured Mesh
      character(6) ::  ManufacturedType
      integer(ip)  ::  ManufacturedNpoinX
      real(rp)     ::  ManufacturedInternalRadius,ManufacturedExternalRadius 
      real(rp)     ::  ManufacturedAxialLength 
      integer(ip)  ::  ManufacturedDivisionsInAngle,ManufacturedNpoinZ
         
      !For Meshes with Hanging Nodes
      logical :: kfl_HangingNodes = .false.
      integer(ip), allocatable :: pHangingList(:),lHangingList(:)
      real(rp), allocatable    :: rHangingList(:)
      integer(ip), allocatable :: seenHanging(:), iwaHanging(:)
      
      !For Adaptive Refinement
      class(AdaptiveRefinerInterface), pointer :: Refiner => NULL()
      logical :: isSetAdaptiveRefiner = .false.
      integer(ip), allocatable :: auxBoundaryMatch1(:)
      integer(ip), allocatable :: auxBoundaryMatch2(:,:)
      integer(ip), allocatable :: BoundaryNewToOld(:)
      integer(ip), allocatable :: iauxBoundaryMatch1(:)
      type(i1p), allocatable   :: iauxBoundaryMatch2(:)
      integer(ip)              :: auxOldNboun

      !For reading elemental data in an efficient manner
      !integer(ip), allocatable :: gLocalNaiveToInitialElementNumbering(:)
      !type(i1p), allocatable   :: gInitialNumberingElementsINeed(:),gInitialNumberingElementsOthersNeed(:)
      
      class(ParallelOrderingInterface), pointer :: ElementLocal2Initial => NULL()
      class(ParallelOrderingInterface), pointer :: ElementNaive2Initial => NULL()
      
      type(i1p), pointer :: gElementInitialSendToProcessor(:) => NULL()
      type(i1p), allocatable :: ElementNaiveSendToProcessor(:)

      character(10),allocatable :: blockName(:)
      
contains  
      
      procedure :: SetParallelLibrary   
      !procedure :: OpenFilesAndReadGeneralInfo => MeshOpenFilesAndReadGeneralInfo
      procedure :: MeshOpenFilesAndReadGeneralInfo
      procedure :: CloseDataFiles => MeshCloseDataFiles
      procedure :: Readat
      procedure :: Initialize
      
      procedure :: countv
      procedure :: ExtnorLpoty
      procedure :: ComputeVmass
      procedure :: ReadatMpi
      procedure :: cshder
      procedure :: cderda
      procedure :: reastr
      procedure :: dombou
      procedure :: blbopo
      procedure :: TestCommunications
      
      procedure :: LocalDimensionsAndIRenumbering
      procedure :: GhostPointsAndLocalOrdering
      procedure :: GhostCoord
      procedure :: LnodsToLocal
      procedure :: LocalGraph
      procedure :: BuildArrayCommunicator
      
      procedure :: Local2InitialScalar
      procedure :: Local2InitialVector      
      generic   :: Local2Initial => Local2InitialScalar, Local2InitialVector
      
      procedure :: Initial2LocalScalar
      procedure :: Initial2LocalVector      
      generic   :: Initial2Local => Initial2LocalScalar, Initial2LocalVector
      
      procedure :: Global2InitialScalar
      procedure :: Global2InitialVector
      generic   :: Global2Initial => Global2InitialScalar, Global2InitialVector
      
      procedure :: Initial2GlobalScalar
      procedure :: Initial2GlobalVector
      generic   :: Initial2Global => Initial2GlobalScalar, Initial2GlobalVector
      
      procedure :: Initial2ProcNumberScalar
      procedure :: Initial2ProcNumberVector      
      generic   :: Initial2ProcNumber => Initial2ProcNumberScalar, Initial2ProcNumberVector
      procedure :: ExportigPointGlobNumber
      
      procedure :: GetNumBlocks
      procedure :: GetNdime
      procedure :: GetGnpoin
      procedure :: GetNpoin
      procedure :: GetNbopo
      procedure :: GetNboun
      procedure :: GetNbody
      procedure :: GetNpoinLocal
      procedure :: GetNpoinGhost
      procedure :: GetMnodb
      procedure :: GetMnode
      procedure :: GetNelem
      procedure :: GetGnelem
      procedure :: GetNelty
      procedure :: GetIelty
      procedure :: GetElementTypeSizes 
      procedure :: Finbou
      procedure :: FindElement
      procedure :: GetTimes => MeshGetTimes
      
      procedure :: GetBoundaryIelem
      procedure :: GetLnodb
      procedure :: GetLnode
      procedure :: GetElemArraySize
      procedure :: GetBounArraySize
      procedure :: GetPointBlock
      procedure :: GetPointCoord
      procedure :: GetVecCoord
      procedure :: GetIbopo
      procedure :: GetBoundaryPoint
      procedure :: GetLelpo
      procedure :: GetPosgp
      procedure :: GetVmass
      procedure :: GetCoord
      procedure :: GetNeighboringPoints
      
      procedure :: ElementSetSizes
      procedure :: ElementLoad
      procedure :: FaceLoad
      procedure :: ElementSetPointers
      procedure :: BoundaryLoad
      procedure :: BoundarySetPointers
      procedure :: ElementAlloc
      procedure :: ElementDealloc

      procedure :: InitializeSystem
      !procedure :: DeallocateSystem
      
      procedure :: DeallocateReadGlobals => msh_DeallocateReadGlobals
      procedure :: DeallocateMeshInfo => MeshDeallocateMeshInfo
      procedure :: DeallocateLocalArrays => MeshDeallocateLocalArrays

      procedure :: Smooth
      procedure :: Project
      procedure :: ComputeMassMatrix
      
      procedure :: Local2GlobalScalar
      procedure :: Local2GlobalVector
      procedure :: Global2LocalScalar
      procedure :: Global2LocalVector
      procedure :: ElementLocal2Global
      
      generic   :: Local2Global => Local2GlobalScalar, Local2GlobalVector
      generic   :: Global2Local => Global2LocalScalar, Global2LocalVector
      
      procedure :: NonLinearElements
      procedure :: IsNonLinearElement
      
      procedure :: Turnof => MeshTurnof
      procedure :: GetLbody
      procedure :: GetWeight
      procedure :: GetIsoRange

      !Body definition
      procedure :: reabody
      
      !Periodic Boundary Conditions
      procedure :: GetPerio
      procedure :: SetPerio
      procedure :: GetArePeriodicDimensions
      procedure :: AssemblyPeriodicBC
      procedure :: AssemblyPeriodicBCToZero
      procedure :: MasterToSlaveInteger
      procedure :: MasterToSlaveReal
      procedure :: MasterToSlaveReal1
      procedure :: MasterToSlaveInteger1
      generic   :: MasterToSlave => MasterToSlaveInteger, MasterToSlaveReal,MasterToSlaveReal1,MasterToSlaveInteger1
      procedure :: Initial2MasterInitial
      procedure :: GetMaster

      procedure :: GetRange => MeshGetRange
      procedure :: GetGRange => MeshGetGRange
      
      !ALE
      procedure :: GetALE
      procedure :: SetALE
      procedure :: SetDisplacements
      procedure :: SetVelocities
      procedure :: DeallocExnorLpoty
      procedure :: DeallocVmass
      procedure :: GetMeshVeloc
      procedure :: GetDeformedCoord
      procedure :: GetDeformedRange
      procedure :: ComputeCheckFoldedALE
      procedure :: GetCheckFoldedALE

      !For initializing from arrays
      procedure :: BuildLnodsFromArray
      procedure :: SetNpoinAndBuildOrderingFromArray
      procedure :: BuildCoordinatesFromArray
      procedure :: BuildBlocksFromArray
      procedure :: BuildExnorFromArray
      procedure :: BuildListOfBodiesFromArray
      
      !For assembly to arrays
      procedure :: AssemblyToArray
      
      !For Hanging nodes
      procedure :: GetHanging
      procedure :: SetHanging
      procedure :: BuildHangingNodesFromArray
      procedure :: BuildHangingGraph
      procedure :: IsHangingNode
      procedure :: IsHangingFace
      
      !For array assembly
      procedure :: AssemblyToArrayHanging
      
      !For system assembly
      procedure :: PrepareHangingMatrices
      procedure :: AssemblyHangingNodesDiag
      procedure :: AssemblyHangingNodesDiagToZero
      procedure :: InterpolateHangingValues
      
      !For Adaptive Refinement
      procedure :: SetAdaptiveRefiner
      procedure :: IsAnAdaptiveMesh
      procedure :: Refine => MeshRefine
      procedure :: InitializeAdaptiveRefinement => MeshInitializeAdaptiveRefinement
      procedure :: GetBoundaryNewToOld
      procedure :: MatchOldBoundariesToRefiner => MeshMatchOldBoundariesToRefiner
      procedure :: GetRebalanceInverseBoundaryMatches => MeshGetRebalanceInverseBoundaryMatches
      procedure :: ComputeRebalancingNumbering => MeshComputeRebalancingNumbering
      
      !Manufactured Mesh
      procedure :: ManufacturedMesh
      procedure :: ManufacturedMeshb
      procedure :: IsManufactured
      
      !For restart
      procedure :: InitialOrdering
      procedure :: DeallocateInitial
      
      

   end type  
   
   interface Global2Local
     module procedure Global2LocalScalar, Global2LocalVector  
   end interface
   
   interface Local2Global
     module procedure Local2GlobalScalar, Local2GlobalVector  
   end interface
   
   interface 
      
      subroutine MeshOpenFilesAndReadGeneralInfo(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine MeshCloseDataFiles(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
   
      subroutine ExtnorLpoty(a)
         import FemMesh
         implicit none
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine ComputeVmass(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine ReadatMPI(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine cshder(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine cderda(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
   
      subroutine countv(a)
         import FemMesh
         class(FemMesh) :: a
      end subroutine
      
      subroutine reastr(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
            
      subroutine reaLnods(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine reacoord(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine reabody(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine      
      
      subroutine Partitionate(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine dombou(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine blbopo(a)
         use typre
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine TestCommunications(a)
         use typre
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
            
      subroutine LocalDimensionsAndIRenumbering(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine

      subroutine GhostPointsAndLocalOrdering(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine BuildArrayCommunicator(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine GhostCoord(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine LnodsToLocal(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine LocalGraph(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine

      subroutine Finbou(a,pnodb,knodb,iboun)
         use typre
         import FemMesh
         implicit none
         
         class(FemMesh) :: a
         integer(ip) :: pnodb,knodb(pnodb)
         integer(ip) :: iboun

       end subroutine
       
       subroutine FindElement(a,pnode,lnode,ielem)
         use typre
         import FemMesh
         implicit none
         
         class(FemMesh) :: a
         integer(ip) :: pnode,lnode(pnode)
         integer(ip) :: ielem

       end subroutine
      
      subroutine ElementSetSizes(a,Element)
         import FemMesh
         import MemoryMan
         import FiniteElement
         implicit none
         class(FemMesh) :: a
         class(FiniteElement) :: Element
      
      
      end subroutine
      
      subroutine msh_DeallocateReadGlobals(a)
         import
         implicit none
         class(FemMesh) :: a
         
      end subroutine

      subroutine ElementSetPointers(a,Element)
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh) :: a
         class(FiniteElement) :: Element
      
      
      end subroutine
      
      subroutine ElementLoad(a,ielem,Element)
         use typre
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh)         :: a
         integer(ip)            :: ielem
         class(FiniteElement)   :: Element
      
      end subroutine
      
      subroutine FaceLoad(a,iface,ielem,Element)
         use typre
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh)         :: a
         integer(ip)            :: iface,ielem
         class(FiniteElement)   :: Element
      
      end subroutine
      
      subroutine ElementAlloc(a,Element,Memor,rule,outstr)
         use Mod_Memor
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh) :: a
         class(FiniteElement), pointer :: Element
         type(MemoryMan)     :: Memor
         character(*)        :: outstr
         character(6)        :: rule
      end subroutine
      
      subroutine ElementDealloc(a,Element,Memor,rule,outstr)
         use Mod_Memor
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh) :: a
         class(FiniteElement), pointer :: Element
         type(MemoryMan)     :: Memor
         character(*)        :: outstr
         character(6)        :: rule
      end subroutine
      
      subroutine BoundarySetPointers(a,Element)
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh) :: a
         class(FiniteElement) :: Element
      
      
      end subroutine
      
      subroutine BoundaryLoad(a,iboun,Element)
         use typre
         import FemMesh
         import FiniteElement
         implicit none
         class(FemMesh)         :: a
         integer(ip)            :: iboun
         class(FiniteElement)    :: Element
      
      end subroutine
      
      subroutine OutputMeshInfo(a)
         import FemMesh
         implicit none
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine AssemblyPeriodicBC(a,ndofn,LinearSystem,Memor)
         use typre
         use Mod_Memor
         use Mod_ParallelSystemInterface
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         class(ParallelSystemInterface) :: LinearSystem
         type(MemoryMan) :: Memor
      end subroutine   
      
      subroutine AssemblyPeriodicBCToZero(a,ndofn,LinearSystem,Memor)
         use typre
         use Mod_Memor
         use Mod_ParallelSystemInterface
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         class(ParallelSystemInterface) :: LinearSystem
         type(MemoryMan) :: Memor
      end subroutine   
      
      subroutine PeriodicBCModify(a)
         import FemMesh
         implicit none
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine Project(a,ndofn,array)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         real(rp)    :: array(ndofn,a%npoin)
      
      end subroutine
      
      subroutine ComputeMassMatrix(a)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         
      end subroutine 
 
      subroutine ComputeRanges(a)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         
      end subroutine  
      
      subroutine reaExnor(a)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         
      end subroutine  

      subroutine AssemblyToArray(a,e,ndofn,elrhs,array)
         use typre
         use Mod_Element
         import FemMesh
         implicit none
         class(FemMesh) :: a
         class(FiniteElement) :: e
         integer(ip) :: ndofn
         real(rp) :: elrhs(ndofn,*),array(ndofn,*)
      end subroutine
      
      subroutine AssemblyToArrayHanging(a,e,ndofn,elrhs,array)
         use typre
         use Mod_Element
         import FemMesh
         implicit none
         class(FemMesh) :: a
         class(FiniteElement) :: e
         integer(ip) :: ndofn
         real(rp) :: elrhs(ndofn,*),array(ndofn,*)
      end subroutine
      
      subroutine BuildHangingGraph(a)
         import FemMesh
         class(FemMesh), intent(inout) :: a
      end subroutine
      
      subroutine InterpolateHangingValues(a,ndofn,array)
         use typre
         use Mod_Element
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         real(rp) :: array(ndofn,*)
      end subroutine
      
      subroutine MeshRefine(a,itask)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         character(6) :: itask
      
      end subroutine
      
      subroutine PrepareHangingMatrices(a,e,elmat,elrhs)
         use typre
         use Mod_Element
         import FemMesh
         implicit none
         class(FemMesh), target :: a
         class(FiniteElement) :: e
         real(rp), allocatable ::  elmat(:,:,:,:),elrhs(:,:)
         integer(ip) :: HangCount
      end subroutine
      
      subroutine AssemblyHangingNodesDiag(a,ndofn,LinearSystem,Memor)
         use typre
         use Mod_Memor
         use Mod_ParallelSystemInterface
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         class(ParallelSystemInterface) :: LinearSystem
         type(MemoryMan) :: Memor
      end subroutine
      
      subroutine AssemblyHangingNodesDiagToZero(a,ndofn,LinearSystem,Memor)
         use typre
         use Mod_Memor
         use Mod_ParallelSystemInterface
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         class(ParallelSystemInterface) :: LinearSystem
         type(MemoryMan) :: Memor
      end subroutine
      
      subroutine MeshInitializeAdaptiveRefinement(a,Refiner)
         use typre
         use Mod_AdaptiveInterface
         import FemMesh
         implicit none
         class(FemMesh) :: a
         class(AdaptiveRefinerInterface) :: Refiner
      end subroutine
      
      subroutine MeshMatchOldBoundariesToRefiner(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine

      subroutine ManufacturedMesh(a)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         
      end subroutine  
      
      subroutine ManufacturedMeshb(a)
         import FemMesh
         implicit none
         class(FemMesh) :: a
         
      end subroutine  
      
      subroutine Smooth(a,ndofn,array)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         real(rp)    :: array(ndofn,*)
         
      end subroutine  
      
      subroutine MeshComputeRebalancingNumbering(a,PointGlobNumber,PointProcNumber)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: PointGlobNumber(:), PointProcNumber(:)
      end subroutine
      
      subroutine MeshTurnof(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine MeshDeallocateMeshInfo(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine MeshDeallocateLocalArrays(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine ComputeCheckFoldedALE(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine GetCheckFoldedALE(a,AreFoldedALE)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
         logical :: AreFoldedALE
      end subroutine
      
      subroutine InitialOrdering(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine DeallocateInitial(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      
      subroutine msh_InitializeElementDataStructures(a)
         use typre
         import FemMesh
         implicit none
         class(FemMesh) :: a
      end subroutine
   
   
   
   end interface   
   
   
   
contains   
  
   subroutine SetParallelLibrary(a,ParallelLibrary)
      implicit none
      class(FemMesh) :: a
      class(ParallelLibraryInterface), target :: ParallelLibrary
      
      a%ParallelLibrary => ParallelLibrary
   end subroutine
   
   subroutine Readat(a)
      implicit none
      class(FemMesh) :: a
      
      real(rp) :: time2,time1
      
      call a%Timer%Total%Tic
      
      !Open files and read info
      !call a%OpenFilesAndReadGeneralInfo

      call a%MeshOpenFilesAndReadGeneralInfo
      
      call a%cderda    !Derivated parameters
      
      call a%cshder    !Compute shape and derivatives

      !Default is read from disk
      if (a%kfl_MeshInputType == 0) then
         !Read lnods and send parts to processes
         call reaLnods(a)
         
         !Partitionate
         call Partitionate(a)
         
         !Renumbering 
         call a%LocalDimensionsAndIRenumbering
         
         !Ghost points and Local Ordering
         call a%GhostPointsAndLocalOrdering
         
         !Read the Coordinates
         call reaCoord(a)
         
         !Read ExternalNormal
         call reaExnor(a)
         
         !Lnods to Local numbering
         call a%LnodsToLocal
         
      !Manufactured Mesh   
      elseif (a%kfl_MeshInputType == 1) then
         call a%ManufacturedMesh
      endif
      
      !Local Graph
      call a%LocalGraph
      
      !Build Array Communicator
      call a%BuildArrayCommunicator
      
      !Build local arrays from Scattered Data
      call a%GhostCoord
      
      !Compute the coordinates Ranges
      call ComputeRanges(a)
      
      call a%Timer%Total%Toc
   end subroutine
  
   subroutine Initialize(a)
      implicit none
      
      class(FemMesh) :: a
      
      logical :: kfl_UseElementDataStructures
      
      call a%Timer%Total%Tic
      
      if (a%kfl_perio == 1) then
         !If element data structures is set, we can not use it here yet
         !We deactivate it and reactivate it later
         kfl_UseElementDataStructures = a%kfl_UseElementDataStructures
         a%kfl_UseElementDataStructures = .false.
         
         call PeriodicBCModify(a)
         
         !Build local arrays from Scattered Data
         call a%GhostCoord
         
         !Reactivate element data structures if it was set  
         a%kfl_UseElementDataStructures = kfl_UseElementDataStructures
         
      endif
      
      !Compute LBOEL 
      call a%Dombou     
      
      if (a%kfl_UseElementDataStructures) then
         call msh_InitializeElementDataStructures(a)
      endif
      
      !Compute Vmass
      call a%ComputeVmass
      
      !Compute EXNOR, LPOTY
      call a%ExtnorLpoty
      
      !Compute List of boundaries for each point
      call a%Blbopo
      
      if (a%kfl_TestCommunications == 1) &
         call a%TestCommunications
      
      !For Manufactured Mesh we do not need it
      if (a%kfl_MeshInputType == 0) then
         !Read boundary bodies
         call reabody(a)      
      elseif (a%kfl_MeshInputType == 1) then
         !Manufactured mesh 
         call ManufacturedMeshb(a)
      endif
         
      call a%CloseDataFiles
      
      call a%Timer%Total%Toc
      
      call OutputMeshInfo(a)
   
   end subroutine
         
   subroutine Initial2LocalScalar(a,InitialArray,LocalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray, LocalArray
      integer(ip) :: ipoin

      if (InitialArray > a%RP_poin0 .and. InitialArray <= a%RP_poin0 + a%RP_npoinLocal) then
         LocalArray = InitialArray - a%RP_poin0
      else
         call a%RP_InitialOrdering%Global2Local(InitialArray,LocalArray)
         if (LocalArray == 0) call runend('Initial2LocalScalar: point does not belong to initial')
         LocalArray = LocalArray+a%RP_npoinLocal
      endif

   end subroutine      
         
   subroutine Initial2LocalVector(a,npoin,InitialArray,LocalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray(npoin), LocalArray(npoin)
      integer(ip) :: ipoin
      
      do ipoin = 1,npoin
         if (InitialArray(ipoin) > a%RP_poin0 .and. InitialArray(ipoin) <= a%RP_poin0 + a%RP_npoinLocal) then
            LocalArray(ipoin) = InitialArray(ipoin) - a%RP_poin0
         else
            call a%RP_InitialOrdering%Global2Local(InitialArray(ipoin),LocalArray(ipoin))
            if (LocalArray(ipoin) == 0) call runend('Initial2LocalScalar: point does not belong to initial')
            LocalArray(ipoin) = LocalArray(ipoin)+a%RP_npoinLocal
         endif
      enddo
      
   end subroutine
         
   subroutine Local2InitialScalar(a,InitialArray,LocalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray, LocalArray
      integer(ip) :: ipoin

      if (LocalArray <= a%RP_npoinLocal) then
         InitialArray = LocalArray + a%RP_poin0 
      else
         LocalArray = LocalArray-a%RP_npoinLocal
         call a%RP_InitialOrdering%Local2Global(LocalArray,InitialArray)
      endif

   end subroutine      
         
   subroutine Local2InitialVector(a,npoin,InitialArray,LocalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray(npoin), LocalArray(npoin)
      integer(ip) :: ipoin
      
      do ipoin = 1,npoin
         if (LocalArray(ipoin) <= a%RP_npoinLocal) then
            InitialArray(ipoin) = LocalArray(ipoin) + a%RP_poin0 
         else
            LocalArray = LocalArray-a%RP_npoinLocal
            call a%RP_InitialOrdering%Local2Global(LocalArray(ipoin),InitialArray(ipoin))
         endif
      enddo
      
   end subroutine
   
   subroutine Initial2GlobalVector(a,npoin,InitialArray,GlobalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray(npoin), GlobalArray(npoin)
      integer(ip) :: ilocal(npoin)
   
      if (a%kfl_ReadType == 0) then
         if (allocated(a%gPointGlobNumber)) then
            GlobalArray = a%gPointGlobNumber(InitialArray)
         else
        		call a%Initial2LocalVector(npoin,InitialArray,InitialArray)
	         GlobalArray = a%lPointGlobNumber(InitialArray)
         end if
      elseif (a%kfl_ReadType == 1) then
        !if local point
        call a%Initial2LocalVector(npoin,InitialArray,GlobalArray)
        call a%RP_LocalOrdering%Local2Global(npoin,GlobalArray,GlobalArray)
      endif

   end subroutine
   
   subroutine Initial2GlobalScalar(a,InitialArray,GlobalArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: InitialArray, GlobalArray
      integer(ip) :: ilocal
   
      if (a%kfl_ReadType == 0) then
         if (allocated(a%gPointGlobNumber)) then
  	         GlobalArray = a%gPointGlobNumber(InitialArray)
         else
        		call a%Initial2LocalScalar(InitialArray,InitialArray)
            GlobalArray = a%lPointGlobNumber(InitialArray)
         end if
      elseif (a%kfl_ReadType == 1) then
        call a%Initial2LocalScalar(InitialArray,GlobalArray)
        call a%RP_LocalOrdering%Local2Global(GlobalArray,GlobalArray)
      endif
   end subroutine
   
   subroutine Global2InitialVector(a,npoin,GlobalArray,InitialArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: InitialArray(npoin), GlobalArray(npoin)
      integer(ip) :: ilocal(npoin)
   
      if (a%kfl_ReadType == 0) then
         if (allocated(a%igPointGlobNumber)) then
            InitialArray = a%igPointGlobNumber(GlobalArray)
         else
            call a%Global2LocalVector(npoin,GlobalArray,GlobalArray)
            InitialArray = a%ilPointGlobNumber(GlobalArray)
         end if
      elseif (a%kfl_ReadType == 1) then
        !if local point
        call a%RP_LocalOrdering%Global2Local(npoin,GlobalArray,GlobalArray)
        call a%Local2InitialVector(npoin,GlobalArray,InitialArray)
      endif

   end subroutine
   
   subroutine Global2InitialScalar(a,GlobalArray,InitialArray)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: InitialArray, GlobalArray
      integer(ip) :: ilocal
   
      if (a%kfl_ReadType == 0) then
         if (allocated(a%igPointGlobNumber)) then
            InitialArray = a%igPointGlobNumber(GlobalArray)
         else
            call a%Global2LocalScalar(GlobalArray,GlobalArray)
            InitialArray = a%ilPointGlobNumber(GlobalArray)
         end if
      elseif (a%kfl_ReadType == 1) then
        call a%RP_LocalOrdering%Global2Local(GlobalArray,GlobalArray)
        call a%Local2InitialScalar(GlobalArray,InitialArray)
      endif
   end subroutine
  
   subroutine Initial2ProcNumberVector(a,npoin,GlobalNumbering,ProcNumber)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: GlobalNumbering(npoin), ProcNumber(npoin)
      integer(ip) :: ilocal(npoin)
      
      if (a%kfl_ReadType == 0) then
         if (allocated(a%gPointProcNumber)) then
            ProcNumber = a%gPointProcNumber(GlobalNumbering)
         else
            call a%Initial2LocalVector(npoin,GlobalNumbering,GlobalNumbering)
            ProcNumber = a%lPointProcNumber(GlobalNumbering)
         end if
         
      elseif (a%kfl_ReadType == 1) then
        !if local point
        call a%Initial2LocalVector(npoin,GlobalNumbering,ProcNumber)
        call a%RP_LocalOrdering%GetProc(npoin,ProcNumber,ProcNumber)
      endif
   end subroutine
   
   subroutine Initial2ProcNumberScalar(a,GlobalNumbering,ProcNumber)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: GlobalNumbering, ProcNumber
      
      if (a%kfl_ReadType == 0) then
         if (allocated(a%gPointProcNumber)) then
            ProcNumber = a%gPointProcNumber(GlobalNumbering)
         else
        		call a%Initial2LocalScalar(GlobalNumbering,GlobalNumbering)
            ProcNumber = a%lPointProcNumber(GlobalNumbering)
         end if
      elseif (a%kfl_ReadType == 1) then
         call a%Initial2LocalScalar(GlobalNumbering,ProcNumber)
         call a%RP_LocalOrdering%GetProc(ProcNumber,ProcNumber)
      endif
   end subroutine

   subroutine GetNumBlocks(a,numB)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: numB
      
      numB = a%numGeoBlocks
   end subroutine
   
   subroutine GetNdime(a,ndime)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndime
      
      ndime = a%ndime
   end subroutine
   
   subroutine GetGnpoin(a,gnpoin)
      class(FemMesh), target :: a
      integer(ip) :: gnpoin
      
      gnpoin = a%gnpoin
   end subroutine
   
   subroutine GetMnodb(a,mnodb)
      class(FemMesh), target :: a
      integer(ip) :: mnodb
      
      mnodb = a%mnodb
   end subroutine
   
   subroutine GetMnode(a,mnode)
      class(FemMesh), target :: a
      integer(ip) :: mnode
      
      mnode = a%mnode
   end subroutine
   
   subroutine GetNpoin(a,npoin)
      class(FemMesh), target :: a
      integer(ip) :: npoin
      
      npoin = a%npoin
   end subroutine
   
   subroutine GetNbopo(a,nbopo)
      class(FemMesh), target :: a
      integer(ip) :: nbopo
      
      nbopo = a%nbopo
   end subroutine
     
   subroutine GetNpoinLocal(a,npoinLocal)
      class(FemMesh), target :: a
      integer(ip) :: npoinLocal
      
      npoinLocal = a%npoinLocal
   end subroutine
   
   subroutine GetNpoinGhost(a,npoinGhost)
      class(FemMesh), target :: a
      integer(ip) :: npoinGhost
      
      npoinGhost = a%npoinGhost
   end subroutine
   
   subroutine GetNboun(a,nboun)
      class(FemMesh), target :: a
      integer(ip) :: nboun
      
      nboun = a%nboun
   end subroutine

   subroutine GetNbody(a,nbody)
      class(FemMesh), target :: a
      integer(ip) :: nbody
      
      nbody = a%nbody
   end subroutine
   
   subroutine GetNelem(a,nelem)
      class(FemMesh) :: a
      integer(ip) :: nelem
      
      nelem = a%nelem
   end subroutine
   
   subroutine GetGnelem(a,gnelem)
      class(FemMesh) :: a
      integer(ip) :: gnelem
      
      gnelem = a%gnelem
   end subroutine
   
   subroutine GetNelty(a,nelty)
      class(FemMesh) :: a
      integer(ip) :: nelty   
      
      nelty = a%nelty   
   end subroutine
   
   subroutine GetIelty(a,ielem,ielty)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ielem,ielty
      
      ielty = 1
      if (a%nelty > 1) then 
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
      endif
   end subroutine
   
   subroutine GetElementTypeSizes(a,nelty,nnode,ngaus,ltopo)
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: nelty
      integer(ip), pointer :: nnode(:),ngaus(:),ltopo(:)      

      nelty = a%nelty
      nnode => a%nnode
      ngaus => a%ngaus
      ltopo => a%ltopo
   end subroutine
   
   subroutine GetBoundaryIelem(a,iboun,ielem)
      implicit none
      class(FemMesh), target  :: a
      integer(ip),intent(in)  :: iboun
      integer(ip),intent(out) :: ielem
      integer(ip) :: pblty,ispos2

      pblty = 1
      ispos2 = (a%nnodb(1)+1)*(iboun-1)
      if (a%nelty > 1) then 
         ispos2 = a%pboel(iboun)-1
         pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      endif
      ielem = a%lboel(ispos2+a%nnodb(pblty)+1)
   end subroutine
   
   subroutine GetLnodb(a,iboun,pnodb,lnodb)
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: iboun,pnodb
      integer(ip) :: lnodb(*)
      integer(ip) :: ielem,pblty,ispos,ispos2,inodb,inode
      
      pblty = 1
      ispos2 = (a%nnodb(1)+1)*(iboun-1)
      if (a%nelty > 1) then 
         ispos2 = a%pboel(iboun)-1
         pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      endif
      ielem = a%lboel(ispos2+a%nnodb(pblty)+1)
      ispos = a%nnode(1)*(ielem-1)
      if (a%nelty > 1) ispos = a%pnods(ielem)-1
      pnodb=a%nnodb(pblty)
      
      do inodb = 1,pnodb
         inode = a%lboel(ispos2+inodb)
         lnodb(inodb) = a%lnods(ispos+inode) 
      enddo 
   end subroutine
   
   subroutine GetLnode(a,ielem,pnode,lnode)
      class(FemMesh), target :: a
      integer(ip) :: ielem,pnode
      integer(ip), pointer :: lnode(:)
      integer(ip) :: ispos
   
      pnode = a%nnode(1)
      ispos = a%nnode(1)*(ielem-1)            
      !------------------------------------------------------------------------------
      if (a%nelty > 1) then  
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         ispos = a%pnods(ielem)-1
      endif
      !------------------------------------------------------------------------------
   
      lnode => a%lnods(ispos+1:ispos+pnode)
   end subroutine
   
   subroutine GetIbopo(a,ipoin,ibopo)
      class(FemMesh), target  :: a
      integer(ip), intent(in) :: ipoin
      integer(ip), intent(out):: ibopo
      
      ibopo = a%Lpoty(ipoin)
   end subroutine
   
   subroutine GetBoundaryPoint(a,ipoin,ibopo,exnor)
      class(FemMesh), target  :: a
      integer(ip), intent(in) :: ipoin
      integer(ip), intent(out):: ibopo
      real(rp), pointer       :: exnor(:,:)
      
      ibopo = a%Lpoty(ipoin)
      if (ibopo > 0) exnor => a%exnor(:,:,ibopo)
   end subroutine
   
   subroutine GetLelpo(a,ipoin,pelpo,lelpo)
      use typre
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: ipoin
      integer(ip), pointer :: lelpo(:)
      integer(ip) :: pelpo
      
      lelpo => a%lelpo(a%pelpo(ipoin):a%pelpo(ipoin+1)-1)
      pelpo = a%pelpo(ipoin+1)-a%pelpo(ipoin)
   end subroutine    
   
   subroutine GetElemArraySize(a,ielem,pnode,pgaus)
      use typre
      class(FemMesh) :: a
      integer(ip) :: ielem,pnode,pgaus
      integer(ip) :: ielty
     
      ielty = 1 
      if (a%nelty > 1) then
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
      endif
      pnode = a%nnode(ielty)
      pgaus = a%ngaus(ielty)
   end subroutine

   subroutine GetBounArraySize(a,iboun,pnodb,pgaub)
      use typre
      class(FemMesh) :: a
      integer(ip) :: iboun,pnodb,pgaub
      integer(ip) :: pblty
      
      pblty = 1
      if (a%nelty > 1) then
         pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      endif
      pnodb=a%nnodb(pblty)
      pgaub = a%ngaub(pblty)
   end subroutine

   subroutine GetPointBlock(a,ipoin,iblock)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ipoin
      real(rp)    :: iblock
      
      iblock = a%geoblock(ipoin)
   end subroutine 
   
   subroutine GetPointCoord(a,ipoin,coord)
      use typre
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: ipoin
      real(rp), pointer :: coord(:)
      
      coord => a%coord(:,ipoin)
   end subroutine 

   subroutine GetVecCoord(a,vec,coord)
      use typre
      implicit none
      class(FemMesh), target :: a
      !real(rp) :: vec(:)
      integer(ip) :: vec(:)
      real(rp) :: coord(a%ndime,size(vec))
      
      coord = a%coord(:,vec)
   end subroutine 
   
   subroutine GetVmass(a,vmass)
      use typre
      implicit none
      class(FemMesh), target :: a
      real(rp), pointer :: vmass(:)
      
      vmass => a%vmass(:)
   end subroutine 
   
   subroutine InitializeSystem(a,ndofn,issym,LinearSystem,Memor,exmod,exfile,lun_solve,opti)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndofn,issym
      class(ParallelSystemInterface), pointer :: LinearSystem
      type(MemoryMan)      :: Memor
      character(150)       :: exmod
      character(150)       :: exfile
      integer(ip)          :: lun_solve
      character(5), optional :: opti
      
      call LinearSystem%Init(a%npoinLocal,a%npoinGhost,a%gnpoin,ndofn,issym,a%ia,a%ja,&
         &a%LocalOrdering,Memor,exmod,exfile,lun_solve,opti)
   end subroutine
   
   subroutine ExportigPointGlobNumber(a)
      use typre
      use def_parame
      use Mod_iofile
      implicit none
      class(FemMesh) :: a
      character(150) :: fil_pdata_dom
      integer(ip)    :: lun_pdata_dom
      
      if (a%MPIrank == a%MPIroot) then
         fil_pdata_dom = 'ignum.m'
         call iofile(zero,lun_pdata_dom,fil_pdata_dom,'ExportigPointGlobNumber')
         write(lun_pdata_dom,*) 'ignumbering = [', a%igPointGlobNumber, '];'
         
         call iofile(two,lun_pdata_dom,fil_pdata_dom,'ExportigPointGlobNumber')
      endif
   end subroutine
   
   subroutine Local2GlobalVector(a,npoin,LocalArray, GlobalArray)
      use typre
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
      
      call a%LocalOrdering%Local2Global(npoin,LocalArray,GlobalArray)
   end subroutine
   
   subroutine Local2GlobalScalar(a,LocalArray, GlobalArray)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: LocalArray, GlobalArray
      
      call a%LocalOrdering%Local2Global(LocalArray,GlobalArray)
   end subroutine
   
   subroutine Global2LocalVector(a,npoin,GlobalArray, LocalArray)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
      
      call a%LocalOrdering%Global2Local(npoin,GlobalArray,LocalArray)
   end subroutine
   
   subroutine Global2LocalScalar(a,GlobalArray, LocalArray)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      call a%LocalOrdering%Global2Local(GlobalArray,LocalArray)
   end subroutine   
      
   subroutine ElementLocal2Global(a,ielem,elemi)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip)    :: ielem, elemi
      
      elemi = ielem + a%elem0
   end subroutine
   
   subroutine NonLinearElements(a,flag)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: flag
      
      integer(ip) :: ielty,linear_nnode
      
      if (a%ndime == 2) then
         linear_nnode = 3
      elseif (a%ndime == 3) then
         linear_nnode = 4
      endif
      
      flag = 0
      do ielty = 1,a%nelty
         if (a%nnode(ielty) /= linear_nnode) flag = 1
      enddo
   end subroutine
   
   subroutine IsNonLinearElement(a,ielem,flag)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ielem,flag
      
      integer(ip) :: ielty,linear_nnode,pnode
      
      if (a%ndime == 2) then
         linear_nnode = 3
      elseif (a%ndime == 3) then
         linear_nnode = 4
      endif
      
      flag = 0
      if (a%nelty == 1) then
         pnode = a%pnods(1)
      else
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
      endif
      if (pnode /= linear_nnode) flag = 1
   end subroutine
         
   subroutine GetPerio(a,kfl_perio)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_perio
      
      kfl_perio = a%kfl_perio
   end subroutine
   
   subroutine SetPerio(a,kfl_perio)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_perio
      
      a%kfl_perio = kfl_perio
   end subroutine
   
   subroutine Initial2MasterInitial(a,ipoin,imaster)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ipoin,imaster
      
      call a%Initial2Global(ipoin,ipoin)
      imaster = a%gMasterSlave(ipoin)
      imaster = a%igPointGlobNumber(imaster)
   end subroutine   
   
   subroutine MasterToSlaveReal1(a,ndofn,array)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndofn
      real(rp)    :: array(*)
      integer(ip) :: ipoin,islave,imaster
      
      call MasterToSlaveReal(a,ndofn,array)
   end subroutine
   
   subroutine MasterToSlaveInteger1(a,ndofn,array)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndofn
      integer(ip) :: array(*)
      integer(ip) :: ipoin,islave,imaster
      
      call MasterToSlaveInteger(a,ndofn,array)
   end subroutine
   
   subroutine MasterToSlaveReal(a,ndofn,array)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndofn
      real(rp)    :: array(ndofn,a%npoin)
      integer(ip) :: ipoin,islave,imaster
      
      do ipoin = 1,a%nslave
         islave = a%SlaveList(ipoin)
         imaster = a%MasterSlave(islave)
         
         array(:,islave) = array(:,imaster)
      enddo
   end subroutine

   subroutine MasterToSlaveInteger(a,ndofn,array)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ndofn
      integer     :: array(ndofn,a%npoin)
      integer(ip) :: ipoin,islave,imaster
      
      do ipoin = 1,a%nslave
         islave = a%SlaveList(ipoin)
         imaster = a%MasterSlave(islave)
         array(:,islave) = array(:,imaster)
      enddo
   end subroutine
   
   subroutine GetMaster(a,islave,iMaster)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: islave,iMaster
      
      imaster = a%MasterSlave(islave)
   end subroutine

   subroutine SetALE(a,aleflag)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: aleflag
      
      a%kfl_alemov = aleflag
   end subroutine

   subroutine GetALE(a,isALE)
      use typre
      implicit none
      class(FemMesh) :: a
      logical :: isALE

      isALE = .FALSE.
      
      if (a%kfl_alemov == 1) then 
         isALE = .TRUE.
      else
      end if
   end subroutine

   subroutine SetDisplacements(a,Displacement)
      use typre
      implicit none
      class(FemMesh)        :: a
      real(rp), target      :: Displacement(:,:,:)

      a%displ=> Displacement
   end subroutine

   subroutine SetVelocities(a,Velocity)
      use typre
      implicit none
      class(FemMesh)        :: a
      real(rp), target      :: Velocity(:,:,:)
 
      a%veloc=>Velocity
   end subroutine
   
   subroutine GetMeshVeloc(a,Velocity)
      use typre
      implicit none
      class(FemMesh), target        :: a
      real(rp), pointer             :: Velocity(:,:,:)
 
      Velocity=>a%veloc
   end subroutine

   subroutine GetCoord(a,coord)
      use typre
      implicit none
      class(FemMesh), target  :: a
      real(rp), pointer       :: coord(:,:)

      coord=>a%coord
   end subroutine
   
   subroutine GetDeformedCoord(a,coord)
      use typre
      implicit none
      class(FemMesh), target  :: a
      real(rp), pointer       :: coord(:,:)

      call runend('deformed coord not ready, careful do not overwrite coord!')
   end subroutine

   subroutine DeallocExnorLpoty(a)
      use typre
      implicit none
      class(FemMesh)                :: a
     
      call a%Memor%dealloc(a%npoin,a%lpoty,'lpoty','exnor')
      call a%Memor%dealloc(a%ndime,a%ndime,a%nbopo,a%exnor,'exnor','exnor')
   end subroutine

   subroutine DeallocVmass(a)
      use typre
      implicit none
      class(FemMesh)        :: a
      
      call a%Memor%dealloc(a%npoin,a%vmass,'VMASS','vmass')
   end subroutine

   subroutine MeshGetTimes(a,cputim)
      implicit none
      class(FemMesh) :: a
      real(rp) :: cputim(25)
      
      call a%Timer%Total%GetValue(cputim(1))
      call a%Timer%reaLnods%GetValue(cputim(2))
      call a%Timer%Partitionate%GetValue(cputim(3))
      call a%Timer%LocalDimensionsAndIRenumbering%GetValue(cputim(4))
      call a%Timer%GhostPoints%GetValue(cputim(5))
      call a%Timer%BuildLocalOrdering%GetValue(cputim(6))
      call a%Timer%reaCoord%GetValue(cputim(7))
      call a%Timer%LnodsToLocal%GetValue(cputim(8))
      call a%Timer%LocalGraph%GetValue(cputim(9))
      call a%Timer%BuildArrayCommunicator%GetValue(cputim(10))
      call a%Timer%GhostCoord%GetValue(cputim(11))
      call a%Timer%PeriodicBoundaryConditions%GetValue(cputim(12))
      call a%Timer%Dombou%GetValue(cputim(13))
      call a%Timer%ExtnorLpoty%GetValue(cputim(14))
      call a%Timer%ComputeVmass%GetValue(cputim(15))
      call a%Timer%Blbopo%GetValue(cputim(16))
      call a%Timer%TestGhostCommunicate%GetValue(cputim(17))
      call a%Timer%TestAllToAll%GetValue(cputim(18))
      call a%Timer%TestReduce%GetValue(cputim(19))
      call a%Timer%TestBcast%GetValue(cputim(20))
      call a%Timer%TestBAndwidth%GetValue(cputim(21))
      call a%Timer%Refine%GetValue(cputim(22))
   end subroutine
   
   subroutine MeshGetRange(a,range)
      implicit none
      class(FemMesh) :: a
      real(rp) :: range(2,*)
      
      range(:,1:a%ndime) = a%range(:,1:a%ndime)
   end subroutine
   
   subroutine GetDeformedRange(a,range)
      implicit none
      class(FemMesh) :: a
      real(rp) :: range(2,*)
      integer(ip) :: idime
      
      !Recalculate local range
      do idime = 1,a%ndime
         a%range(1,idime) = minval(a%coord(idime,:) + a%displ(idime,:,1))
         a%range(2,idime) = maxval(a%coord(idime,:) + a%displ(idime,:,1))
      enddo
      range(:,1:a%ndime) = a%range(:,1:a%ndime)
   end subroutine
   
   subroutine MeshGetGRange(a,grange)
      implicit none
      class(FemMesh) :: a
      real(rp) :: grange(2,*)

      grange(:,1:a%ndime) = a%grange(:,1:a%ndime)
   end subroutine
   
   subroutine GetPosgp(a,posgp)
      implicit none
      class(FemMesh), target :: a
      real(rp), pointer :: posgp(:,:,:)
      
      posgp => a%posgp
   end subroutine

   subroutine GetLbody(a,lbody,iboun)
      implicit none
      class(FemMesh), target   :: a
      integer(ip), pointer     :: lbody
      integer(ip)              :: iboun
      
      lbody => a%lbody(iboun)
   end subroutine

   subroutine GetWeight(a,weight,iboun)
      implicit none
      class(FemMesh), target   :: a
      real (rp), pointer       :: weight
      integer(ip)              :: iboun
      
      weight => a%weight(iboun)
   end subroutine
   
   subroutine BuildLnodsFromArray(a,nelem,pnods,lnods)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: nelem, pnods(*), lnods(*)
      integer(ip) :: pnode, lnodssize
      
      a%nelem = nelem
      if (a%nelty == 1) then
         call a%Memor%alloc(1,a%pnods,'pnods','BuildLnodsFromArray')
         pnode = pnods(2)-pnods(1)
         a%pnods(1) = pnode
         lnodssize = pnode*a%nelem
      else
         call a%Memor%alloc(a%nelem+1,a%pnods,'pnods','BuildLnodsFromArray')
         a%pnods = pnods(1:a%nelem+1)
         lnodssize = a%pnods(a%nelem+1)-1
      endif
      
      call a%Memor%alloc(lnodssize,a%lnods,'lnods','BuildLnodsFromArray')
      a%lnods(1:lnodssize) = lnods(1:lnodssize)
   end subroutine
   
   subroutine SetNpoinAndBuildOrderingFromArray(a,gnpoin,npoin,npoinLocal,LocalToGlobal,ProcessorList,ndofn)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: gnpoin,npoin, npoinLocal, LocalToGlobal(*), ProcessorList(*),ndofn
      integer(ip) :: poin1
      
      a%gnpoin = gnpoin
      a%npoin = npoin
      a%npoinLocal = npoinLocal
      a%npoinGhost = npoin-npoinLocal
      
      !Create and Initialize the LocalOrdering
      call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      if (a%npoinLocal /= 0) then
         call a%LocalOrdering%Local2Global(1_ip,poin1)
         a%poin0 = poin1-1
      else
         a%poin0 = 0
      endif   
   end subroutine
   
   subroutine BuildCoordinatesFromArray(a,newcoord)
      implicit none
      class(FemMesh) :: a
      real(rp) :: newcoord(a%ndime,a%npoin)
      
      call a%Memor%alloc(a%ndime,a%npoin,a%coord,'coord','BuildCoordinatesFromArray')
      a%coord = newcoord
   end subroutine

   subroutine BuildBlocksFromArray(a,newblock)
      implicit none
      class(FemMesh) :: a
      real(rp) :: newblock(a%npoin)
      
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','BuildBlocksFromArray')
      a%geoblock= newblock
   end subroutine
   
   subroutine BuildExnorFromArray(a,isExnor,ExternalNormal)
      implicit none
      class(FemMesh) :: a
      logical     :: isExnor(:)
      type(r1p)   :: ExternalNormal(*)
      integer(ip) :: alloc_count,ipoin
      
      call a%Memor%alloc(a%npoin,a%isExnor,'isExnor','BuildExtnorFromArray')
      call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal', 'BuildExtnorFromArray')
      a%isExnor = isExnor
      alloc_count = 0
      do ipoin = 1,a%npoin
         if (a%isExnor(ipoin) .eqv. .true.) then
            alloc_count = alloc_count + rp*a%ndime
            allocate(a%ExternalNormal(ipoin)%a(a%ndime))
            a%ExternalNormal(ipoin)%a=ExternalNormal(ipoin)%a
         endif
      enddo
      call a%Memor%allocObj(0,'ExternalNormal','BuildExnorFromArray',alloc_count)
   end subroutine
   
   subroutine BuildListOfBodiesFromArray(a,lbody)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: lbody(*)
      
      call a%Memor%alloc(a%nboun,a%lbody,'lbody','BuildListOfBodiesFromArray')
      a%lbody = lbody(1:a%nboun)
   end subroutine
   
   subroutine BuildHangingNodesFromArray(a,pHangingList,lHangingList,rHangingList)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: pHangingList(*), lHangingList(*) 
      real(rp)    :: rHangingList(*)
      integer(ip) :: nz
      
      a%kfl_HangingNodes = .true.
      call a%Memor%alloc(a%npoin+1,a%pHangingList,'pHangingList','BuildHangingNodesFromArray')
      nz = pHangingList(a%npoin+1)-1
      call a%Memor%alloc(nz,a%lHangingList,'lHangingList','BuildHangingNodesFromArray')
      call a%Memor%alloc(nz,a%rHangingList,'rHangingList','BuildHangingNodesFromArray')
      
      a%pHangingList = pHangingList(1:a%npoin+1)
      a%lHangingList = lHangingList(1:nz)
      a%rHangingList = rHangingList(1:nz)
      
      call a%Memor%alloc(a%npoin,a%iwaHanging,'iwaHanging','BuildHangingNodesFromArray')
      call a%Memor%alloc(a%npoin,a%seenHanging,'seenHanging','BuildHangingNodesFromArray')
   end subroutine
   
   subroutine IsHangingNode(a,ipoin,isHanging)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ipoin, isHanging
      
      if (a%kfl_HangingNodes) then
         if (a%pHangingList(ipoin+1)-a%pHangingList(ipoin) > 0) then
            isHanging = 1
         else
            isHanging = 0
         endif
      else
         isHanging = 0
      endif
   end subroutine   
   
   subroutine IsHangingFace(a,iface,isHanging)
      implicit none
      class(FemMesh) :: a
      integer(ip) :: iface, isHanging
      
      if (a%kfl_HangingNodes) then
         if (a%auxBoundaryMatch2(1,iface) > 0) then
            isHanging = 1
         else
            isHanging = 0
         endif
      else
         isHanging = 0
      endif
   end subroutine
   
   subroutine IsAnAdaptiveMesh(a,IsAdaptive)
      implicit none
      class(FemMesh) :: a
      logical :: IsAdaptive 
      
      isAdaptive = a%IsSetAdaptiveRefiner
   end subroutine
   
   subroutine SetAdaptiveRefiner(a,Refiner)
      use typre
      implicit none
      class(FemMesh) :: a
      class(AdaptiveRefinerInterface), target :: Refiner
      
      a%IsSetAdaptiveRefiner = .true.
      a%Refiner => Refiner
   end subroutine
   
   subroutine GetHanging(a,kfl_HangingNodes)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_HangingNodes
      
      if (a%kfl_HangingNodes .eqv. .true.) then
         kfl_HangingNodes = 1
      else
         kfl_HangingNodes = -1
      endif
   end subroutine
   
   subroutine SetHanging(a,kfl_HangingNodes)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_HangingNodes
      
      if (kfl_HangingNodes == 1) then
         a%kfl_HangingNodes = .true.
      else
         a%kfl_HangingNodes = .false.
      endif
   end subroutine
   
   subroutine GetBoundaryNewToOld(a,BoundaryNewToOld)
      use typre
      implicit none
      class(FemMesh), target :: a
      integer(ip), pointer :: BoundaryNewToOld(:)
      
      BoundaryNewToOld => a%BoundaryNewToOld
   end subroutine
   
   subroutine IsManufactured(a,kfl_IsManufactured)
      use typre
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_IsManufactured
      
      kfl_IsManufactured = 0
      if (a%kfl_MeshInputType == 1) then
         kfl_IsManufactured = 1
      endif
   end subroutine
   
   subroutine GetIsoRange(a,pnode,isorange)
      implicit none
      class(FemMesh), target   :: a
      real(rp)                 :: isorange(2)
      integer(ip)              :: ndime,pnode

      if (a%ndime==2) then
         if (pnode==3 .or. pnode==6 .or. pnode==10) then
            isorange(1) = 0.0
            isorange(2) = 1.0
         elseif (pnode==4 .or. pnode==9 .or. pnode==16) then
            isorange(1) = -1.0
            isorange(2) = 1.0
         endif
      else
         if (pnode==4 .or. pnode==6 .or. pnode==10 .or. pnode==20) then
            isorange(1) = 0.0
            isorange(2) = 1.0
         elseif (pnode==8 .or. pnode==27 .or. pnode==64) then
            isorange(1) = -1.0
            isorange(2) = 1.0
         endif
      endif
   end subroutine
   
   subroutine MeshGetRebalanceInverseBoundaryMatches(a,iauxBoundaryMatch1,iauxBoundaryMatch2)
      implicit none
      class(FemMesh), target :: a
      integer(ip), pointer :: iauxBoundaryMatch1(:) 
      type(i1p), pointer   :: iauxBoundaryMatch2(:)
      
      iauxBoundaryMatch1 => a%iauxBoundaryMatch1
      iauxBoundaryMatch2 => a%iauxBoundaryMatch2
   end subroutine
   
   subroutine GetNeighboringPoints(a,ipoin,neigh,numneigh)
      implicit none
      class(FemMesh), target :: a
      integer(ip) :: numneigh, ipoin
      integer(ip), pointer :: neigh(:)

      numneigh = a%ia(ipoin+1)-a%ia(ipoin)
      neigh => a%ja(a%ia(ipoin):(a%ia(ipoin+1)-1))
   end subroutine
   
   subroutine GetArePeriodicDimensions(a,PerDime)
      implicit none
      class(FemMesh) :: a
      logical :: PerDime(3)
      
      PerDime = a%PerDime
   end subroutine

end module Mod_mesh

module Mod_IsotropicAdaptiveRefiner
   use Mod_AdaptiveInterface
   use typre
   use Mod_Memor
   use Mod_MyParallelOrdering, only : MyParallelOrdering
   use Mod_MyParallelOrdering, only : MyParallelOrdering
   !PETScCommunicator does not work, because it only communicates the first layer of ghost nodes it makes use of ia and ja!
   use Mod_PetscCommunicator
   use Mod_BuildGraph
   use Mod_GlobalElementIdentifier
   use Mod_RebalanceCommunicator
   use Mod_LinkedList
   use Mod_AdaptiveElementDependentSubroutines
   use Mod_Timer
   implicit none
   private
   public IsotropicAdaptiveRefiner, IsotropicAdaptiveRefiner_Const
   
   type :: ElemTypeList
      integer(1) :: ielty = -1
      integer(1) :: ivariation = -1
   end type
   
   type :: CSRList
      integer(ip), allocatable :: p(:)
      integer(ip), allocatable :: l(:)
      real(rp), allocatable    :: r(:)
   end type
   integer(ip), parameter :: MaxNodesPerElement = 25_ip
   
   type :: WeirdFaceList
      integer(ip), allocatable :: Faces(:)
      type(GlobalElementIdentifier), allocatable :: GEI(:)
      integer(ip), allocatable :: iface(:)
   end type
   
   type :: WeirdFaceList2
      integer(ip), allocatable :: Elements(:)
      type(GlobalElementIdentifier), allocatable :: GEI(:)
      integer(ip), allocatable :: iface(:)
   end type
   
   type :: AdaptiveTimers
      type(Timer) :: &
         Initialize, &
         PrepareOldExternalBoundariesInfoToKeep,&
         Refine,     &
            AmendMarkel,&
            Do21Balancing,&
            RefineSerial,&
            NumberingToParallel,&
            UpdateArraysToNewParallelNumbering,&
            HangingNodes,&
            SelectComponentsToExport,&
            MarkDeadElements,&
         LoadRebalance, &
            RebalanceRenumbering,&
            RebalanceElemProcs,&
            PrepareRebalanceFaces_a,&
            RebalanceSendNodesAndElements,&
            RebuildRefinementStructure,&
            PrepareRebalanceFaces_b,&
            PrepareRebalanceCommunicator
   end type
      
   
   type, extends(AdaptiveRefinerInterface) :: IsotropicAdaptiveRefiner
      
      integer(ip) :: ndime
      integer(ip) :: npoin, nelem      !Current number of points and elements
      
      type(ElemTypeList), allocatable :: ElementTypeList(:)
      
      integer(ip), allocatable :: pPointerToChildren(:)  !Pointer to the first child element
      integer(ip), allocatable :: lPointerToChildren(:)
      integer(ip), allocatable :: PointerToParent(:), OldPointerToParent(:),OldChildNumber(:)
      integer(ip), allocatable :: pedger(:), ledger(:,:)   !graph of the edge refinement
                                                           !ledger(1,ispos) = same level node of the edge
                                                           !ledger(2,ispos) = lower level node between ipoin and ledger(1,ispos)
                                                        !For higher order elements this will need to be enhanced.
      integer(ip) :: tnewpo, tnunrefpo
      integer(ip) :: tnewel, tnunref
      integer(ip) :: oldnpoin,oldnelem
      
      
      integer(ip), allocatable :: level(:)
      integer(ip), allocatable :: renpo(:), irenpo(:) !point renumbering from old to new mesh !only existing points in the
      
      !List of elements
      integer(ip), allocatable :: pnods(:),lnods(:)
      integer(ip), allocatable :: markel(:),oldmarkel(:)
      integer(ip), allocatable :: renel(:), irenel(:)
      
      !List of faces (points to the neighbour element)
      integer(ip), allocatable :: pfacs(:), lfacs(:,:)
      
      !Maximum dimensions
      integer(ip) :: mnode,mnodb,mface,mchild,mchildface
      
      !For parallel rearranging
      integer(ip) :: gnpoin
      integer(ip) :: npoinLocal,npoinGhost
      integer(ip) :: oldgnpoin
      integer(ip) :: oldnpoinLocal, oldnpoinGhost
      integer(ip), allocatable :: oldpnods(:), oldlnods(:)
      
      integer(ip) :: poin0
      integer(ip), allocatable :: ParallelLocalReordering(:)
      
      !Ordering and communicator
      type(MyParallelOrdering), pointer :: LocalOrdering => NULL(), OldLocalOrdering => NULL()
      type(MyParallelOrdering)      :: TestLocalOrdering
      type(PetscCommunicator) :: Communicator

!       integer(ip) :: npoinTop
!       integer(ip), allocatable :: EdgerEntryPoints(:)
      
      !Global Element Identifier
      type(GlobalElementIdentifier), allocatable :: GEI_ElementList(:), OldGEI_ElementList(:)
      type(MyParallelOrdering) :: ElementLocalOrdering
      type(MyParallelOrdering)      :: TestElementLocalOrdering
      
      !For HangingNodes
      integer(ip), allocatable :: pHangingList(:)
      integer(ip), allocatable :: lHangingList(:)
      real(rp), allocatable    :: rHangingList(:)
      type(CSRList) :: HangingListForElements
      integer(ip), allocatable :: isHangingNode(:)
      
      !For exporting to mesh
      integer(ip) :: ExportNelem,ExportLnodsSizes,ExportNpoin,ExportHangingSizes
      logical, allocatable :: ExportElement(:)
      integer(ip), allocatable :: ExportPointNumbering(:),ExportElementNumbering(:)
      integer(ip), allocatable :: iRenelExport(:)
      
      integer(ip)              :: OldExportNpoin,OldExportNelem
      integer(ip), allocatable :: OldExportPointNumbering(:),OldExportElementNumbering(:)
      type(PetscCommunicator) :: ExportCommunicator
      
      !For boundaries
      integer(ip), allocatable :: OldExternalBoundaries(:,:)
      type(LinkedListHeader), allocatable :: HeaderOldElementBoundaries(:)
      type(LinkedListStorage) :: StorageOldElementBoundaries
      
      
      !ExtraPoints
      integer(ip) :: ExtraNpoin,ExtraNelem,ExtraHangingSizes,ExtraLnodsSizes
      integer(ip), allocatable :: ExtraPoints(:,:)
      type(CSRList) :: ExtraHangingList, ExtraElements
      type(i1p), allocatable :: ExtraElementsListOthersNeed(:)
      integer(ip), allocatable :: ExtraNelemINeed(:), ExtraNelemOthersNeed(:)
      
      
      !Load Rebalancing
      integer(ip) :: kfl_RebalanceAlgorithm = 1
      
      integer(ip) :: RebalanceOldNelem
      integer(ip), allocatable :: RebalancePointNumbering(:), RebalanceProcessorList(:)
      type(i1p), allocatable   :: RebalanceElementProcessors(:)
      integer(ip) :: RebalanceOldNpoinLocal
      type(RebalanceCommunicator) :: RebalanceComm

      type(i1p), allocatable :: RebalanceListElementsINeed(:)
      type(i1p), allocatable :: RebalanceListElementsOthersNeed(:)
      integer(ip), allocatable :: RebalanceNelemINeed(:), RebalanceNelemOthersNeed(:)
      integer(ip), allocatable :: RebalanceNfaceINeed(:), RebalanceNfaceOthersNeed(:)
      type(WeirdFaceList), allocatable  :: RebalanceListFacesOthersNeed(:)
      type(WeirdFacelist2), allocatable :: RebalanceListFacesINeed(:)
      
      !Two-1-Balancing
      integer(ip) :: kfl_21Balancing = 0
      
      !Coordinates for deciding tetra variation
      real(rp), pointer :: coord(:,:) => NULL()
      
      !List of recently refined elements
      integer(ip) :: nJustRefinedElements
      integer(ip), allocatable :: lJustRefinedElements(:)
      type(AdaptiveTimers) :: Timers
      
      !For element information transfering
      integer(ip) :: ARET_ielty,ARET_ivariation,ARET_ichild,ARET_InterpolationTask
                   
      !ForDeallocation
      logical :: NotRefined = .true.
      
contains

      procedure :: Initialize
      procedure :: SetCoordArray
      procedure :: Refine

      procedure :: InitSerial
      procedure :: InitParallel
      procedure :: InitGlobalElementNumbering
      procedure :: AmendMarkel
      procedure :: RefineSerial
      procedure :: NumberingToParallel
      procedure :: UpdateArraysToNewParallelNumbering
      procedure :: HangingNodes
      procedure :: SelectComponentsToExport
      procedure :: AdditionalCommunicateHangingElements
      procedure :: MarkDeadElements
      procedure :: LocalElementFromGEI
      procedure :: CheckElementResponsibility
      procedure :: GetLnodsFromElement
      
      !For Transfering to Mesh
      procedure :: GetHangingListDimensions
      procedure :: GetHangingList
      procedure :: GetHangingFacesList
      procedure :: GetPointDimensions
      procedure :: GetElementDimensions
      procedure :: GetLocalOrdering
      procedure :: GetLnods
      procedure :: GetLevel
      !Info on boundaries
      procedure :: PrepareOldExternalBoundariesInfoToKeep
      procedure :: MatchOldLboelToOldExternalBoundaries
      procedure :: MatchNewLboelToNewElementsAndFaces
      procedure :: MatchOldBoundariesToNewBoundaries
      
      procedure :: UpdateVariableReal
      procedure :: UpdateVariableInteger
      procedure :: UpdateVariableLogical
      procedure :: UpdateVariableReal1
      procedure :: UpdateVariableInteger1
      procedure :: UpdateVariableLogical1
      
      !Load Rebalancing
      procedure :: LoadRebalance
      procedure :: RebalanceRenumbering
      procedure :: RebalanceElemProcs
      procedure :: RebalanceSendNodesAndElements
      procedure :: RebuildRefinementStructure
      procedure :: PrepareRebalanceCommunicator
      
      procedure :: RebalanceVariableReal
      procedure :: RebalanceVariableInteger
      procedure :: RebalanceVariableLogical
      procedure :: RebalanceVariableReal1
      procedure :: RebalanceVariableInteger1
      procedure :: RebalanceVariableLogical1
      !InfoOnBoundaries
      procedure :: PrepareRebalanceFaces_a
      procedure :: PrepareRebalanceFaces_b
      procedure :: GetiBoundaryMatch
      procedure :: GetiBoundaryToElementsAndFacesMatch
      procedure :: RebalanceVariableIntegerBoundary1
      procedure :: RebalanceVariableIntegerBoundary
      procedure :: RebalanceVariableRealBoundary1
      procedure :: RebalanceVariableRealBoundary
      procedure :: RebalanceVariableTyper1pBoundary
      procedure :: Set_2_1_Balancing
      procedure :: Do21Balancing
      procedure :: SetRebalancingNumbering
      
      !Turnof
      procedure :: WriteTimes
      procedure :: Dealloc
      
      procedure :: ElementalRefine
      procedure :: ElementalRebalance
      procedure :: InterpolateElement => ISO_InterpolateElement
      
      

   end type  
   
   type :: EdgerCriteriaTester
      integer(ip) :: parents(2)
      integer(ip) :: levels(2)
      integer(ip) :: orderedparents(2)
   end type
      
   
   interface
      subroutine BuildLfacs(a)
         import IsotropicAdaptiveRefiner
         implicit none
         class(IsotropicAdaptiveRefiner) :: a
      end subroutine
   end interface
   
   interface IsotropicAdaptiveRefiner_Const
      procedure constructor
   end interface IsotropicAdaptiveRefiner_Const

contains
   
   function constructor()
       class(IsotropicAdaptiveRefiner), pointer :: constructor

       allocate(constructor)

   end function constructor

   subroutine Initialize(a,InitialData)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerInitialData) :: InitialData
      
      call a%Timers%Initialize%Tic
      
      call a%InitSerial(InitialData)
      call a%InitParallel(InitialData)
      call a%InitGlobalElementNumbering
      
      call a%Timers%Initialize%Toc
      
   end subroutine

   subroutine SetCoordArray(a,coord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      real(rp), target :: coord(:,:)
      
      a%coord => coord
   end subroutine
   
   subroutine Refine(a,markelLastLevel)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: markelLastLevel(*)

      a%NotRefined = .false.
      
      call a%Timers%Refine%Tic
      
      call a%AmendMarkel(markelLastLevel)
      call a%Do21Balancing
      call a%RefineSerial
      call a%NumberingToParallel
      call a%UpdateArraysToNewParallelNumbering
      call a%HangingNodes  
      call a%SelectComponentsToExport
      call a%MarkDeadElements
      
      call a%Timers%refine%Toc
      
   end subroutine
   
   subroutine InitSerial(a,InitialData)
      use typre
      use Mod_Mesh
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerInitialData) :: InitialData
      
      integer(ip) :: npoin,nelem,ndime,nchild,nchildface,nface,pnodb
      integer(ip) :: ielem,ispos, lnodsSize,pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: ielty,ivariation
      
      a%ndime = InitialData%ndime
      a%npoin = InitialData%npoin
      a%nelem = InitialData%nelem
      
      a%oldnpoin = a%npoin
      a%oldnelem = a%nelem
      
      !Allocate some arrays
      !No children nor parents so far
      allocate(a%pPointerToChildren(a%nelem+1))
      allocate(a%lPointerToChildren(0))
      allocate(a%PointerToParent(a%nelem))
      allocate(a%OldPointerToParent(0))
      allocate(a%OldChildNumber(0))
      a%pPointerToChildren = 1
      a%PointerToParent = 0
      
      allocate(a%pedger(a%npoin+1), a%ledger(2,0))
      a%pedger = 1
      
!       allocate(a%pfacer(a%nelem+1), a%facer(2,0))
!       a%pfacer = 1
      
      !Build pnods and lnods for the mesh refiner
      !First we count
      lnodsSize = 0
      !We also set some dimensions
      a%mnode = 0
      a%mnodb = 0
      a%mface = 0
      a%mchild = 0
      a%mchildface = 0
      do ielem = 1,a%nelem
         call InitialData%GetLnode(ielem,pnode,lnode)
         lnodsSize = lnodsSize + pnode
         call AdaptiveGetElementType(a%ndime,pnode,ielty)
         call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
         a%mnode = max(a%mnode,pnode)
         a%mnodb = max(a%mnodb,pnodb)
         a%mchild = max(a%mchild,nchild)
         a%mface = max(a%mface,nface)
         a%mchildface = max(a%mchildface,nchildface)
      enddo
      
      allocate(a%pnods(a%nelem+1))
      allocate(a%lnods(lnodsSize))
      
      !Also ElementTypeList
      allocate(a%ElementTypeList(a%nelem))
      
      !Second we fill pnods and lnods
      a%pnods(1) = 1
      do ielem = 1,a%nelem
         call InitialData%GetLnode(ielem,pnode,lnode)
         a%pnods(ielem+1) = a%pnods(ielem)+pnode
         ispos = a%pnods(ielem)-1
         a%lnods(ispos+1:ispos+pnode) = lnode
         
         !ElementTypeList
         call AdaptiveGetElementType(a%ndime,pnode,ielty)
         a%ElementTypeList(ielem)%ielty = ielty
      enddo
      
      
      !Level: we initialize it to 0
      allocate(a%level(a%npoin))
      a%level(:) = 0
      
      !we initialize markel to 0
      allocate(a%markel(a%nelem))
      a%markel = 0
      
      !Build the list of faces
      call BuildLfacs(a)
      
      !HangingList
      allocate(a%pHangingList(a%npoin+1))
      allocate(a%lHangingList(0))
      allocate(a%rHangingList(0))
      a%pHangingList(:) = 1
      allocate(a%isHangingNode(a%npoin))
      a%isHangingNode = 0
      
   end subroutine
   
   subroutine InitParallel(a,InitialData)
      use typre
      use Mod_Mesh
      use MPI
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerInitialData) :: InitialData
      integer(ip), allocatable :: LocalToGlobal(:),ProcessorList(:)
      integer(ip), allocatable :: pelpo(:),lelpo(:),ia(:),ja(:)
      integer(ip) :: ipoin,auxnnode(2),ielem
      
      integer(ip) :: auxdimensions(10), auxdimensions2(10)
      integer(ip) :: ierr
      
      !Scatter maximum dimensions values to all processes
      auxdimensions(1) = a%mnode
      auxdimensions(2) = a%mnodb
      auxdimensions(3) = a%mchild
      auxdimensions(4) = a%mface
      auxdimensions(5) = a%mchildface
      call  MPI_AllREDUCE( auxdimensions, auxdimensions2, 5, MPI_INTEGER4, MPI_MAX,a%MPIcomm, ierr )
      a%mnode = auxdimensions2(1)
      a%mnodb = auxdimensions2(2) 
      a%mchild = auxdimensions2(3)
      a%mface = auxdimensions2(4)
      a%mchildface = auxdimensions2(5)
      
      
      !For parallel rearrangement
      a%npoinLocal = InitialData%npoinLocal
      a%npoinGhost = a%npoin- a%npoinLocal
      
      a%gnpoin = InitialData%gnpoin
      
      allocate(LocalToGlobal(a%npoin))
      allocate(ProcessorList(a%npoin))
      do ipoin = 1,a%npoin
         LocalToGlobal(ipoin) = ipoin
         ProcessorList(ipoin) = ipoin
      enddo
      call InitialData%Local2Global(a%npoin,LocalToGlobal,LocalToGlobal)
      call InitialData%GetProcessor(a%npoin,ProcessorList,ProcessorList)
      
      allocate(a%LocalOrdering)
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,LocalToGlobal,ProcessorList,a%Memor)
      
      !Build a communicator
      !call BuildElpo(pelpo,lelpo,a%npoin,a%nelem,a%pnods,a%lnods,2_ip,auxnnode,a%Memor)     !Compute glelpo and gpelpo
      !call BuildGraph(ia,ja,pelpo,lelpo,a%npoin,a%nelem,a%pnods,a%lnods,2_ip,auxnnode,a%Memor)
      allocate(ia(1), ja(1),lelpo(1),pelpo(1))
      call a%Communicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,a%npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)
      !For export communicator
      call a%ExportCommunicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,a%npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)
      
      deallocate(ia,ja,lelpo,pelpo)
      
      deallocate(LocalToGlobal, ProcessorList)
      
!       !EdgerEntryPoints ALWAYS contains the array of top level nodes
!       allocate(a%EdgerEntryPoints(a%npoin))
!       do ipoin = 1,a%npoin
!          a%EdgerEntryPoints(ipoin) = ipoin
!       enddo
!       a%npoinTop = a%npoin
      
      !For Exporting to Mesh
      !Elements
      allocate(a%ExportElement(a%nelem))
      a%ExportElement = .true.
      a%ExportNelem = a%nelem
      allocate(a%ExportElementNumbering(a%nelem))
      do ielem = 1,a%nelem
         a%ExportElementNumbering(ielem) = ielem
      enddo
      allocate(a%OldExportElementNumbering(0_ip))
      
      !Points
      allocate(a%ExportPointNumbering(a%npoin))
      do ipoin = 1,a%npoin
         a%ExportPointNumbering(ipoin) = ipoin
      enddo
      a%ExportNpoin = a%npoin
      !OldPoints
      allocate(a%OldExportPointNumbering(0_ip))
      a%OldExportNpoin = 0_ip
      !inverse export renumbering
      allocate(a%iRenelExport(a%ExportNelem))
      do ielem = 1,a%nelem
         a%irenelExport(ielem) = ielem
      enddo
      
      !ExtraPoints
      a%ExtraNpoin = 0
      allocate(a%ExtraPoints(1,1))
      allocate(a%ExtraHangingList%p(1))
      allocate(a%ExtraHangingList%l(1))
      allocate(a%ExtraHangingList%r(1))
      allocate(a%ExtraElements%p(1))
      allocate(a%ExtraElements%l(1))
   end subroutine   
   
   subroutine InitGlobalElementNumbering(a)
      use typre
      use Mod_HeapSortExternal
      use Mod_LinkedList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip) :: lprocs(MaxNodesPerElement), GlobalLnods(MaxNodesPerElement)
      
      type(LinkedListStorage) :: SToSend
      type(LinkedListHeader), allocatable :: HToSend(:)
      
      type(LinkedListStorage) :: SToReceive
      type(LinkedListHeader), allocatable :: HToReceive(:)
      
      type(i1p), allocatable :: ListToSend(:), ListToReceive(:),ArrayToSend(:),ArrayToReceive(:)
      
      logical, allocatable :: NeedToSend(:), NeedToReceive(:)
      integer(ip), allocatable :: Elem0(:),auxElem0(:),auxElementSend(:)
      
      integer(ip) :: minpoin,maxpoin,LocalElementCount,ielem,ispos,pnode
      integer(ip) :: iminloc,iminval,imaxloc,imaxval,iminvalLocal,iminrank
      integer(ip) :: iaux1,iaux2,ListSize,inode

      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      integer(ip) :: irank,jrank
      
      !I need to assign each element to a processor
      !Elements with all nodes local are mine
      !If some nodes are ghost, the element belongs to the processor with the lowest rank
      !This means that it belongs to the processor to which the node with the lowest global number belongs
      
      !Allocate the structure
      allocate(a%GEI_ElementList(a%nelem))
      
      !First I count how many elements are mine
      if (a%npoinLocal /= 0) then
         call a%LocalOrdering%Local2Global(1_ip,minpoin)
         call a%LocalOrdering%Local2Global(a%npoinLocal,maxpoin)
      else
         minpoin = -1
         maxpoin = -1
      endif
      
      allocate(HToReceive(a%MPIsize))
      allocate(HToSend(a%MPIsize))
      
      call SToReceive%Initialize(a%npoinLocal,'NonRemov','Repeatable',a%Memor)
      call SToSend%Initialize(a%npoinLocal,'NonRemov','Repeatable',a%Memor)
      
      !Count and mark local elements
      !Also count to whom I should ask the element
      !Also count to whom I should send the element
      LocalElementCount = 0
      
      allocate(auxElementSend(a%MPIsize))
      auxElementSend = 0_ip
      
      do ielem = 1,a%nelem
         ispos = a%pnods(ielem)-1
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         
         call a%LocalOrdering%Local2Global(pnode,a%lnods(ispos+1:ispos+pnode),GlobalLnods(1:pnode))
         iminloc = minloc(GlobalLnods(1:pnode),1)
         iminval = GlobalLnods(iminloc)
         imaxloc = maxloc(GlobalLnods(1:pnode),1)
         imaxval = GlobalLnods(imaxloc)
         !If the element is mine
         if (iminval >= minpoin .and. iminval <= maxpoin) then
            a%GEI_ElementList(ielem)%GlobalOriginalElement = 1
            LocalElementCount = LocalElementCount+1
            
            !If I need to send it, add it to the list to send
            if (imaxval > maxpoin) then
               call a%LocalOrdering%GetProc(pnode,a%lnods(ispos+1:ispos+pnode),lprocs)
               do inode = 1,pnode
                  if (lprocs(inode) /= a%MPIrank .and. auxElementSend(lprocs(inode)+1) /= ielem) then
                     !Add it to the list to send of this processor
                     call SToSend%Add(HToSend(lprocs(inode)+1),ielem)
                     
                     !This is to make sure that I do not send the same element twice
                     !to the same processor
                     auxElementSend(lprocs(inode)+1) = ielem
                  endif
               enddo
            endif
         
         !Else add it to the list to receive
         else
            iminvalLocal = a%lnods(ispos+iminloc)
            call a%LocalOrdering%GetProc(iminvalLocal,iminrank)
            call SToReceive%Add(HToReceive(iminrank+1),ielem)
         endif
      enddo
      
      !Allocate and fill the list of elements to be sent and received
      allocate(ListToReceive(a%MPIsize))
      allocate(ListToSend(a%MPIsize))
      
      allocate(NeedToSend(a%MPIsize))
      NeedToSend = .false.
      allocate(NeedToReceive(a%MPIsize))
      NeedToReceive = .false.
      
      !Transform lists to i1p
      call LinkedListToList(HToSend,SToSend,ListToSend,NeedToSend)
      call LinkedListToList(HToReceive,SToReceive,ListToReceive,NeedToReceive)
      
      !Build GlobalElementNumbering
      allocate(Elem0(a%MPIsize))
      allocate(auxElem0(a%MPIsize))
      auxElem0(1:a%MPIsize) = LocalElementCount
      call MPI_AllToAll(auxElem0, 1, MPI_INTEGER4, Elem0, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      deallocate(auxElem0)
      
      iaux1 = 0
      do irank = 0,a%MPIsize-1
         iaux2 = Elem0(irank+1)
         Elem0(irank+1) = iaux1
         iaux1 = iaux1+iaux2
      enddo
      
      !Set the global element numbering for local elements
      LocalElementCount = 0
      do ielem = 1,a%nelem
         if (a%GEI_ElementList(ielem)%GlobalOriginalElement == 1) then
            LocalElementCount = LocalElementCount+1
            
            a%GEI_ElementList(ielem)%GlobalOriginalElement = LocalElementCount+Elem0(a%MPIrank+1)
         endif
      enddo
      
      !Sort the element lists prior to communications (so that the order is the same in both sides)
      do irank = 0,a%MPIsize-1
         if (NeedToSend(irank+1) .eqv. .true.) then
            ListSize = size(ListToSend(irank+1)%l)
            call HeapSortExternal(ListSize,ListToSend(irank+1)%l,ElementIsLowerThan)
         endif
         
         if (NeedToReceive(irank+1) .eqv. .true.) then
            ListSize = size(ListToReceive(irank+1)%l)
            call HeapSortExternal(ListSize,ListToReceive(irank+1)%l,ElementIsLowerThan)
         endif
      enddo
      
      !Allocate actual arrays to send and receive
      !Fill actual arrays to send
      allocate(ArrayToSend(a%MPIsize))
      allocate(ArrayToReceive(a%MPIsize))
      do irank = 0,a%MPIsize-1
         !Send
         if (NeedToSend(irank+1) .eqv. .true.) then
            ListSize = size(ListToSend(irank+1)%l)
            allocate(ArrayToSend(irank+1)%l(ListSize))
            
            !The elements to send are local and I can fill their global numbering already
            ArrayToSend(irank+1)%l = a%GEI_ElementList(ListToSend(irank+1)%l)%GlobalOriginalElement
         endif
         
         !Receive
         if (NeedToReceive(irank+1) .eqv. .true.) then
            ListSize = size(ListToReceive(irank+1)%l)
            allocate(ArrayToReceive(irank+1)%l(ListSize))
         endif
      enddo
      
      !Allocate and Initialize
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      
      !Communicate the element global numbering
      do irank = 0,a%MPIsize-1
         if (NeedToSend(irank+1) .eqv. .true.) then
            ListSize = size(ArrayToSend(irank+1)%l)
            call MPI_ISEND(ArrayToSend(irank+1)%l,ListSize,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NeedToReceive(irank+1) .eqv. .true.) then
            ListSize = size(ArrayToReceive(irank+1)%l)
            call MPI_IRECV(ArrayToReceive(irank+1)%l,ListSize,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
         
      enddo

      !Wait for communications to finish
      do irank = 0,a%MPIsize-1
         if (NeedToSend(irank+1) .eqv. .true.)  call MPI_WAIT(irequest01(irank+1), status, ierr)
         if (NeedToReceive(irank+1) .eqv. .true.) call MPI_WAIT(irequest02(irank+1), status, ierr)
      enddo
      
      !Finally store the global element numbers
      do irank = 0,a%MPIsize-1
         if (NeedToReceive(irank+1) .eqv. .true.) then
            a%GEI_ElementList(ListToReceive(irank+1)%l)%GlobalOriginalElement = ArrayToReceive(irank+1)%l
         endif
      enddo
            
      !Deallocates      
      do irank = 0,a%MPIsize-1
         if (NeedToSend(irank+1) .eqv. .true.) then
            deallocate(ListToSend(irank+1)%l,ArrayToSend(irank+1)%l)
         endif
         if (NeedToReceive(irank+1) .eqv. .true.) then
            deallocate(ListToReceive(irank+1)%l,ArrayToReceive(irank+1)%l)
         endif
      enddo
      deallocate(ListToSend,ArrayToSend,ListToReceive,ArrayToReceive)
      deallocate(Elem0,auxElementSend)
      call SToSend%Dealloc
      call SToReceive%Dealloc
      deallocate(HToReceive,HToSend)
      deallocate(irequest01,irequest02)
      
      call BuildElementLocalOrdering(a%MPIcomm,a%nelem,a%GEI_ElementList,a%ElementLocalOrdering,a%Memor)

      
contains
      subroutine LinkedListToList(Header,Storage,List,NeedToCommunicate)
         use typre
         use Mod_LinkedList
         implicit none
         type(LinkedListHeader) :: Header(*)
         type(LinkedListStorage) :: Storage
         type(i1p) :: List(*)
         logical :: NeedToCommunicate(a%MPIsize)
         
         integer(ip) :: inelem,irank,ielem,icountlist
         
         NeedToCommunicate = .false.
         !Fill ListToSend
         do irank = 0,a%MPIsize-1
            call Header(irank+1)%GetNelem(inelem)
            
            if (inelem /= 0) then
               NeedToCommunicate(irank+1) = .true.
               
               allocate(List(irank+1)%l(inelem))
               
               call Storage%GoToFirstOfList(Header(irank+1))
               
               call Storage%GetNext(ielem)
               icountlist = 0
               do while (ielem /= -1) 
                  icountlist = icountlist+1
                  List(irank+1)%l(icountlist) = ielem
               
                  call Storage%GetNext(ielem)
               end do
            endif
         enddo
      end subroutine
      
      logical function ElementIsLowerThan(elem1,elem2)
         use typre
         implicit none
         integer(ip) :: elem1, elem2
         
         integer(ip) :: lnods1(MaxNodesPerElement),lnods2(MaxNodesPerElement)
         integer(ip) :: ispos1,ispos2,pnode1,pnode2
         
         ispos1 = a%pnods(elem1)-1
         pnode1 = a%pnods(elem1+1)-a%pnods(elem1)
         
         ispos2 = a%pnods(elem2)-1
         pnode2 = a%pnods(elem2+1)-a%pnods(elem2)
         

         if (pnode1 < pnode2) then
            ElementIsLowerThan = .true.
            return 
         elseif (pnode2 < pnode1) then
            ElementIsLowerThan = .false.
            return 
         else
            call a%LocalOrdering%Local2Global(pnode1,a%lnods(ispos1+1:ispos1+pnode1),lnods1)
            call a%LocalOrdering%Local2Global(pnode2,a%lnods(ispos2+1:ispos2+pnode2),lnods2)
            
            do inode = 1,pnode1
               if (lnods1(inode) < lnods2(inode)) then 
                  ElementIsLowerThan = .true.
                  return 
               elseif (lnods1(inode) > lnods2(inode)) then
                  ElementIsLowerThan = .false.
                  return 
               endif
            enddo
         endif
         
         ElementIsLowerThan = .false.
         return 
      end function
         
   end subroutine
   
   subroutine AmendMarkel(a,markelLastLevel)
      use typre
      use Mod_LinkedList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: markelLastLevel(*)
      
      integer(ip) :: ElementCount
      logical, allocatable :: ElementToCommunicate(:)
      
      integer(ip) :: nUnrefElements
      integer(ip), allocatable :: ListUnrefElements(:)
      
      type(ListGEI), allocatable :: GEI_NelemINeed(:),GEI_NelemOthersNeed(:)
      integer(ip) :: GEI_Size
      
      type(LinkedListHeader), allocatable :: HeaderListToAsk(:)
      type(LinkedListStorage) :: StorageListToAsk
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: irank,jrank
      logical, allocatable :: MightCommunicate(:)
      integer(ip), allocatable :: NelemINeed(:), NelemOthersNeed(:)
      type(i1p), allocatable :: MarkelINeed(:), MarkelOthersNeed(:)
      
      integer(ip), allocatable :: auxLocal2GlobalElementList(:,:)
      
      !Aux
      integer(ip), allocatable :: auxGhostNodes(:)
      
      integer(ip) :: nGhostChildren,nUnrefChildren,nLocalChildren
      integer(ip) :: ielem,elemi,iparent,ippos,pchild,jchild,jelem,jspos,jpoin
      integer(ip) :: ipoin,icount,itopParentLocal,iTopParentGlobal,ichild
      integer(ip) :: ifoundElement,ilevel,TopLevelNelem
      type(GlobalElementIdentifier) :: iGEI
      type(GlobalElementIdentifier) :: GEI_Sizer
      
      call a%Timers%AmendMarkel%Tic

      !First of all we need to go from markelLastLevel to markel
      !For later
      allocate(ListUnrefElements(a%nelem))
      nUnrefElements = 0
      
      ElementCount = 0
      do ielem = 1,a%nelem
         if (a%ExportElement(ielem) .eqv. .true.) then
            ElementCount = ElementCount+1
            if (markelLastLevel(ElementCount) == 1) then
               a%markel(ielem) = 2
            elseif (markelLastLevel(ElementCount) == 0) then
               a%markel(ielem) = 0
            elseif (markelLastLevel(ElementCount) == -1) then
               a%markel(ielem) = -1
            endif   
            
            !We also prepare a list of unrefined elements for later use
            if (a%markel(ielem) == -1) then
               nUnrefElements = nUnrefElements+1
               ListUnrefElements(nUnrefElements) = ielem
            endif
            
         endif
      enddo
      
      !If I have an unrefined element, I need to know if its siblings
      !are also unrefined
      !If the siblings are ghost-ghost elements (not computed locally by 
      !the mesh and the physical problem) then I need to communicate markel
      !for these siblings
      !But this is only if all the local siblings are marked as -1
      
      !Auxiliary
      allocate(ElementToCommunicate(a%nelem))
      ElementToCommunicate = .false.
      
      allocate(HeaderListToAsk(a%MPIsize))
      call StorageListToAsk%Initialize(2*a%npoinGhost,'NonRemovable','Repeatable',a%Memor)
      
      !Loop through unrefined elements
      do elemi = 1,nUnrefElements
         ielem = ListUnrefElements(elemi)
            
         iparent = a%PointerToParent(ielem)
         if (iparent /= 0) then
            ippos = a%pPointerToChildren(iparent)-1
            pchild = a%pPointerToChildren(iparent+1)-a%pPointerToChildren(iparent)
            
            nGhostChildren = 0
            nUnrefChildren = 0
            nLocalChildren = 0
            do jchild = 1,pchild
               jelem = a%lPointerToChildren(ippos+jchild)
               
               if (a%ExportElement(jelem) .eqv. .true.) then
                  nLocalChildren = nLocalChildren+1
                  if (a%markel(jelem) == -1) then
                     nUnrefChildren = nUnrefChildren+1
                  endif
               else
                  !Last level but not exported, this means GHOST!
                  if (a%markel(jelem) == 0) then
                     nGhostChildren = nGhostChildren+1
                  endif   
               endif
            enddo
            
            !If all local elements are marked as unrefine
            !then I need to communicate all ghost elements
            if ((nLocalChildren == nUnrefChildren) .and. (nGhostChildren /= 0)) then
               do jchild = 1,pchild
                  jelem = a%lPointerToChildren(ippos+jchild)
                  if ( (a%markel(jelem) == 0) .and. (a%ExportElement(jelem) .eqv. .false.) .and. (ElementToCommunicate(jelem) .eqv. .false.)) then
                     ElementToCommunicate(jelem) = .true.
                                       
                     !Decide to which processor should I ask each element
                     !Since all the nodes of the element are ghost, I will
                     !simply ask it to the processor of the first node
                     jspos = a%pnods(jelem)-1
                     jpoin = a%lnods(jspos+1)
                     call a%LocalOrdering%GetProc(jpoin,jrank)
                     
                     !Add to list
                     call StorageListToAsk%Add(HeaderListToAsk(jrank+1),jelem)
                  endif
               enddo
            endif
         endif
      enddo
      deallocate(ListUnrefElements,ElementToCommunicate)
      
      !At most, I will communicate with those processors with whom I have a ghost node
!       allocate(MightCommunicate(a%MPIsize))
!       MightCommunicate = .false.
!       allocate(auxGhostNodes(a%npoinGhost))
!       do ipoin = 1,a%npoinGhost
!          auxGhostNodes(ipoin) = a%npoinLocal+ipoin
!       enddo
!       call a%LocalOrdering%GetProc(a%npoinGhost,auxGhostNodes,auxGhostNodes)
!       do ipoin = 1,a%npoinGhost
!          MightCommunicate(auxGhostNodes(ipoin)+1) = .true.
!       enddo
!       deallocate(auxGhostNodes)
         
      
      !Allocate the communications structure
      allocate(NelemINeed(a%MPIsize),GEI_NelemINeed(a%MPIsize),MarkelINeed(a%MPIsize))
      allocate(NelemOthersNeed(a%MPIsize),GEI_NelemOthersNeed(a%MPIsize),MarkelOthersNeed(a%MPIsize))
      
      
      !Allocate communications arrays
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      allocate(irequest03(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      
      !For knowing the sizes to send
      GEI_Size = sizeof(GEI_Sizer)
      
      !Communicate to my neighbours how many of their elements I need
      NelemINeed = 0
      NelemOthersNeed = 0
      do irank = 0,a%MPIsize-1
         call HeaderListToAsk(irank+1)%GetNelem(NelemINeed(irank+1))
      enddo
      
      call MPI_AllToAll(NelemINeed, 1, MPI_INTEGER4, NelemOthersNeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);      
      
      !Communicate to my neighbours how many of their elements I need
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) then
            allocate(GEI_NelemINeed(irank+1)%a(NelemINeed(irank+1)))
            call StorageListToAsk%GoToFirstOfList(HeaderListToAsk(irank+1))
            
            !Load from linked list to transfer list
            call StorageListToAsk%GetNext(ielem)
            icount = 1
            do while (ielem /= -1)
               GEI_NelemINeed(irank+1)%a(icount) = a%GEI_ElementList(ielem)
            
               icount = icount+1
               call StorageListToAsk%GetNext(ielem)
            enddo
            
            !Send GEI_NelemINeed
            call MPI_ISEND(GEI_NelemINeed(irank+1)%a,NelemINeed(irank+1)*GEI_Size,MPI_BYTE, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
         endif
        
      enddo
      
      !Receive the elements others need from me
      do irank = 0,a%MPIsize-1
         !if (MightCommunicate(irank+1) .eqv. .true.) then
            if (NelemOthersNeed(irank+1) /= 0) then
               allocate(GEI_NelemOthersNeed(irank+1)%a(NelemOthersNeed(irank+1)))
               call MPI_IRECV(GEI_NelemOthersNeed(irank+1)%a,NelemOthersNeed(irank+1)*GEI_Size,MPI_byte,irank,mtag2,a%MPIcomm,irequest03(irank+1),ierr) 
            endif
         !endif
      enddo
      !deallocate(MightCommunicate)
      
      !Wait for communications to finish
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      
      !Allocate to send arrays
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) then
            allocate(MarkelINeed(irank+1)%l(NelemINeed(irank+1)))
         endif
         if (NelemOthersNeed(irank+1) /= 0) then
            allocate(MarkelOthersNeed(irank+1)%l(NelemOthersNeed(irank+1)))
         endif
      enddo
      
      !Find the elements in my processor so that I can mark markel
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1)/= 0) then
            do ielem= 1,NelemOthersNeed(irank+1)
               iGEI =  GEI_NelemOthersNeed(irank+1)%a(ielem)
               call a%LocalElementFromGEI(iGEI,ifoundElement)
               
               !Store the markel value
               MarkelOthersNeed(irank+1)%l(ielem) = a%markel(ifoundElement)
            enddo
         endif
      enddo
      
      !Send and receive MarkelOthersNeed and MarkelINeed
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      
      !Send MarkelOthersNeed, Receive MarkelINeed
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(MarkelOthersNeed(irank+1)%l,NelemOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NelemINeed(irank+1) /= 0) then
            call MPI_IRECV(MarkelINeed(irank+1)%l,NelemINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
         endif
      enddo
      !Wait for communications to finish
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Finally loop through the list I need and update markel
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) then
            call StorageListToAsk%GoToFirstOfList(HeaderListToAsk(irank+1))
               
            !Load from linked list
            call StorageListToAsk%GetNext(ielem)
            icount = 1
            do while (ielem /= -1)
               !I only overwrite if it is -1, otherwise I leave it as it is
               if (MarkelINeed(irank+1)%l(icount) == -1_ip) then
                  a%markel(ielem) = -1_ip
               endif
               
               !Next element
               call StorageListToAsk%GetNext(ielem)
               icount = icount+1
            enddo
         endif
      enddo
      
      !Now we deallocate everything
      deallocate(HeaderListToAsk)
      call StorageListToAsk%Dealloc
      
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) deallocate(MarkelINeed(irank+1)%l,GEI_NelemINeed(irank+1)%a)
         if (NelemOthersNeed(irank+1) /= 0) deallocate(MarkelOthersNeed(irank+1)%l,GEI_NelemOthersNeed(irank+1)%a)
      enddo
      deallocate(NelemINeed,NelemOthersNeed,MarkelINeed,MarkelOthersNeed,GEI_NelemINeed,GEI_NelemOthersNeed)
      deallocate(irequest01,irequest02,irequest03)
      
      call a%Timers%AmendMarkel%Toc
      
   end subroutine   
   
   subroutine Do21Balancing(a)
      use typre
      use Mod_LinkedList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip), pointer :: lnods(:) => NULL()
      
      integer(ip) :: nfacs,ielem,ifpos,jface,jelem,inode,pnode,ipoin,ichild,maxlevel,minlevel,level,leveldiff,ispos,parent,pchild
      integer(ip), allocatable :: MaxMinLevels(:,:)
      
      call a%Timers%Do21Balancing%Tic
      
      if (a%kfl_21Balancing == 1) then
         allocate(MaxMinLevels(2,a%npoin))
         MaxMinLevels(1,:) = 1e6
         MaxMinLevels(2,:) = -1
      
         !Maximum and minimum level of the elements to which a node belongs
         do ielem = 1,a%nelem
            if (a%ExportElement(ielem) .eqv. .true.) then
               call a%GetLnodsFromElement(ielem,pnode,lnods)
               level = a%GEI_ElementList(ielem)%level
               
               !testing
               if (a%markel(ielem) == 2) level = level+1
               
               do inode = 1,pnode
                  ipoin = lnods(inode)
                  MaxMinLevels(1,ipoin) = min(MaxMinLevels(1,ipoin),level)
                  MaxMinLevels(2,ipoin) = max(MaxMinLevels(2,ipoin),level)
               enddo
            endif
         enddo
         call a%Communicator%GhostCommunicate(2_ip,MaxMinLevels)
               
         do ielem = 1,a%nelem
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            
            !If hanging nodes and trying to refine, do not refine
            do inode = 1,pnode
               ipoin = lnods(inode)
               if (a%isHangingNode(ipoin) /= 0) then
                  if (a%markel(ielem) == 2) then
                     a%markel(ielem) = 0
                  endif
               endif
            enddo
            
            !do not allow 3 level jumps between elements which share a node
            maxlevel = maxval(MaxMinLevels(2,lnods))
            minlevel = minval(MaxMinLevels(1,lnods))
            level = a%GEI_ElementList(ielem)%level
            leveldiff = maxlevel-minlevel
            if (a%markel(ielem) == -1 .and. level < maxlevel) then
               a%markel(ielem) = 0
            elseif (a%markel(ielem) == 2 .and. level > minlevel) then
               a%markel(ielem) = 0
            endif
         enddo
         
         deallocate(MaxMinLevels)
      endif
      
      call a%Timers%Do21Balancing%Toc
   end subroutine   
   
   subroutine RefineSerial(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      !Working arrays
      integer(ip), allocatable :: eliwa(:)
      integer(ip), allocatable :: pelpo(:), lelpo(:)
      
      integer(ip) :: oldnpoin, newel, newLnodsSize, oldLnodsSize
      integer(ip) :: newLChildSize,oldLChildSize,newLfacsSize,oldLFacsSize
      integer(ip) :: newpo,nunref,nunrefpo,tnunref,nz,tnpoin,tnze
      
      integer(ip), allocatable :: itemp(:), itemp2(:)
      
      integer(ip), allocatable :: seen(:), iwa(:,:)
      integer(ip), allocatable :: itemppPointerToChildren(:),itemplPointerToChildren(:),itempPointerToParent(:)
      integer(ip), allocatable :: itempPfacs(:), itempLfacs(:,:)
      integer(ip), allocatable :: tpedger(:), tledger(:,:)
      
      
      integer(ip) :: ielem,ispos,elemi,ichild,childi,isposnew,inode,ipoin,iface, ifpos
      integer(ip) :: ifposnew, ispos2,iparent,icposnew,icpos,ichild2
      integer(ip) :: iz
      integer(ip) :: pnode,pnodb,nchild,nface,nchildface
      
      integer(ip) :: jelem,elemj,jspos,jface,jnchild,jnchildface,jnode,jnface,jpnodb,jpnode,jpoin
      integer(ip) :: jv
      
      integer(ip) :: kelem,elemk,kelem2
      integer(ip) :: lelem,lspos,eleml,lnchild,lnface,lnchildface,lpnodb,lpnode
      integer(ip) :: lchild,lface

      integer(ip) :: ElementCount,iaux,auxnnode(2)
      
      type(EdgerCriteriaTester) :: ECT
      type(GlobalElementIdentifier), allocatable :: itempGEI(:)
      
      type(ElemTypeList), allocatable :: itempETL(:)
      
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: newnpointop
      integer(ip) :: tnewinteriorpoints
      
      integer(ip), allocatable :: touched(:)
      
      integer(ip) :: ielty, jelty,lelty,jvariation,lvariation
      real(rp) :: elcod(a%ndime,25)
      integer(ip) :: lnode(25)
      integer(ip) :: ivariation
      
      integer(ip) :: jnnode,jnodelist(25),nodei,nodej, testipoin,iadd,ninteriornodes,pointlevel
      integer(ip) :: kchild,knode,kpnode,kpoin,nInTheFaceNodes
      integer(ip), pointer :: klnods(:) => NULL()
      integer(ip) :: tnewInTheFaceNodes
      
      call a%Timers%RefineSerial%Tic
      
      a%oldnpoin = a%npoin
      a%oldnelem = a%nelem   
      
      newel = 0
      newLnodsSize = 0
      oldLnodsSize = size(a%lnods)
      newLfacsSize = 0
      oldLfacsSize = size(a%lfacs,2)
      newLChildSize = 0
      oldLChildSize = size(a%lPointerToChildren)
      
      allocate(eliwa(a%nelem))
      eliwa(:) = 0
      
      !Refined Element counter
      a%nJustRefinedElements = 0
      !1rst, loops to count dimensions
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 2) then  !refine
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            ielty = a%ElementTypeList(ielem)%ielty
            call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
            
            !Since I am going to refine, I need to set the variation (for tetras)
            !First from internal to external numbering (coord is in external)
            lnode(1:pnode) = a%lnods(ispos+1:ispos+pnode)
            lnode(1:pnode) = a%ExportPointNumbering(lnode(1:pnode))
            elcod(:,1:pnode) = a%coord(:,lnode(1:pnode))
            call AdaptiveGetElementVariation(ielty,elcod,ivariation)
            a%ElementTypeList(ielem)%ivariation = ivariation
            
            !Counting and dimensions
            a%nJustRefinedElements = a%nJustRefinedElements + 1
            newel = newel+nchild
            newLnodsSize = newLnodsSize + nchild*pnode
            newLfacsSize = newLfacsSize + nchild*nface
            newLChildSize = newLChildSize + nchild
         elseif (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then !unrefine
            if (a%PointerToParent(ielem) /= 0) then !the element has a parent 
               iparent = a%PointerToParent(ielem)
               eliwa(iparent) = eliwa(iparent) + 1 !we will only unrefine if all the child are marked as unrefined
            else
               a%markel(ielem) = 0 !this element is top level, can not be unrefined
            endif     
         endif
      enddo
      
      !Initialize list of just refined elements
      if (allocated(a%lJustRefinedElements)) then
         deallocate(a%lJustRefinedElements)
      endif
      allocate(a%lJustRefinedElements(a%nJustRefinedElements))
      
      a%nJustRefinedElements = 0
      !2nd loop to count dimensions, it knows truly unrefined elements
      tnunref = 0
      do ielem = 1,a%nelem
         nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
         pnode = a%pnods(ielem+1) - a%pnods(ielem)
         
         if (a%markel(ielem) == 1 .and. eliwa(ielem) == nchild) then !this is a parent with all of its child marked as -1, unrefine
            !Count dimensions update
            ispos = a%pPointerToChildren(ielem)-1
            do ichild = 1,nchild
               jelem = a%lPointerToChildren(ispos + ichild)
               jspos = a%pnods(jelem)-1
               jpnode = a%pnods(jelem+1) - a%pnods(jelem)
               jelty = a%ElementTypeList(jelem)%ielty
               call ElementGetDimensions(jelty,pnodb,nchild,nface,nchildface)
               
               !For measuring sizes
               tnunref = tnunref +1
               newLnodsSize = newLnodsSize-jpnode
               newLfacsSize = newLfacsSize-nface
            enddo
            newLChildSize = newLChildSize - nchild

            a%markel(ielem) = 0 !the element is not refined anymore
         elseif (a%markel(ielem) == 1 .and. eliwa(ielem) > 0) then !not all the child are -1, do not unrefine
            ispos = a%pPointerToChildren(ielem)-1
            do ichild = 1,nchild
               jelem = a%lPointerToChildren(ispos + ichild)
               if (a%markel(jelem) == -1) then
                  a%markel(jelem) = 0
               endif
            enddo    
         elseif (a%markel(ielem) == -2) then
!             tnunref = tnunref + 1
!             !update dimensions taking into account unrefinement
!             call ElementGetDimensions(a%ndime,pnode,pnodb,nchild,nface,nchildface)
!             newLnodsSize = newLnodsSize-pnode
!             newLfacsSize = newLfacsSize-nface
         elseif (a%markel(ielem) == 2) then
            !Just refined elements
            a%nJustRefinedElements = a%nJustRefinedElements + 1
            a%lJustRefinedElements(a%nJustRefinedElements) = ielem
         endif   
      enddo
      
      !new element numbers (unref 2)
      nunref = 0
      do ielem = 1,a%nelem
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then
            eliwa(ielem) = 0
            nunref = nunref +1
         else
            eliwa(ielem) = ielem - nunref  !eliwa has now the new element numbers
         endif    
      enddo
      
      !Update the List of JustRefinedElements to new element numbering
      a%lJustRefinedElements = eliwa(a%lJustRefinedElements)
       
      !reallocate pnods and lnods
      allocate(itemp(a%nelem+newel-nunref+1))
      allocate(itemp2(oldLnodsSize+newLnodsSize))
      
      !reallocate ElementList (GlobalElementIdentifier)
      allocate(itempGEI(a%nelem+newel-nunref))
      
      !reallocate ElementTypeList
      allocate(itempETL(a%nelem+newel-nunref))
      
      !--------------------------------------------------
      !Store new number of elements
      a%tnewel = newel
      a%tnunref = tnunref
      !--------------------------------------------------
      
      nunref = 0
      itemp(1) = 1
      do ielem = 1,a%nelem
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then
            nunref = nunref+1
         else
            !Update Lnods
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            isposnew = itemp(ielem-nunref)-1
            itemp2(isposnew+1:isposnew+pnode) = a%lnods(ispos+1:ispos+pnode)
            itemp(ielem-nunref+1) = itemp(ielem-nunref)+pnode
            
            !Update ElementList
            itempGEI(ielem-nunref) = a%GEI_ElementList(ielem)
            
            !Update ElementList
            itempETL(ielem-nunref) = a%ElementTypeList(ielem)
         endif
      enddo
      call move_alloc(a%pnods,a%oldpnods)
      call move_alloc(a%lnods,a%oldlnods)
      
      call move_alloc(itemp,a%pnods)
      call move_alloc(itemp2,a%lnods)
      
      !GlobalElementIdentifier LIst
      if (allocated(a%OldGEI_ElementList)) deallocate(a%OldGEI_ElementList)
      
      call move_alloc(a%GEI_ElementList,a%OldGEI_ElementList)
      call move_alloc(itempGEI,a%GEI_ElementList)
      
      !ElementTypeList
      call move_alloc(itempETL,a%ElementTypeList)
      
      !reallocate lfacs
      !we read in the old numbering (ielem)
      !and we write in the new numbering (elemi)
      !pnods and lnods are already in the new numbering (elemi)
      !PointerToParent and PointerToChildren are still in the old numbering (ielem)
      allocate(itempPfacs(a%nelem+newel-nunref+1))
      allocate(itempLfacs(2,oldLFacsSize+newLfacsSize))
      itempPfacs(1) = 1
      tnunref = nunref
      nunref = 0
      do ielem = 1,a%nelem
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then
               nunref = nunref+1
         else
            elemi = eliwa(ielem)
            ispos = a%pnods(elemi)-1
            pnode = a%pnods(elemi+1)-a%pnods(elemi)
            ielty = a%ElementTypeList(elemi)%ielty
            call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
            
            itempPfacs(elemi+1) = itempPfacs(elemi)+nface
            ifposnew = itempPfacs(elemi)-1
               
            ifpos = a%pfacs(ielem)-1   
            do iface = 1,nface
               jelem = a%lfacs(1,ifpos+iface)
               jface = a%lfacs(2,ifpos+iface)
               
               
               if (jelem > 0 .and. jelem <= a%nelem) then !this is a normal face connection
                  elemj = eliwa((jelem))
                  
                  if (elemj == 0) then !the neighbour is gone
                     
                     !The element and the parents are dead
                     if (a%markel(jelem) == -2) then
                        itempLfacs(1,ifposnew+iface) = 0
                        itempLfacs(2,ifposnew+iface) = 0
                        
                     !The neighbour is gone, natural unrefinement, look up
                     else !if(a%markel(jelem) == -1) then
                  

                        !I need their info!!
                        lelem= a%PointerToParent(jelem)
                        eleml = eliwa(lelem)
                        lpnode = a%pnods(eleml+1)-a%pnods(eleml)
                        lelty = a%ElementTypeList(eleml)%ielty
                        lvariation = a%ElementTypeList(eleml)%ivariation
                        call ElementGetDimensions(lelty,lpnodb,lnchild,lnface,lnchildface)
                        
                        !Get the child number (need to loop through all the children)
                        lnchild = a%pPointerToChildren(lelem+1)-a%pPointerToChildren(lelem)
                        lspos = a%pPointerToChildren(lelem)-1
                        do kelem = 1,lnchild
                           kelem2 = a%lPointerToChildren(lspos+kelem)
                           if (kelem2 == jelem) lchild = kelem
                        enddo
                        !Get the face of its parent to which I connect
                        lface = jface  !entry face for seeking up faces
                        call up_faces(lelty,lvariation,lnface,lchild,lface)
                        
                        itempLfacs(1,ifposnew+iface) = -eleml
                        itempLfacs(2,ifposnew+iface) = lface
                     endif
                  else
                     !just change to new element numbering, same face
                     itempLfacs(1,ifposnew+iface) = elemj
                     itempLfacs(2,ifposnew+iface) = jface
                  endif 
               elseif (jelem > a%nelem) then !this is a boundary exterior face (no estic segur de què fa aquí)
                  itempLfacs(1,ifposnew+iface) = jelem + newel - tnunref
                  itempLfacs(2,ifposnew+iface) = jface

               elseif (jelem < 0) then           !this is a Discontinuous-Galerkin face (different levels)   
                  elemj = eliwa(abs(jelem))
                  if (elemj == 0) then !the neighbour is gone
                     !The element and their parents are dead
                     if (a%markel(abs(jelem)) == -2) then
                        itempLfacs(1,ifposnew+iface) = 0
                        itempLfacs(2,ifposnew+iface) = 0
                     
                     !Normal unrefinement, look up
                     else !if (a%markel(abs(jelem)) == -1) then
                        lelem= a%PointerToParent(-jelem)
                        eleml = eliwa(lelem)
                        lpnode = a%pnods(eleml+1)-a%pnods(eleml)
                        lelty = a%ElementTypeList(eleml)%ielty
                        lvariation = a%ElementTypeList(eleml)%ivariation
                        call ElementGetDimensions(lelty,lpnodb,lnchild,lnface,lnchildface)
                     
                        !Get the child number (need to loop and test)
                        lnchild = a%pPointerToChildren(lelem+1)-a%pPointerToChildren(lelem)
                        lspos = a%pPointerToChildren(lelem)-1
                        do kelem = 1,lnchild
                           kelem2 = a%lPointerToChildren(lspos+kelem)
                           if (kelem2 == -jelem) lchild = kelem
                        enddo
                        !Get the face of its parent to which I connect
                        lface = jface !entry face for up faces
                        call up_faces(lelty,lvariation,lnface,lchild,lface)
                     
                        itempLfacs(1,ifposnew+iface) = -eleml
                        itempLfacs(2,ifposnew+iface) = lface
                     endif
                  else
                     itempLfacs(1,ifposnew+iface) = -elemj
                     itempLfacs(2,ifposnew+iface) = jface
                  endif   
               elseif (jelem == 0) then
                  itempLfacs(1,ifposnew+iface) = 0
                  itempLfacs(2,ifposnew+iface) = 0
               endif    
            enddo   
         endif
      enddo
      call move_alloc(itempPfacs,a%pfacs)
      call move_alloc(itempLFacs,a%lfacs)    
      
      !Build renel and irenel
      !Build renel and irenel
      if (allocated(a%renel)) deallocate(a%renel)
      if (allocated(a%irenel)) deallocate(a%irenel)
      allocate(a%renel(a%nelem))
      a%renel = 0_ip
      allocate(a%irenel(a%nelem+newel-nunref))
      a%irenel = 0_ip
      nunref = 0
      do ielem = 1,a%nelem
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then
            nunref = nunref+1
            
            !We store the Pointer to the parent in the OLD numbering
            a%renel(ielem) = -a%PointerToParent(ielem)
         else
            !renel and irenel
            a%renel(ielem) = ielem - nunref
            a%irenel(ielem-nunref) = ielem
         endif
      enddo
      
      
      
      !Reallocation of pchild and pparent
      allocate(itemppPointerToChildren(a%nelem+newel-nunref+1))
      allocate(itempLPointerToChildren(oldLChildSize+newLChildSize))
      allocate(itempPointerToParent(a%nelem+newel-nunref))
      
      itemppPointerToChildren(1) = 1
      nunref = 0
      newel = 0
      do ielem = 1,a%nelem
         !unrefined element
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then 
               nunref = nunref+1
         else
            elemi = eliwa(ielem)
            !Who is my parent
            if (a%PointerToParent(ielem) /= 0) then
               itempPointerToParent(elemi) = eliwa(a%PointerToParent(ielem))
            else
               itempPointerToParent(elemi) = 0
            endif
            
            !Now for the children
            !Element without children
            if (a%markel(ielem) == 0) then
               itemppPointerToChildren(elemi+1) = itemppPointerToChildren(elemi)
            
            !Element with old children
            !element with same status as in the previous iteration, just change element numbering if required
            elseif (a%markel(ielem) == 1) then
               ispos = a%pPointerToChildren(ielem)-1
               nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
               
               elemi = eliwa(ielem)
               icposnew = itemppPointerToChildren(elemi)-1
               
               !starting index, next element
               itemppPointerToChildren(elemi+1) = itemppPointerToChildren(elemi)+nchild
               !Update my children to my new numbering
               do ichild = 1,nchild
                  ichild2 = a%lPointerToChildren(ispos+ichild)
                  itemplPointerToChildren(icposnew+ichild) = eliwa(ichild2)
               enddo
                     
           
            !element with new children
             elseif(a%markel(ielem) == 2) then      
               !I add the new elements in new numbering:
               elemi = eliwa(ielem)
               ispos = a%pnods(elemi)-1
               pnode = a%pnods(elemi+1)-a%pnods(elemi)
               ielty = a%ElementTypeList(elemi)%ielty
               call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
               
               !Set the element children
               itemppPointerToChildren(ielem-nunref+1) = itemppPointerToChildren(ielem-nunref)+nchild
               icposnew = itemppPointerToChildren(elemi)-1
               do ichild = 1,nchild
                  itemplPointerToChildren(icposnew+ichild) = a%nelem-tnunref+newel+ichild
               enddo
               
               !Set the element as the parent of its children
               do ichild = 1,nchild
                  itempPointerToParent(a%nelem-tnunref+newel+ichild) = elemi
               enddo
               
               !Set the Global Element Identifier (GEI)
               do ichild = 1,nchild
                  !I also set the Global Element Identifier for the new children
                  call GEI_SetupFromParent(a%GEI_ElementList(a%nelem-tnunref+newel+ichild),a%GEI_ElementList(elemi),ichild)
               enddo
               
               !Set the ElementTypeList (ETL)
               !I cannot set the variation, in order to do that I need the coordinates of the nodes, which I do not have
               !I will ask for the variation when I am going to refine
               do ichild = 1,nchild
                  a%ElementTypeList(a%nelem-tnunref+newel+ichild)%ielty = a%ElementTypeList(elemi)%ielty
               enddo    
               
               !Update the current number of new elements
               newel = newel+nchild
            endif   
               
         endif
      enddo
      
      !Create OldChildNumber which is required for transfering data between elements later
      deallocate(a%OldChildNumber)
      allocate(a%OldChildNumber(a%nelem))
      do ielem = 1,a%nelem
         nchild = a%pPointerToChildren(ielem+1)- a%pPointerToChildren(ielem)
         ispos = a%pPointerToChildren(ielem)-1
         do ichild = 1,nchild
            childi = a%lPointerToChildren(ispos+ichild)
            a%OldChildNumber(childi) = ichild
         enddo
      enddo   
      
      !Moveallocs
      itemppPointerToChildren(a%nelem-tnunref+1:a%nelem-tnunref+newel+1) = itemppPointerToChildren(a%nelem-tnunref+1)
      call move_alloc(a%PointerToParent,a%OldPointerToParent)
      call move_alloc(itempPointerToParent,a%PointerToParent)
      call move_alloc(itemppPointerToChildren,a%pPointerToChildren)
      call move_alloc(itemplPointerToChildren,a%lPointerToChildren)
      

      
      !reallocate markel and mark touched
      allocate(touched(a%oldnpoin))
      touched = 0
      allocate(itemp2(a%nelem+newel-nunref))
      itemp2 = 0
      nunref = 0
      do ielem = 1,a%nelem
         if (a%markel(ielem) == -1 .or. a%markel(ielem) == -2) then
            nunref = nunref+1
         else
            itemp2(ielem-nunref) = a%markel(ielem)
            call a%GetLnodsFromElement(eliwa(ielem),pnode,lnods)
            touched(lnods) = 1
         endif
      enddo
      call move_alloc(a%markel,a%oldmarkel)
      call move_alloc(itemp2,a%markel)
      
      
      
!------------------------------------------------------------------      
       !THE NUMBER OF ELEMENTS IS CHANGED DUE TO UNREFINEMENT!
       !Element numbering is the new one for all arrays after this point
       a%nelem = a%nelem-tnunref
!------------------------------------------------------------------       
      
      !We need to do the same for points and go to the new numbering!
      !Now we want to update edger, ledger and build the remaining (refined elements) parts of lnods
      
      !First, Build pelpo and lelpo only with refined elements
      allocate(pelpo(a%npoin+1))
      pelpo(:) = 0
      pelpo(1) = 1
      
      !First to count
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         ispos = a%pnods(ielem)-1
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         ielty = a%ElementTypeList(ielem)%ielty
         call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
   
         do inode = 1,pnode
            ipoin = a%lnods(ispos+inode)
            pelpo(ipoin+1) = pelpo(ipoin+1)+1
         enddo
      enddo
         
      !Compress pelpo
      do ipoin=1,a%npoin
         pelpo(ipoin+1)=pelpo(ipoin+1)+pelpo(ipoin)
      end do

      !2nd to fill
      !Create lelpo 
      allocate(lelpo(pelpo(a%npoin+1)))   
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         ispos = a%pnods(ielem)-1
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         ielty = a%ElementTypeList(ielem)%ielty
         call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
         
         !Fill lelpo
         do inode = 1,pnode
            ipoin = a%lnods(ispos+inode)
            
            ispos2 = pelpo(ipoin)
            lelpo(ispos2) = ielem
            pelpo(ipoin) = pelpo(ipoin)+1
         enddo
      enddo
      
      !Correct pelpo
      do ipoin = a%npoin,2,-1
         pelpo(ipoin) = pelpo(ipoin-1)
      enddo
      pelpo(1) = 1
      
      
      !Let us see which are the new node numbers
      allocate(seen(a%npoin))
      allocate(iwa(2,a%npoin))
      seen = 0
      iwa = 0
      
      
      
      
      !1rst loop to determine the new size of edger
      newpo = 0
      do ipoin = 1,a%npoin
         !Decompress already existing refined edges
         nz = 0
         do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
            !decompress
            jpoin = a%ledger(1,jv)
            if (touched(a%ledger(2,jv)) /= 0) then
               nz = nz+1
               seen(jpoin) = nz
               iwa(1,nz) = jpoin
               iwa(2,nz) = a%ledger(2,jv)
            endif
         enddo
         
         !add new edges from the current refinement
         do jv = pelpo(ipoin),pelpo(ipoin+1)-1
            jelem = lelpo(jv)
            jspos = a%pnods(jelem)-1
            pnode = a%pnods(jelem+1)-a%pnods(jelem)
            jelty = a%ElementTypeList(jelem)%ielty
            call ElementGetDimensions(jelty,pnodb,nchild,nface,nchildface)
                           
            !Seek for inode
            nodeido: do nodei = 1,pnode
               testipoin = a%lnods(jspos+nodei)
               if (testipoin == ipoin) then
                  inode = nodei
                  exit nodeido
               endif
            enddo nodeido
            !Get the connected edge nodes
            call Adaptive_GetEdgeNodesForInode(jelty,inode,jnnode,jnodelist)
            
            do nodej = 1,jnnode
               jnode = jnodelist(nodej)
               jpoin = a%lnods(jspos+jnode)
               
               if (ipoin /= jpoin) then
                  ECT%parents(1) = ipoin
                  ECT%parents(2) = jpoin
                  ECT%levels(1) = a%level(ipoin)
                  ECT%levels(2) = a%level(jpoin)
                  call EdgerCriteria(ECT)
                  
                  if (ipoin == ECT%orderedparents(1)) then
                     if (seen(jpoin)/= 0) then 
                        
                     else
                        nz = nz +1
                        newpo = newpo+1
                        
!                         auxcounter(1,newpo) = ipoin
!                         auxcounter(2,newpo) = jpoin
                        
                        seen(jpoin) = nz
                        iwa(1,nz) = jpoin
                        iwa(2,nz) = newpo+a%npoin
                     endif
                  endif
               endif
            enddo 
         enddo
         !clean seen
         do iz = 1,nz
            seen(iwa(1,iz)) = 0
         enddo
      enddo
      
      tnewinteriorpoints = 0
      !Loop through elements and add interior nodes which do not go to edges nor faces (quads, hexas)
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         call Adaptive_GetNumberOfInteriorNodesToAdd(int(a%ElementTypeList(ielem)%ielty,ip),nInteriorNodes)
         newpo = newpo+nInteriorNodes
         tnewinteriorpoints = tnewinteriorpoints + nInteriorNodes
      enddo
      
      tnewInTheFaceNodes = 0
      !Loop through elements and add in the face nodes which do not go to edges nor faces (quads, hexas)
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         ielty = a%ElementTypeList(ielem)%ielty
         call Adaptive_GetNumberOfInTheFaceNodesToAdd(ielty,nInTheFaceNodes)
         do inode = 1,nInTheFaceNodes
            call Adaptive_GetFaceForInTheFaceNode(ielty,inode,iface)
            
            ispos = a%pfacs(ielem)-1
            jelem = a%lfacs(1,ispos+iface)
            jface = a%lfacs(2,ispos+iface)
            
            !Same level neighbour
            if (jelem > 0 .and. jelem <= a%nelem) then
               !If it is already refined in previous steps, do nothing
               if (a%markel(jelem) == 1) then
                  
               !If it is going to be refined now, count only once   
               elseif (a%markel(jelem) == 2) then
                  if (jelem > ielem) then
                     tnewInTheFaceNodes = tnewInTheFaceNodes + 1
                     newpo = newpo + 1
                  endif
               
               !Same level neighbour, it is not going to be refined, I count it
               elseif (a%markel(jelem) == 0) then
                     tnewInTheFaceNodes = tnewInTheFaceNodes + 1
                     newpo = newpo + 1
               endif
            
            !neighbour is higher level, I count it
            else
               tnewInTheFaceNodes = tnewInTheFaceNodes + 1
               newpo = newpo + 1
            endif
         enddo
      enddo
      
      
      
      
      !new point numbering due to unrefinement
      if (allocated(a%renpo)) deallocate(a%renpo)
      allocate(a%renpo(a%npoin))
      nunrefpo = 0
      a%renpo(:) = 0
      do ipoin = 1,a%npoin
         if (touched(ipoin) == 0) then
            nunrefpo = nunrefpo +1
         else
            a%renpo(ipoin) = ipoin - nunrefpo
         endif
      enddo

      !--------------------------------------
      !Store dimensions
      a%tnewpo = newpo
      a%tnunrefpo = nunrefpo
      !--------------------------------------
      
      !new dimensions
      tnpoin = a%npoin + newpo - nunrefpo  
      
      !reallocate level (before touched!!)
      nunrefpo = 0
      newnpointop = 0
      allocate(itemp2(tnpoin))
      do ipoin = 1,a%npoin
         if (touched(ipoin) == 0) then
            nunrefpo = nunrefpo + 1
         else
            itemp2(a%renpo(ipoin)) = a%level(ipoin)
            if (a%level(ipoin) == 0) newnpoinTop = newnpoinTop+1
         endif
      enddo
      itemp2(a%npoin-nunrefpo+1:tnpoin) = 0   
      call move_alloc(itemp2,a%level)
      
      !Number of edges
      tnze = tnpoin - newnpointop -tnewinteriorpoints -tnewInTheFaceNodes
      
      !deallocate touched
      deallocate(touched)
      
      !modify lnods to the updated point numbering
      do ispos = 1,a%pnods(a%nelem+1)-1
         iaux = a%renpo(a%lnods(ispos))
         a%lnods(ispos) = iaux
      enddo
      
      !modify pnods to the updated elements
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 2) then
            pnode  =a%pnods(ielem+1)-a%pnods(ielem)
            icpos = a%pPointerToChildren(ielem)-1
            nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
            do ichild = 1,nchild
               childi = a%lPointerToChildren(icpos+ichild)
               a%pnods(childi+1) = a%pnods(childi)+pnode
            enddo
         endif
      enddo
      
      !2nd loop to fill
      allocate(tpedger(tnpoin+1),tledger(2,tnze))  
      tpedger(1) = 1
      newpo = 0
      seen(:) = 0
      do ipoin = 1,a%npoin
         if (a%renpo(ipoin) > 0) then
               !Decompress already existing refined edges
               !only points which are not in the lower-most level will have 
               !something in the graph, do not worry about unrefine
               nz = 0
               do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1 
                  jpoin = a%renpo(a%ledger(1,jv))  !take into account new point numbering
                  if (a%renpo(a%ledger(2,jv)) /= 0) then
                  
                     nz = nz+1
                     !decompress
                     seen(jpoin) = nz
                     iwa(1,nz) = jpoin
                     iwa(2,nz) = (a%renpo(a%ledger(2,jv)))  !take into account new point numbering
                  endif    
               enddo
            
               !add new edges from the current refinement
               do jv = pelpo(ipoin),pelpo(ipoin+1)-1
                  jelem = lelpo(jv)
                  jspos = a%pnods(jelem)-1
                  pnode = a%pnods(jelem+1)-a%pnods(jelem)
                  jelty = a%ElementTypeList(jelem)%ielty
                  jvariation = a%ElementTypeList(jelem)%ivariation
                  call ElementGetDimensions(jelty,pnodb,nchild,nface,nchildface)
                  
                  !first add the ipoin vertex to the new elements
                  do jnode = 1,pnode
                     jpoin = a%lnods(jspos+jnode)
                     if (jpoin == a%renpo(ipoin)) then
                        inode = jnode
                        
                        lspos = a%pPointerToChildren(jelem)-1
                        call add_node_to_child(jelty,jvariation,nchild,inode,jnode,a%pnods,a%lnods,a%renpo(ipoin),a%lPointerToChildren(lspos+1))
                     endif
                  enddo
               
                  !secondly add the vertices on the edges (refined) to the new elements and the refined edge list
                  !Get the connected edge nodes
                  call Adaptive_GetEdgeNodesForInode(jelty,inode,jnnode,jnodelist)
                  do nodej = 1,jnnode
                     jnode = jnodelist(nodej)
                     jpoin = a%lnods(jspos+jnode)
                     
                     if (a%renpo(ipoin) /= jpoin) then
                        ECT%parents(1) = a%renpo(ipoin)
                        ECT%parents(2) = jpoin
                        ECT%levels(1) = a%level(a%renpo(ipoin))
                        ECT%levels(2) = a%level(jpoin)
                        call EdgerCriteria(ECT)
                        
                        if (a%renpo(ipoin) == ECT%orderedparents(1)) then
                           if (seen(jpoin) /= 0) then  !already existing refined edge
                              
                           else !new edge refinement
                              nz = nz +1
                              newpo = newpo+1
                              seen(jpoin) = nz
                              iwa(1,nz) = jpoin
                              iwa(2,nz) = newpo+a%npoin-nunrefpo
                              a%level(iwa(2,nz)) = max(a%level(jpoin),a%level(a%renpo(ipoin)))+1
                           endif
                           
                           lspos = a%pPointerToChildren(jelem)-1
                           call add_node_to_child(jelty,jvariation,nchild,inode,jnode,a%pnods,a%lnods,iwa(2,seen(jpoin)),a%lPointerToChildren(lspos+1))
                        endif
                     endif
                  enddo 
               enddo
               !clean seen
               do iz = 1,nz
                  seen(iwa(1,iz)) = 0
               enddo
               !compress and store in tedger (temporary) !with the new renumbering!!!
               tpedger(a%renpo(ipoin)+1) = tpedger(a%renpo(ipoin))+nz
               tledger(1:2,tpedger(a%renpo(ipoin)):tpedger(a%renpo(ipoin)+1)-1) = iwa(1:2,1:nz)
                   
         endif
      enddo
      tpedger(a%npoin+2-nunrefpo:tnpoin+1) = tpedger(a%npoin+1-nunrefpo)
      call move_alloc(tpedger,a%pedger)
      call move_alloc(tledger,a%ledger)
      
      !Loop through elements and add nodes which do not go to edges nor faces (quads, hexas)
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         call a%GetLnodsFromElement(ielem,pnode,lnods)
         pointlevel = maxval(a%level(lnods))+1
         
         call Adaptive_GetNumberOfInteriorNodesToAdd(int(a%ElementTypeList(ielem)%ielty,ip),nInteriorNodes)
         ielty = a%ElementTypeList(ielem)%ielty
         ivariation = a%ElementTypeList(ielem)%ivariation
         do iadd = 1,nInteriorNodes
            newpo = newpo + 1
            ipoin = newpo+a%npoin-nunrefpo
            call Adaptive_AddNodeToChildInterior(ielty,ivariation,iadd,ipoin,a%lPointerToChildren(a%pPointerToChildren(ielem)),a%pnods,a%lnods)
            a%level(ipoin) = pointlevel
         enddo
      enddo
      
      tnewInTheFaceNodes = 0
      !Loop through elements and add in the face nodes which do not go to edges nor faces (quads, hexas)
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         call a%GetLnodsFromElement(ielem,pnode,lnods)
         pointlevel = maxval(a%level(lnods))+1
         
         ielty = a%ElementTypeList(ielem)%ielty
         ivariation = a%ElementTypeList(ielem)%ivariation
         call Adaptive_GetNumberOfInTheFaceNodesToAdd(ielty,nInTheFaceNodes)
         do inode = 1,nInTheFaceNodes
            call Adaptive_GetFaceForInTheFaceNode(ielty,inode,iface)
            
            ispos = a%pfacs(ielem)-1
            jelem = a%lfacs(1,ispos+iface)
            jface = a%lfacs(2,ispos+iface)
            
            !Same level neighbour
            if (jelem > 0 .and. jelem <= a%nelem) then
               jelty = a%ElementTypeList(jelem)%ielty
               jvariation = a%ElementTypeList(jelem)%ivariation
               
               !If it is already refined in previous steps, just ask the neighbour for the node number
               if (a%markel(jelem) == 1) then
                  !Ask to neighbour
                  
                  call Adaptive_GetInTheFacePointFromParent(jelty,jface,kchild,knode)
                  jspos = a%pPointerToChildren(jelem)-1
                  kelem = a%lPointerToChildren(jspos+kchild)
                  call a%GetLnodsFromElement(kelem,kpnode,klnods)
                  kpoin = klnods(knode)
                  
                  ispos2 = a%pPointerToChildren(ielem)-1
                  !Set in ALL of my children
                  call Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,kpoin,a%lPointerToChildren(ispos2+1),a%pnods,a%lnods)
                  
                  
               !If it is going to be refined now, count only once   
               elseif (a%markel(jelem) == 2) then
                  if (jelem > ielem) then
                     tnewInTheFaceNodes = tnewInTheFaceNodes + 1
                     newpo = newpo + 1
                     
                     !Set in ALL of my children
                     ispos2 = a%pPointerToChildren(ielem)-1
                     ipoin = newpo+a%npoin-nunrefpo
                     call Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,ipoin,a%lPointerToChildren(ispos2+1),a%pnods,a%lnods)
                     a%level(ipoin) = pointlevel
                  else
                     !Ask to neighbour
                     call Adaptive_GetInTheFacePointFromParent(jelty,jface,kchild,knode)
                     jspos = a%pPointerToChildren(jelem)-1
                     kelem = a%lPointerToChildren(jspos+kchild)
                     call a%GetLnodsFromElement(kelem,kpnode,klnods)
                     kpoin = klnods(knode)
                     
                     !Set in ALL of my children
                     ispos2 = a%pPointerToChildren(ielem)-1
                     call Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,kpoin,a%lPointerToChildren(ispos2+1),a%pnods,a%lnods)
                 endif    
               
               !Same level neighbour, neighbour is not going to be refined, I count it
               elseif (a%markel(jelem) == 0) then
                  tnewInTheFaceNodes = tnewInTheFaceNodes + 1
                  newpo = newpo + 1
                  
                  !Set in ALL of my children
                  ispos2 = a%pPointerToChildren(ielem)-1
                  ipoin = newpo+a%npoin-nunrefpo
                  call Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,ipoin,a%lPointerToChildren(ispos2+1),a%pnods,a%lnods)
                  a%level(ipoin) = pointlevel
               endif
            
            !neighbour is higher level, I count it
            else
               tnewInTheFaceNodes = tnewInTheFaceNodes + 1
               newpo = newpo + 1
               
               !Set in ALL of my children
               ispos2 = a%pPointerToChildren(ielem)-1
               ipoin = newpo+a%npoin-nunrefpo
               call Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,ipoin,a%lPointerToChildren(ispos2+1),a%pnods,a%lnods)
               a%level(ipoin) = pointlevel
            endif
         enddo
      enddo
      
      a%nelem = a%nelem+newel
      a%npoin = tnpoin
      a%tnewpo = newpo

      deallocate(pelpo,lelpo)
      deallocate(eliwa)
      deallocate(seen,iwa)
      
      !update the list of faces, new elements
      !first we set pfacs
      !modify pfacs to the updated elements
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 2) then
            !first we set pfacs for all the children
            pnode  =a%pnods(ielem+1)-a%pnods(ielem)
            ielty = a%ElementTypeList(ielem)%ielty
            call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
            icpos = a%pPointerToChildren(ielem)-1
            do ichild = 1,nchild
               childi = a%lPointerToChildren(icpos+ichild)
               a%pfacs(childi+1) = a%pfacs(childi)+nface
            enddo
         endif
      enddo
      
      !second we fill lfacs
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 2) then
            call upd_facs_newel(a,ielem)
            
            !we are done with the refinement of this element
            a%markel(ielem) = 1
         endif
      enddo  
   call a%Timers%RefineSerial%Toc
      
   end subroutine
   
   subroutine EdgerCriteria(ECT)
      implicit none
      type(EdgerCriteriaTester) :: ECT
      
      if (ECT%parents(1) == ECT%parents(2) ) then
         ECT%orderedParents = -1
      else
         
         if (ECT%levels(1) == ECT%levels(2)) then
            ECT%orderedparents(1) = minval(ECT%parents)
            ECT%orderedparents(2) = maxval(ECT%parents)
         elseif (ECT%levels(1) > ECT%levels(2)) then
            ECT%orderedparents(1) = ECT%parents(1)
            ECT%orderedparents(2) = ECT%parents(2)
         else
            ECT%orderedparents(1) = ECT%parents(2)
            ECT%orderedparents(2) = ECT%parents(1)
         endif
      endif
   end subroutine
   
   subroutine upd_facs_newel(a,ielem)
      !This subroutine updates the faces of the elements when a new element is added
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      integer(ip):: ielem,jchild,iface,jface,jface2,jelem,jfacnod(a%mnodb),nchild
      integer(ip) :: nface,pnode
      integer(ip) :: facchild(2,a%mchildface),kfacchild,kface,kelem,isequal,kfacnod(a%mnodb),jfpos(a%mchild)
      integer(ip) :: ifpos,icpos,jcpos,jnchild,jnchildface,jnface,jpnodb,jpnode,jspos
      integer(ip) :: knchild,knchildface,knface,kpnodb,kpnode,kspos
      integer(ip) :: lcpos,lelem,lface,lnchild,lnchildface,lnface,lpnodb,lpnode,lspos
      integer(ip) :: lfacchild(2,a%mchildface)
      integer(ip), pointer :: son(:) => NULL()
      integer(ip) :: seek_value(2), replace_value(2),lelty,jelty,kelty,ielty,ivariation,lvariation
      integer(ip) :: childj,facej,inchild,jjchild(8),jjnface(8),jjface(6,8),kfpos

      !index i represents the parent
      !index j represents the children
      !index l represents the parent's neighbour
      !index k represents the children of the parent's neighbour

      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      icpos = a%pPointerToChildren(ielem)-1
      nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
      ifpos = a%pfacs(ielem)-1

      son => a%lPointerToChildren(icpos+1:icpos+nchild)
      do jchild = 1,nchild
         jfpos(jchild) = a%pfacs(son(jchild))-1
      enddo

      ielty = a%ElementTypeList(ielem)%ielty
      ivariation = a%ElementTypeList(ielem)%ivariation
      call AdaptiveConnectInteriorFaces(ielty,ivariation,jfpos,son,a%lfacs)
      
      !Get the list of faces which need to be checked against other elements
      call AdaptiveGetChildrenExternalFacesList(ielty,ivariation,inchild,jjchild,jjnface,jjface)
      
      !Loop through all the external faces I need to check
      do childj = 1,inchild
         jchild = jjchild(childj)
         do facej = 1,jjnface(childj)
            jface = jjface(facej,childj)
         
            iface = jface
            call up_faces(ielty,ivariation,nface,jchild,iface) !get my parent face
            if (a%lfacs(1,ifpos+iface) < 0) then  ! parent is connected to upper
                  a%lfacs(1,jfpos(jchild)+jface) = a%lfacs(1,ifpos+iface)
                  a%lfacs(2,jfpos(jchild)+jface) = a%lfacs(2,ifpos+iface)
            elseif (a%lfacs(1,ifpos+iface) == 0) then
                  a%lfacs(1,jfpos(jchild)+jface) = a%nelem+ielem !parent face is in the boundary (level0 element)
                  a%lfacs(2,jfpos(jchild)+jface) = iface         !we mark it with element numbering past nelem
            elseif (a%lfacs(1,ifpos+iface) > a%nelem) then
                  a%lfacs(1,jfpos(jchild)+jface) = a%lfacs(1,ifpos+iface)
                  a%lfacs(2,jfpos(jchild)+jface) = a%lfacs(2,ifpos+iface)
            else  !my parent is connected to his neighbour
                  lelem = a%lfacs(1,ifpos+iface)
                  lface = a%lfacs(2,ifpos+iface)
                  lcpos = a%pPointerToChildren(lelem)-1
                  lspos = a%pnods(lelem)-1
                  lpnode = a%pnods(lelem+1)-a%pnods(lelem)
                  lelty = a%ElementTypeList(lelem)%ielty
                  lvariation = a%ElementTypeList(lelem)%ivariation
                  call ElementGetDimensions(lelty,lpnodb,lnchild,lnface,lnchildface)
                  
                  lnchild = a%pPointerToChildren(lelem+1)- a%pPointerToChildren(lelem)
                  if (lnchild == 0) then !the neighbour of my parent has no child
                     a%lfacs(1,jfpos(jchild)+jface) = -lelem
                     a%lfacs(2,jfpos(jchild)+jface) = lface
                  else !the difficult case, move along the tree and update   
                     jelem = son(jchild)
                     jspos = a%pnods(jelem)-1
                     jpnode = a%pnods(jelem+1)-a%pnods(jelem)
                     jelty = a%ElementTypeList(jelem)%ielty
                     call ElementGetDimensions(jelty,jpnodb,jnchild,jnface,jnchildface)
                     jcpos = a%pPointerToChildren(jelem)-1
                     
                     call face_from_el(jelty,jnface,jpnodb,a%lnods(jspos+1:jspos+jpnode),jface,1,jfacnod)
                  
                     if (lvariation < 1 .or. lvariation > 3) then
                        !write(*,*) 'uops'
                     endif
                     
                     call down_faces(lelty,lvariation,lface,lfacchild)
                     do kfacchild = 1,lnchildface
                        kface = lfacchild(2,kfacchild)
                        kelem = lfacchild(1,kfacchild)  !neighbour son number
                        kelem = a%lPointerToChildren(lcpos+kelem)
                        kspos = a%pnods(kelem)-1
                        kpnode = a%pnods(kelem+1)-a%pnods(kelem)
                        kelty = a%ElementTypeList(kelem)%ielty
                        call ElementGetDimensions(kelty,kpnodb,knchild,knface,knchildface)
                        kfpos = a%pfacs(kelem)-1
                        
                        call face_from_el(kelty,knface,kpnodb,a%lnods(kspos+1),kface,1,kfacnod)
                        call face_compare(jfacnod,kfacnod,jpnodb,isequal)
                        if (isequal==1) then
                           a%lfacs(1,jfpos(jchild)+jface) = kelem
                           a%lfacs(2,jfpos(jchild)+jface) = kface
                           
                           a%lfacs(1,kfpos+kface) = jelem
                           a%lfacs(2,kfpos+kface) = jface

                           seek_value(1) = -ielem
                           seek_value(2) = iface
                           
                           replace_value(1) = -jelem
                           replace_value(2) = jface
                           call rectree(a,kelem,kface,seek_value,replace_value)
                        endif
                     enddo
                  endif
            endif
         enddo
      enddo

   end subroutine
   
   recursive subroutine rectree(a,parent,parent_face,seek_value,replace_value)
      !This subroutine moves through all the descendents of a given element seeking for a neighbour face and replacing it
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip):: parent,parent_face,seek_value(2),replace_value(2)

      integer(ip) :: iface,jelem,facchild(2,a%mchildface),jface
      integer(ip) :: icpos,jfpos,nchild,nface,nchildface,pnode,pnodb,ielty,ivariation

      !If I have children
      nchild = a%pPointerToChildren(parent+1)-a%pPointerToChildren(parent)
      icpos = a%pPointerToChildren(parent)-1
      if (nchild /= 0) then
         !Dimensions of the element
         pnode = a%pnods(parent+1)-a%pnods(parent)
         ielty = a%ElementTypeList(parent)%ielty
         ivariation = a%ElementTypeList(parent)%ivariation
         call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
         
         !Go
         call down_faces(ielty,ivariation,parent_face,facchild)
         do iface = 1,nchildface
            jelem = a%lPointerToChildren(icpos+facchild(1,iface))
            jface = facchild(2,iface)
            jfpos = a%pfacs(jelem)-1
            if (a%lfacs(1,jfpos+jface) == seek_value(1) .and. a%lfacs(2,jfpos+jface) == seek_value(2)) then
               a%lfacs(1,jfpos+jface) = replace_value(1)
               a%lfacs(2,jfpos+jface) = replace_value(2)
               
               call rectree(a,jelem,jface,seek_value,replace_value)

            endif
         enddo
      endif
   end subroutine
   
   subroutine NumberingToParallel(a)
   !This subroutine computes:
   !newpoin0
   !ParallelLocalReordering
   !LocalOrdering
      use Mod_Mesh
      use MPI 
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a

      integer(ip), allocatable :: irenpo(:)
      
      integer(ip), allocatable :: newPoin0(:)
      integer(ip), allocatable :: auxNewPoin0(:)
      
      integer(ip), allocatable :: NGhostsINeed(:), NGhostsOthersNeed(:)
      type(i1p), allocatable :: GhostsINeed(:), GhostsOthersNeed(:)
      type(i2p), allocatable :: GhostsINeed2(:), GhostsOthersNeed2(:)
      type(i1p), allocatable :: auxGhostsINeed(:)
      
      integer(ip), allocatable :: itemp(:)
      integer(ip), allocatable :: itempPedger(:), itempLedger(:,:)
      
      integer(ip) :: newnpoin
      integer(ip) :: newnpoinGhost,newnpoinLocal,newRemainingNpoinGhost
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      integer(ip) :: irank,jrank
      
      !Auxiliary variables
      integer(ip) :: finalPointLocal, finalPointGlobal
      integer(ip) :: intermediatePointLocal
      integer(ip) :: ifinalPointLocal, jfinalPointLocal, jfinalSonLocal
      integer(ip) :: iOldLocalNumbering,iOldGlobalNumbering
      integer(ip) :: jOldLocalNumbering,jOldGlobalNumbering
      
      integer(ip) :: oldiparent,oldjparent,goldiparent,goldjparent
      integer(ip) :: OldLocalParent1,OldLocalParent2 
      integer(ip) :: OldGlobalParent1, OldGlobalParent2
      integer(ip) :: IntermediateParent1, IntermediateParent2
      integer(ip) :: parent1, parent2
      
      integer(ip) :: iaux1,iaux2
      
      integer(ip) :: ielem,ipoin,inode,pnode,ispos,ripoin
      integer(ip) :: jpoin,jv,json,jparent,jspos
      integer(ip) :: tnze
      
      logical :: kfl_isLocal
      integer(ip) :: nedge
      
      logical, allocatable :: NeedToCommunicate(:)
      integer(ip), allocatable :: NewLocalToGlobal(:),NewProcessorList(:)
      
      
      integer(ip) :: ElementCount, GhostCount
      logical     :: isGhostElement
      
      integer(ip), allocatable :: pelpo(:), lelpo(:), ia(:), ja(:)
      integer(ip) :: auxnnode(2)
      
      type(EdgerCriteriaTester) :: ECT
      
      integer(ip), allocatable :: auxParallelLocalReordering(:)
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: icount, newNpoinGhostToExport, ExportNelem, ExportLnodsSizes,jnode,poini,kpoin
      integer(ip) :: finalPointRank
      
      integer(ip) :: elemi
      integer(ip) :: oldlnods(25),goldlnods(25),goldprocs(25)
      
      integer(ip), pointer :: jlnods(:) => NULL()
      integer(ip) :: minnode,ielty,iinteriornode,ninteriornodes,jchild
      integer(ip) :: jelem,jpnode,minpoin
      
      type(LinkedListHeader), allocatable :: HeaderInteriorGhostsINeed(:)
      type(LinkedListStorage)             :: StorageInteriorGhostsINeed
      
      integer(ip), allocatable :: auxInteriorGhostsINeed(:,:)
      integer(ip) :: nInteriorGhostsINeed,minrank
      
      type :: InteriorGhostsIdentifier
         type (GlobalElementIdentifier) :: GEI
         integer(ip)                    :: iinteriornode    
      end type
      type :: InteriorGhostsList
         type(InteriorGhostsIdentifier), allocatable :: a(:)
      end type
      type(InteriorGhostsIdentifier) :: IGI_Sizer
      type(InteriorGhostsList), allocatable :: InteriorGhostsINeed(:), InteriorGhostsOthersNeed(:)
      integer(ip) :: IGI_Size 
      
      integer(ip) :: newnpoinGhostEdger,newnpoinGhostInterior,next,facnod(25),facproc(25),iface,iinthefacenode
      integer(ip) :: kchild,kelem,kpnode,nchild,nchildface,newnpoinGhostInTheFace,nInTheFaceGhostsINeed,nInTheFaceNodes
      integer(ip) :: knode,nface,pnodb
      integer(ip), pointer :: klnods(:) => NULL()
      
      call a%Timers%NumberingToParallel%Tic
      
      a%oldgnpoin = a%gnpoin
      a%oldnpoinLocal = a%npoinLocal
      a%oldnpoinGhost = a%npoinGhost

      newnpoin = a%npoin
      
      !I will need an inverse a%renpo which lets me go from (intermediate) numbering to the old one
      allocate(irenpo(a%npoin))
      irenpo = 0
      do ipoin = 1,a%oldnpoin
         if (a%renpo(ipoin) /= 0) then
            irenpo(a%renpo(ipoin)) = ipoin
         endif
      enddo
      
      !1rst, redo my local numbering (newlocal before old ghosts and new ghosts)
      allocate(a%ParallelLocalReordering(a%npoin))
      a%ParallelLocalReordering = 0
      
      newnpoinLocal = 0
      !Select which are local and which are ghost
      do ipoin = 1,a%oldnpoinLocal
         ripoin = a%renpo(ipoin)
         if (ripoin /= 0) then
            !This is an old local point which still exists, we fill it
            newnpoinLocal = newnpoinLocal+1
            a%ParallelLocalReordering(ripoin) = newnpoinLocal !ripoin will indeed be equal to newnpoinLocal
         endif
      enddo

      !I start filling the old ghosts
      newnpoinGhost = 0
      do ipoin = a%oldnpoinLocal+1,a%oldnpoinLocal+a%oldnpoinGhost
         ripoin = a%renpo(ipoin)
         if (ripoin /= 0) then
            !This is an old ghost point which still exists, we fill it
            newnpoinGhost = newnpoinGhost+1
            a%ParallelLocalReordering(ripoin) = -newnpoinGhost
         endif
      enddo
      newRemainingNpoinGhost = newnpoinGhost
      
      !-----------------------------------------------------------------------------------
      !2nd, obtain the new global numbering for both:
      !locals: easy, once I know the number of locals of all processors
         
      !new locals
            
      !ghosts: this can be difficult.
         !Old ghosts: send them with old global numbering and wait for answer
         
         !New ghosts: These are identified by their parents
            !Send how many I will ask
            !For each I will ask, send both parents in old global numbering
            !The other processor needs to loop through pedger and ledger and find them
      !-----------------------------------------------------------------------------------

      !Loop through the new nodes and decide which are local and which are ghosts
      !First I loop through all the candidate parents
      newnpoinghostEdger = 0
      do ipoin = 1,newnpoin-a%tnewpo
         nedge = a%pedger(ipoin+1)-a%pedger(ipoin)
         if (nedge > 0) then
            oldiparent = irenpo(ipoin)
            call a%LocalOrdering%Local2Global(oldiparent,goldiparent)   
         endif
         do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1 
            jparent = a%ledger(1,jv)
            json = a%ledger(2,jv)
            !Check if the son is new
            if (json > newnpoin - a%tnewpo) then
               oldjparent = irenpo(jparent)
               
               kfl_isLocal = .false.
               if (oldjparent <= a%oldnpoinLocal .and. oldiparent <= a%oldnpoinLocal) then
                  kfl_isLocal = .true.
               else
                  !now I know who my parent was (pre to renumbering 1)
                  call a%LocalOrdering%Local2Global(oldjparent,goldjparent)
                  if (goldiparent < goldjparent) then
                     !if node can belong to 2 processors (1 parent from each processor)
                     !then it belongs to the processor with a lower number
                     if (oldiparent <= a%oldnpoinLocal) then
                        kfl_isLocal = .true.
                     endif
                  elseif (goldjparent < goldiparent) then
                     if (oldjparent <= a%oldnpoinLocal) then
                        kfl_isLocal = .true.
                     endif
                  endif
               endif
               
               if (kfl_isLocal .eqv. .true.) then
                  newnpoinLocal = newnpoinLocal+1
                  a%ParallelLocalReordering(json) = newnpoinLocal
               else
                  newnpoinghostEdger = newnpoinGhostEdger +1
                  newnpoinGhost = newnpoinGhost +1 
                  a%ParallelLocalReordering(json) = -newnpoinGhost
               endif
            endif
         enddo
      enddo
      
      
      !Loop through new elements and do the same for new interior nodes
      newnpoinGhostInterior = 0
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         ielty = a%ElementTypeList(ielem)%ielty
         !Check if there is an interior node
         call Adaptive_GetNumberOfInteriorNodesToAdd(ielty,nInteriorNodes)
         
         if (nInteriorNodes /= 0) then
            !Check if I am responsible for the element
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            oldlnods(1:pnode) = irenpo(lnods)
            call a%LocalOrdering%Local2Global(pnode,oldlnods,goldlnods)
            minnode = minloc(goldlnods(1:pnode),1)
            minpoin = oldlnods(minnode)
            
            do iinteriornode = 1,nInteriorNodes
               call Adaptive_GetInteriorPointFromParent(ielty,iinteriornode,jchild,jnode)
               ispos = a%pPointerToChildren(ielem)-1
               jelem = a%lPointerToChildren(ispos+jchild)
               call a%GetLnodsFromElement(jelem,jpnode,jlnods)
               jpoin = jlnods(jnode)
               
               !If I am responsible the interior node will be mine
               if (minpoin <= a%oldnpoinLocal) then
                  !new point will be local
                  newnpoinLocal = newnpoinLocal+1
                  a%ParallelLocalReordering(jpoin) = newnpoinLocal
               else
                  newnpoinGhost = newnpoinGhost+1
                  newnpoinGhostInterior = newnpoinGhostInterior+1
                  a%ParallelLocalReordering(jpoin) = -newnpoinGhost
                  
               endif   
            enddo
         endif
      enddo
      
      
      !Loop through new elements and do the same for new InTheFace nodes
      newnpoinGhostInTheFace = 0
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         ielty = a%ElementTypeList(ielem)%ielty
         call Adaptive_GetNumberOfInTheFaceNodesToAdd(ielty,nInTheFaceNodes)
         do inode = 1,nInTheFaceNodes
            call Adaptive_GetFaceForInTheFaceNode(ielty,inode,iface)
            call Adaptive_GetInTheFacePointFromParent(ielty,iface,kchild,knode)
            ispos = a%pPointerToChildren(ielem)-1
            kelem = a%lPointerToChildren(ispos+kchild)
            call a%GetLnodsFromElement(kelem,kpnode,klnods)
            kpoin = klnods(knode)
            
            !If the point has still not been classified, then we classify it
            if (a%ParallelLocalReordering(kpoin) == 0) then
               call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
               call a%GetLnodsFromElement(ielem,pnode,lnods)
               call face_from_el(ielty,nface,pnodb,lnods,iface,0_ip,facnod)
               
               call a%LocalOrdering%GetProc(pnodb,irenpo(facnod(1:pnodb)),facproc)
               minrank = minval(facproc(1:pnodb))
               if (minrank == a%MPIrank) then
                  newnpoinLocal = newnpoinLocal+1
                  a%ParallelLocalReordering(kpoin) = newnpoinLocal
               else
                  newnpoinGhost = newnpoinGhost+1
                  newnpoinGhostInTheFace = newnpoinGhostInTheFace+1
                  a%ParallelLocalReordering(kpoin) = -newnpoinGhost
               endif
            endif   
         enddo
      enddo
      
      
 
      !Count locals and ghosts
      !newbornNPoinGhost = newNpoinGhost - newRemainingNpoinGhost;
      a%npoinGhost = newNpoinGhost
      a%npoinLocal = newnpoinLocal
      
      
      !write(*,*) 'Previous loops: A better choice (balanced) for the ownership of the new conflict nodes must be possible WITHOUT communications'
      !write(*,*) 'Be careful because this ownership is assumed in the asking loops below, so if changed everything should be doublechecked'
      
      !Modify ghosts
      do ipoin = 1,newnpoin
         !ghosts are numbered negatively
         if (a%ParallelLocalReordering(ipoin) < 0) then
            a%ParallelLocalReordering(ipoin) = newnpoinLocal + (-a%ParallelLocalReordering(ipoin)) !this is an addition, ghosts are numbered negatively
         endif
      enddo
      
      !Next step is, build the new globalNumbering for my ghostPoints
      !First everybody says how many local nodes they have
      allocate(newPoin0(a%MPIsize))
      allocate(auxnewPoin0(a%MPIsize))
      auxnewPoin0(1:a%MPIsize) = newnpoinLocal
      call MPI_AllToAll(auxnewPoin0, 1, MPI_INTEGER4, newPoin0, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      deallocate(auxnewPoin0)
      
      !GlobalNpoin
      a%gnpoin = sum(newPoin0)
      !New Poin0 for each processor
      iaux1 = 0
      do irank = 0,a%MPIsize-1
         iaux2 = newPoin0(irank+1)
         newPoin0(irank+1) = iaux1
         iaux1 = iaux1+iaux2
      enddo
      a%poin0 = newPoin0(a%MPIrank+1)
      deallocate(newPoin0)
      
      !NewLocalToGlobal and Processor LIst
      allocate(NewLocalToGlobal(a%npoin))
      allocate(NewProcessorList(a%npoin))
      NewLocalToGlobal = 0
      
      !First we fill the local nodes
      do ipoin = 1,newnpoinLocal
         NewLocalToGlobal(ipoin) = ipoin + a%poin0
         NewProcessorList(ipoin) = a%MPIrank
      enddo
            
      !Now I need to ask the neighbour for the new numbering of my ghost points
      !There are several types of ghosts:
      !Old Ghosts, of which I already know they global numbering (a%npoinLocal+1:a%npoinLocal+newRemainingNpoinGhost)
      ! and new Edger ghosts, of which I only know who their parents are (a%npoinLocal+newRemainingNpoinGhost+1:a%npoinLocal+newRemainingNpoinGhost+newnpoinGhostEdger)
      ! new interior ghosts, of which I only know element, son node (a%npoinLocal+newRemainingNpoinGhost+newnpoinGhostEdger+1: a%npoinLocal+newRemainingNpoinGhost+newnpoinGhostEdger+newnpoinGhostInterior)
      
      !Let us start with the old ghosts
      !I will only need to communicate with those processors with which I have a ghost.
      !They will also have ghosts from me, so they will know that I need to communicate with them
      !Let us count how many of them I am connected to
      
      allocate(NGhostsINeed(a%MPIsize))
      NGhostsINeed = 0
      allocate(NGhostsOthersNeed(a%MPIsize))
      NGhostsOthersNeed = 0
      
      do ipoin = a%oldnpoinLocal+1,a%oldnpoin
         !If the old ghost point continues to exist
         if (a%renpo(ipoin) /= 0) then
            call a%LocalOrdering%GetProc(ipoin,irank)
            NGhostsINeed(irank+1) = NGhostsINeed(irank+1)+1
         endif
      enddo
      call MPI_AllToAll(NGhostsINeed, 1, MPI_INTEGER4, NGhostsOthersNeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      !I need to store this, because in the new ghost nodes information exchange might not be symmetric
      allocate(NeedToCommunicate(a%MPIsize))
      NeedToCommunicate = .false.
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0 .or. NGhostsOthersNeed(irank+1) /= 0) NeedToCommunicate(irank+1) = .true.
      enddo
      
      !Store the ghostsINeed
      allocate(GhostsINeed(a%MPIsize))
      allocate(GhostsOthersNeed(a%MPIsize))
      
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then
            allocate(GhostsINeed(irank+1)%l(NGhostsINeed(irank+1)))
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            allocate(GhostsOthersNeed(irank+1)%l(NGhostsOthersNeed(irank+1)))
         endif
      enddo
      
      !Fill the ghosts I need
      NGhostsINeed = 0
      do ipoin = a%oldnpoinLocal+1,a%oldnpoin
         !If the old ghost point continues to exist
         if (a%renpo(ipoin) /= 0) then
            call a%LocalOrdering%GetProc(ipoin,irank)
            NGhostsINeed(irank+1) = NGhostsINeed(irank+1)+1
            call a%LocalOrdering%Local2Global(ipoin,GhostsINeed(irank+1)%l(NGhostsINeed(irank+1)))
         endif
      enddo
      
      
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      
      !Send and receive the ghostsINeed /GhostsOthersNeed
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then
            call MPI_ISEND(GhostsINeed(irank+1)%l,NGhostsINeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_IRECV(GhostsOthersNeed(irank+1)%l,NGhostsOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
         endif
      enddo
      !Wait for communications to finish
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then
            call MPI_WAIT(irequest01(irank+1), status, ierr)
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_WAIT(irequest02(irank+1), status, ierr)
         endif
      enddo
      
      !Transform from Old Global Numbering to NewLocalToGlobal
      !It is not sufficient to send the point number, it is also necessary to send the processor number
      !This is due to load rebalancing: a child of two local nodes does not need to be local anymore
      do irank = 0,a%MPIsize-1
         do ipoin = 1,NGhostsOthersNeed(irank+1)
            iOldGlobalNumbering = GhostsOthersNeed(irank+1)%l(ipoin)
            call a%LocalOrdering%Global2Local(iOldGlobalNumbering,iOldLocalNumbering)
            intermediatePointLocal = a%renpo(iOldLocalNumbering)
            if (intermediatePointLocal == 0) then
               write(*,*) 'problemo'
            endif
            finalPointLocal = a%ParallelLocalReordering(intermediatePointLocal)
            finalPointGlobal = finalPointLocal + a%poin0
            GhostsOthersNeed(irank+1)%l(ipoin) = finalPointGlobal
         enddo
      enddo
      
      !send/receive the new global numbering for the old ghost points
      !Send and receive the ghostsINeed /GhostsOthersNeed
      do irank = 0,a%MPIsize-1
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(GhostsOthersNeed(irank+1)%l,NGhostsOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NGhostsINeed(irank+1) /= 0) then
            call MPI_IRECV(GhostsINeed(irank+1)%l,NGhostsINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
         endif
      enddo
      !Wait for communications to finish
      do irank = 0,a%MPIsize-1
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_WAIT(irequest01(irank+1), status, ierr)
         endif
         if (NGhostsINeed(irank+1) /= 0) then
            call MPI_WAIT(irequest02(irank+1), status, ierr)
         endif
      enddo
      
      !Now I can fill the new global numbering
      NGhostsINeed = 0
      do ipoin = a%oldnpoinLocal+1,a%oldnpoin
         !If the old ghost point continues to exist
         if (a%renpo(ipoin) /= 0) then
            call a%LocalOrdering%GetProc(ipoin,irank)
            NGhostsINeed(irank+1) = NGhostsINeed(irank+1)+1
            
            finalPointGlobal = GhostsINeed(irank+1)%l(NGhostsINeed(irank+1))
            finalPointLocal = a%ParallelLocalReordering(a%renpo(ipoin))
            NewLocalToGlobal(finalPointLocal) = finalPointGlobal
            NewProcessorList(finalPointLocal) = irank
         endif
      enddo
      
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then
            deallocate(GhostsINeed(irank+1)%l)
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            deallocate(GhostsOthersNeed(irank+1)%l)
         endif
      enddo
      deallocate(GhostsINeed,GhostsOthersNeed)
      
      !---------------------------------------------------------------------------------------------
      !Second part, the new  edger ghosts, the only thing I know about them is who their parents are
      !Again, I first count how many of each processor I need
      !The information about them is in ledger, numbering is intermediate
      NGhostsINeed = 0
      do ipoin = 1,newnpoin
         iOldLocalNumbering = irenpo(ipoin)
         !ipoin did already exist
         if (iOldLocalNumbering /= 0) then
            call a%LocalOrdering%Local2Global(iOldLocalNumbering,iOldGlobalNumbering)
            call a%LocalOrdering%GetProc(iOldLocalNumbering,irank)
            if (irank == a%MPIrank) irank = a%MPIsize+1
            
            do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
               jpoin = a%ledger(1,jv)
               json = a%ledger(2,jv)
               !If it is a newborn ghost
               !We identify it because these nodes remain to be numbered
               if (NewLocalToGlobal(a%ParallelLocalReordering(json)) == 0) then
                  !Now I need to decide from whom I need it
                  !Both my parents might be ghosts, I need to select the one with the lowest global number to ask its processor
                  
                  jOldLocalNumbering = irenpo(jpoin)
                  call a%LocalOrdering%Local2Global(jOldLocalNumbering,jOldGlobalNumbering)
                  call a%LocalOrdering%GetProc(jOldLocalNumbering,jrank)
                  if (jrank == a%MPIrank) jrank = a%MPIsize+1
                  
                  jrank = min(irank,jrank)
                  !-------------------------------------------------------------------------------------------------------
                  !Add one to the counter for this processor
                  NGhostsINeed(jrank+1) = NGhostsINeed(jrank+1)+1
                  !-------------------------------------------------------------------------------------------------------
               endif
            enddo
         endif
      enddo
      
      !Allocates
      allocate(GhostsINeed2(a%MPIsize),auxGhostsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(GhostsINeed2(irank+1)%l(2,NGhostsINeed(irank+1)))
         allocate(auxGhostsINeed(irank+1)%l(NGhostsINeed(irank+1)))
      enddo
      
      
      !Second to fill
      NGhostsINeed = 0
      do ipoin = 1,newnpoin
         iOldLocalNumbering = irenpo(ipoin)
         !ipoin did already exist
         if (iOldLocalNumbering /= 0) then
            call a%LocalOrdering%Local2Global(iOldLocalNumbering,iOldGlobalNumbering)
            call a%LocalOrdering%GetProc(iOldLocalNumbering,irank)
            if (irank == a%MPIrank) irank = a%MPIsize+1
            
            do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
               jpoin = a%ledger(1,jv)
               json = a%ledger(2,jv)
               !If it is a newborn ghost
               !We identify them because these nodes remain to be numbered
               if (NewLocalToGlobal(a%ParallelLocalReordering(json)) == 0) then
                  !Now I need to decide from whom I need it
                  !Both my parents might be ghosts, I need to select the one with the lowest global number to ask its processor
                  
                  jOldLocalNumbering = irenpo(jpoin)
                  call a%LocalOrdering%Local2Global(jOldLocalNumbering,jOldGlobalNumbering)
                  call a%LocalOrdering%GetProc(jOldLocalNumbering,jrank)
                  if (jrank == a%MPIrank) jrank = a%MPIsize+1
                  
                  jrank = min(irank,jrank)
                  
                  !-------------------------------------------------------------------------------------------------------
                  !Fill the list to ask for
                  NGhostsINeed(jrank+1) = NGhostsINeed(jrank+1)+1
                  GhostsINeed2(jrank+1)%l(1,NGhostsINeed(jrank+1)) = iOldGlobalNumbering
                  GhostsINeed2(jrank+1)%l(2,NGhostsINeed(jrank+1)) = jOldGlobalNumbering
                  
                  auxGhostsINeed(jrank+1)%l(NGhostsINeed(jrank+1)) = a%ParallelLocalReordering(json)
                  !-------------------------------------------------------------------------------------------------------
               endif
            enddo
         endif
      enddo
            
      !Send sizes      
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      !Send and receive the number of ghostsINeed /GhostsOthersNeed
      do irank = 0,a%MPIsize-1
         if (NeedToCommunicate(irank+1) .eqv. .true.) then   
            call MPI_ISEND(NghostsINeed(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_IRECV(NGhostsOthersNeed(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
         endif
      enddo
      !Wait for communications to finish
      do irank = 0,a%MPIsize-1
         if (NeedToCommunicate(irank+1) .eqv. .true.) then
            call MPI_WAIT(irequest01(irank+1), status, ierr)
            call MPI_WAIT(irequest02(irank+1), status, ierr)
         endif
      enddo
      
      allocate(GhostsOthersNeed2(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(GhostsOthersNeed2(irank+1)%l(2,NGhostsOthersNeed(irank+1)))
      enddo
      
      !Send and receive the ghostsINeed /GhostsOthersNeed
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then
            call MPI_ISEND(GhostsINeed2(irank+1)%l,2*NGhostsINeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_IRECV(GhostsOthersNeed2(irank+1)%l,2*NGhostsOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
         endif
      enddo
      !Wait for communications to finish
      do irank = 0,a%MPIsize-1
         if (NGhostsINeed(irank+1) /= 0) then   
            call MPI_WAIT(irequest01(irank+1), status, ierr)
         endif
         if (NGhostsOthersNeed(irank+1) /= 0) then
            call MPI_WAIT(irequest02(irank+1), status, ierr)
         endif
      enddo
      
      !Now I need to find which the global number of each son is
      !First I translate the oldGlobalNumbering to intermediate numbering (so that I can use edger)
      do irank = 0,a%MPIsize-1
         do ipoin = 1,NGhostsOthersNeed(irank+1)
            OldGlobalParent1 = GhostsOthersNeed2(irank+1)%l(1,ipoin)
            OldGlobalParent2 = GhostsOthersNeed2(irank+1)%l(2,ipoin)
            
            call a%LocalOrdering%Global2Local(OldGlobalParent1,OldLocalParent1)
            call a%LocalOrdering%Global2Local(OldGlobalParent2,OldLocalParent2)
            
            IntermediateParent1 = a%renpo(OldLocalParent1)
            IntermediateParent2 = a%renpo(OldLocalParent2)
            
            !Who is the parent?
            ECT%parents(1) = IntermediateParent1
            ECT%parents(2) = IntermediateParent2
            ECT%levels(1) = a%level(IntermediateParent1)
            ECT%levels(2) = a%level(IntermediateParent2)
            call EdgerCriteria(ECT)
                  
            parent1 = ECT%orderedparents(1)
            parent2 = ECT%orderedparents(2)
            
            pedgerloop: do jv = a%pedger(parent1),a%pedger(parent1+1)-1
               jpoin = a%ledger(1,jv)
               if (jpoin == parent2) then
                  json = a%ledger(2,jv)
                  finalPointLocal = a%ParallelLocalReordering(json)
                  finalPointGlobal = NewLocalToGlobal(finalPointLocal)
                  finalPointRank   = NewProcessorList(finalPointLocal)
                  GhostsOthersNeed2(irank+1)%l(1,ipoin) = finalPointGlobal
                  GhostsOthersNeed2(irank+1)%l(2,ipoin) = finalPointRank
                  exit pedgerloop
               endif
            enddo pedgerloop
         enddo
      enddo
      
      !This is repeated so I put it in a subroutine
      call SendAndReceiveGhostsWeNeed2_AndDeallocate
      

      
    
      
      !--------------------------------------------------------------------
      !Third part, the new interior ghosts, I know the parent element and 
      !its interior node number
      allocate(HeaderInteriorGhostsINeed(a%MPIsize))
      call StorageInteriorGhostsINeed%Initialize(newnpoinGhostInterior,'NonRemovable','Repeatable',a%Memor)
      allocate(auxInteriorGhostsINeed(2,newnpoinGhostInterior))
      
      allocate(GhostsINeed2(a%MPIsize),GhostsOthersNeed2(a%MPIsize),auxGhostsINeed(a%MPIsize))
      
      nInteriorGhostsINeed = 0
      NghostsINeed = 0
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         ielty = a%ElementTypeList(ielem)%ielty
         !Check if there is an interior node
         call Adaptive_GetNumberOfInteriorNodesToAdd(ielty,nInteriorNodes)
         
         if (nInteriorNodes /= 0) then
            !Check if I am responsible for the element
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            oldlnods(1:pnode) = irenpo(lnods)
            call a%LocalOrdering%Local2Global(pnode,oldlnods,goldlnods)
            call a%LocalOrdering%GetProc(pnode,oldlnods,goldprocs)
            minnode = minloc(goldlnods(1:pnode),1)
            minpoin = oldlnods(minnode)
            minrank = goldprocs(minnode)
            
            do iinteriornode = 1,nInteriorNodes
               call Adaptive_GetInteriorPointFromParent(ielty,iinteriornode,jchild,jnode)
               ispos = a%pPointerToChildren(ielem)-1
               jelem = a%lPointerToChildren(ispos+jchild)
               call a%GetLnodsFromElement(jelem,jpnode,jlnods)
               jpoin = jlnods(jnode)
               
               if (NewLocalToGlobal(a%ParallelLocalReordering(jpoin)) == 0) then
                  nInteriorGhostsINeed = nInteriorGhostsINeed + 1
                  auxInteriorGhostsINeed(1,nInteriorGhostsINeed) = ielem
                  auxInteriorGhostsINeed(2,nInteriorGhostsINeed) = iinteriornode
                  call StorageInteriorGhostsINeed%Add(HeaderInteriorGhostsINeed(minrank+1),nInteriorGhostsINeed)
                  nGhostsINeed(minrank+1) = nGhostsINeed(minrank+1)+1
               endif   
            enddo
         endif
      enddo
      
      allocate(InteriorGhostsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(auxGhostsINeed(irank+1)%l(nGhostsINeed(irank+1)))
         allocate(GhostsINeed2(irank+1)%l(2,nGhostsINeed(irank+1)))
         allocate(InteriorGhostsINeed(irank+1)%a(nGhostsINeed(irank+1)))
         call StorageInteriorGhostsINeed%GoToFirstOfList(HeaderInteriorGhostsINeed(irank+1))
         call StorageInteriorGhostsINeed%GetNext(next)
         icount = 0
         do while (next /= -1)
            icount = icount + 1 
            ielem = auxInteriorGhostsINeed(1,next) 
            iinteriornode = auxInteriorGhostsINeed(2,next) 
            InteriorGhostsINeed(irank+1)%a(icount)%GEI = a%GEI_ElementList(ielem)
            InteriorGhostsINeed(irank+1)%a(icount)%iinteriornode = iinteriornode
            
            !Store the point number here so that we can use it when we get the results back
            ielty = a%ElementTypeList(ielem)%ielty
            call Adaptive_GetInteriorPointFromParent(ielty,iinteriornode,jchild,jnode)
            ispos = a%pPointerToChildren(ielem)-1
            jelem = a%lPointerToChildren(ispos+jchild)
            call a%GetLnodsFromElement(jelem,jpnode,jlnods)
            jpoin = jlnods(jnode)
            auxGhostsINeed(irank+1)%l(icount) = a%ParallelLocalReordering(jpoin)
            
            
            call StorageInteriorGhostsINeed%GetNext(next)
         enddo
      enddo
      
      call SendAndReceiveInteriorAndFaceNodes
      
      do irank = 0,a%MPIsize-1
         do ipoin = 1,NGhostsOthersNeed(irank+1)
            call a%LocalElementFromGEI(InteriorGhostsOthersNeed(irank+1)%a(ipoin)%GEI,ielem)
            iinteriornode = InteriorGhostsOthersNeed(irank+1)%a(ipoin)%iinteriornode
            
            ielty = a%ElementTypeList(ielem)%ielty
            call Adaptive_GetInteriorPointFromParent(ielty,iinteriornode,jchild,jnode)
            ispos = a%pPointerToChildren(ielem)-1
            jelem = a%lPointerToChildren(ispos+jchild)
            call a%GetLnodsFromElement(jelem,jpnode,jlnods)
            jpoin = jlnods(jnode)
            
            finalPointLocal = a%ParallelLocalReordering(jpoin)
            finalPointGlobal = NewLocalToGlobal(finalPointLocal)
            finalPointRank   = NewProcessorList(finalPointLocal)
            GhostsOthersNeed2(irank+1)%l(1,ipoin) = finalPointGlobal
            GhostsOthersNeed2(irank+1)%l(2,ipoin) = finalPointRank
         enddo
      enddo
      
      call SendAndReceiveGhostsWeNeed2_AndDeallocate
      
      !More deallocates
      deallocate(HeaderInteriorGhostsINeed)
      call StorageInteriorGhostsINeed%Dealloc
      deallocate(auxInteriorGhostsINeed)
      do irank = 0,a%MPIsize-1
         deallocate(InteriorGhostsINeed(irank+1)%a)
         deallocate(InteriorGhostsOthersNeed(irank+1)%a)
      enddo
      deallocate(InteriorGhostsINeed,InteriorGhostsOthersNeed)
      
      
      !--------------------------------------------------------------------------------
      !Fourth part, the new in the face nodes
      !Loop through new elements and do the same for new InTheFace nodes
      !I know the parent and the intheface node
      
      allocate(HeaderInteriorGhostsINeed(a%MPIsize))
      call StorageInteriorGhostsINeed%Initialize(newnpoinGhostInTheFace,'NonRemovable','Repeatable',a%Memor)
      allocate(auxInteriorGhostsINeed(2,newnpoinGhostInTheFace))
      
      allocate(GhostsINeed2(a%MPIsize),GhostsOthersNeed2(a%MPIsize),auxGhostsINeed(a%MPIsize))
      
      nInTheFaceGhostsINeed = 0
      NGhostsINeed = 0
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         ielty = a%ElementTypeList(ielem)%ielty
         call Adaptive_GetNumberOfInTheFaceNodesToAdd(ielty,nInTheFaceNodes)
         do inode = 1,nInTheFaceNodes
            call Adaptive_GetFaceForInTheFaceNode(ielty,inode,iface)
            call Adaptive_GetInTheFacePointFromParent(ielty,iface,kchild,knode)
            ispos = a%pPointerToChildren(ielem)-1
            kelem = a%lPointerToChildren(ispos+kchild)
            call a%GetLnodsFromElement(kelem,kpnode,klnods)
            kpoin = klnods(knode)
            
            !If the point has still not been classified, then we classify it
            if (NewLocalToGlobal(a%ParallelLocalReordering(kpoin)) == 0) then
               !This is to make sure that we do the node only once
               NewLocalToGlobal(a%ParallelLocalReordering(kpoin)) = -1
               
               call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
               call a%GetLnodsFromElement(ielem,pnode,lnods)
               call face_from_el(ielty,nface,pnodb,lnods,iface,0_ip,facnod)
               
               call a%LocalOrdering%GetProc(pnodb,irenpo(facnod(1:pnodb)),facproc)
               minrank = minval(facproc(1:pnodb))
               
               nInTheFaceGhostsINeed = nInTheFaceGhostsINeed + 1
               auxInteriorGhostsINeed(1,nInTheFaceGhostsINeed) = ielem
               auxInteriorGhostsINeed(2,nInTheFaceGhostsINeed) = inode
               call StorageInteriorGhostsINeed%Add(HeaderInteriorGhostsINeed(minrank+1),nInTheFaceGhostsINeed)
               nGhostsINeed(minrank+1) = nGhostsINeed(minrank+1)+1
            endif   
         enddo
      enddo
      
      
      
      allocate(InteriorGhostsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(auxGhostsINeed(irank+1)%l(nGhostsINeed(irank+1)))
         allocate(GhostsINeed2(irank+1)%l(2,nGhostsINeed(irank+1)))
         allocate(InteriorGhostsINeed(irank+1)%a(nGhostsINeed(irank+1)))
         call StorageInteriorGhostsINeed%GoToFirstOfList(HeaderInteriorGhostsINeed(irank+1))
         call StorageInteriorGhostsINeed%GetNext(next)
         icount = 0
         do while (next /= -1)
            icount = icount + 1 
            ielem = auxInteriorGhostsINeed(1,next) 
            iinthefacenode = auxInteriorGhostsINeed(2,next) 
            InteriorGhostsINeed(irank+1)%a(icount)%GEI = a%GEI_ElementList(ielem)
            InteriorGhostsINeed(irank+1)%a(icount)%iinteriornode = iinthefacenode
            
            !Store the point number here so that we can use it when we get the results back
            ielty = a%ElementTypeList(ielem)%ielty
            call Adaptive_GetFaceForInTheFaceNode(ielty,iinthefacenode,iface)
            call Adaptive_GetInTheFacePointFromParent(ielty,iface,kchild,knode)
            ispos = a%pPointerToChildren(ielem)-1
            kelem = a%lPointerToChildren(ispos+kchild)
            call a%GetLnodsFromElement(kelem,kpnode,klnods)
            kpoin = klnods(knode)
            auxGhostsINeed(irank+1)%l(icount) = a%ParallelLocalReordering(kpoin)
            
            
            call StorageInteriorGhostsINeed%GetNext(next)
         enddo
      enddo
      
      call SendAndReceiveInteriorAndFaceNodes
      
      do irank = 0,a%MPIsize-1
         do ipoin = 1,NGhostsOthersNeed(irank+1)
            call a%LocalElementFromGEI(InteriorGhostsOthersNeed(irank+1)%a(ipoin)%GEI,ielem)
            iinthefacenode = InteriorGhostsOthersNeed(irank+1)%a(ipoin)%iinteriornode
            
            ielty = a%ElementTypeList(ielem)%ielty
            call Adaptive_GetFaceForInTheFaceNode(ielty,iinthefacenode,iface)
            call Adaptive_GetInTheFacePointFromParent(ielty,iface,kchild,knode)
            ispos = a%pPointerToChildren(ielem)-1
            kelem = a%lPointerToChildren(ispos+kchild)
            call a%GetLnodsFromElement(kelem,kpnode,klnods)
            kpoin = klnods(knode)
            
            finalPointLocal = a%ParallelLocalReordering(kpoin)
            finalPointGlobal = NewLocalToGlobal(finalPointLocal)
            finalPointRank   = NewProcessorList(finalPointLocal)
            GhostsOthersNeed2(irank+1)%l(1,ipoin) = finalPointGlobal
            GhostsOthersNeed2(irank+1)%l(2,ipoin) = finalPointRank
         enddo
      enddo
      
      call SendAndReceiveGhostsWeNeed2_AndDeallocate
      
      deallocate(HeaderInteriorGhostsINeed)
      call StorageInteriorGhostsINeed%Dealloc
      
      
      !Deallocates
      deallocate(NeedToCommunicate)
      deallocate(NGhostsINeed,NGhostsOthersNeed)
      
      
      !We reinitialize the LocalOrdering
      if (associated(a%OldLocalOrdering)) then 
         call a%OldLocalOrdering%Deallocate(a%Memor)
         deallocate(a%OldLocalOrdering)
      endif
      a%OldLocalOrdering => a%LocalOrdering
      allocate(a%LocalOrdering)
      
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,NewLocalToGlobal,NewProcessorList,a%Memor)
      deallocate(NewLocalToGlobal)
      deallocate(NewProcessorList)
      
      !We reinitialize the communicator
      call a%Communicator%Deallocate
      
      !this needs to be seriously improved (even the required data (ia, ja) for building the communicator
      call a%Communicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,a%npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)
      
      !We also need a local to global ordering for top level elements
      call a%ElementLocalOrdering%Deallocate(a%Memor)
      call BuildElementLocalOrdering(a%MPIcomm,a%nelem,a%GEI_ElementList,a%ElementLocalOrdering,a%Memor)

      call a%Timers%NumberingToParallel%Toc
contains


      subroutine SendAndReceiveGhostsWeNeed2_AndDeallocate
         implicit none
         !Send and receive the ghostsINeed /GhostsOthersNeed
         do irank = 0,a%MPIsize-1
            if (NghostsOthersNeed(irank+1) /= 0) then
               call MPI_ISEND(GhostsOthersNeed2(irank+1)%l,2*NGhostsOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            endif
            if (NGhostsINeed(irank+1) /= 0) then
               call MPI_IRECV(GhostsINeed2(irank+1)%l,2*NghostsINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
            endif
         enddo
         !Wait for communications to finish
         do irank = 0,a%MPIsize-1
            if (NGhostsOthersNeed(irank+1) /= 0) then
               call MPI_WAIT(irequest01(irank+1), status, ierr)
            endif
            if (NGhostsINeed(irank+1) /= 0) then
               call MPI_WAIT(irequest02(irank+1), status, ierr)
            endif
         enddo
         
         !Finally I store the new node numbering, and processor
         do irank = 0,a%MPIsize-1
            do ipoin = 1,NGhostsINeed(irank+1)
               finalPointLocal = auxGhostsINeed(irank+1)%l(ipoin)
               finalPointGlobal = GhostsINeed2(irank+1)%l(1,ipoin)
               finalPointRank   = GhostsINeed2(irank+1)%l(2,ipoin)
               NewLocalToGlobal(finalPointLocal) = finalPointGlobal
               NewProcessorList(finalPointLocal) = finalPointRank
            enddo
         enddo
         
         !Deallocates
         do irank = 0,a%MPIsize-1
            deallocate( GhostsINeed2(irank+1)%l, auxGhostsINeed(irank+1)%l)
            deallocate( GhostsOthersNeed2(irank+1)%l)
         enddo
         deallocate(GhostsINeed2,auxGhostsINeed)
         deallocate(GhostsOthersNeed2)
      end subroutine
      
      
      subroutine SendAndReceiveInteriorAndFaceNodes
         implicit none
            
         !Send sizes      
         irequest01 = MPI_REQUEST_NULL
         irequest02 = MPI_REQUEST_NULL
         !Send and receive the number of ghostsINeed /GhostsOthersNeed
         do irank = 0,a%MPIsize-1
            if (NeedToCommunicate(irank+1) .eqv. .true.) then   
               call MPI_ISEND(NghostsINeed(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
               call MPI_IRECV(NGhostsOthersNeed(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
            endif
         enddo
         !Wait for communications to finish
         do irank = 0,a%MPIsize-1
            if (NeedToCommunicate(irank+1) .eqv. .true.) then
               call MPI_WAIT(irequest01(irank+1), status, ierr)
               call MPI_WAIT(irequest02(irank+1), status, ierr)
            endif
         enddo
         
         allocate(InteriorGhostsOthersNeed(a%MPIsize))
         do irank = 0,a%MPIsize-1
            allocate(InteriorGhostsOthersNeed(irank+1)%a(nGhostsOthersNeed(irank+1)))
            allocate(GhostsOthersNeed2(irank+1)%l(2,nGhostsOthersNeed(irank+1)))
         enddo
         
         !Send and receive the ghostsINeed /GhostsOthersNeed
         IGI_Size = sizeof(IGI_Sizer)
         do irank = 0,a%MPIsize-1
            if (NGhostsINeed(irank+1) /= 0) then
               call MPI_ISEND(InteriorGhostsINeed(irank+1)%a,IGI_Size*NGhostsINeed(irank+1),MPI_BYTE, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            endif
            if (NGhostsOthersNeed(irank+1) /= 0) then
               call MPI_IRECV(InteriorGhostsOthersNeed(irank+1)%a,IGI_Size*NGhostsOthersNeed(irank+1),MPI_BYTE,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr) 
            endif
         enddo
         !Wait for communications to finish
         do irank = 0,a%MPIsize-1
            if (NGhostsINeed(irank+1) /= 0) then   
               call MPI_WAIT(irequest01(irank+1), status, ierr)
            endif
            if (NGhostsOthersNeed(irank+1) /= 0) then
               call MPI_WAIT(irequest02(irank+1), status, ierr)
            endif
         enddo
      end subroutine
      
      
   end subroutine
   
   subroutine BuildElementLocalOrdering(MPIcomm,nelem,GEI_ElementList,ElementLocalOrdering,Memor)
      use typre
      use Mod_ParallelOrderingInterface
      implicit none
      integer(ip) :: nelem,MPIcomm
      type(GlobalElementIdentifier) :: GEI_ElementList(nelem)
      class(ParallelOrderingInterface) :: ElementLocalOrdering
      type(MemoryMan) :: Memor
      
      integer(ip), allocatable :: auxLocal2GlobalElementList(:)
      integer(ip) :: TopLevelNelem,ielem
      
      !We need to fill MarkelOthersNeed
      !In order to look for them, I need some kind of local to global of topwise elements (level 0) elements
      allocate(auxLocal2GlobalElementList(nelem))
      auxLocal2GlobalElementList = 0
      TopLevelNelem = 0
      do ielem = 1,nelem
         if (GEI_ElementList(ielem)%level == 0) then
            auxLocal2GlobalElementList(ielem) = GEI_ElementList(ielem)%GlobalOriginalElement
            TopLevelNelem = TopLevelNelem+1
         endif
      enddo
      !The second argument is dummy (ProcessorList not required)
      call ElementLocalOrdering%Init(MPIcomm,TopLevelNelem,auxLocal2GlobalElementList,auxLocal2GlobalElementList,Memor)
      
      deallocate(auxLocal2GlobalElementList)
      
      
   
   end subroutine
   
   subroutine UpdateArraysToNewParallelNumbering(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip), allocatable :: itemp(:)
      integer(ip) :: IntermediatePointLocal,finalPointLocal
      integer(ip), allocatable :: itempPedger(:), itempLedger(:,:)
      
      integer(ip) :: iaux1,iaux2,ielem
      
      integer(ip) :: ifinalPointLocal,inode,ipoin,jfinalPointLocal,ispos,jfinalSonLocal,jpoin,json,jv,pnode,tnze
      integer(ip) :: ientryPoint
      
      type(EdgerCriteriaTester) :: ECT
      integer(ip) :: newnpoinTop
      
      call a%Timers%UpdateArraysToNewParallelNumbering%Tic
!-------------------------------------------------------------------------------      
      !After this point I already have the new numbering  
      !I need to modify:
      !pnods and lnods
      !EdgerEntryPoints
      !tedger, pedger
      !a%renpo
      !touched
     
      !Also in the mesh lots of things need to be changed and adapted
      !I need to be able to go from the prerefinement to the after parallel setting
      
      !we start with pnods, lnods in the mesh refiner
      do ielem = 1,a%nelem
         ispos = a%pnods(ielem)-1
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         do inode = 1,pnode
            IntermediatePointLocal = a%lnods(ispos+inode)
            finalPointLocal = a%ParallelLocalReordering(IntermediatePointLocal)
            a%lnods(ispos+inode) = finalPointLocal
         enddo
      enddo
      
      !level
      allocate(itemp(a%npoin))
      do ipoin = 1,a%npoin
         finalPointLocal = a%ParallelLocalReordering(ipoin)
         itemp(finalPointLocal) = a%level(ipoin)
      enddo
      call move_alloc(itemp,a%level)
      
      !Pedger and Ledger
      tnze = a%pedger(a%npoin+1)-1
      allocate(itempPedger(a%npoin+1),itempLedger(2,tnze))
      itempPedger = 0
      !first we count how many of them for each point
      do ipoin = 1,a%npoin
         ifinalPointLocal = a%ParallelLocalReordering(ipoin)
         do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
            jpoin = a%ledger(1,jv)
            jfinalPointLocal = a%ParallelLocalReordering(jpoin)
            
            !Who is the parent?
            ECT%parents(1) = ifinalPointLocal
            ECT%parents(2) = jfinalPointLocal
            ECT%levels(1) = a%level(ifinalPointLocal)
            ECT%levels(2) = a%level(jfinalPointLocal)
            call EdgerCriteria(ECT)
            
            
            if (ifinalPointLocal == ECT%orderedParents(1)) then
               itempPedger(ifinalPointLocal) = itempPedger(ifinalPointLocal)+1
            else
               itempPedger(jfinalPointLocal) = itempPedger(jfinalPointLocal)+1
            endif
         enddo
      enddo
      !We compress itempPedger
      iaux1 = 1
      do ipoin = 1,a%npoin
         iaux2 = itempPedger(ipoin)
         itempPedger(ipoin) = iaux1
         iaux1 = iaux1+iaux2
      enddo
      itempPedger(a%npoin+1) = iaux1
      !Second to fill
      do ipoin = 1,a%npoin
         ifinalPointLocal = a%ParallelLocalReordering(ipoin)
         do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
            jpoin = a%ledger(1,jv)
            jfinalPointLocal = a%ParallelLocalReordering(jpoin)
            
            json = a%ledger(2,jv)
            jfinalSonLocal = a%ParallelLocalReordering(json)
            
            !Who is the parent?
            ECT%parents(1) = ifinalPointLocal
            ECT%parents(2) = jfinalPointLocal
            ECT%levels(1) = a%level(ifinalPointLocal)
            ECT%levels(2) = a%level(jfinalPointLocal)
            call EdgerCriteria(ECT)
            
            if (ifinalPointLocal == ECT%orderedParents(1)) then
               ispos = itempPedger(ifinalPointLocal)
               
               itempLedger(1,ispos) = jfinalPointLocal
               itempLedger(2,ispos) = jfinalSonLocal
               
               itempPedger(ifinalPointLocal) = itempPedger(ifinalPointLocal)+1
            else
               ispos = itempPedger(jfinalPointLocal)
               
               itempLedger(1,ispos) = ifinalPointLocal
               itempLedger(2,ispos) = jfinalSonLocal
            
               itempPedger(jfinalPointLocal) = itempPedger(jfinalPointLocal)+1
            endif
         enddo
      enddo
      !Fix itempPedger
      do ipoin = a%npoin,2,-1
         itempPedger(ipoin) = itempPedger(ipoin-1)
      enddo
      itempPedger(1) = 1
      !Moveallocs
      call move_alloc(itempPedger,a%pedger)
      call move_alloc(itempLedger,a%ledger)
      
!       !EdgerEntryPoints
!       !Some of the old entry points might be dead
!       !be careful
!       newnpoinTop = 0
!       do ipoin = 1,a%npoinTop
!          ientryPoint = a%EdgerEntryPoints(ipoin)
!          intermediatePointLocal = a%renpo(ientryPoint)
!          if (intermediatePointLocal /= 0) then
!             newnpoinTop = newnpoinTop+1
!             ifinalPointLocal = a%ParallelLocalReordering(intermediatePointLocal)
!             a%EdgerEntryPoints(newnpoinTop) = ifinalPointLocal
!          endif
!       enddo
!       if (newnpoinTop /= a%npoinTop) then
!          a%npoinTop = newnpoinTop
!          call a%Memor%realloc(a%npoinTop,a%EdgerEntryPoints,'EdgerEntryPoints','UpdateArraysToNewParallelNumbering')
!       endif
      
      !a%renpo and a%irenpo
      if (allocated(a%irenpo)) deallocate(a%irenpo)
      allocate(a%irenpo(a%npoin))
      a%irenpo = 0_ip
      do ipoin = 1,a%oldnpoin
         if (a%renpo(ipoin) /= 0) then
            a%renpo(ipoin) = a%ParallelLocalReordering(a%renpo(ipoin))
            a%irenpo(a%renpo(ipoin)) = ipoin
         endif
      enddo

      !Deallocate lots of things
      deallocate(a%ParallelLocalReordering)
      
      call a%Timers%UpdateArraysToNewParallelNumbering%Toc
   end subroutine
   
   subroutine HangingNodes(a)
      use typre
      use Mod_LinkedList
      use Mod_AdaptiveInterpolationList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip)              :: TestResult
      integer(ip) :: ielem,ispos,ipnode,ipnodb,inchild,inface,inchildface,iface
      integer(ip) :: ifpos,jelem,jface,jspos,jpnode,jpnodb,jnchild,jnface,jnchildface
      integer(ip) :: ifacnod(a%mnodb),jfacnod(a%mnodb)
      integer(ip) :: inodb,ipoin,parents(2),jnodb,jpoin,jv
      
      !integer(ip) :: first,last
      !integer(ip), allocatable :: next(:)
      
      type(LinkedListHeader), allocatable :: HeaderPointLevelList(:)
      type(LinkedListStorage)             :: StoragePointLevelList
      integer(ip) :: ilevel,npoinlevel
      
      type(LinkedListHeader), allocatable :: HeaderElementLevelList(:)
      type(LinkedListStorage)             :: StorageElementLevelList
      
      integer(ip), allocatable :: seen(:,:),iwa(:)
      real(rp), allocatable :: rwa(:)
      integer(ip) :: poini
      integer(ip) :: jparent, json,nz,lparent,lpoin,kpoin,poink
      real(rp) :: weightk
      
      type :: HangingNodesList
         integer(ip) :: nz = 0
         integer(ip), allocatable :: AssemblyPoints(:)
         real(rp), allocatable ::    AssemblyWeights(:)
      end type
      type (HangingNodesList), allocatable :: auxHangingList(:)
      
      integer(ip) :: totalnz, totalnzElementList
      
      !For list of hangings for each element
      type(LinkedListHeader), allocatable :: HeaderElementHanging(:)
      type(LinkedListStorage)             :: StorageElementHanging
      integer(ip), allocatable :: eseen(:)
      integer(ip) :: nhang, realnhang
      integer(ip), allocatable :: auxHangingListL(:)
      integer(ip) :: ielty,jelty
      
      type(CSRList) :: ElementEdgesList
      
      type(LinkedListHeader), allocatable :: HeaderEdgeElementsList(:)
      type(LinkedListStorage)             :: StorageEdgeElementsList
      integer(ip), allocatable            :: EdgePointerToParent(:)
      
      type(EdgerCriteriaTester) :: ECT
      integer(ip) :: npoinToAdd, PointsToAdd(3)
      
      integer(ip) :: edgej, edgek,ichild,iedge,inedge,iepos
      integer(ip) :: iparent,iparentedge,ipointoAdd,jedge,jnedge,parent1,parent2
      integer(ip) :: parentedgei, parentepos,totalnedges,lnodedge(2)
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: testparent
      
      type(InterpolationList) :: InterpList
      integer(ip) :: inode,ivariation,jchild,jnode,pnode,lnode
      integer(ip), pointer :: jlnods(:) => NULL()
      

      call a%Timers%HangingNodes%Tic
      
      !Hanging nodes
      allocate(HeaderElementHanging(a%nelem))
      call StorageElementHanging%Initialize(a%nelem,'NonRemovable','Repeatable',a%Memor)
      
      if (allocated(a%isHangingNode)) deallocate(a%isHangingNode)
      
      allocate(a%isHangingNode(a%npoin))
      a%isHangingNode = 0
      
      totalnzElementList = 0_ip
      !Identify hanging nodes
      do ielem = 1,a%nelem
         ispos = a%pnods(ielem)-1
         ipnode = a%pnods(ielem+1)-a%pnods(ielem)
         ielty = a%ElementTypeList(ielem)%ielty
         call ElementGetDimensions(ielty,ipnodb,inchild,inface,inchildface)
         
         ifpos = a%pfacs(ielem)-1
         do iface = 1,inface
            jelem = a%lfacs(1,ifpos+iface) 
            !If this is a hanging face
            if (jelem < 0) then
               jelem = -jelem
               jface = a%lfacs(2,ifpos+iface)
               jspos = a%pnods(jelem)-1
               jpnode = a%pnods(jelem+1)-a%pnods(jelem)
               jelty = a%ElementTypeList(jelem)%ielty
               call ElementGetDimensions(jelty,jpnodb,jnchild,jnface,jnchildface)
               
               !Check all the nodes of the face against the nodes in the other face
               !If they are not in the other face, then I am a hanging node
               call face_from_el(ielty,inface,ipnodb,a%lnods(ispos+1),iface,0,ifacnod)
               call face_from_el(jelty,jnface,jpnodb,a%lnods(jspos+1),jface,0,jfacnod)
               
               !Check if any of my nodes is hanging
               do inodb = 1,ipnodb
                  ipoin = ifacnod(inodb)
                  !If  the node has still not been classified as hanging, 
                  !then test it
                  if (a%isHangingNode(ipoin) == 0) then
                     TestResult = 1_ip
                     jpnodbloop: do jnodb = 1,jpnodb
                        jpoin = jfacnod(jnodb)
                        if (ipoin == jpoin) then
                           TestResult = 0_ip
                           exit jpnodbloop
                        endif
                     enddo jpnodbloop
                     a%isHangingNode(ipoin) = TestResult
                  endif
                  
                  !If a%isHangingNode, add it to the list of hanging nodes for the high level element
                  !I can only be sure that local nodes will be properly added to the list
                  !but I think this is enough
                  if (a%isHangingNode(ipoin) == 1_ip) then
                     call StorageElementHanging%Add(HeaderElementHanging(jelem),ipoin)
                     totalnzElementList = totalnzElementList + 1
                  endif
               enddo
            endif
         enddo
      enddo
      
      !We ghost communicate to make sure that hanging nodes are properly 
      !identified through domain boundaries
      call a%Communicator%GhostCommunicate(1_ip,a%isHangingNode)
      
      !-----------------------------------------------------------------
      !Only for 3D, additional to HangingListForElements
      !Elements which do not have a hanging face but
      !have a "hanging edge" do also have to be added!!
      !This is what we do here using StorageElementHanging
      
      if (a%ndime == 3) then
         !To mark already done edges
         totalnedges = a%pedger(a%npoin+1)-1
         if (totalnedges /= 0) then
            allocate(EdgePointerToParent(totalnedges))
            EdgePointerToParent = -1
            
            allocate(HeaderEdgeElementsList(totalnedges))
            call StorageEdgeElementsList%Initialize(4*totalnedges,'NonRemovable','Repeatable',a%Memor)
            
            allocate(ElementEdgesList%p(a%nelem+1))
            allocate(ElementEdgesList%l(12*a%nelem))
            ElementEdgesList%p(1) = 1
            do ielem = 1,a%nelem
               call a%GetLnodsFromElement(ielem,ipnode,lnods)
               !ipnode = a%pnods(ielem+1)-a%pnods(ielem)
               !ispos = a%pnods(ielem)-1
               ielty = a%ElementTypeList(ielem)%ielty
               call AdaptiveGetNedge(ielty,inedge)
               
               !List of edges for each element
               ElementEdgesList%p(ielem+1) = ElementEdgesList%p(ielem)+inedge
               iepos = ElementEdgesList%p(ielem)-1
               
               do iedge = 1,inedge
                  call AdaptiveGetEdges(ielty,iedge,lnodedge)
                  !ElementLocal to Local
                  lnodedge = lnods(lnodedge)
                  
                  !Who is the parent?
                  ECT%parents(1) = lnodedge(1)
                  ECT%parents(2) = lnodedge(2)
                  ECT%levels(1) = a%level(lnodedge(1))
                  ECT%levels(2) = a%level(lnodedge(2))
                  call EdgerCriteria(ECT)
                        
                  parent1 = ECT%orderedparents(1)
                  parent2 = ECT%orderedparents(2)
                  
                  !I do only need edges which are refined, 
                  !so I am sure they will be in pedger, edger
                  !Seek for the edge
                  jspos = a%pedger(parent1)-1
                  jnedge= a%pedger(parent1+1)-a%pedger(parent1)
                  
                  PointEdges_Loop: do jedge = 1,jnedge  
                     edgej = jspos+jedge
                     testparent = a%ledger(1,edgej)
                     if (parent2 == testparent) then
                        !add the element to the edge list 
                        !only if the element is last level (markel == 0) (to be exported, etc)
                        !Otherwise, the element is already added by its son
                        if (a%markel(ielem) == 0) then
                           call StorageEdgeElementsList%Add(HeaderEdgeElementsList(edgej),ielem) 
                        endif  
                        ElementEdgesList%l(iepos+iedge) = edgej
                        
                        !I take the oportunity to note down the parent edge if this has not been done
                        if (EdgePointerToParent(edgej) == -1) then
                           iparent = a%PointerToParent(ielem)
                           if (iparent /= 0) then
                              
                              call GEI_GetPositionInElementForLevel(a%GEI_ElementList(ielem),int(a%GEI_ElementList(ielem)%level,ip),ichild)
                              call Adaptive_up_edges(ielty,ichild,iedge,iparentedge)
                              if (iparentedge /= -1) then
                                 parentepos = ElementEdgesList%p(iparent)-1
                                 parentedgei = ElementEdgesList%l(parentepos+iparentedge)
                                 EdgePointerToParent(edgej) = parentedgei
                                 
                              else
                                 EdgePointerToParent(edgej) = 0
                              endif
                           else
                              EdgePointerToParent(edgej) = 0
                           endif
                        endif
                        exit PointEdges_Loop
                     endif
                  enddo PointEdges_Loop
               
               enddo
            enddo
            
            !Now, loop through edges and move up the hierarchy adding hanging nodes to elements
            !We are indeed looping through parent edges (pedger, ledger)
            do ipoin = 1,a%npoin
               ispos = a%pedger(ipoin)-1
               inedge = a%pedger(ipoin+1)-a%pedger(ipoin)
               
               do jedge = 1,inedge
                  edgej = ispos + jedge
                  
                  jpoin = a%ledger(1,edgej)
                  kpoin = a%ledger(2,edgej)
                  
                  !Add to elements in the current (parent) edge
                  !if (npoinToAdd /= 0) then
                  if (a%isHangingNode(kpoin) == 1) then
                     
                     !First my edge
                     call StorageEdgeElementsList%GoToFirstOfList(HeaderEdgeElementsList(edgej))
                     call StorageEdgeElementsList%GetNext(ielem)
                     do while (ielem /= -1)
                        call StorageElementHanging%Add(HeaderElementHanging(ielem),kpoin)
                        totalnzElementList = totalnzElementList+1
                        
                        call StorageEdgeElementsList%GetNext(ielem)
                     enddo
                     
                     !Second go up through parent edges
                     edgek = EdgePointerToParent(edgej)
                     do while (edgek > 0)
                        !Add to elements in the current (parent) edge
                        call StorageEdgeElementsList%GoToFirstOfList(HeaderEdgeElementsList(edgek))
                        call StorageEdgeElementsList%GetNext(ielem)
                        do while (ielem /= -1)
                           call StorageElementHanging%Add(HeaderElementHanging(ielem),kpoin)
                           totalnzElementList = totalnzElementList+1
                           
                           call StorageEdgeElementsList%GetNext(ielem)
                        enddo
                     
                        edgek = EdgePointerToparent(edgek)
                     enddo
                     
                     
                  endif
                  
                  
               enddo
            enddo
         
            call StorageEdgeElementsList%Dealloc
            deallocate(HeaderEdgeElementsList)
            deallocate(ElementEdgesList%p,ElementEdgesList%l)
            deallocate(EdgePointerToParent)
         endif
      endif
      !------------------------------------------------------------------
      
      
      !List of hanging nodes for each (higher level) element
      !Now I decompress StorageElementHanging, make sure nodes are not repeated
      !then compress again
      if (allocated(a%HangingListForElements%p)) then
         deallocate(a%HangingListForElements%p,a%HangingListForElements%l)
      endif
      allocate(a%HangingListForElements%p(a%nelem+1))
      allocate(a%HangingListForElements%l(totalnzElementList))
      
      
      allocate(eseen(a%npoin))
      eseen = 0
      a%HangingListForElements%p(1) = 1
      totalnzElementList = 0
      do ielem = 1,a%nelem
         call HeaderElementHanging(ielem)%GetNelem(nhang)
         realnhang = nhang
         if (nhang /= 0) then
            realnhang = 0
            call StorageElementHanging%GoToFirstOfList(HeaderElementHanging(ielem))
            call StorageElementHanging%GetNext(ipoin)
            do while (ipoin /= -1)
               if (eseen(ipoin) /= ielem) then
                  eseen(ipoin) = ielem
                  realnhang = realnhang + 1
                  totalnzElementList = totalnzElementList + 1
                  a%HangingListForElements%l(totalnzElementList) = ipoin
               endif
               call StorageElementHanging%GetNext(ipoin)
            enddo
         endif
         a%HangingListForElements%p(ielem+1) = a%HangingListForElements%p(ielem)+ realnhang
      enddo
      deallocate(HeaderElementHanging)
      call StorageElementHanging%Dealloc
      deallocate(eseen)
      !Reallocate if there were repetitions and sizes do not coincide
      if (totalnzElementList /= size(a%HangingListForElements%l)) then
         allocate(auxHangingListL(totalnzElementList))
         auxHangingListL(1:totalnzElementList) = a%HangingListForElements%l(1:totalnzElementList)
         call move_alloc(auxHangingListL,a%HangingListForElements%l)
         !call a%Memor%realloc(totalnzElementList,a%HangingListForElements%l,'HangingListForElements','HangingNodes')
      endif
      
      !I Need to do all first level nodes, then all second level nodes, etc...
      !Entry points cannot be used because for quads and hexas not all refined points are in edges
      allocate(HeaderPointLevelList(GEI_ElementMaxLevels+1))
      call StoragePointLevelList%Initialize(a%npoin,'NonRemov','NonRepeat',a%Memor)
      
      allocate(HeaderElementLevelList(GEI_ElementMaxLevels+1))
      call StorageElementLevelList%Initialize(a%nelem,'NonRemov','NonRepeat',a%Memor)
      
      do ipoin = 1,a%npoin
         call StoragePointLevelList%Add(HeaderPointLevelList(a%level(ipoin)+1),ipoin)
      enddo
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 1) then
            call StorageElementLevelList%Add(HeaderElementLevelList(a%GEI_ElementList(ielem)%level+1),ielem)
         endif
      enddo
      
      !Loop through the list of hierarchycally ordered points   
      !working arrays
      allocate(seen(2,a%npoin))
      seen = 0
      allocate(iwa(a%npoin))
      allocate(rwa(a%npoin))
      
      
      allocate(auxHangingList(a%npoin))
      
      totalnz = 0
      do ilevel = 1,GEI_ElementMaxLevels+1
         call StoragePointLevelList%GoToFirstOfList(HeaderPointLevelList(ilevel))
         call StoragePointLevelList%GetNext(ipoin)
         do while (ipoin /= -1) 
            !For all my children
            do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
               jparent = a%ledger(1,jv)
               json = a%ledger(2,jv)
               
               if (a%isHangingNode(json) == 1) then
                  nz = 0
                  
                  !I need to add to the assembly list the ipoin and jpoin
                  !own assembly lists multiplied times 0.5_rp
                  parents(1) = ipoin
                  parents(2) = jparent
                  
                  do lparent = 1,2
                     lpoin = parents(lparent)
                     
                     !If it is not a hanging node, simply add it
                     if (a%isHangingNode(lpoin) == 0) then
                        kpoin = lpoin
                        weightk = 0.5_rp
                        
                        call add_to_working_arrays
                     !Parent is also hanging, need to add all of them assemblypoints to the list
                     else
                        do poink = 1,auxHangingList(lpoin)%nz
                           weightk = auxHangingList(lpoin)%AssemblyWeights(poink)*0.5_rp
                           kpoin = auxHangingList(lpoin)%AssemblyPoints(poink)
                           
                           call add_to_working_arrays
                        enddo   
                     endif
                  enddo
                  
                  !Now I need to compress the results
                  auxHangingList(json)%nz = nz
                  allocate(auxHangingList(json)%AssemblyPoints(nz))
                  allocate(auxHangingList(json)%AssemblyWeights(nz))
                  
                  auxHangingList(json)%AssemblyPoints = iwa(1:nz)
                  auxHangingList(json)%AssemblyWeights = rwa(1:nz)
                  totalnz = totalnz + nz
                  
               endif
               
               !Add the point to the linked list
               !next(last) = json
               !last = json
            enddo
            call StoragePointLevelList%GetNext(ipoin)
         enddo
         
         call StorageElementLevelList%GoToFirstOfList(HeaderElementLevelList(ilevel))
         call StorageElementLevelList%GetNext(ielem)
         do while (ielem /= -1)
            ielty = a%ElementTypeList(ielem)%ielty
            ivariation = a%ElementTypeList(ielem)%ivariation
            call Adaptive_GetInterpolationList(ielty,ivariation,InterpList)
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            do inode = 1,InterpList%n
               jchild = InterpList%childElementList(inode)
               jnode  = InterpList%childObjectiveNodeList(inode)
               ispos = a%pPointerToChildren(ielem)-1
               jelem = a%lPointerToChildren(ispos+jchild)
               call a%GetLnodsFromElement(jelem,jpnode,jlnods)
               json = jlnods(jnode)
               
               if (a%isHangingNode(json) == 1) then
                  nz = 0
                  
                  do lparent = 1,InterpList%nInterpoLatingNodes(inode)
                     lnode = InterpList%InterpolatingNodesList(lparent,inode)
                     lpoin = lnods(lnode)
                     
                     !If it is not a hanging node, simply add it
                     if (a%isHangingNode(lpoin) == 0) then
                        kpoin = lpoin
                        weightk = InterpList%InterpolatingCoefficientsList(lparent,inode)
                        
                        call add_to_working_arrays
                     !Parent is also hanging, need to add all of them assemblypoints to the list
                     else
                        do poink = 1,auxHangingList(lpoin)%nz
                           weightk = auxHangingList(lpoin)%AssemblyWeights(poink)*InterpList%InterpolatingCoefficientsList(lparent,inode)
                           kpoin = auxHangingList(lpoin)%AssemblyPoints(poink)
                           
                           call add_to_working_arrays
                        enddo   
                     endif
                  enddo
                  
                  !Now I need to compress the results
                  auxHangingList(json)%nz = nz
                  allocate(auxHangingList(json)%AssemblyPoints(nz))
                  allocate(auxHangingList(json)%AssemblyWeights(nz))
                  
                  auxHangingList(json)%AssemblyPoints = iwa(1:nz)
                  auxHangingList(json)%AssemblyWeights = rwa(1:nz)
                  totalnz = totalnz + nz
               endif
            
            enddo
            call StorageElementLevelList%GetNext(ielem)
         enddo
         
      enddo
      
      deallocate(HeaderPointLevelList)
      call StoragePointLevelList%Dealloc
      call StorageElementLevelList%Dealloc
      
      if (allocated(a%pHangingList)) then
         deallocate(a%pHangingList)
         deallocate(a%lHangingList)
         deallocate(a%rHangingLIst)
      endif
      
      allocate(a%pHangingList(a%npoin+1))
      allocate(a%lHangingList(totalnz))
      allocate(a%rHangingList(totalnz))
      
      a%pHangingList(1) = 1
      
      do ipoin = 1,a%npoin
         nz = auxHangingList(ipoin)%nz
         a%pHangingList(ipoin+1) = nz + a%pHangingList(ipoin)
         ispos = a%pHangingList(ipoin)-1
         if (nz /= 0) then
            a%lHangingList(ispos+1:ispos+nz) = auxHangingList(ipoin)%AssemblyPoints
            a%rHangingList(ispos+1:ispos+nz) = auxHangingList(ipoin)%AssemblyWeights
            
            deallocate(auxHangingList(ipoin)%AssemblyPoints)
            deallocate(auxHangingList(ipoin)%AssemblyWeights)
         endif
      enddo
         
      deallocate(auxHangingList)
   
      call a%Timers%HangingNodes%Toc
   
contains
      
      !This internal subroutine adds the assembly point and weights to the working arrays
      subroutine add_to_working_arrays
         implicit none
         
         integer(ip) :: kz
         
         if (seen(1,kpoin) /= json) then
            nz = nz+1
            seen(1,kpoin) = json
            seen(2,kpoin) = nz
            iwa(nz) = kpoin
            rwa(nz) = weightk
         else
            kz = seen(2,kpoin)
            rwa(kz) = rwa(kz)+weightk
         endif
      end subroutine
   end subroutine
   
   subroutine SelectComponentsToExport(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      
      integer(ip) :: ielem,ispos,pnode,ipoin,icount,jnode,jpoin,phang,jhpos,knode,kpoin
      integer(ip), pointer :: lnods(:) => NULL(), lhangs(:) => NULL()
      
      type(MyParallelOrdering) :: ExportLocalOrdering
      type(MyParallelOrdering) :: TestExportLocalOrdering
      integer(ip), allocatable :: ExportLocalToGlobal(:),ExportProcessorList(:)
      integer(ip), allocatable :: ia(:),ja(:),pelpo(:),lelpo(:)
      integer(ip) :: ExportIpoin, ExportNpoinGhost
      logical :: doexport
      integer(ip) :: ispos2,pnode2,ipoin2
      integer(ip), pointer :: lnods2(:) => NULL()
      integer(ip) :: iExportElement
      !I want to export only last level elements with at least one local point
      !I want to export only points belonging to at least one exported element
      
      call a%Timers%SelectComponentsToExport%Tic
      
      !Allocations
      deallocate(a%ExportElement)
      deallocate(a%OldExportPointNumbering)
      
      
      !I need the old ExportPointNumbering
      call move_alloc(a%ExportPointNumbering,a%OldExportPointNumbering)
      a%OldExportNpoin = a%ExportNpoin
      
      call move_alloc(a%ExportElementNumbering,a%OldExportElementNumbering)
      a%OldExportNelem = a%ExportNelem
      
      allocate(a%ExportElement(a%nelem))
      allocate(a%ExportPointNumbering(a%npoin))
      allocate(a%ExportElementNumbering(a%nelem))
      a%ExportElementNumbering = 0
      
      
      !Point Initializations
      do ipoin = 1,a%npoinLocal
         a%ExportPointNumbering(ipoin) = ipoin
      enddo
      a%ExportPointNumbering(a%npoinLocal+1:a%npoin) = 0_ip
      a%ExportNpoin = a%npoinLocal
      
      !Element Initializations
      a%ExportElement = .false.
      a%ExportNelem = 0
      a%ExportLnodsSizes = 0
      
      !Choose points and elements
      do ielem = 1,a%nelem
         if (a%markel(ielem) == 0) then
            doexport = .false.
            
            !First criteria, my nodes
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            lnods => a%lnods(ispos+1:ispos+pnode)
            ipoin = minval(lnods)
            if (ipoin <= a%npoinLocal) then
               doexport = .true.
            endif
            
            !Second criteria
            !If I am the higher level hanging element of a hanging node
            !only evaluate if doexport is still false
            if (doexport .eqv. .false.) then
               pnode2 = a%HangingListForElements%p(ielem+1) - a%HangingListForElements%p(ielem)
               if (pnode2 /= 0) then
                  ispos2 = a%HangingListForElements%p(ielem)-1
                  lnods2 => a%HangingListForElements%l(ispos2+1:ispos2+pnode2)
                  ipoin2 = minval(lnods2)
                  if (ipoin2 <= a%npoinLocal) then
                     doexport = .true.
                  endif
               endif
            endif
            
            if (doexport .eqv. .true.) then
            !if (.true.) then
                
            
               !Element
               a%ExportElement(ielem) = .true.
               a%ExportNelem = a%ExportNelem+1
               a%ExportLnodsSizes = a%ExportLnodsSizes + pnode
               
               a%ExportElementNumbering(ielem) = a%ExportNelem
               
               !Nodes
               do jnode = 1,pnode
                  jpoin = lnods(jnode)
                  if (a%ExportPointNumbering(jpoin) == 0) then
                     a%ExportPointNumbering(jpoin) = 1_ip
                     
                     !If these are hanging nodes, mark the parents also
                     phang = a%pHangingList(jpoin+1)-a%pHangingList(jpoin)
                     if (phang /= 0) then
                        jhpos = a%pHangingList(jpoin)-1
                        do knode = 1,phang
                           kpoin = a%lHangingList(jhpos+knode)
                           if (a%ExportPointNumbering(kpoin) == 0) then
                              a%ExportPointNumbering(kpoin) = 1_ip
                           endif
                        enddo
                     endif
                  endif
               enddo
            endif
         endif
      enddo
      
      !Export Point Numbering
      a%ExportNpoin = a%npoinLocal
      do ipoin = a%npoinLocal+1,a%npoin
         if (a%ExportPointNumbering(ipoin) /= 0) then
            a%ExportNpoin = a%ExportNpoin+1
            a%ExportPointNumbering(ipoin) = a%ExportNpoin
         endif
      enddo
      
      !For going from exported elements to internal elements
      deallocate(a%iRenelExport)
      allocate(a%iRenelExport(a%ExportNelem))
      iexportElement = 0
      do ielem = 1,a%nelem
         if (a%ExportElement(ielem) .eqv. .true.) then
            iexportElement = iexportElement + 1
            a%irenelExport(iexportElement) = ielem
         endif
      enddo
      
      !Stablish dimensions for the HangingList to be exported
      a%ExportHangingSizes = 0
      do ipoin = 1,a%npoin
         if (a%ExportPointNumbering(ipoin) /= 0) then
            phang = a%pHangingList(ipoin+1)- a%pHangingList(ipoin)
            a%ExportHangingSizes = a%ExportHangingSizes + phang
         endif
      enddo
      !Communicate elements with hanging nodes from other processes
      !wich assembly to a node in this process but which I CAN'T SEE FROM HERE
      call a%AdditionalCommunicateHangingElements
      
      !BuildAnExportCommunicator
      allocate(ExportLocalToGlobal(a%ExportNpoin+a%ExtraNpoin))
      allocate(ExportProcessorList(a%ExportNpoin+a%ExtraNpoin))
      do ipoin = 1,a%npoin
         ExportIpoin = a%ExportPointNumbering(ipoin)
         if (ExportIpoin /= 0) then
            ExportLocalToGlobal(ExportIpoin) = ipoin
            ExportProcessorList(ExportIpoin) = ipoin
         endif
      enddo
      !ExportPointsToGlobal
      call a%LocalOrdering%Local2Global(a%ExportNpoin,ExportLocalToGlobal,ExportLocalToGlobal)
      call a%LocalOrdering%GetProc(a%ExportNpoin,ExportProcessorList,ExportProcessorList)
      !ExtraPoints
      do ipoin = 1,a%ExtraNpoin
         ExportLocalToGlobal(a%ExportNpoin+ipoin) = a%ExtraPoints(1,ipoin)
         ExportProcessorList(a%ExportNpoin+ipoin) = a%ExtraPoints(2,ipoin)
      enddo
      call ExportLocalOrdering%Init(a%MPIcomm,a%ExportNpoin+a%ExtraNpoin,ExportLocalToGlobal,ExportProcessorList,a%Memor)
      deallocate(ExportLocalToGlobal)
      deallocate(ExportProcessorList)
      
      !We reinitialize the communicator
      call a%ExportCommunicator%Deallocate
      

      ExportNpoinGhost = a%ExportNpoin + a%ExtraNpoin - a%npoinLocal
      call a%ExportCommunicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,ExportNpoinGhost,a%gnpoin,ExportLocalOrdering,a%Memor)
      
      call ExportLocalOrdering%Deallocate(a%Memor)
      
      call a%Timers%SelectComponentsToExport%Toc
      
   end subroutine
   
   subroutine AdditionalCommunicateHangingElements(a)
      use typre
      use Mod_LinkedList
      use Mod_HeapSortExternal
      use Mod_ToCSR
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      !After this subroutine is run, we have:
      !nodes I Need, hanging list for NodesINeed, ElementsINeed
      !This subroutine communicates the additional elements whith hanging nodes
      !which assembly in a parent node in another subdomain.
      !And they were called EXTRA.
      
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: minpoin,maxpoin,ielem,glnods(MaxNodesPerElement),lranks(MaxNodesPerElement),ispos,pnode
      integer(ip) :: iminval,jnode,jpoin,phang,jhpos,knode,kpoin,krank,ispos2,ippos,ippos2
      integer(ip) :: ElemCount, ExtraNpoinAux ,inode,ipoin,ipoinLocal,ipoinGlobal,jhang
      integer(ip) :: sumNpoinINeed,sumHangingListSizeINeed
      
      
      type(LinkedListHeader), allocatable :: HeaderElementsOthersNeed(:)
      type(LinkedListStorage) :: StorageElementsOthersNeed
      
      type(CSRList), allocatable :: ElementsOthersNeed(:),ElementsINeed(:),HangingNodesOthersNeed(:),HangingNodesINeed(:)
      type :: NodeList
         integer(ip), allocatable :: l(:,:)
      end type
      type(NodeList), allocatable :: NodesOthersNeed(:),NodesINeed(:)
      integer(ip), allocatable :: NpoinINeed(:),NpoinOthersNeed(:)
      integer(ip), allocatable :: HangingListSizeINeed(:),HangingListSizeOthersNeed(:),lnodsSizeINeed(:),lnodsSizeOthersNeed(:)
      
      integer(ip), allocatable :: iwa(:),seen(:),auxSizesINeed(:,:), auxSizesOthersNeed(:,:)
      integer(ip), allocatable :: seenProc(:)
      logical, allocatable :: needToAdd(:)
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3,mtag4 = 4, mtag5 = 0, mtag6 = 6
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:),irequest05(:),irequest06(:)
      integer, allocatable :: irequest07(:),irequest08(:),irequest09(:),irequest10(:),irequest11(:),irequest12(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: irank,jrank
      integer(ip), allocatable :: auxGlobalNpoinINeed(:), auxGlobalNpoinINeed2(:)
      integer(ip) :: binarysearchvalue,npoinCount,npoinCount2,prepoint
      logical :: testneedtoadd
      
      integer(ip), allocatable :: auxNeedToAddEnsureNoRankRepeat(:),auxBinarySearchForRank(:)
      
      integer(ip) :: iaux1,iaux2,maxNpoinINeed
      type(i1p), allocatable :: auxi1pHanging(:)
      type(r1p), allocatable :: auxr1pHanging(:)
      integer(ip) :: auxcountHanging
      
      
      !We start by deallocating the already existing structure
      a%ExtraNpoin = 0
      a%ExtraNelem = 0
      a%ExtraLnodsSizes = 0
      a%ExtraHangingSizes = 0
      deallocate(a%ExtraPoints)
      deallocate(a%ExtraHangingList%p)
      deallocate(a%ExtraHangingList%l)
      deallocate(a%ExtraHangingList%r)
      deallocate(a%ExtraElements%p)
      deallocate(a%ExtraElements%l)
      
      if (allocated(a%ExtraElementsListOthersNeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%ExtraElementsListOthersNeed(irank+1)%l)
         enddo
      else
         allocate(a%ExtraElementsListOthersNeed(a%MPIsize))
      endif
      
      
      if (a%npoinLocal > 0) then
         call a%LocalOrdering%Local2Global(1_ip,minpoin)
         call a%LocalOrdering%Local2Global(a%npoinLocal,maxpoin)
      else
         maxpoin = 0
         minpoin = 0
      endif
      
      !Auxiliary array to know if the element has a node from this rank
      allocate(seenProc(a%MPIsize))
      seenProc = 0
      
      !allocate linked list of elements
      allocate(HeaderElementsOthersNeed(a%MPIsize))
      call StorageElementsOthersNeed%Initialize(a%npoinGhost,'NonRemov','Repeatable',a%Memor)
      
      if (allocated(a%ExtraNelemOthersNeed)) deallocate(a%ExtraNelemOthersNeed)
      allocate(a%ExtraNelemOthersNeed(a%MPIsize))
      a%ExtraNelemOthersNeed = 0
      allocate(lnodsSizeOthersNeed(a%MPIsize))
      lnodsSizeOthersNeed = 0
      
      !First to select and fill the elements to send
      do ielem = 1,a%nelem
         !only last level elements
         if (a%markel(ielem) == 0) then
            !only elements of which I am responsible
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            lnods => a%lnods(ispos+1:ispos+pnode)
            call a%LocalOrdering%Local2Global(pnode,lnods,glnods)
            iminval = minval(glnods(1:pnode))
            !If I am responsible for the element
            if ((iminval >= minpoin) .and. (iminval <= maxpoin)) then
               call a%LocalOrdering%GetProc(pnode,lnods,lranks)
               
               !mark the rank as they have a node in the element (no need to send)
               seenProc(lranks(1:pnode)+1) = ielem
               
               do jnode = 1,pnode
                  jpoin = lnods(jnode)
                  !If these are hanging nodes, mark the parents also
                  phang = a%pHangingList(jpoin+1)-a%pHangingList(jpoin)
                  if (phang /= 0) then
                     jhpos = a%pHangingList(jpoin)-1
                     do knode = 1,phang
                        kpoin = a%lHangingList(jhpos+knode)
                        call a%LocalOrdering%GetProc(kpoin,krank)
                        
                        if (seenProc(krank+1) /= ielem) then
                           seenProc(krank+1) = ielem
                           
                           !Add the element to the linked list of elements to send
                           call StorageElementsOthersNeed%Add(HeaderElementsOthersNeed(krank+1),ielem)
                           !increase the sizes to store the list of elements for this processor
                           a%ExtraNelemOthersNeed(krank+1) = a%ExtraNelemOthersNeed(krank+1)+1
                           lnodsSizeOthersNeed(krank+1) = lnodsSizeOthersNeed(krank+1)+pnode
                        endif
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
      
      !Second to select and fill the nodes to send to each rank
      !allocate linked list of elements
      allocate(iwa(a%npoin))
      allocate(seen(a%npoin))
      seen = 0
      
      allocate(HangingListSizeOthersNeed(a%MPIsize))
      HangingListSizeOthersNeed = 0
      
      allocate(ElementsOthersNeed(a%MPIsize))
      allocate(HangingNodesOthersNeed(a%MPIsize))
      allocate(NodesOthersNeed(a%MPIsize))
      allocate(npoinOthersNeed(a%MPIsize))
      NpoinOthersNeed = 0
      
      do irank = 0,a%MPIsize-1
         npoinOthersNeed(irank+1) = 0
         
         !Allocate the elements to send list
         allocate(ElementsOthersNeed(irank+1)%p(a%ExtraNelemOthersNeed(irank+1)+1))
         allocate(ElementsOthersNeed(irank+1)%l(lnodsSizeOthersNeed(irank+1)))
         ElementsOthersNeed(irank+1)%p(1) = 1
         ElemCount = 0
         
         allocate(a%ExtraElementsListOthersNeed(irank+1)%l(a%ExtraNelemOthersNeed(irank+1)))
         
         
         !Start looping list
         call StorageElementsOthersNeed%GoToFirstOfList(HeaderElementsOthersNeed(irank+1))
         call StorageElementsOthersNeed%GetNext(ielem)
         do while (ielem /= -1)
            
         
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            lnods => a%lnods(ispos+1:ispos+pnode)
            do inode = 1,pnode
               ipoin = lnods(inode)
               if (seen(ipoin) == 0) then
                  npoinOthersNeed(irank+1) = npoinOthersNeed(irank+1)+1
                 
                  !Fill seen and iwa
                  seen(ipoin) = npoinOthersNeed(irank+1)
                  iwa(npoinOthersNeed(irank+1)) = ipoin
                  
!                   call a%LocalOrdering%Local2Global(ipoin,iwa(1,npoinOthersNeed(irank+1)))
!                   call a%LocalOrdering%GetProc(ipoin,iwa(2,npoinOthersNeed(irank+1)))
                  
                  !add to sizes for hanging list to send for this rank
                  phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
                  HangingListSizeOthersNeed(irank+1) = HangingListSizeOthersNeed(irank+1)+phang
                  
                  !if hanging node also add its parents
                  ippos = a%pHangingList(ipoin)-1
                  do jhang = 1,phang
                     jpoin = a%lHangingList(ippos+jhang)
                     
                     if (seen(jpoin) == 0) then
                        npoinOthersNeed(irank+1) = npoinOthersNeed(irank+1)+1
                     
                        !Fill seen and iwa
                        seen(jpoin) = npoinOthersNeed(irank+1)
                        iwa(npoinOthersNeed(irank+1)) = jpoin
!                         call a%LocalOrdering%Local2Global(jpoin,iwa(1,npoinOthersNeed(irank+1)))
!                         call a%LocalOrdering%GetProc(jpoin,iwa(2,npoinOthersNeed(irank+1)))   
                     endif
                  enddo
                     
                 
               endif
            enddo
            
            !Fill pnods and lnods to send to this rank
            ElemCount = ElemCount+1
            ispos2 = ElementsOthersNeed(irank+1)%p(ElemCount)-1
            ElementsOthersNeed(irank+1)%p(ElemCount+1) = ispos2 + 1 + pnode 
            !Fill lnods with the send local numbering
            ElementsOthersNeed(irank+1)%l(ispos2+1:ispos2+pnode) = seen(lnods)
            
            a%ExtraElementsListOthersNeed(irank+1)%l(ElemCount) = ielem
            
            call StorageElementsOthersNeed%GetNext(ielem)
         enddo
         
         !Allocate and fill the list of nodes
         allocate(NodesOthersNeed(irank+1)%l(NpoinOthersNeed(irank+1),2))
         call a%LocalOrdering%Local2Global(NpoinOthersNeed(irank+1),iwa,NodesOthersNeed(irank+1)%l(:,1))
         call a%LocalOrdering%GetProc(NpoinOthersNeed(irank+1),iwa,NodesOthersNeed(irank+1)%l(:,2))
         !NodesOthersNeed(irank+1)%l = iwa(:,1:NpoinOthersNeed(irank+1))
         
         !Allocate the list of hanging nodes to send
         allocate(HangingNodesOthersNeed(irank+1)%p(NpoinOthersNeed(irank+1)+1))
         allocate(HangingNodesOthersNeed(irank+1)%l(HangingListSizeOthersNeed(irank+1)))
         allocate(HangingNodesOthersNeed(irank+1)%r(HangingListSizeOthersNeed(irank+1)))
         
         !Fill the list of hanging nodes to send
         !iwa has the nodes to send in local numbering
         HangingNodesOthersNeed(irank+1)%p(1) = 1
         do inode = 1,NpoinOthersNeed(irank+1)
            ipoin = iwa(inode)
            
            phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
            ippos = a%pHangingList(ipoin)-1
            
            ippos2 = HangingNodesOthersNeed(irank+1)%p(inode)-1
            
            !I send the hanging nodes in the local numbering of the send array
            HangingNodesOthersNeed(irank+1)%l(ippos2+1:ippos2+phang) = seen(a%lHangingList(ippos+1:ippos+phang))
            HangingNodesOthersNeed(irank+1)%r(ippos2+1:ippos2+phang) = a%rHangingList(ippos+1:ippos+phang)
            !call a%LocalOrdering%Local2Global(phang,a%lHangingList(ippos+1:ippos+phang),HangingNodesOthersNeed(irank+1)%l(ippos2+1:ippos2+phang)) 
            HangingNodesOthersNeed(irank+1)%p(inode+1) = HangingNodesOthersNeed(irank+1)%p(inode) + phang
         enddo
         
         !clean seen
         do inode = 1,NpoinOthersNeed(irank+1)
            ipoin = iwa(inode)
            seen(ipoin) = 0
         enddo
      enddo
      
      call StorageElementsOthersNeed%Dealloc
      
      !Communicate with whom I need to communicate
      !maybe this could be only with neighbours, not sure after load rebalancing
      if (allocated(a%ExtraNelemINeed)) deallocate(a%ExtraNelemINeed)
      allocate(a%ExtraNelemINeed(a%MPIsize))
      call MPI_AllToAll(a%ExtraNelemOthersNeed, 1, MPI_INTEGER4, a%ExtraNelemINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
     
      !Allocate all MPI irequests
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      allocate(irequest05(a%MPIsize))
      allocate(irequest06(a%MPIsize))
      allocate(irequest07(a%MPIsize))
      allocate(irequest08(a%MPIsize))
      allocate(irequest09(a%MPIsize))
      allocate(irequest10(a%MPIsize))
      allocate(irequest11(a%MPIsize))
      allocate(irequest12(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL
      irequest05 = MPI_REQUEST_NULL
      irequest06 = MPI_REQUEST_NULL
      irequest07 = MPI_REQUEST_NULL
      irequest08 = MPI_REQUEST_NULL
      irequest09 = MPI_REQUEST_NULL
      irequest10 = MPI_REQUEST_NULL
      irequest11 = MPI_REQUEST_NULL
      irequest12 = MPI_REQUEST_NULL
      
      !Send and receive
      allocate(auxSizesOthersNeed(3,a%MPIsize))
      allocate(auxSizesINeed(3,a%MPIsize))
      do irank = 0,a%MPIsize-1
         !Send
         if (a%ExtraNelemOthersNeed(irank+1) /= 0) then
            auxSizesOthersNeed(1,irank+1) = lnodsSizeOthersNeed(irank+1)
            auxSizesOthersNeed(2,irank+1) = NpoinOthersNeed(irank+1)
            auxSizesOthersNeed(3,irank+1) = HangingListSizeOthersNeed(irank+1)
            
            call MPI_ISEND(auxSizesOthersNeed(:,irank+1),3,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         !Receive
         if (a%ExtraNelemINeed(irank+1) /= 0) then 
            call MPI_IRECV(auxSizesINeed(:,irank+1),3,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      !Wait for communications to finish
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)   
      
      !Allocate with the received sizes
      allocate(HangingListSizeINeed(a%MPIsize))
      allocate(lnodsSizeINeed(a%MPIsize))
      allocate(npoinINeed(a%MPIsize))
      HangingListSizeINeed = 0
      lnodsSizeINeed = 0
      NpoinINeed = 0
      
      allocate(ElementsINeed(a%MPIsize))
      allocate(NodesINeed(a%MPIsize))
      allocate(HangingNodesINeed(a%MPIsize))
      
      do irank = 0,a%MPIsize-1
         if (a%ExtraNelemINeed(irank+1) /= 0) then
            lnodsSizeINeed(irank+1) = auxSizesINeed(1,irank+1)
            NpoinINeed(irank+1) = auxSizesINeed(2,irank+1)
            HangingListSizeINeed(irank+1) = auxSizesINeed(3,irank+1)
            
            allocate(ElementsINeed(irank+1)%p(a%ExtraNelemINeed(irank+1)+1))
            allocate(ElementsINeed(irank+1)%l(lnodsSizeINeed(irank+1)))
            allocate(NodesINeed(irank+1)%l(NpoinINeed(irank+1),2))
            allocate(HangingNodesINeed(irank+1)%p(NpoinINeed(irank+1)+1))
            allocate(HangingNodesINeed(irank+1)%l(HangingListSizeINeed(irank+1)))
            allocate(HangingNodesINeed(irank+1)%r(HangingListSizeINeed(irank+1)))
         endif
      enddo
      
      !Send and receive the large amount of information
      !Send
      do irank = 0,a%MPIsize-1
         if (a%ExtraNelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(ElementsOthersNeed(irank+1)%p,a%ExtraNelemOthersNeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(ElementsOthersNeed(irank+1)%l,lnodsSizeOthersNeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
            call MPI_ISEND(NodesOthersNeed(irank+1)%l,2*NpoinOthersNeed(irank+1),MPI_INTEGER4, irank, mtag3, a%MPIcomm,irequest03(irank+1), ierr)
            call MPI_ISEND(HangingNodesOthersNeed(irank+1)%p,NpoinOthersNeed(irank+1)+1,MPI_INTEGER4, irank, mtag4, a%MPIcomm,irequest04(irank+1), ierr)
            call MPI_ISEND(HangingNodesOthersNeed(irank+1)%l,HangingListSizeOthersNeed(irank+1),MPI_INTEGER4, irank, mtag5, a%MPIcomm,irequest05(irank+1), ierr)
            call MPI_ISEND(HangingNodesOthersNeed(irank+1)%r,HangingListSizeOthersNeed(irank+1),MPI_REAL8, irank, mtag6, a%MPIcomm,irequest11(irank+1), ierr)
         endif
      enddo
      !Receive
      do irank = 0,a%MPIsize-1
         if (a%ExtraNelemINeed(irank+1) /= 0) then
            call MPI_IRECV(ElementsINeed(irank+1)%p,a%ExtraNelemINeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest06(irank+1), ierr)
            call MPI_IRECV(ElementsINeed(irank+1)%l,lnodsSizeINeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest07(irank+1), ierr)
            call MPI_IRECV(NodesINeed(irank+1)%l,2*NpoinINeed(irank+1),MPI_INTEGER4, irank, mtag3, a%MPIcomm,irequest08(irank+1), ierr)
            call MPI_IRECV(HangingNodesINeed(irank+1)%p,NpoinINeed(irank+1)+1,MPI_INTEGER4, irank, mtag4, a%MPIcomm,irequest09(irank+1), ierr)
            call MPI_IRECV(HangingNodesINeed(irank+1)%l,HangingListSizeINeed(irank+1),MPI_INTEGER4, irank, mtag5, a%MPIcomm,irequest10(irank+1), ierr)
            call MPI_IRECV(HangingNodesINeed(irank+1)%r,HangingListSizeINeed(irank+1),MPI_REAL8, irank, mtag6, a%MPIcomm,irequest12(irank+1), ierr)
         endif
      enddo
      !Wait
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest05, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest06, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest07, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest08, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest09, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest10, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest11, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest12, MPI_STATUSES_IGNORE, ierr)
      
      !Now I need to merge the NodesINeed, HangingNodesINeed and ElementsINeed into a single array
      !We start with the nodes
      sumNpoinINeed = sum(NpoinINeed)
      allocate(a%ExtraPoints(2,sumNpoinINeed))
      a%ExtraNpoin = 0
      
      !ExtraHanging
      allocate(a%ExtraHangingList%p(sumNpoinINeed+1))
      sumHangingListSizeINeed = sum(HangingListSizeINeed)
      allocate(a%ExtraHangingList%l(sumHangingListSizeINeed))
      allocate(a%ExtraHangingList%r(sumHangingListSizeINeed))
      a%ExtraHangingList%p(1) = 1
      
      !ExtraElements
      allocate(a%ExtraElements%p(sum(a%ExtraNelemINeed)+1))
      allocate(a%ExtraElements%l(sum(lnodsSizeINeed)))
      a%ExtraElements%p(1) = 1
      a%ExtraNelem = 0
      
      !This is necessary to check which nodes have already been added and avoid repeating them
      deallocate(seen)
!       allocate(seen(a%gnpoin))
!       seen = 0
      maxNpoinINeed = maxval(NpoinINeed)
      allocate(needToAdd(maxNpoinINeed))
      
      !here
      !--------------------------------------------------------
      allocate(auxGlobalNpoinINeed(sumNpoinINeed))
      npoinCount = 0
      do irank = 0,a%MPIsize-1
         !Only if I Need points from them
         do inode = 1,NpoinINeed(irank+1)
            ipoinGlobal = NodesINeed(irank+1)%l(inode,1)
            testneedtoadd = .false.
            
            !If not an existing local or ghost point
            call a%LocalOrdering%Global2Local(ipoinGlobal,ipoinLocal)
            if (ipoinLocal == 0) then
               testneedtoadd = .true. 
            elseif (a%ExportPointNumbering(ipoinLocal) == 0) then
               testneedtoadd = .true.
            endif
            
            if (testneedtoadd .eqv. .true.) then
               npoinCount = npoinCount+1
               auxGlobalNpoinINeed(npoinCount) = ipoinGlobal
            endif
         enddo
      enddo
      call HeapSort(npoinCount,auxGlobalNpoinINeed)
      
      !Eliminate repeated entries
!       allocate(auxGlobalNpoinINeed2(sumNpoinINeed))
!       prepoint = -1
!       npoinCount2 = 0
!       do ipoin = 1,npoinCount
!          ipoinGlobal = auxGlobalNpoinINeed(ipoin)
!          if (ipoinGlobal /= prepoint) then
!             npoinCount2 = npoinCount2+1
!             auxGlobalNpoinINeed2(npoinCount2) = ipoinGlobal
!             prepoint = ipoinGlobal
! 
!             a%ExtraNpoin = a%ExtraNpoin+1
!          endif
!       enddo
!       call move_alloc(auxGlobalNpoinINeed2,auxGlobalNpoinINeed)   
!       npoinCount = npoinCount2

      call EliminateRepeatedEntries(npoinCount,auxGlobalNpoinINeed,npoinCount2)
      a%ExtraNpoin = a%ExtraNpoin+npoinCount2
      npoinCount = npoinCount2
      
      !We want the need to add to be only once per point,
      !even if the point is sent by more than one domain
      allocate(auxNeedToAddEnsureNoRankRepeat(npoinCount))
      auxNeedToAddEnsureNoRankRepeat = 0
      allocate(auxBinarySearchForRank(maxNpoinINeed))
      allocate(auxi1pHanging(npoinCount),auxr1pHanging(npoinCount))
      auxcountHanging = 0
      
      ExtraNpoinAux = 0
      do irank = 0,a%MPIsize-1
         auxBinarySearchForRank = 0
         
         !Only if I Need points from them
         do inode = 1,NpoinINeed(irank+1)
            ipoinGlobal = NodesINeed(irank+1)%l(inode,1)
            
            needToAdd(inode) = .false.
            !If not an existing local or ghost point
            call a%LocalOrdering%Global2Local(ipoinGlobal,ipoinLocal)
            if (ipoinLocal == 0) then
               needToAdd(inode) = .true. 
            elseif (a%ExportPointNumbering(ipoinLocal) == 0) then
               needToAdd(inode) = .true.
            endif
            
            !Store the ExportPointNumbering in NodesINeed
            if (needToAdd(inode) .eqv. .false.) then
               NodesINeed(irank+1)%l(inode,1) = a%ExportPointNumbering(ipoinLocal)
            
            elseif (needToAdd(inode) .eqv. .true.) then
               !Check if it was already added
               binarysearchvalue = binarySearch_I(auxGlobalNpoinINeed(1:npoinCount),ipoinGlobal)
               auxBinarySearchForRank(inode) = binarysearchvalue
               if (auxNeedToAddEnsureNoRankRepeat(binarysearchvalue) == 0) then
                  a%ExtraPoints(1,binarysearchvalue) = NodesINeed(irank+1)%l(inode,1)
                  a%ExtraPoints(2,binarysearchvalue) = NodesINeed(irank+1)%l(inode,2)
                  auxNeedToAddEnsureNoRankRepeat(binarysearchvalue) = binarysearchvalue
               else
                  needToAdd(inode) = .false.
               endif
               
               !Put in Nodes I Need the export position
               NodesINeed(irank+1)%l(inode,1) = binarysearchvalue + a%ExportNpoin
            endif
         enddo
         
      !--------------------------------------------------------------------
      
!          !Only if I Need points from them
!          do inode = 1,NpoinINeed(irank+1)
!             ipoinGlobal = NodesINeed(irank+1)%l(inode,1)
!             
!             needToAdd(inode) = .false.
!             !If not an existing local or ghost point
!             call a%LocalOrdering%Global2Local(ipoinGlobal,ipoinLocal)
!             if (ipoinLocal == 0) then
!                needToAdd(inode) = .true. 
!             elseif (a%ExportPointNumbering(ipoinLocal) == 0) then
!                needToAdd(inode) = .true.
!             endif
!             
!             !Store the ExportPointNumbering in NodesINeed
!             if (needToAdd(inode) .eqv. .false.) then
!                NodesINeed(irank+1)%l(inode,1) = a%ExportPointNumbering(ipoinLocal)
!             
!             elseif (needToAdd(inode) .eqv. .true.) then
!                !Check if it was already added
!                if (seen(ipoinGlobal) == 0) then
!                   a%ExtraNpoin = a%ExtraNpoin+1
!                   seen(ipoinGlobal) = a%ExtraNpoin
!                   
!                   a%ExtraPoints(1,a%ExtraNpoin) = ipoinGlobal
!                   a%ExtraPoints(2,a%ExtraNpoin) = NodesINeed(irank+1)%l(inode,2)
!                   
!                   !Put in seen the export position
!                   seen(ipoinGlobal) = a%ExportNpoin + a%ExtraNpoin
!                else
!                   needToAdd(inode) = .false.
!                
!                endif
!                !Put in Nodes I Need the export position
!                NodesINeed(irank+1)%l(inode,1) = seen(ipoinGlobal)
!                
!             endif
!          enddo
         
         !Now I need to write extra nodes and extra hanging ALREADY IN EXPORT NUMBERING
         !We start with extra hanging
         !First we count the sizes
!          do inode = 1,NpoinINeed(irank+1)
!             if (NeedToAdd(inode) .eqv. .true.) then
!                ipoin = auxBinarySearchForRank(inode)
!                
!                phang = HangingNodesINeed(irank+1)%p(inode+1)-HangingNodesINeed(irank+1)%p(inode)
!                a%ExtraHangingList%p(ipoin) = phang
!             endif
!          enddo
!          
!          iaux1 = 1
!          do ipoin = 1,npoinCount
!             iaux2 = a%ExtraHangingList%p(ipoin)
!             a%ExtraHangingList%p(ipoin) = iaux1
!             iaux1 = iaux1+iaux2
!          enddo
!          a%ExtraHangingList%p(npoinCount+1) = iaux1
         
         do inode = 1,NpoinINeed(irank+1)
            if (NeedToAdd(inode) .eqv. .true.) then
               ipoin = auxBinarySearchForRank(inode)
               
               phang = HangingNodesINeed(irank+1)%p(inode+1)-HangingNodesINeed(irank+1)%p(inode)
               ippos = HangingNodesINeed(irank+1)%p(inode)-1
               
               allocate(auxi1pHanging(ipoin)%l(phang))
               allocate(auxr1pHanging(ipoin)%a(phang))
               auxcountHanging = auxcountHanging+phang
               
               do jhang = 1,phang
                  jnode = HangingNodesINeed(irank+1)%l(ippos+jhang)
                  jpoin = NodesINeed(irank+1)%l(jnode,1)
                  
                  !Fill the list
                  auxi1pHanging(ipoin)%l(jhang) = jpoin
                  auxr1pHanging(ipoin)%a(jhang) = HangingNodesINeed(irank+1)%r(ippos+jhang)
                  
               enddo
               
            endif
         enddo
         
         
         !Next we need to export the ExtraElements
         do ielem = 1,a%ExtraNelemINeed(irank+1)
            a%ExtraNelem = a%ExtraNelem+1
            
            pnode = ElementsINeed(irank+1)%p(ielem+1)-ElementsINeed(irank+1)%p(ielem)
            ispos = ElementsINeed(irank+1)%p(ielem)-1
            
            ispos2 = a%ExtraElements%p(a%ExtraNelem)-1
            !store, but from Sending to Export numbering
            a%ExtraElements%l(ispos2+1:ispos2+pnode) = NodesINeed(irank+1)%l(ElementsINeed(irank+1)%l(ispos+1:ispos+pnode),1)
            
            a%ExtraElements%p(a%ExtraNelem+1) = a%ExtraElements%p(a%ExtraNelem)+pnode
         enddo
      enddo
      
      call a%Memor%allocObj(0,'a%ExtraHangingList','additionalCommunicate',auxcountHanging*(rp))
      call a%Memor%allocObj(0,'a%ExtraHangingList','additionalCommunicate',auxcountHanging*(ip))
      call i1p2CSR(npoinCount,auxi1pHanging,a%ExtraHangingList%p,a%ExtraHangingList%l,a%Memor,'a%ExtraHangingList')
      call r1p2CSR(npoinCount,auxr1pHanging,a%ExtraHangingList%p,a%ExtraHangingList%r,a%Memor,'a%ExtraHangingList')
      
      
      deallocate(auxi1pHanging,auxr1pHanging)
      
      deallocate(auxNeedToAddEnsureNoRankRepeat,auxBinarySearchForRank)
      a%ExtraHangingSizes = a%ExtraHangingList%p(a%ExtraNpoin+1)-1
      a%ExtraLnodsSizes = a%ExtraElements%p(a%ExtraNelem+1)-1
      
   end subroutine
   
   subroutine EliminateRepeatedEntries(n0,iwa,nend)
      use typre
      implicit none
      integer(ip) :: n0,  nend
      integer(ip), allocatable :: iwa(:)
      
      integer(ip), allocatable :: iwa2(:)
      integer(ip) :: prepoint,ipoinGlobal
      integer(ip) :: ipoin
      
      allocate(iwa2(n0))
       !Eliminate repeated entries
      prepoint = -1
      nend = 0
      do ipoin = 1,n0
         ipoinGlobal = iwa(ipoin)
         if (ipoinGlobal /= prepoint) then
            nend = nend+1
            iwa2(nend) = ipoinGlobal
            prepoint = ipoinGlobal

            !a%ExtraNpoin = a%ExtraNpoin+1
         endif
      enddo
      deallocate(iwa)
      allocate(iwa(nend))
      iwa = iwa2(1:nend)
      deallocate(iwa2)
!       npoinCount = nend
   end subroutine
   
   subroutine LocalElementFromGEI(a,GEI,ielem)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      type(GlobalElementIdentifier) :: GEI
      integer(ip) :: ielem
      
      integer(ip) :: iTopParentGlobal, itopParentLocal,ilevel,ippos,ichild
      
      iTopParentGlobal = GEI%GlobalOriginalElement
               
      !change twice because auxLocal2GlobalElementList is not of nelem size, but of TopNelemSize
      call a%ElementLocalOrdering%Global2Local(iTopParentGlobal,itopParentLocal)
      
      !Recursively move through the list of elements till I find the element
      ielem = itopParentLocal
      ilevel = 1
      do while (ilevel <= GEI%level)
         call GEI_GetPositionInElementForLevel(GEI,ilevel,ichild)
         ippos = a%pPointerToChildren(ielem)-1
         ielem = a%lPointerToChildren(ippos+ichild)   
         
         ilevel = ilevel+1
      enddo
   end subroutine
   
   
   
!For Transfering to Mesh   
   subroutine GetPointDimensions(a,newnpoinLocal,newnpoin,newgnpoin)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: NewGnpoin,NewNpoin,NewNpoinLocal,newnelem
      
      newnpoinLocal = a%npoinLocal
      newnpoin = a%ExportNpoin+a%ExtraNpoin
      newgnpoin = a%gnpoin
   end subroutine
   
   subroutine GetLocalOrdering(a,LocalToGlobal,ProcessorList)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: LocalToGlobal(*),ProcessorList(*)
      
      integer(ip) :: ipoin, poini
         
         ProcessorList(1:a%npoinLocal) = a%MPIrank
         do ipoin = 1,a%npoinLocal
            LocalToGlobal(ipoin) = ipoin
         enddo
         call a%LocalOrdering%Local2Global(a%npoinLocal,LocalToGlobal(1:a%npoinLocal),LocalToGlobal(1:a%npoinLocal))
         
         do ipoin = a%npoinLocal+1,a%npoin
            poini = a%ExportPointNumbering(ipoin)
            if (poini /= 0) then
               call a%LocalOrdering%Local2Global(ipoin,LocalToGlobal(poini))
               call a%LocalOrdering%GetProc(ipoin,ProcessorList(poini))
            endif
         enddo
         !ExtraPoints
         do ipoin = 1,a%ExtraNpoin
            LocalToGlobal(a%ExportNpoin+ipoin) = a%ExtraPoints(1,ipoin)
            ProcessorList(a%ExportNpoin+ipoin) = a%ExtraPoints(2,ipoin)
         enddo
   end subroutine
   
   subroutine GetElementDimensions(a,newnelem,LnodsSize)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: newnelem,LnodsSize
      
      integer(ip) :: ielem, pnode
      
      !Account for extra elements
      newnelem = a%ExportNelem+a%ExtraNelem
      LnodsSize = a%ExportLnodsSizes+a%ExtraLnodsSizes
   end subroutine    
   
   subroutine GetLnods(a,pnods,lnods)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: NewNelem
      integer(ip) :: pnods(*), lnods(*)
      
      integer(ip) :: ielem,ElementCount, LnodsCount,ispos,pnode,isposnew
      
      !I only export elements which are marked as for export
      pnods(1) = 1
      ElementCount = 0
      do ielem = 1,a%nelem
         !only last level elements to export
         if (a%ExportElement(ielem) .eqv. .true.) then
            ElementCount = ElementCount + 1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            
            isposnew = pnods(ElementCount) -1
            ispos = a%pnods(ielem)-1
            
            pnods(ElementCount+1) = pnods(ElementCount)+pnode
            lnods(isposnew+1:isposnew+pnode) = a%ExportPointNumbering(a%lnods(ispos+1:ispos+pnode))
         endif
      enddo
      !Also export additional elements
      do ielem = 1,a%ExtraNelem
         ElementCount = ElementCount+1
         pnode = a%ExtraElements%p(ielem+1)-a%ExtraElements%p(ielem)
         
         isposnew = pnods(ElementCount)-1
         ispos = a%ExtraElements%p(ielem)-1
         
         pnods(ElementCount+1) = pnods(ElementCount)+pnode
         lnods(isposnew+1:isposnew+pnode) = a%ExtraElements%l(ispos+1:ispos+pnode)
      enddo
      
   end subroutine
   
   subroutine GetHangingListDimensions(a,HangingListDimensions)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: HangingListDimensions
      
      HangingListDimensions = a%ExportHangingSizes + a%ExtraHangingSizes
   end subroutine
      
   subroutine GetHangingList(a,pHangingList,lHangingList,rHangingList)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: pHangingList(*),lHangingList(*)
      real(rp) :: rHangingList(*)
      
      integer(ip) :: ihpos,ipoin,phang,poini
      integer(ip) :: LocalHangingSizes,ExportHangingSizes
      
      !Local nodes (all of them)
      LocalHangingSizes = a%pHangingList(a%npoinlocal+1)-1
      pHangingList(1:a%npoinLocal+1) = a%pHangingList(1:a%npoinLocal+1)

      lHangingList(1:LocalHangingSizes) = a%ExportPointNumbering(a%lHangingList(1:LocalHangingSizes))
      rHangingList(1:LocalHangingSizes) = a%rHangingList(1:LocalHangingSizes)
      
      !Ghost nodes (only those selected)
      ExportHangingSizes = LocalHangingSizes
      do ipoin = a%npoinLocal+1,a%npoin
         poini = a%ExportPointNumbering(ipoin)
         if (poini /= 0) then
            ihpos = a%pHangingList(ipoin)-1
            phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
            pHangingList(poini+1) = ExportHangingSizes+1+phang
            
            lHangingList(ExportHangingSizes+1:ExportHangingSizes+phang) = a%ExportPointNumbering(a%lHangingList(ihpos+1:ihpos+phang))
            rHangingList(ExportHangingSizes+1:ExportHangingSizes+phang) = a%rHangingList(ihpos+1:ihpos+phang)
            
            !Update next starting position
            ExportHangingSizes = ExportHangingSizes+phang
         endif
      enddo
      !ExtraNodes
      do ipoin = 1,a%ExtraNpoin
         ihpos = a%ExtraHangingList%p(ipoin)-1
         phang = a%ExtraHangingList%p(ipoin+1)-a%ExtraHangingList%p(ipoin) 
         
         pHangingList(a%ExportNpoin+ipoin+1) = pHangingList(a%ExportNpoin+ipoin)+phang
         
         lHangingList(ExportHangingSizes+1:ExportHangingSizes+phang) = a%ExtraHangingList%l(ihpos+1:ihpos+phang)
         rHangingList(ExportHangingSizes+1:ExportHangingSizes+phang) = a%ExtraHangingList%r(ihpos+1:ihpos+phang)
         
         !Update next starting position
         ExportHangingSizes = ExportHangingSizes+phang
      enddo   
      
   end subroutine
   
   subroutine GetHangingFacesList(a,pHangingList,lHangingList)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: pHangingList(*),lHangingList(2,*)
      
      integer(ip) :: ielem,iexport, nfacs, istart, iend, fstart,fend
      
      pHangingList(1) = 1
      do ielem = 1,a%Nelem
         iexport = a%ExportElementNumbering(ielem)
         if (iexport > 0) then
            nfacs = a%pfacs(ielem+1)-a%pfacs(ielem)
            istart = pHangingList(iexport)
            iend = istart + nfacs
            PHangingList(iexport+1) = iend
            
            fstart = a%pfacs(ielem)
            fend   = a%pfacs(ielem+1)
            lHangingList(2,istart:iend) = a%lfacs(2,fstart:fend)
            lHangingList(1,istart:iend) = a%ExportElementNumbering(abs(a%lfacs(1,fstart:fend)))
         
         
         endif
      enddo
      
      call runend('testing pending')
   
   
   
   end subroutine
      
   subroutine UpdateVariableReal(a,ndime,coord,newcoord)
      use typre
      use Mod_AdaptiveInterpolationList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      real(rp) :: coord(ndime,*)
      real(rp) :: newcoord(ndime,*)
      
      logical, allocatable :: isSet(:)

      integer(ip) :: nunrefpo, norigpo,ipoin,jv,jpoin,newpoin,json,oldipoin
      integer(ip) :: oldExportIpoin, ExportIpoin,ExportJpoin,ExportJson,oldExportJpoin,oldJpoin
      
      type(InterpolationList) :: InterpList
      integer(ip), pointer :: jlnode(:) => NULL(), ilnode(:) => NULL()
      integer(ip) :: elemi,ielem,ielty,ivariation,inode,ipnode,jchild,jelem,jelemspos
      integer(ip) :: jnode,jpnode,knode,kpoin,nodek,OldExportKpoin,oldKpoin
      
      !reallocate coord taking into account unrefined edges
      !allocate(newcoord(ndime,a%MeshNpoin))
      allocate(isSet(a%Exportnpoin))
      isSet = .false.
      
      do oldipoin = 1,a%oldnpoin
         oldExportIpoin = a%OldExportPointNumbering(oldipoin)
         if (oldExportIpoin /= 0) then
            ipoin = a%renpo(oldipoin)
            if (ipoin /= 0) then
               ExportIpoin = a%ExportPointNumbering(ipoin)
               if (ExportIpoin /= 0) then
                  newcoord(1:ndime,ExportIpoin) = coord(:,oldExportIpoin)
                  isSet(Exportipoin) = .true.
               endif
            endif
         endif
      enddo

      
      !I loop through all the edges and fill rtemp
      do ipoin = 1,a%npoin
         oldIpoin = a%irenpo(ipoin)
         if (oldIpoin /= 0) then
            OldExportIpoin = a%OldExportPointNumbering(oldIpoin)
            if (oldExportIpoin /= 0) then
               do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
                  jpoin = a%ledger(1,jv)
                  oldJpoin = a%irenpo(jpoin)
                  if (oldJpoin /= 0) then
                     OldExportJpoin = a%OldExportPointNumbering(OldJpoin)
                     if (OldExportJpoin /= 0) then
                        json = a%ledger(2,jv)
                        ExportJson= a%ExportPointNumbering(json)
                        if (ExportJson /= 0) then
                           if (isSet(Exportjson) .eqv. .false.) then
                              isSet(Exportjson) = .true.
                              newcoord(:,ExportJson) = (coord(:,OldExportIpoin)+coord(:,OldExportJpoin))*0.5
                           endif   
                        endif
                     endif
                  endif
               enddo
            endif
         endif
      enddo
      
      !Interpolations done through the elements
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         ielty = a%ElementTypeList(ielem)%ielty
         ivariation = a%ElementTypeList(ielem)%ivariation
         call Adaptive_GetInterpolationList(ielty,ivariation,InterpList)
         !Loop through all points to interpolate
         do inode = 1,InterpList%n
            jchild = InterpList%childElementList(inode)
            jelemspos = a%pPointerToChildren(ielem)-1
            jelem = a%lPointerToChildren(jelemspos+jchild)
            
            call a%GetLnodsFromElement(jelem,jpnode,jlnode)
            jnode = InterpList%ChildObjectiveNodeList(inode)
            jpoin = jlnode(jnode)
            
            ExportJpoin = a%ExportPointNumbering(jpoin)
            if (ExportJpoin /= 0) then
               if (isSet(ExportJpoin) .eqv. .false.) then
                  !The node to interpolate is still not set
                  !Check that I have the info of all the nodes required to interpolate
                  
                  !mark it, we will unmark it later if the element cannot interpolate
                  isSet(ExportJpoin) = .true.
                  
                  newcoord(:,ExportJpoin) = 0.0_rp
                  call a%GetLnodsFromElement(ielem,ipnode,ilnode)
                  knodeloop: do knode = 1,InterpList%nInterpolatingNodes(inode)
                     nodek = InterpList%InterpolatingNodesList(knode,inode)
                     kpoin = ilnode(nodek)
                     
                     oldKpoin = a%irenpo(kpoin)
                     if (oldKpoin /= 0) then
                        OldExportKpoin = a%OldExportPointNumbering(oldKpoin)
                        if (oldExportKpoin /= 0) then
                           newcoord(:,ExportJpoin) = newcoord(:,ExportJpoin) + coord(:,OldExportKpoin)*InterpList%InterpolatingCoefficientsList(knode,inode)
                        else
                           !we were not ready for interpolating
                           
                           isSet(ExportJpoin) = .false.
                           exit knodeloop
                        endif
                     endif
                  enddo knodeloop
               endif
            endif
         enddo
      enddo
      
      !Ghost communications, only required for certain unrefinement cases
      !but I don't know if the other subdomains will need it, so I have to do it always
      
      !Also, always required for extra fields
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:a%ExportNpoin))
      
      
      deallocate(isSet)
   end subroutine
   
   subroutine UpdateVariableReal1(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      real(rp), target :: coord(*)
      real(rp), target :: newcoord(*)
      
      real(rp), pointer :: auxcoord(:,:) => NULL()
      real(rp), pointer :: auxnewcoord(:,:) => NULL()
      
      auxcoord(1:ndime,1:a%oldnpoin) => coord(1:ndime*a%oldnpoin)
      auxnewcoord(1:ndime,1:a%npoin) => newcoord(1:ndime*a%npoin)
      
      call a%UpdateVariableReal(ndime,auxcoord,auxnewcoord)
   end subroutine
   
   subroutine UpdateVariableInteger(a,ndime,coord,newcoord,maxmin)
      use typre
      use Mod_AdaptiveInterpolationList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      integer(ip) :: coord(ndime,*)
      integer(ip) :: newcoord(ndime,*)
      character(5) :: maxmin
      
      integer(ip) :: multiplier
      
      logical, allocatable :: isSet(:)

      integer(ip) :: nunrefpo, norigpo,ipoin,jv,jpoin,newpoin,json,oldipoin
      integer(ip) :: oldExportIpoin, ExportIpoin,ExportJpoin,ExportJson,oldExportJpoin,oldJpoin
      
      type(InterpolationList) :: InterpList
      integer(ip), pointer :: jlnode(:) => NULL(), ilnode(:) => NULL()
      integer(ip) :: elemi,ielem,ielty,ivariation,inode,ipnode,jchild,jelem,jelemspos
      integer(ip) :: jnode,jpnode,knode,kpoin,nodek,OldExportKpoin,oldKpoin,idime
      
      !reallocate coord taking into account unrefined edges
      !allocate(newcoord(ndime,a%MeshNpoin))
      allocate(isSet(a%ExportNpoin))
      isSet = .false.
      
      do oldipoin = 1,a%oldnpoin
         oldExportIpoin = a%OldExportPointNumbering(oldipoin)
         if (oldExportIpoin /= 0) then
            ipoin = a%renpo(oldipoin)
            if (ipoin /= 0) then
               ExportIpoin = a%ExportPointNumbering(ipoin)
               if (ExportIpoin /= 0) then
                  newcoord(1:ndime,ExportIpoin) = coord(:,oldExportIpoin)
                  isSet(Exportipoin) = .true.
               endif
            endif
         endif
      enddo
      
      if (maxmin == 'maxim') then
         multiplier = 1_ip
      elseif (maxmin == 'minim') then
         multiplier = -1_ip
      elseif (maxmin == 'minab') then
         multiplier = -2
      elseif (maxmin == 'maxab') then
         multiplier = 2
      endif
      
      !I loop through all the edges and fill rtemp
      do ipoin = 1,a%npoin
         oldIpoin = a%irenpo(ipoin)
         if (oldIpoin /= 0) then
            OldExportIpoin = a%OldExportPointNumbering(oldIpoin)
            if (oldExportIpoin /= 0) then
               do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
                  jpoin = a%ledger(1,jv)
                  oldJpoin = a%irenpo(jpoin)
                  if (oldJpoin /= 0) then
                     OldExportJpoin = a%OldExportPointNumbering(OldJpoin)
                     if (OldExportJpoin /= 0) then
                        json = a%ledger(2,jv)
                        ExportJson= a%ExportPointNumbering(json)
                        if (ExportJson /= 0) then   
                           if (isSet(Exportjson) .eqv. .false.) then
                              isSet(Exportjson) = .true.
                              if (multiplier == 1 .or. multiplier == -1) then
                                 newcoord(:,ExportJson) = max(coord(:,OldExportIpoin)*multiplier,coord(:,OldExportJpoin)*multiplier)*multiplier
                              elseif (multiplier == -2 .or. multiplier == 2) then
                                 do idime = 1,ndime
                                    if (abs(coord(idime,OldExportIpoin))*sign(1,multiplier) > abs(coord(idime,OldExportJpoin))*sign(1,multiplier)) then
                                       newcoord(idime,ExportJson) = coord(idime,OldExportIpoin)
                                    else
                                       newcoord(idime,ExportJson) = coord(idime,OldExportJpoin)
                                    endif
                                 enddo
                              endif
                                 
                           endif   
                        endif
                     endif
                  endif
               enddo
            endif
         endif
      enddo
      
      !Interpolations done through the elements
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         ielty = a%ElementTypeList(ielem)%ielty
         ivariation = a%ElementTypeList(ielem)%ivariation
         call Adaptive_GetInterpolationList(ielty,ivariation,InterpList)
         !Loop through all points to interpolate
         do inode = 1,InterpList%n
            jchild = InterpList%childElementList(inode)
            jelemspos = a%pPointerToChildren(ielem)-1
            jelem = a%lPointerToChildren(jelemspos+jchild)
            
            call a%GetLnodsFromElement(jelem,jpnode,jlnode)
            jnode = InterpList%ChildObjectiveNodeList(inode)
            jpoin = jlnode(jnode)
            
            ExportJpoin = a%ExportPointNumbering(jpoin)
            if (ExportJpoin /= 0) then
               if (isSet(ExportJpoin) .eqv. .false.) then
                  !The node to interpolate is still not set
                  !Check that I have the info of all the nodes required to interpolate
                  
                  !mark it, we will unmark it later if the element cannot interpolate
                  isSet(ExportJpoin) = .true.
                  
                  newcoord(:,ExportJpoin) = -multiplier*1e9
                  call a%GetLnodsFromElement(ielem,ipnode,ilnode)
                  knodeloop: do knode = 1,InterpList%nInterpolatingNodes(inode)
                     nodek = InterpList%InterpolatingNodesList(knode,inode)
                     kpoin = ilnode(nodek)
                     
                     oldKpoin = a%irenpo(kpoin)
                     if (oldKpoin /= 0) then
                        OldExportKpoin = a%OldExportPointNumbering(oldKpoin)
                        if (oldExportKpoin /= 0) then
                           newcoord(:,ExportJpoin) = max(newcoord(:,ExportJpoin)*multiplier,coord(:,OldExportKpoin)*multiplier)*multiplier
                        else
                           !we were not ready for interpolating
                           isSet(ExportJpoin) = .false.
                           exit knodeloop
                        endif
                     endif
                  enddo knodeloop
               endif
            endif
         enddo
      enddo
      
      !Ghost communications, only required for certain unrefinement cases
      !but I don't know if the other subdomains will need it, so I have to do it always
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:a%ExportNpoin))
      
      deallocate(isSet)
   end subroutine
   
   subroutine UpdateVariableInteger1(a,ndime,coord,newcoord,maxmin)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      integer(ip), target :: coord(*)
      integer(ip), target :: newcoord(*)
      character(5) :: maxmin
      
      integer(ip), pointer :: auxcoord(:,:) => NULL()
      integer(ip), pointer :: auxnewcoord(:,:) => NULL()
      
      auxcoord(1:ndime,1:a%oldnpoin) => coord(1:ndime*a%oldnpoin)
      auxnewcoord(1:ndime,1:a%npoin) => newcoord(1:ndime*a%npoin)
      
      call a%UpdateVariableInteger(ndime,auxcoord,auxnewcoord,maxmin)
   end subroutine
   
   subroutine UpdateVariableLogical(a,ndime,coord,newcoord,maxmin)
      use typre
      use Mod_AdaptiveInterpolationList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      logical :: coord(ndime,*)
      logical :: newcoord(ndime,*)
      character(5) :: maxmin
      
      integer(ip) :: multiplier
      
      logical, allocatable :: isSet(:)

      integer(ip) :: nunrefpo, norigpo,ipoin,jv,jpoin,newpoin,json,oldipoin
      integer(ip) :: oldExportIpoin, ExportIpoin,ExportJpoin,ExportJson,oldExportJpoin,oldJpoin
      integer(ip) :: IntCoordIpoin(ndime), IntCoordJpoin(ndime), IntCoordJson(ndime),IntCoordKpoin(ndime)
      
      type(InterpolationList) :: InterpList
      integer(ip), pointer :: jlnode(:) => NULL(), ilnode(:) => NULL()
      integer(ip) :: elemi,ielem,ielty,ivariation,inode,ipnode,jchild,jelem,jelemspos
      integer(ip) :: jnode,jpnode,knode,kpoin,nodek,OldExportKpoin,oldKpoin
      
      !reallocate coord taking into account unrefined edges
      !allocate(newcoord(ndime,a%MeshNpoin))
      allocate(isSet(a%ExportNpoin))
      isSet = .false.
      
      do oldipoin = 1,a%oldnpoin
         oldExportIpoin = a%OldExportPointNumbering(oldipoin)
         if (oldExportIpoin /= 0) then
            ipoin = a%renpo(oldipoin)
            if (ipoin /= 0) then
               ExportIpoin = a%ExportPointNumbering(ipoin)
               if (ExportIpoin /= 0) then
                  newcoord(1:ndime,ExportIpoin) = coord(:,oldExportIpoin)
                  isSet(Exportipoin) = .true.
               endif
            endif
         endif
      enddo
      
      if (maxmin == 'maxim') then
         multiplier = 1_ip
      elseif (maxmin == 'minim') then
         multiplier = -1_ip
      endif
      
      !I loop through all the edges and fill rtemp
      do ipoin = 1,a%npoin
         oldIpoin = a%irenpo(ipoin)
         if (oldIpoin /= 0) then
            OldExportIpoin = a%OldExportPointNumbering(oldIpoin)
            if (oldExportIpoin /= 0) then
               do jv = a%pedger(ipoin),a%pedger(ipoin+1)-1
                  jpoin = a%ledger(1,jv)
                  oldJpoin = a%irenpo(jpoin)
                  if (oldJpoin /= 0) then
                     OldExportJpoin = a%OldExportPointNumbering(OldJpoin)
                     if (OldExportJpoin /= 0) then
                        json = a%ledger(2,jv)
                        ExportJson= a%ExportPointNumbering(json)
                        if (ExportJson /= 0) then   
                           if (isSet(Exportjson) .eqv. .false.) then
                              isSet(Exportjson) = .true.
                              call Logical2Integer(ndime,coord(:,OldExportIpoin),IntCoordIpoin)
                              call Logical2Integer(ndime,coord(:,OldExportJpoin),IntCoordJpoin)
                              
                               IntCoordJson= max(IntCoordIpoin*multiplier,IntCoordJpoin*multiplier)*multiplier
                               call Integer2Logical(ndime,IntCoordJson,newcoord(:,ExportJson))
                           endif   
                        endif
                     endif
                  endif
               enddo
            endif
         endif
      enddo
      
      !Interpolations done through the elements
      do elemi = 1,a%nJustRefinedElements
         ielem = a%lJustRefinedElements(elemi)
         
         ielty = a%ElementTypeList(ielem)%ielty
         ivariation = a%ElementTypeList(ielem)%ivariation
         call Adaptive_GetInterpolationList(ielty,ivariation,InterpList)
         !Loop through all points to interpolate
         do inode = 1,InterpList%n
            jchild = InterpList%childElementList(inode)
            jelemspos = a%pPointerToChildren(ielem)-1
            jelem = a%lPointerToChildren(jelemspos+jchild)
            
            call a%GetLnodsFromElement(jelem,jpnode,jlnode)
            jnode = InterpList%ChildObjectiveNodeList(inode)
            jpoin = jlnode(jnode)
            
            ExportJpoin = a%ExportPointNumbering(jpoin)
            if (ExportJpoin /= 0) then
               if (isSet(ExportJpoin) .eqv. .false.) then
                  !The node to interpolate is still not set
                  !Check that I have the info of all the nodes required to interpolate
                  
                  !mark it, we will unmark it later if the element cannot interpolate
                  isSet(ExportJpoin) = .true.
                  
                  IntCoordJpoin = -multiplier
                  call a%GetLnodsFromElement(ielem,ipnode,ilnode)
                  knodeloop: do knode = 1,InterpList%nInterpolatingNodes(inode)
                     nodek = InterpList%InterpolatingNodesList(knode,inode)
                     kpoin = ilnode(nodek)
                     
                     oldKpoin = a%irenpo(kpoin)
                     if (oldKpoin /= 0) then
                        OldExportKpoin = a%OldExportPointNumbering(oldKpoin)
                        if (oldExportKpoin /= 0) then
                           call Logical2Integer(ndime,coord(:,OldExportKpoin),IntCoordKpoin)
                           IntCoordJpoin = max(IntCoordJpoin*multiplier,IntCoordKpoin*multiplier)*multiplier
                        else
                           !we were not ready for interpolating
                           isSet(ExportJpoin) = .false.
                           exit knodeloop
                        endif
                     endif
                  enddo knodeloop
                  call Integer2Logical(ndime,IntCoordJpoin,newcoord(:,ExportJpoin))
               endif
            endif
         enddo
      enddo
      
      !Ghost communications, only required for certain unrefinement cases
      !but I don't know if the other subdomains will need it, so I have to do it always
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:a%ExportNpoin))
      
      deallocate(isSet)
   end subroutine
   
   subroutine UpdateVariableLogical1(a,ndime,coord,newcoord,maxmin)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      logical, target :: coord(*)
      logical, target :: newcoord(*)
      character(5) :: maxmin
      
      logical, pointer :: auxcoord(:,:) => NULL()
      logical, pointer :: auxnewcoord(:,:) => NULL()
      
      auxcoord(1:ndime,1:a%oldnpoin) => coord(1:ndime*a%oldnpoin)
      auxnewcoord(1:ndime,1:a%npoin) => newcoord(1:ndime*a%npoin)
      
      call a%UpdateVariableLogical(ndime,auxcoord,auxnewcoord,maxmin)
   end subroutine
      
   subroutine GetLevel(a,level)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: level(*)
      
      integer(ip) :: ipoin, poini
         
         do ipoin = 1,a%npoinLocal
            level(ipoin) = a%level(ipoin)
         enddo
         
         do ipoin = a%npoinLocal+1,a%npoin
            poini = a%ExportPointNumbering(ipoin)
            if (poini /= 0) then
               level(poini) = a%level(ipoin)
            endif
         enddo
   end subroutine
      
!For load rebalancing

   subroutine LoadRebalance(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      call a%Timers%LoadRebalance%Tic
      
      call a%RebalanceRenumbering
      call a%RebalanceElemProcs
      call a%PrepareRebalanceFaces_a
      call a%RebalanceSendNodesAndElements
      call a%RebuildRefinementStructure
      call a%PrepareRebalanceFaces_b
      call a%PrepareRebalanceCommunicator
      
      call a%Timers%LoadRebalance%Toc
      
   end subroutine

   subroutine RebalanceRenumbering(a)
      use typre
      use Mod_NaivePartition
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      type(NaivePartition) :: NaivePartitioner
      integer(ip) :: ipoin,gipoin,gipoinmax,i,irank
      
      call a%Timers%RebalanceRenumbering%Tic
      
      !In this subroutine we simply decide the new numbering and processor for each node
      
      !First, store the old number of npoinLocal
      a%RebalanceOldNpoinLocal = a%npoinLocal
      a%RebalanceOldNelem = a%nelem
      
      
      if (a%kfl_RebalanceAlgorithm == 0) then
         if (allocated(a%RebalancePointNumbering)) then
            deallocate(a%RebalancePointNumbering)
            deallocate(a%RebalanceProcessorList)
         endif
         
         allocate(a%RebalancePointNumbering(a%npoin))
         allocate(a%RebalanceProcessorList(a%npoin))
         
         !No options now, just a naive partitioner. 
         !This is a general case, so if I work with this one 
         !other cases should also work
         !Not so general: global point numbering does not change
         do ipoin = 1,a%npoinLocal
            a%RebalancePointNumbering(ipoin) = ipoin
         enddo
         call a%LocalOrdering%Local2Global(a%npoinLocal,a%RebalancePointNumbering(1:a%npoinLocal),a%RebalancePointNumbering(1:a%npoinLocal))
         a%RebalanceProcessorList = -20
         a%RebalanceProcessorList(1:a%npoinLocal) = a%MPIrank
         a%RebalanceProcessorList(a%npoinLocal+1:a%npoin) = -1
         
         !Assign a global number and processor through a partitioning process to each node
      
         call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
         
         !I Invert the global point numbering to make it general
         if (a%npoinLocal > 0) then
            gipoinmax = a%RebalancePointNumbering(a%npoinLocal)
         else
            gipoinmax = 0
         endif
         do ipoin = 1,a%npoinLocal
            a%RebalancePointNumbering(ipoin) = gipoinmax - ipoin +1
         enddo
         
         do ipoin = 1,a%npoinLocal
            gipoin = a%RebalancePointNumbering(ipoin)
            call NaivePartitioner%GetPointRank(gipoin,a%RebalanceProcessorList(ipoin))
         enddo
         
         !GhostCommunicate new global numbering and processors for all ghost nodes
         call a%Communicator%GhostCommunicate(1_ip,a%RebalancePointNumbering)
         call a%Communicator%GhostCommunicate(1_ip,a%RebalanceProcessorList)
         
      elseif (a%kfl_RebalanceAlgorithm == 1) then
         !Do nothing, rebalancing Numbering has been set externally
         
      endif
      
      call a%Timers%RebalanceRenumbering%Toc
      
   end subroutine
   
   subroutine RebalanceElemProcs(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      call a%Timers%RebalanceElemProcs%Tic
      
      call RebalanceElemProcs_a(a)
      call RebalanceElemProcs_b(a)
      
      call a%Timers%RebalanceElemProcs%Toc
      
   end subroutine   
   
   subroutine RebalanceElemProcs_a(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip), allocatable :: ichild(:), seen(:,:), iwa(:,:), nz(:),ielem(:)
      integer(ip), allocatable :: nchild(:),iepos(:)
      
      integer(ip) :: elemi,ilevel,ispos,pnode,inode,ipoin,irank,ispos2,pnode2
      integer(ip) :: siblingi, isibling,iz
      logical :: kfl_continue, kfl_lastsibling
      
      !In this subroutine we choose to which processors every element needs to be sent
      !we need to guarantee that if an element is sent, its siblings and ancestors
      !(and the ancestors' siblings) also need to be sent

      
      !We know where each node goes, ghost communicates will do the rest
      
      !For elements it is more involved, we need the whole hierarchy!!
      
      !Fer primer sense hierarchy, nomes d'acord amb nodes
      !Nomes per a elements I am responsible for
      !pero marcar d'alguna manera els ranks dels nodes de cada element encara 
      !que no en sigui responsable
      !per a així saber almenys si no cal enviar-los.
      
      !Depart from elements I am responsible and move up till I find someone
      !with the rank I want to transmit. In the meanwhile mark everybody for the rank
      
      
      !This is extra expensive
      !I need to decompress at each element check if the rank is already there then add
      
      !After I know which element goes where
      !Send the list of ranks to the responsible of the element
      
      !Then the responsible sends the element to the merged list of ranks
      
      !Can be done better in a recursive way?
      !move through the tree
      !have a list for each level (affordable in terms of memory: 32 levels, 10000 processors: 1.2 MBytes)
      !iwa and seen for each level
      !Fins i tot si un element es 100% local, els seus fills no han de ser-ho per força
      !(despres de multiples load rebalances)
      !Per tant cal enviar tots els ranks excepte els dels seus nodes segur al processador responsable
      !Ell s'encarregarà de fer el merge i enviar l'element a tots els processadors que corresponguin
      
      !un cop enviats els elements puc construir els ghosts nodes
      
      !Després cal reconstruir-ho tot a partir del GEI
      if (allocated(a%RebalanceElementProcessors)) then
         do irank = 0,size(a%RebalanceElementProcessors)-1
            deallocate(a%RebalanceElementProcessors(irank+1)%l)
         enddo
         deallocate(a%RebalanceElementProcessors)
      endif
      allocate(a%RebalanceElementProcessors(a%nelem))
      
      allocate(ichild(GEI_ElementMaxLevels))
      allocate(nchild(GEI_ElementMaxLevels))
      ichild = 0
      allocate(seen(a%MPIsize,GEI_ElementMaxLevels))
      allocate(iwa(a%MPIsize,GEI_ElementMaxLevels))
      allocate(nz(GEI_ElementMaxLevels))
      seen = 0
      nz = 0
      
      allocate(ielem(GEI_ElementMaxLevels))
      allocate(iepos(GEI_ElementMaxLevels))

      do elemi = 1,a%nelem
         ielem(1) = elemi
         !Enter through top level elements, move through all of them recursively
         if (a%GEI_ElementList(ielem(1))%level == 0) then
            ilevel = 0
            nchild(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1)+1)-a%pPointerToChildren(ielem(ilevel+1))
            iepos(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1))-1
            ichild(ilevel+1) = 0
            kfl_continue = .true.
                        
            do while (kfl_continue .eqv. .true.)
               
               if (ichild(ilevel+1) < nchild(ilevel+1)) then
                  !Move down to the next children
                  ichild(ilevel+1) = ichild(ilevel+1)+1
                  !Move down
                  ilevel = ilevel + 1
                  ielem(ilevel+1) = a%lPointerToChildren(iepos(ilevel)+ichild(ilevel))
                  
                  iepos(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1))-1
                  nchild(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1)+1)-a%pPointerToChildren(ielem(ilevel+1))
                  ichild(ilevel+1) = 0
                  
               else
                  !Mark the new ranks of the nodes for the level
                  ispos = a%pnods(ielem(ilevel+1))-1
                  pnode = a%pnods(ielem(ilevel+1)+1)-a%pnods(ielem(ilevel+1))
                  do inode = 1,pnode
                     ipoin = a%lnods(ispos+inode)
                     irank = a%RebalanceProcessorList(ipoin)
                     
                     if (seen(irank+1,ilevel+1) == 0) then
                        nz(ilevel+1) = nz(ilevel+1)+1
                        seen(irank+1,ilevel+1) = 1
                        iwa(nz(ilevel+1),ilevel+1) = irank
                     endif
                  enddo
                  
                  !Do the same for the hanging nodes of which I am high level element
                  pnode2 = a%HangingListForElements%p(ielem(ilevel+1)+1)-a%HangingListForElements%p(ielem(ilevel+1))
                  if (pnode2 /= 0) then
                     ispos2 = a%HangingListForElements%p(ielem(ilevel+1))-1
                     do inode = 1,pnode2
                        ipoin = a%HangingListForElements%l(ispos2+inode)
                        irank = a%RebalanceProcessorList(ipoin)
                        
                        if (seen(irank+1,ilevel+1) == 0) then
                           nz(ilevel+1) = nz(ilevel+1)+1
                           seen(irank+1,ilevel+1) = 1
                           iwa(nz(ilevel+1),ilevel+1) = irank
                        endif
                     enddo
                  endif
                  
                  !Top level, just copy and exit
                  if (ilevel == 0) then
                     !Copy data
                     allocate(a%RebalanceElementProcessors(ielem(ilevel+1))%l(nz(ilevel+1)))
                     a%RebalanceElementProcessors(ielem(ilevel+1))%l = iwa(1:nz(ilevel+1),ilevel+1)
                     
                     !Clean the level
                     do iz = 1,nz(ilevel+1)
                        irank = iwa(iz,ilevel+1)
                        
                        !Clean
                        seen(irank+1,ilevel+1) = 0
                     enddo
                     !Clean the current level
                     nz(ilevel+1) = 0
                     
                     !We are done
                     ilevel = ilevel - 1
                     kfl_continue = .false.
                  !Lower levels
                  else
                     !If this is the last sibling, then do all the siblings of the level
                     if (ichild(ilevel) == nchild(ilevel)) then
                        
                        !Loop through my parents children and fill them             
                        do siblingi = 1,nchild(ilevel)
                           isibling = a%lPointerToChildren(iepos(ilevel)+siblingi)
                           
                           allocate(a%RebalanceElementProcessors(isibling)%l(nz(ilevel+1)))
                           a%RebalanceElementProcessors(isibling)%l = iwa(1:nz(ilevel+1),ilevel+1)
                        enddo
                     
                        !Copy info to the upper level and clean the current level
                        do iz = 1,nz(ilevel+1)
                           irank = iwa(iz,ilevel+1)
                           if (seen(irank+1,ilevel) == 0) then
                              nz(ilevel) = nz(ilevel)+1
                              seen(irank+1,ilevel) = 1
                              iwa(nz(ilevel),ilevel) = irank
                           endif
                           
                           !Clean
                           seen(irank+1,ilevel+1) = 0
                        enddo
                        !Clean the current level
                        nz(ilevel+1) = 0
                     endif
                     
                     !And move up
                     ilevel = ilevel - 1
                  endif
               endif
            
            enddo
            
         endif
      enddo
      
      !After this I should know where every element needs to be sent
      !But only locally, I need to send the info about the elements I am not responsible of 
      !to their processors
   end subroutine
   
   subroutine RebalanceElemProcs_b(a)
      use typre
      use Mod_LinkedList
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      
      integer(ip) :: glnods(MaxNodesPerElement), glprocs(MaxNodesPerElement)
      integer(ip), pointer :: lnods(:) => NULL()
      
      integer(ip), allocatable :: NelemOthersNeed(:), NelemINeed(:),SizesElementsOthersNeed(:),SizesElementsINeed(:)
      type(CSRList), allocatable :: ElementsOthersNeed(:), ElementsINeed(:)
      
      type(ListGEI), allocatable :: GEI_ElementsINeed(:),GEI_ElementsOthersNeed(:)
      integer(ip) :: GEI_Size
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:),irequest05(:),irequest06(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      integer(ip) :: irank,jrank
      
      type(LinkedListHeader), allocatable :: HeaderElementRanks(:)
      type(LinkedListStorage)             :: StorageElementRanks
      
      integer(ip), allocatable :: seen(:), iwa(:)
      integer(ip) :: elemi,ielem,irpos,ispos,iz,minrank,nz,pnode,prank
      logical :: IsResponsible
      type(GlobalElementIdentifier) :: GEI_Sizer
      
      !Now I know locally where every element needs to go, I need to merge the info from all processors
      
      !Allocate and Initialize some MPI
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      allocate(irequest05(a%MPIsize))
      allocate(irequest06(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL      
      irequest05 = MPI_REQUEST_NULL
      irequest06 = MPI_REQUEST_NULL      

      
      !For counting sizes
      allocate(NelemOthersNeed(a%MPIsize))
      allocate(SizesElementsOthersNeed(a%MPIsize))
      NelemOthersNeed = 0
      SizesElementsOthersNeed = 0
            
      !First to count
      do ielem = 1,a%nelem
         !Refined elements only
         !if (a%markel(ielem) == 1) then
            !If I am not responsible for the element, send the info
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            lnods => a%lnods(ispos+1:ispos+pnode)
            call a%LocalOrdering%GetProc(pnode,lnods,glprocs)
            
            !I am not responsible, list must be increased
            minrank = minval(glprocs(1:pnode))
            if (minrank /= a%MPIrank) then
               NelemOthersNeed(minrank+1) = NelemOthersNeed(minrank+1) + 1
               SizesElementsOthersNeed(minrank+1) = SizesElementsOthersNeed(minrank+1) + size(a%RebalanceElementProcessors(ielem)%l)
            endif
         !endif
      enddo
      !write(*,*) 'not refined elements (markel = 0) should not be sent'
      !write(*,*) 'elements I am not responsible, but of which I do not own any node, should not be sent either'
      !write(*,*) 'the info should be rebuilt locally from their siblings info'
      !write(*,*) 'needs to be refactored if it affects performance'
      
      !Allocate Communication Structures
      allocate(ElementsOthersNeed(a%MPIsize))
      allocate(GEI_ElementsOthersNeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1)/= 0) then
            allocate(GEI_ElementsOthersNeed(irank+1)%a(NelemOthersNeed(irank+1)))
            allocate(ElementsOthersNeed(irank+1)%p(NelemOthersNeed(irank+1)+1))
            ElementsOthersNeed(irank+1)%p(1) = 1
            allocate(ElementsOthersNeed(irank+1)%l(SizesElementsOthersNeed(irank+1)))
         endif
      enddo
      
      NelemOthersNeed = 0
      !Second to fill
      do ielem = 1,a%nelem
         !Refined elements only
         !if (a%markel(ielem) == 1) then
            !If I am not responsible for the element, send the info
            ispos = a%pnods(ielem)-1
            pnode = a%pnods(ielem+1)-a%pnods(ielem)
            lnods => a%lnods(ispos+1:ispos+pnode)
            call a%LocalOrdering%GetProc(pnode,lnods,glprocs)
            
            !I am not responsible, list must be increased
            minrank = minval(glprocs(1:pnode))
            if (minrank /= a%MPIrank) then
               NelemOthersNeed(minrank+1) = NelemOthersNeed(minrank+1) + 1
               iz = NelemOthersNeed(minrank+1)
               irpos = ElementsOthersNeed(minrank+1)%p(iz)-1
               prank = size(a%RebalanceElementProcessors(ielem)%l)
               ElementsOthersNeed(minrank+1)%p(iz + 1) = ElementsOthersNeed(minrank+1)%p(iz) + prank
               ElementsOthersNeed(minrank+1)%l(irpos+1:irpos+prank) = a%RebalanceElementProcessors(ielem)%l
               
               GEI_ElementsOthersNeed(minrank+1)%a(iz) = a%GEI_ElementList(ielem)
            endif
         !endif
      enddo
      
      !Communicate how many elements I will send 
      allocate(NelemINeed(a%MPIsize))
      call MPI_AllToAll(NelemOthersNeed, 1, MPI_INTEGER4, NelemINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      !Send and receive SizesElementsOthersNeed, but only the required ones
      allocate(SizesElementsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(SizesElementsOthersNeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (NelemINeed(irank+1) /= 0) then
            call MPI_IRECV(SizesElementsINeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Allocate Communication Structures
      allocate(ElementsINeed(a%MPIsize))
      allocate(GEI_ElementsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) then
            allocate(GEI_ElementsINeed(irank+1)%a(NelemINeed(irank+1)))
            allocate(ElementsINeed(irank+1)%p(NelemINeed(irank+1)+1))
            ElementsINeed(irank+1)%p(1) = 1
            allocate(ElementsINeed(irank+1)%l(SizesElementsINeed(irank+1)))
         endif
      enddo
      
      !For knowing the sizes to send
      GEI_Size = sizeof(GEI_Sizer)
      
      !Send all the remaining info
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(ElementsOthersNeed(irank+1)%p,NelemOthersNeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(ElementsOthersNeed(irank+1)%l,SizesElementsOthersNeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
            call MPI_ISEND(GEI_ElementsOthersNeed(irank+1)%a,NelemOthersNeed(irank+1)*GEI_Size,MPI_byte,irank,mtag3,a%MPIcomm,irequest03(irank+1),ierr) 
         endif
      enddo
      !Receive all the remaining info
      do irank = 0,a%MPIsize-1
         if (NelemINeed(irank+1) /= 0) then
            call MPI_IRECV(ElementsINeed(irank+1)%p,NelemINeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest04(irank+1), ierr)
            call MPI_IRECV(ElementsINeed(irank+1)%l,SizesElementsINeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest05(irank+1), ierr)
            call MPI_IRECV(GEI_ElementsINeed(irank+1)%a,NelemINeed(irank+1)*GEI_Size,MPI_byte,irank,mtag3,a%MPIcomm,irequest06(irank+1),ierr) 
         endif
      enddo
      !Wait
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest05, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest06, MPI_STATUSES_IGNORE, ierr)
      
      !Next I need to merge received info for each element with already available info
      !decompress, seen and iwa?
      !or fill a repeated linked list, then decompress, seen, iwa and compress without repeating
      allocate(HeaderElementRanks(a%nelem))
      call StorageElementRanks%Initialize(2*a%nelem,'NonRemovable','Repeatable',a%Memor)
      
      !First add what I already had locally
      do ielem = 1,a%nelem
         nz = size(a%RebalanceElementProcessors(ielem)%l)
         do iz = 1,nz
            irank = a%RebalanceElementProcessors(ielem)%l(iz)
            call StorageElementRanks%Add(HeaderElementRanks(ielem),irank)
         enddo
      enddo
      
      !Now add what I have received
      do irank = 0,a%MPIsize-1
         do elemi = 1, NelemINeed(irank+1)
            !From GEI to local
            call a%LocalElementFromGEI(GEI_ElementsINeed(irank+1)%a(elemi),ielem)
            
            prank = ElementsINeed(irank+1)%p(elemi+1) - ElementsINeed(irank+1)%p(elemi)
            irpos = ElementsINeed(irank+1)%p(elemi)-1
            do iz = 1,prank
               jrank = ElementsINeed(irank+1)%l(irpos+iz)
               call StorageElementRanks%Add(HeaderElementRanks(ielem),jrank)
            enddo
         enddo
      enddo
      
      !now loop elements and merge lists without repeating
      allocate(seen(a%MPIsize))
      allocate(iwa(a%MPIsize))
      seen = 0
      
      do ielem = 1,a%nelem
         nz = 0
         call StorageElementRanks%GoToFirstOfList(HeaderElementRanks(ielem))
         call StorageElementRanks%GetNext(irank)
         
         do while (irank /= -1)
            if (seen(irank+1) /= ielem) then
               seen(irank+1) = ielem
               nz = nz+1
               iwa(nz) = irank
            endif
            call StorageElementRanks%GetNext(irank)
         enddo
      
         deallocate(a%RebalanceElementProcessors(ielem)%l)
         allocate(a%RebalanceElementProcessors(ielem)%l(nz))
         
         a%RebalanceElementProcessors(ielem)%l(1:nz) = iwa(1:nz)
      enddo
      call StorageElementRanks%Dealloc
      
      
      !We are now ready to send nodes and elements
      
   end subroutine
   
   subroutine RebalanceSendNodesAndElements(a)
      use typre
      use Mod_LinkedList
      use Mod_PointToProc
      use Mod_HeapSortExternal
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      type(LinkedListHeader), allocatable :: HeaderNodesOthersNeed(:)
      type(LinkedListStorage)             :: StorageNodesOthersNeed
      
      integer(ip), allocatable :: NpoinOthersNeed(:),NpoinINeed(:),NpoinLocals(:)
      type(CSRList), allocatable, target :: ElementsOthersNeed(:), ElementsINeed(:)
      type :: ListGEI
         type(GlobalElementIdentifier), allocatable :: a(:)
      end type
      
      
      type :: ListETL
         type(ElemTypeList), allocatable :: a(:)
      end type
      type(ListETL), allocatable :: ETL_ElementsINeed(:), ETL_ElementsOthersNeed(:)
      integer(ip) :: ETL_Size
      type(ElemTypeList) :: ETL_Sizer
      
      integer(ip), allocatable :: SizesElementsINeed(:), SizesElementsOthersNeed(:)
      integer(ip), allocatable :: auxRebalanceListINeed(:,:)
      
      type(PointToProc) :: PointToProcessor
      integer(ip) :: ipoin,irank,ielem,pnode,ranki,ispos,ispos2,elemi,inode
      integer(ip), pointer :: lnods(:) => NULL()
      logical:: IsResponsible
      integer(ip) :: GEI_Size
      
      integer(ip), allocatable :: seen(:),iwaGhosts(:),auxiwaGhosts(:)
      integer(ip), allocatable :: auxpnods(:),auxlnods(:)
      type(GlobalElementIdentifier), allocatable :: auxGEI_ElementList(:)
      
      type(ListGEI), allocatable :: GEI_ElementsINeed(:),GEI_ElementsOthersNeed(:)
      
      type(LinkedListHeader), allocatable :: HeaderElementLevels(:)
      type(LinkedListStorage) :: StorageElementLevels
      integer(ip) :: ilevel,elemcount
      
      integer(ip), allocatable :: NewLocalToGlobal(:),NewProcessorList(:)
      integer(ip), allocatable :: ia(:),ja(:)
      integer(ip) :: isize1
      
      integer(ip) :: auxminghost, auxmaxghost,gipoin,ighost,nauxGhostNodes,nauxGhostNodes2
      integer(ip), allocatable :: auxGhostNodes(:)
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3,mtag4 = 4
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:),irequest05(:),irequest06(:),irequest07(:),irequest08(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      type(ElemTypeList), allocatable :: auxETL_ElementList(:)
      type(GlobalElementIdentifier) :: GEI_Sizer
      integer(ip) :: i
      integer(ip), allocatable :: auxReorDeringElements(:)
      
      call a%Timers%RebalanceSendNodesAndElements%Tic
      
      !allocate(HeaderNodesOthersNeed(a%MPIsize))
      !call StorageNodesOthersNeed%Initialize(a%npoinLocal,'NonRepeatable','NonRemovable',a%Memor)
      
      allocate(NpoinOthersNeed(a%MPIsize))
      allocate(NpoinINeed(a%MPIsize))
      allocate(NpoinLocals(a%MPIsize))
      NpoinOthersNeed = 0
      
      do ipoin = 1,a%npoinLocal
         irank = a%RebalanceProcessorList(ipoin)
         NpoinOthersNeed(irank+1) = NpoinOthersNeed(irank+1) + 1
      enddo   
      
      !I Need to: 
      !prepare a local numbering
      !Compute the number of local points
      call MPI_AllToAll(NpoinOthersNeed, 1, MPI_INTEGER4, NpoinINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
!--------------------------------------------------------------
      !npoinLocal to Rebalanced from now on
      
      !Number of local points for my processor
      a%npoinLocal = sum(NpoinINeed)
      
      NpoinOthersNeed = a%npoinLocal
      !NpoinINeed now stores the number of local points of each processor
      call MPI_AllToAll(NpoinOthersNeed, 1, MPI_INTEGER4, NpoinLocals, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
!-------------------------------------------------------------
      !poin0 to Rebalanced from now on
      !Compute my poin0
      a%poin0 = 0
      do irank = 0,a%MPIrank-1
         a%poin0 = a%poin0 + NpoinLocals(irank+1)
      enddo
      
      !Let us take the oportunity to start the PointToProc object, avoiding one AllToAll
      call PointToProcessor%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call PointToProcessor%SetMemor(a%Memor)
      call PointToProcessor%InitializeFromList(NpoinLocals)
      
      !deallocate if necessary
      if (allocated(a%RebalanceNelemINeed)) then
         deallocate(a%RebalanceNelemINeed)
         deallocate(a%RebalanceNelemOthersNeed)
      endif
      
      !Send the elements in global numbering, send also the GEI
      allocate(a%RebalanceNelemOthersNeed(a%MPIsize))
      allocate(SizesElementsOthersNeed(a%MPIsize))
      a%RebalanceNelemOthersNeed = 0
      SizesElementsOthersNeed = 0
      
      !First to count
      do ielem = 1,a%nelem
         !Check if I am responsible for the element
         call a%GetLnodsFromElement(ielem,pnode,lnods)
         call a%CheckElementResponsibility(pnode,lnods,IsResponsible)
         
         if (IsResponsible) then
            !Add the element to the list of elements to send for a processor
            do ranki = 1,size(a%RebalanceElementProcessors(ielem)%l)
               irank = a%RebalanceElementProcessors(ielem)%l(ranki)
               a%RebalanceNelemOthersNeed(irank+1) = a%RebalanceNelemOthersNeed(irank+1) + 1
               SizesElementsOthersNeed(irank+1) = SizesElementsOthersNeed(irank+1) + pnode
            enddo
         endif
      enddo
      
      !Allocate
      !List to be used later
      if (allocated(a%RebalanceListElementsOthersNeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%RebalanceListElementsOthersNeed(irank+1)%l)
         enddo
         deallocate(a%RebalanceListElementsOthersNeed)
      endif
      
      
      allocate(a%RebalanceListElementsOthersNeed(a%MPIsize))
      
      allocate(ElementsOthersNeed(a%MPIsize))
      allocate(GEI_ElementsOthersNeed(a%MPIsize))
      allocate(ETL_ElementsOthersNeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(ElementsOthersNeed(irank+1)%p(a%RebalanceNelemOthersNeed(irank+1)+1))
         ElementsOthersNeed(irank+1)%p(1) = 1
         allocate(GEI_ElementsOthersNeed(irank+1)%a(a%RebalanceNelemOthersNeed(irank+1)))
         allocate(ETL_ElementsOthersNeed(irank+1)%a(a%RebalanceNelemOthersNeed(irank+1)))
         allocate(ElementsOthersNeed(irank+1)%l(SizesElementsOthersNeed(irank+1)))
         
         !to be used later
         allocate(a%RebalanceListElementsOthersNeed(irank+1)%l(a%RebalanceNelemOthersNeed(irank+1)))
      enddo
      
      !Second to fill
      a%RebalanceNelemOthersNeed = 0
      do ielem = 1,a%nelem
         !Check if I am responsible for the element
         call a%GetLnodsFromElement(ielem,pnode,lnods)
         call a%CheckElementResponsibility(pnode,lnods,IsResponsible)
         
         if (IsResponsible) then
            !Add the element to the list of elements to send for a processor
            do ranki = 1,size(a%RebalanceElementProcessors(ielem)%l)
               irank = a%RebalanceElementProcessors(ielem)%l(ranki)
               a%RebalanceNelemOthersNeed(irank+1) = a%RebalanceNelemOthersNeed(irank+1) + 1
               
               ElementsOthersNeed(irank+1)%p(a%RebalanceNelemOthersNeed(irank+1)+1) = ElementsOthersNeed(irank+1)%p(a%RebalanceNelemOthersNeed(irank+1)) + pnode
               ispos = ElementsOthersNeed(irank+1)%p(a%RebalanceNelemOthersNeed(irank+1))-1
               ElementsOthersNeed(irank+1)%l(ispos+1:ispos+pnode) = a%RebalancePointNumbering(lnods)
               
               GEI_ElementsOthersNeed(irank+1)%a(a%RebalanceNelemOthersNeed(irank+1)) = a%GEI_ElementList(ielem)
               ETL_ElementsOthersNeed(irank+1)%a(a%RebalanceNelemOthersNeed(irank+1)) = a%ElementTypeList(ielem)
               
               !to be used later
               a%RebalanceListElementsOthersNeed(irank+1)%l(a%RebalanceNelemOthersNeed(irank+1)) = ielem
            enddo
         endif
      enddo
      
      !Communicate how many elements I will send 
      allocate(a%RebalanceNelemINeed(a%MPIsize))
      call MPI_AllToAll(a%RebalanceNelemOthersNeed, 1, MPI_INTEGER4, a%RebalanceNelemINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      !Send and receive SizesElementsOthersNeed, but only the required ones
      !Allocate and Initialize some MPI
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      allocate(SizesElementsINeed(a%MPIsize))
      SizesElementsINeed = 0
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(SizesElementsOthersNeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%RebalanceNelemINeed(irank+1) /= 0) then
            call MPI_IRECV(SizesElementsINeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Allocate Communication Structures
      allocate(ElementsINeed(a%MPIsize))
      allocate(GEI_ElementsINeed(a%MPIsize))
      allocate(ETL_ElementsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNelemINeed(irank+1) /= 0) then
            allocate(GEI_ElementsINeed(irank+1)%a(a%RebalanceNelemINeed(irank+1)))
            allocate(ETL_ElementsINeed(irank+1)%a(a%RebalanceNelemINeed(irank+1)))
            allocate(ElementsINeed(irank+1)%p(a%RebalanceNelemINeed(irank+1)+1))
            ElementsINeed(irank+1)%p(1) = 1
            allocate(ElementsINeed(irank+1)%l(SizesElementsINeed(irank+1)))
         endif
      enddo
      
      !For knowing the sizes to send
      GEI_Size = sizeof(GEI_Sizer)
      ETL_Size = sizeof(ETL_Sizer)
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      allocate(irequest05(a%MPIsize))
      allocate(irequest06(a%MPIsize))
      allocate(irequest07(a%MPIsize))
      allocate(irequest08(a%MPIsize))
      
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL      
      irequest05 = MPI_REQUEST_NULL
      irequest06 = MPI_REQUEST_NULL      
      irequest07 = MPI_REQUEST_NULL
      irequest08 = MPI_REQUEST_NULL
      
      !Send all the remaining info
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNelemOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(ElementsOthersNeed(irank+1)%p,a%RebalanceNelemOthersNeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(ElementsOthersNeed(irank+1)%l,SizesElementsOthersNeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
            call MPI_ISEND(GEI_ElementsOthersNeed(irank+1)%a,a%RebalanceNelemOthersNeed(irank+1)*GEI_Size,MPI_byte,irank,mtag3,a%MPIcomm,irequest03(irank+1),ierr)
            call MPI_ISEND(ETL_ElementsOthersNeed(irank+1)%a,a%RebalanceNelemOthersNeed(irank+1)*ETL_Size,MPI_byte,irank,mtag4,a%MPIcomm,irequest07(irank+1),ierr) 
         endif
      enddo
      !Receive all the remaining info
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNelemINeed(irank+1) /= 0) then
            call MPI_IRECV(ElementsINeed(irank+1)%p,a%RebalanceNelemINeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest04(irank+1), ierr)
            call MPI_IRECV(ElementsINeed(irank+1)%l,SizesElementsINeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest05(irank+1), ierr)
            call MPI_IRECV(GEI_ElementsINeed(irank+1)%a,a%RebalanceNelemINeed(irank+1)*GEI_Size,MPI_byte,irank,mtag3,a%MPIcomm,irequest06(irank+1),ierr) 
            call MPI_IRECV(ETL_ElementsINeed(irank+1)%a,a%RebalanceNelemINeed(irank+1)*ETL_Size,MPI_byte,irank,mtag4,a%MPIcomm,irequest08(irank+1),ierr) 
         endif
      enddo
      !Wait
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest05, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest06, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest07, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest08, MPI_STATUSES_IGNORE, ierr)
      
!----------------------------------------------------------------------------------------------------------------
      !nelem to Rebalanced from now on
      
      !Build the new pnods,lnods from the received ones
      a%nelem = sum(a%RebalanceNelemINeed)
      
!-----------------------------------------------------------------------------------------------------------------
      !pnods, lnods, GEI_ElementList, npoinGhost to Rebalanced from now on
      deallocate(a%pnods,a%lnods,a%GEI_ElementList,a%ElementTypeList)
      allocate(a%pnods(a%nelem+1))
      allocate(a%lnods(sum(SizesElementsINeed)))
      allocate(a%GEI_ElementList(a%nelem))
      allocate(a%ElementTypeList(a%nelem))
      
      
      
      
      !Sorted list of ghost points 
      allocate(auxGhostNodes(size(a%lnods)))
      nauxGhostNodes = 0
      
      auxminghost = a%poin0
      auxmaxghost = a%poin0+a%npoinLocal+1
      
      do irank = 0, a%MPIsize-1
         do elemi = 1,a%RebalanceNelemINeed(irank+1)
            
            ispos = ElementsINeed(irank+1)%p(elemi)-1
            pnode = ElementsINeed(irank+1)%p(elemi+1)-ElementsINeed(irank+1)%p(elemi)
            lnods => ElementsINeed(irank+1)%l(ispos+1:ispos+pnode)
            
            do inode = 1,pnode
               ipoin = lnods(inode)
               if (ipoin <= auxminghost .or. ipoin >= auxmaxghost) then
                  nauxGhostNodes = nauxGhostNodes+1
                  auxGhostNodes(nauxGhostNodes) = ipoin
               endif
            enddo
            
           
            
         enddo
      enddo
      call HeapSort(nauxGhostNodes,auxGhostNodes)
      call EliminateRepeatedEntries(nauxGhostNodes,auxGhostNodes,nauxGhostNodes2)
      a%npoinGhost = nauxGhostNodes2
      
      !Copy elements, find ghost nodes
      a%pnods(1) = 1
      ielem = 0
      
      auxminghost = a%poin0
      auxmaxghost = a%poin0+a%npoinLocal+1
      
      if (allocated(a%RebalanceListElementsINeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%RebalanceListElementsINeed(irank+1)%l)
         enddo
         deallocate(a%RebalanceListElementsINeed)
      endif 
      
      allocate(a%RebalanceListElementsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(a%RebalanceListElementsINeed(irank+1)%l(a%RebalanceNelemINeed(irank+1)))
      enddo
      
      do irank = 0, a%MPIsize-1
         do elemi = 1,a%RebalanceNelemINeed(irank+1)
            ispos = ElementsINeed(irank+1)%p(elemi)-1
            pnode = ElementsINeed(irank+1)%p(elemi+1)-ElementsINeed(irank+1)%p(elemi)
            lnods => ElementsINeed(irank+1)%l(ispos+1:ispos+pnode)
            
            !Fill the list of elements in local numbering
            ielem = ielem+1
            ispos2 = a%pnods(ielem)-1
            a%pnods(ielem+1) = a%pnods(ielem) + pnode
            do inode = 1,pnode
               gipoin = lnods(inode)
               if (gipoin > auxminghost .and. gipoin < auxmaxghost) then
                  ipoin = gipoin - a%poin0
               else
                  ighost = binarySearch_I(auxGhostNodes(1:a%npoinGhost),gipoin)
                  if (ighost == 0) then
                     call runend('ghost not found')
                  endif
                  ipoin = a%npoinLocal+ighost
               endif
               a%lnods(ispos2+inode) = ipoin
            enddo
            !a%lnods(ispos2+1:ispos2+pnode) = seen(lnods)
            a%GEI_ElementList(ielem) = GEI_ElementsINeed(irank+1)%a(elemi)
            a%ElementTypeList(ielem) = ETL_ElementsINeed(irank+1)%a(elemi)
            
            a%RebalanceListElementsINeed(irank+1)%l(elemi) = ielem
         enddo
      enddo
      deallocate(ETL_ElementsINeed)
      
      do i = 1,size(GEI_ElementsINeed)
         if (allocated(GEI_ElementsINeed(i)%a) )deallocate(GEI_ElementsINeed(i)%a)
      enddo
      deallocate(GEI_ElementsINeed)
      
      do i = 1,size(GEI_ElementsOthersNeed)
         if (allocated(GEI_ElementsOthersNeed(i)%a)) deallocate(GEI_ElementsOthersNeed(i)%a)
      enddo
      deallocate(GEI_ElementsOthersNeed)
      
      
      !Reorder elements so that high level elements appear first
      allocate(auxReorderingElements(a%nelem))
      
      allocate(auxpnods(a%nelem+1))
      allocate(auxlnods(size(a%lnods)))
      allocate(auxGEI_ElementList(a%nelem))
      allocate(auxETL_ElementList(a%nelem))
      !Linked list of elements for each level
      allocate(HeaderElementLevels(GEI_ElementMaxLevels+1))
      call StorageElementLevels%Initialize(a%nelem,'NonRemov','NonRepeatable',a%Memor)
      do ielem = 1,a%nelem
         ilevel = a%GEI_ElementList(ielem)%Level
         call StorageElementLevels%Add(HeaderElementLevels(ilevel+1),ielem)
      enddo
      auxpnods(1) = 1
      elemcount = 0
      do ilevel = 0,GEI_ElementMaxLevels
         !mark the level0 nodes
         call StorageElementLevels%GoToFirstOfList(HeaderElementLevels(ilevel+1))
         call StorageElementLevels%GetNext(ielem)
         do while (ielem /= -1) 
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            
            elemcount = elemcount+1
            ispos = auxpnods(elemcount)-1
            auxpnods(elemcount+1) = auxpnods(elemcount)+pnode
            auxlnods(ispos+1:ispos+pnode) = lnods
            auxGEI_ElementList(elemcount) = a%GEI_ElementList(ielem)
            auxETL_ElementList(elemcount) = a%ElementTypeList(ielem)

            auxReorderingElements(ielem) = elemcount

            call StorageElementLevels%GetNext(ielem)
            
            
         enddo
      enddo
      call StorageElementLevels%Dealloc
      call move_alloc(auxpnods,a%pnods)
      call move_alloc(auxlnods,a%lnods)
      call move_alloc(auxGEI_ElementList,a%GEI_ElementList)
      call move_alloc(auxETL_ElementList,a%ElementTypeList)
      
      
      do irank = 0, a%MPIsize-1
         do elemi = 1,a%RebalanceNelemINeed(irank+1)
            a%RebalanceListElementsINeed(irank+1)%l(elemi) = auxReorderingElements(a%RebalanceListElementsINeed(irank+1)%l(elemi))
         enddo
      enddo
      
      deallocate(auxReorderingElements)
      
      
      !Build the new local ordering
      a%npoin = a%npoinLocal + a%npoinGhost
      
      allocate(NewLocalToGlobal(a%npoin))
      allocate(NewProcessorList(a%npoin))
      
      do ipoin = 1,a%npoinLocal
         NewLocalToGlobal(ipoin) = a%poin0 + ipoin
         NewProcessorList(ipoin) = a%MPIrank
      enddo
      
      do ipoin = 1,a%npoinGhost
         NewLocalToGlobal(a%npoinLocal+ipoin) = auxGhostNodes(ipoin)
         call PointToProcessor%GetProc(auxGhostNodes(ipoin),NewProcessorList(a%npoinLocal+ipoin))
      enddo
      
      call a%LocalOrdering%Deallocate(a%Memor)
      call a%LocalOrdering%Init(a%MPIcomm,a%npoin,NewLocalToGlobal,NewProcessorList,a%Memor)

      deallocate(NewLocalToGlobal, NewProcessorList)
      
      !We reinitialize the communicator
      call a%Communicator%Deallocate
      allocate(ia(1), ja(1))
      call a%Communicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,a%npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)
      deallocate(ia,ja)
      !I "only" need to rebuild the refinement structure
      
      
      !We also need a local to global ordering for top level elements
      call a%ElementLocalOrdering%Deallocate(a%Memor)
      call BuildElementLocalOrdering(a%MPIcomm,a%nelem,a%GEI_ElementList,a%ElementLocalOrdering,a%Memor)
      
      call PointToProcessor%Finalize
      
      call a%Timers%RebalanceSendNodesAndElements%Toc
      
   end subroutine
      
   subroutine RebuildRefinementStructure(a)
      use typre
      use Mod_LinkedList
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      type(LinkedListHeader), allocatable :: HeaderChildren(:)
      type(LinkedListStorage) :: StorageChildren
      
      integer(ip) :: npoinLevel0,ielem,ilevel,pnode
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: inode,ipoin,iparent
      type(GlobalElementIdentifier) :: GEI_Parent
      
      integer(ip) :: nelemLevel0,nchild,ispos,nedge,ledger(3,MaxNodesPerElement),nz
      integer(ip) :: iedge,jpoin,json,iz,icpos
      integer(ip), allocatable :: iwa(:,:)
      
      type(LinkedListHeader), allocatable :: HeaderEdger(:)
      type(LinkedListStorage) :: StorageEdger
      
      type(EdgerCriteriaTester) :: ECT
      
      integer(ip) :: ilevel0,ielem0,ielem00,jelem,sizeGEIList
      logical :: continuetest
      integer(ip), allocatable :: touched(:)
      
      integer(ip) :: auxChildrenList(50),kelem,kchild,childk
      integer(ip) :: ielty
      
      call a%Timers%RebuildRefinementStructure%Tic
      
      !We might start by rebuilding PointerToParent and PointerToChildren
      !linked list of elements for each level
      deallocate(a%PointerToParent)
      deallocate(a%OldPointerToParent)
      allocate(a%PointerToParent(a%nelem),a%OldPointerToParent(0))
      a%PointerToParent = 0_ip
      
      !markel
      deallocate(a%markel)
      allocate(a%markel(a%nelem))
      a%markel = 0
      
      !level
      deallocate(a%level)
      allocate(a%level(a%npoin))
      a%level = -1
      
      !Linked list of children for each element
      allocate(HeaderChildren(a%nelem))
      call StorageChildren%Initialize(a%nelem,'NonRemov','NonRepeatable',a%Memor)
      
      !Count number of top level elements
      !Make sure there are no problems if all elements are top or if nelem == 0
      nelemLevel0 = 0
      continuetest = .true.
      sizeGEIList = size(a%GEI_ElementList)
      if (nelemLevel0 + 1 > sizeGEIList) continuetest = .false.
      do while (continuetest .eqv. .true.) 
         if (a%GEI_ElementList(nelemLevel0+1)%level == 0) then
            nelemLevel0 = nelemLevel0+1
         else
            continuetest = .false.
         endif
         if (nelemLevel0 + 1 > sizeGEIList) continuetest = .false.
      enddo

      !Compress PointerToChildren
      deallocate(a%lPointerToChildren,a%pPointerToChildren)
      allocate(a%lPointerToChildren(a%nelem-nelemLevel0))
      allocate(a%pPointerToChildren(a%nelem+1))
      a%pPointerToChildren(1) = 1
      
      npoinLevel0 = 0
      nelemLevel0 = 0
      ilevel0 = 0
      ielem0 = 0
      ielem00 = 0
      !Level, touched, PointerToParent, PointerToChildren, markel
      !Elements are already ordered by levels (level 0 first)
      do ielem = 1,a%nelem
         !level of points, touched
         call a%GetLnodsFromElement(ielem,pnode,lnods)
         ilevel = a%GEI_ElementList(ielem)%Level
         
         !On change of level, prepare pointer to children of previous level
         if (ilevel /= ilevel0 ) then
            if (ilevel >= 2) then
               do jelem = ielem00+1,ielem0
                  call HeaderChildren(jelem)%GetNelem(nchild)
                  a%pPointerToChildren(jelem+1)  = a%pPointerToChildren(jelem)+nchild
                  
                  if (nchild /= 0) then
                     ispos = a%pPointerToChildren(jelem)-1
                     
                     !Children need to be added in the proper order
                     !we first load them, then see which son they are, then add them
                     call StorageChildren%ListToArray(HeaderChildren(jelem),auxChildrenList)
                     do childk = 1,nchild
                        kelem = auxChildrenList(childk)
                        call GEI_GetPositionInElementForLevel(a%GEI_ElementList(kelem),ilevel0,kchild)
                        a%lPointerToChildren(ispos+kchild) = kelem
                     enddo
                  endif
               enddo
            endif
            
            ilevel0 = ilevel
            ielem00 = ielem0
            ielem0 = ielem-1
         endif
         
         if (ilevel == 0) then
            nelemLevel0 = nelemLevel0 + 1
            do inode = 1,pnode
               ipoin = lnods(inode)
               if (a%level(ipoin) == -1) then
                  a%level(ipoin) = ilevel
                  npoinLevel0 = npoinLevel0+1
               endif
               
            enddo
         else
            do inode = 1,pnode
               ipoin = lnods(inode)
               if (a%level(ipoin) == -1) then
                  a%level(ipoin) = ilevel
               endif
            enddo
            
            
            GEI_Parent = a%GEI_ElementList(ielem)
            GEI_Parent%level = GEI_Parent%level-1
            call a%LocalElementFromGEI(GEI_Parent,iparent)
           
            a%PointerToParent(ielem) = iparent
            call StorageChildren%Add(HeaderChildren(iparent),ielem)
            
            !markel: parent is refined
            a%markel(iparent) = 1_ip
         endif
      enddo
      
      !Complete pointer to children for last level elements
      do jelem = ielem00+1,a%nelem
         call HeaderChildren(jelem)%GetNelem(nchild)
         a%pPointerToChildren(jelem+1)  = a%pPointerToChildren(jelem)+nchild
         
         if (nchild /= 0) then
            ispos = a%pPointerToChildren(jelem)-1
            
            !Children need to be added in the proper order
            !we first load them, then see which son they are, then add them
            call StorageChildren%ListToArray(HeaderChildren(jelem),auxChildrenList)
            do childk = 1,nchild
               kelem = auxChildrenList(childk)
               call GEI_GetPositionInElementForLevel(a%GEI_ElementList(kelem),ilevel0,kchild)
               a%lPointerToChildren(ispos+kchild) = kelem
            enddo
            
         endif
      enddo
      call StorageChildren%Dealloc
      
      
      !Pedger, ledger, touched
      deallocate(a%pedger,a%ledger)
      allocate(a%pedger(a%npoin+1))
      allocate(a%ledger(2,a%npoin-npoinLevel0))
      
      allocate(iwa(2,a%npoin-npoinLevel0))
      nz = 0
      
      allocate(HeaderEdger(a%npoin))
      call StorageEdger%Initialize(a%npoin,'NonRemov','NonRepeat',a%Memor)

      !touched
      allocate(touched(a%npoin))
      touched = 0
      
      do ielem = 1,a%nelem
         nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
         
         if (nchild /= 0) then
            icpos = a%pPointerToChildren(ielem)-1
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            ielty = a%ElementTypeList(ielem)%ielty
            call get_edges_from_child(ielty,a%pnods,a%lnods,nchild,a%lPointerToChildren(icpos+1:icpos+nchild),nedge,ledger)
            
            do iedge = 1,nedge
               !need to add to edger and pedger, but just once per node
               if (touched(ledger(3,iedge)) == 0) then
                  ipoin = ledger(1,iedge)
                  jpoin = ledger(2,iedge)
                  json = ledger(3,iedge)
                  
                  !Proper ordering of edge parents
                  ECT%parents(1) = ipoin
                  ECT%parents(2) = jpoin
                  ECT%levels(1) = a%level(ipoin)
                  ECT%levels(2) = a%level(jpoin)
                  call EdgerCriteria(ECT)   
                  ipoin = ECT%orderedParents(1)
                  jpoin = ECT%orderedParents(2)
                  
                  nz = nz+1
                  iwa(1,nz) = jpoin
                  iwa(2,nz) = json
                  
                  call StorageEdger%Add(HeaderEdger(ipoin),nz)
               endif
               !touched
               touched(ledger(3,iedge)) = 1
            enddo
         endif
      enddo
      deallocate(touched)
      
      !Compress pedger and ledger
      a%pedger(1) = 1
      nz = 0
      do ipoin = 1,a%npoin
         call StorageEdger%GoToFirstOfList(HeaderEdger(ipoin))
         call StorageEdger%GetNext(iz)
         do while (iz /= -1)
            nz = nz+1
            
            a%ledger(:,nz) = iwa(:,iz)
            
            call StorageEdger%GetNext(iz)
         enddo
         a%pedger(ipoin+1) = nz+1
      enddo
      call StorageEdger%Dealloc
            
      !Pfacer, lfacer
      deallocate(a%pfacs,a%lfacs)
      !This allocates, sets pfacer, and does everything except for hanging faces
      call BuildLfacs(a)
      !this fills lfacs for hanging faces
      do ielem = 1,a%nelem
         nchild = a%pPointerToChildren(ielem+1)-a%pPointerToChildren(ielem)
         if (nchild /= 0) then
            call upd_facs_newel(a,ielem)
         endif
      enddo  
      
!       !npoinTop, EdgerEntryPoints
!       deallocate(a%EdgerEntryPoints)
!       allocate(a%EdgerEntryPoints(npoinLevel0))
!       a%npoinTop = 0
!       do ipoin = 1,a%npoin
!          if (a%level(ipoin) == 0) then
!             a%npoinTop = a%npoinTop+1
!             a%EdgerEntryPoints(a%npoinTop) = ipoin
!          endif
!       enddo
      
      !HangingList, Export, Extra
      call a%HangingNodes  
      call a%SelectComponentsToExport
      
      call a%Timers%RebuildRefinementStructure%Toc
      
   end subroutine
      
   subroutine PrepareRebalanceCommunicator(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      logical :: kfl_isallocated
      
      call a%Timers%PrepareRebalanceCommunicator%Tic
      
      call a%RebalanceComm%IsAllocated(kfl_isallocated)
      if (kfl_isallocated .eqv. .true.) then
         call a%RebalanceComm%Deallocate
      endif
      call a%RebalanceComm%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%RebalanceComm%SetOutputFiles(a%lun_memo,a%lun_outpu)
      call a%RebalanceComm%Initialize(a%RebalanceOldNpoinLocal,a%npoinLocal,a%RebalancePointNumbering,a%RebalanceProcessorList)
      
      call a%Timers%PrepareRebalanceCommunicator%Toc
   end subroutine
      
   subroutine GetLnodsFromElement(a,ielem,pnode,lnods)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      integer(ip) :: ielem,pnode
      integer(ip), pointer :: lnods(:)
      
      integer(ip) :: ispos
      
      ispos = a%pnods(ielem)-1
      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      
      lnods => a%lnods(ispos+1:ispos+pnode)
   end subroutine
   
   subroutine CheckElementResponsibility(a,pnode,lnods,IsResponsible)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: pnode,lnods(pnode)
      logical :: IsResponsible
      
      integer(ip) :: lranks(pnode)
      
      IsResponsible = .false.
      call a%LocalOrdering%GetProc(pnode,lnods,lranks)
      
      if (minval(lranks) == a%MPIrank) IsResponsible = .true.
   end subroutine
      
   !ToMesh
   subroutine RebalanceVariableReal(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      real(rp) :: coord(ndime,*)
      real(rp) :: newcoord(ndime,*)
      
      call a%RebalanceComm%Communicate(ndime,coord(:,1:1),newcoord(:,1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:1))
   end subroutine
   
   subroutine RebalanceVariableReal1(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      real(rp), target :: coord(*)
      real(rp), target :: newcoord(*)
      
      
      call a%RebalanceComm%Communicate(ndime,coord(1:1),newcoord(1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(1:1))
   end subroutine
   
   subroutine RebalanceVariableInteger(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      integer(ip) :: coord(ndime,*)
      integer(ip) :: newcoord(ndime,*)

      
      call a%RebalanceComm%Communicate(ndime,coord(:,1:1),newcoord(:,1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:1))
   end subroutine
   
   subroutine RebalanceVariableInteger1(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      integer(ip), target :: coord(*)
      integer(ip), target :: newcoord(*)
      
      call a%RebalanceComm%Communicate(ndime,coord(1:1),newcoord(1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(1:1))
   end subroutine
   
   subroutine RebalanceVariableLogical(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      logical :: coord(ndime,*)
      logical :: newcoord(ndime,*)

      
      call a%RebalanceComm%Communicate(ndime,coord(:,1:1),newcoord(:,1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(:,1:1))
   end subroutine
   
   subroutine RebalanceVariableLogical1(a,ndime,coord,newcoord)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndime !array dimension
      logical, target :: coord(*)
      logical, target :: newcoord(*)
      
      call a%RebalanceComm%Communicate(ndime,coord(1:1),newcoord(1:1))
      call a%ExportCommunicator%GhostCommunicate(ndime,newcoord(1:1))
   end subroutine
   
   subroutine MarkDeadElements(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip), allocatable :: ichild(:),ielem(:)
      integer(ip), allocatable :: nchild(:),iepos(:)
      
      logical, allocatable :: isAlive(:), LevelIsAlive(:)
      
      integer(ip) :: elemi,ilevel,ispos,pnode,inode,ipoin,irank,ispos2,pnode2,jelem
      integer(ip) :: siblingi, isibling,iz
      logical :: kfl_continue, kfl_lastsibling
      
      call a%Timers%MarkDeadElements%Tic
      
      !I need to look for elements which are dead
      !This only happens when there has been load rebalancing
      !They will be erased in the next round
      
      !Dead Elements are identified because:
         ! They do not have a local node
         ! Their siblings do not have a local node
         ! Their children do not have a local node
         ! They are not high level element of a hanging node
         ! Their siblings are not the high level element of a hanging node
         ! Tehir children are not the high level element of a hanging node
      
      !Move recursively through the element structure
      allocate(ichild(GEI_ElementMaxLevels))
      allocate(nchild(GEI_ElementMaxLevels))
      ichild = 0
      allocate(ielem(GEI_ElementMaxLevels))
      allocate(iepos(GEI_ElementMaxLevels))
      
      allocate(IsAlive(a%nelem))
      isAlive = .false.
      allocate(LevelIsAlive(GEI_ElementMaxLevels))
      LevelIsALive = .false.

      do elemi = 1,a%nelem
         ielem(1) = elemi
         !Enter through top level elements, move through all of them recursively
         if (a%GEI_ElementList(ielem(1))%level == 0) then
            ilevel = 0
            nchild(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1)+1)-a%pPointerToChildren(ielem(ilevel+1))
            iepos(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1))-1
            ichild(ilevel+1) = 0
            kfl_continue = .true.
                        
            do while (kfl_continue .eqv. .true.)
               
               if (ichild(ilevel+1) < nchild(ilevel+1)) then
                  !Move down to the next children
                  ichild(ilevel+1) = ichild(ilevel+1)+1
                  !Move down
                  ilevel = ilevel + 1
                  ielem(ilevel+1) = a%lPointerToChildren(iepos(ilevel)+ichild(ilevel))
                  
                  iepos(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1))-1
                  nchild(ilevel+1) = a%pPointerToChildren(ielem(ilevel+1)+1)-a%pPointerToChildren(ielem(ilevel+1))
                  ichild(ilevel+1) = 0
                  
               else
                  !Mark the new ranks of the nodes for the level
                  ispos = a%pnods(ielem(ilevel+1))-1
                  pnode = a%pnods(ielem(ilevel+1)+1)-a%pnods(ielem(ilevel+1))
                  do inode = 1,pnode
                     ipoin = a%lnods(ispos+inode)
                     if (ipoin <= a%npoinLocal) then
                        LevelIsAlive(ilevel+1) = .true.
                     endif
                  enddo
                  
                  !Do the same for the hanging nodes of which I am high level element
                  pnode2 = a%HangingListForElements%p(ielem(ilevel+1)+1)-a%HangingListForElements%p(ielem(ilevel+1))
                  if (pnode2 /= 0) then
                     ispos2 = a%HangingListForElements%p(ielem(ilevel+1))-1
                     do inode = 1,pnode2
                        ipoin = a%HangingListForElements%l(ispos2+inode)
                        if (ipoin <= a%npoinLocal) then
                           LevelIsAlive(ilevel+1) = .true.
                        endif
                     enddo
                  endif
                  
                  !Top level, just copy and exit
                  if (ilevel == 0) then
                     !Copy data
                     isAlive(ielem(ilevel+1)) = LevelIsAlive(ilevel+1)
                     LevelIsALive(ilevel+1) = .false.
                     
                     !We are done
                     ilevel = ilevel - 1
                     kfl_continue = .false.
                  !Lower levels
                  else
                     !If this is the last sibling, then do all the siblings of the level
                     if (ichild(ilevel) == nchild(ilevel)) then
                        
                        !Loop through my parents children and fill them             
                        do siblingi = 1,nchild(ilevel)
                           isibling = a%lPointerToChildren(iepos(ilevel)+siblingi)
                           
                           isAlive(isibling) = LevelIsAlive(ilevel+1)
                        enddo
                     
                        !Transfer info to the upper level
                        if (LevelIsAlive(ilevel+1)) then
                           LevelIsALive(ilevel) = .true.
                        endif
                        LevelIsAlive(ilevel+1) = .false.
                     endif
                     
                     !And move up
                     ilevel = ilevel - 1
                  endif
               endif
            
            enddo
         endif
      enddo
      
      !transfer the info to markel
      do jelem = 1,a%nelem
         if (isAlive(jelem) .eqv. .false.) then
            a%markel(jelem) = -2
         endif
      enddo
      
      call a%Timers%MarkDeadElements%Toc
      
      
   end subroutine
   
   !For boundary info
   subroutine PrepareOldExternalBoundariesInfoToKeep(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip) :: ielem,iboun,pface,ifpos,iface,iparent,ichild,iparentface,nchild,nchildface,nface,pnodb,pnode
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: facnod(a%mnodb),ielty,ivariation
      
      call a%Timers%PrepareOldExternalBoundariesInfoToKeep%Tic
      
      !First count sizes
      iboun = 0
      do ielem = 1,a%nelem
         if (a%ExportElement(ielem) .eqv. .true.) then
            pface = a%pfacs(ielem+1)-a%pfacs(ielem)
            ifpos = a%pfacs(ielem)-1
         
            do iface = 1,pface
                  if (a%lfacs(1,ifpos+iface) == 0 .or. a%lfacs(1,ifpos+iface) > a%nelem) then
                  iboun = iboun+1
               endif
            enddo
         endif
      enddo
      
      !Second fill
      if (allocated(a%OldExternalBoundaries)) then
         deallocate(a%OldExternalBoundaries)
         deallocate(a%HeaderOldElementBoundaries)
         call a%StorageOldElementBoundaries%Dealloc
      endif
      
      
      allocate(a%OldExternalBoundaries(4,iboun))
      a%OldExternalBoundaries = 0
      
      allocate(a%HeaderOldElementBoundaries(a%nelem))
      call a%StorageOldElementBoundaries%Initialize(iboun,'NonRemov','Repeatable',a%Memor)
      
      iboun = 0
      do ielem = 1,a%nelem
         if (a%ExportElement(ielem) .eqv. .true.) then
            pface = a%pfacs(ielem+1)-a%pfacs(ielem)
            ifpos = a%pfacs(ielem)-1
         
            do iface = 1,pface
               if (a%lfacs(1,ifpos+iface) == 0 .or. a%lfacs(1,ifpos+iface) > a%nelem) then
                  iboun = iboun+1
                  
                  a%OldExternalBoundaries(1,iboun) = ielem
                  a%OldExternalBoundaries(2,iboun) = iface
                  
                  call a%StorageOldElementBoundaries%Add(a%HeaderOldElementBoundaries(ielem),iboun)
                  
                  !My external parent
                  if (a%lfacs(1,ifpos+iface) > a%nelem) then
                     iparent = a%PointerToParent(ielem)
                     call GEI_GetPositionInElementForLevel(a%GEI_ElementList(ielem),int(a%GEI_ElementList(ielem)%level,4),ichild)
                     call a%GetLnodsFromElement(ielem,pnode,lnods)
                     ielty = a%ElementTypeList(iparent)%ielty
                     ivariation = a%ElementTypeList(iparent)%ivariation
                     call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
                     iparentface = iface  !entry face for seeking up faces
                     call up_faces(ielty,ivariation,nface,ichild,iparentface)
                     
                     !With load rebalancing, faces with all nodes ghost will not have the boundary info
                     !but the parent of an all-ghost face might be local!
                     !I can only give the info to my parent if at least one of my nodes is local
                     !This ensures that I have the info
                     call face_from_el(ielty,nface,pnodb,lnods,iface,0,facnod)
                     if (minval(facnod(1:pnodb)) <= a%npoinLocal) then
                        a%OldExternalBoundaries(3,iboun) = iparent
                        a%OldExternalBoundaries(4,iboun) = iparentface
                        call a%StorageOldElementBoundaries%Add(a%HeaderOldElementBoundaries(iparent),iboun)
                     endif
                  endif
               endif
            enddo
            
         endif   
      enddo   
      
      call a%Timers%PrepareOldExternalBoundariesInfoToKeep%Toc
      
   end subroutine
   
   !For the old faces of the mesh, say to which external boundary they belong
   subroutine MatchOldLboelToOldExternalBoundaries(a,nboun,pboel,lboel,BoundaryMatch,ExportOldNboun)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: nboun,pboel(nboun+1),lboel(pboel(nboun+1)-1),BoundaryMatch(nboun)
      
      integer(ip) :: iboun,ibpos,pnodb,ielemExport,ielem,pnode,iface,testboun,ExportOldNboun
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: ielty
      
      BoundaryMatch = -1
      ExportOldNboun = size(a%OldExternalBoundaries,2)
      
      do iboun = 1,nboun
         ibpos = pboel(iboun)-1
         pnodb = pboel(iboun+1)-pboel(iboun)-1
         ielemExport = lboel(ibpos+pnodb+1)
         
         if (ielemExport <= a%ExportNelem) then
            ielem = a%iRenelExport(ielemExport)
         else 
            ielem = 0
         endif
         !Not for extra elements
         if (ielem /= 0) then
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            ielty = a%ElementTypeList(ielem)%ielty
            call Nodes_To_Face(ielty,pnodb,lboel(ibpos+1),iface)
            
            call a%StorageOldElementBoundaries%GoToFirstOfList(a%HeaderOldElementBoundaries(ielem))
            call a%StorageOldElementBoundaries%GetNext(testboun)
            do while (testboun /= -1)
               if (a%OldExternalBoundaries(1,testboun) == ielem .and. a%OldExternalBoundaries(2,testboun) == iface) then
                  BoundaryMatch(iboun) = testboun
                  
                  !exit 
                  testboun = -1
               else
                  call a%StorageOldElementBoundaries%GetNext(testboun)
               endif
            enddo
         else
            
         endif
      enddo
   end subroutine
   
   !For the new faces of the mesh, match them to element and face
   subroutine MatchNewLboelToNewElementsAndFaces(a,nboun,pboel,lboel,BoundaryMatch)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: nboun,pboel(nboun+1),lboel(pboel(nboun+1)-1),BoundaryMatch(2,nboun)
      
      integer(ip) :: iboun,ibpos,pnodb,ielemExport,ielem,iface,pnode,ifpos,pface
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: ielty
      
      BoundaryMatch(1,:) = -1
      
      do iboun = 1,nboun
         ibpos = pboel(iboun)-1
         pnodb = pboel(iboun+1)-pboel(iboun)-1
         ielemExport = lboel(ibpos+pnodb+1)
         
         if (ielemExport <= a%ExportNelem) then
            ielem = a%IrenelExport(ielemExport)
         else 
            ielem = 0
         endif
         
         !Not for extra elements
         if (ielem /= 0) then
            call a%GetLnodsFromElement(ielem,pnode,lnods)
            ielty = a%ElementTypeList(ielem)%ielty
            call Nodes_To_Face(ielty,pnodb,lboel(ibpos+1),iface)
            
            !only for external faces (external lboel might not know properly which are external)
            pface = a%pfacs(ielem+1)-a%pfacs(ielem)
            ifpos = a%pfacs(ielem)-1
            if (a%lfacs(1,ifpos+iface) == 0 .or. a%lfacs(1,ifpos+iface) > a%nelem) then
               BoundaryMatch(1,iboun) = ielem
               BoundaryMatch(2,iboun) = iface
            endif
         else 

         endif
      enddo
   end subroutine
         
   subroutine MatchOldBoundariesToNewBoundaries(a,oldnboun,oldBoundaryMatch,newnboun,newBoundaryMatch,BoundaryNewToOld)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun, oldBoundaryMatch(oldnboun),newnboun,newBoundaryMatch(2,newnboun),BoundaryNewToOld(newnboun)
      
      integer(ip), allocatable :: ioldBoundaryMatch(:)
      integer(ip) :: iboun,ielem,iface,oldielem,oldiboun,iparent,oldiparent,pnode,pnodb,nchild,nface,nchildface,iparentface
      integer(ip) :: ichild,matchboun
      integer(ip), pointer :: lnods(:) => NULL()
      integer(ip) :: ielty,ivariation
      
      !We need the inverse for the old one
      allocate(ioldBoundaryMatch(size(a%OldExternalBoundaries,2)))
      ioldBoundaryMatch = -1
      do iboun = 1,oldnboun
         matchboun = oldBoundaryMatch(iboun)
         if (matchboun > 0) then
            ioldBoundaryMatch(oldBoundaryMatch(iboun)) = iboun
         endif
      enddo
      
      do iboun = 1,newnboun
         ielem = newBoundaryMatch(1,iboun)
         !Not for extra elements
         if (ielem /= -1) then
            iface = newBoundaryMatch(2,iboun)
            
            oldielem = a%irenel(ielem)
            !This conditions is fulfilled if I was and I am
            !or if my children were unrefined
            if (oldielem /= 0) then
               !I was in the old one and I am in the new one, but not sure if i was exported in the previous one
               !There are two possibilities:
               !1. I was markel 0 in the previous iteration: I will find myself in positions 1:2
               !2. I was markel 1 in the previous iteration: I will find myself in positions 3:4
               !Loop through boundaries in the list and see if we match
               call a%StorageOldElementBoundaries%GoToFirstOfList(a%HeaderOldElementBoundaries(oldielem))

               call a%StorageOldElementBoundaries%GetNext(oldiboun)
               do while (oldiboun /= -1)
                  !I was markel0 in the previous iteration or I was markel 1
                  if ((a%OldExternalBoundaries(1,oldiboun) == oldielem .and. &
                     a%OldExternalBoundaries(2,oldiboun) == iface) .or. &
                     (a%OldExternalBoundaries(3,oldiboun) == oldielem .and. &
                     a%OldExternalBoundaries(4,oldiboun) == iface)) then
                     
                     BoundaryNewToOld(iboun) = ioldBoundaryMatch(oldiboun)
                     !Exit
                     oldiboun = -1
                  
                  !If not found continue
                  else
                     call a%StorageOldElementBoundaries%GetNext(oldiboun)
                  endif
               enddo
            else
               !I am a new element, seek for my parent, who was refined
               !1. My parent was refined, I need to inherit from him
               
               iparent = a%PointerToParent(ielem)
               oldiparent = a%irenel(iparent)
               call GEI_GetPositionInElementForLevel(a%GEI_ElementList(ielem),int(a%GEI_ElementList(ielem)%level,4),ichild)
               
               ielty = a%ElementTypeList(iparent)%ielty
               ivariation = a%ElementTypeList(iparent)%ivariation
               call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
               iparentface = iface  !entry face for seeking up faces
               call up_faces(ielty,ivariation,nface,ichild,iparentface)
               
               !Now I need to seek for oldiparent, lface
               call a%StorageOldElementBoundaries%GoToFirstOfList(a%HeaderOldElementBoundaries(oldiparent))
               call a%StorageOldElementBoundaries%GetNext(oldiboun)
               do while (oldiboun /= -1)
                  !I was markel0 in the previous iteration or I was markel 1
                  if (a%OldExternalBoundaries(1,oldiboun) == oldiparent .and. &
                     a%OldExternalBoundaries(2,oldiboun) == iparentface) then
                     BoundaryNewToOld(iboun) = ioldBoundaryMatch(oldiboun)
                     !Exit
                     oldiboun = -1
                  
                  !If not found continue
                  else
                     call a%StorageOldElementBoundaries%GetNext(oldiboun)
                  endif
               enddo
            endif
         endif
      enddo   
   end subroutine   
   
   subroutine PrepareRebalanceFaces_a(a)
      !The first part of the subroutine needs to be done when the numbering and elements are still the old ones
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      type(LinkedListHeader), allocatable :: HeaderFacesOthersNeed(:)
      type(LinkedListStorage) :: StorageFacesOthersNeed
      
      integer(ip) :: iboun,ielem,ranki,irank,iexternalBoundary,iface
      integer(ip) :: GEI_Size
      type(GlobalElementIdentifier) :: iGEI
      integer(ip) :: ifoundElement
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: facnod(a%mnodb),ispos,nchild,nchildface,nface,pnodb,pnode,inodb,ipoin
      logical :: IsResponsible
      integer(ip), allocatable :: IsSentFaceToRank(:)
      integer(ip) :: ielty
      type(GlobalElementIdentifier) :: GEI_Sizer
      
      call a%Timers%PrepareRebalanceFaces_a%Tic
      
      
      !Deallocate if necessary
      if (allocated(a%RebalanceNFaceINeed)) then
         deallocate(a%RebalanceNFaceINeed)
         deallocate(a%RebalanceNFaceOthersNeed)
         deallocate(a%RebalanceListFacesINeed)
         deallocate(a%RebalanceListFacesOthersNeed)
      endif
      
      !For each face I am responsible of, I need to send it
      allocate(a%RebalanceNFaceINeed(a%MPIsize))
      allocate(a%RebalanceNFaceOthersNeed(a%MPIsize))
      a%RebalanceNFaceINeed(1:a%MPIsize) = 0
      a%RebalanceNfaceOthersNeed(1:a%MPIsize) = 0
      
     
      
      allocate(a%RebalanceListFacesOthersNeed(a%MPIsize))
      
      allocate(HeaderFacesOthersNeed(a%MPIsize))
      call StorageFacesOthersNeed%Initialize(2*size(a%OldExternalBoundaries,2),'NonRemov','Repeatable',a%Memor)
      
      !With whom do we need to communicate:
      !This are the same processors we needed to send elements to
      !Recover the info
      !Loop through external boundaries, see to which element they belong
      
      allocate(IsSentFaceToRank(a%MPIsize))
      IsSentFaceToRank = 0
      do iboun = 1,size(a%OldExternalBoundaries,2)
         !Check if I am responsible for the BOUNDARY
         ielem = a%OldExternalBoundaries(1,iboun)
         iface = a%OldExternalBoundaries(2,iboun)
         
         
         ispos = a%pnods(ielem)-1
         pnode = a%pnods(ielem+1)-a%pnods(ielem)
         ielty = a%ElementTypeList(ielem)%ielty
         call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
         call face_from_el(ielty,nface,pnodb,a%lnods(ispos+1),iface,0,facnod)
!          !Check if I am responsible for the boundary
!          !NO: call a%CheckElementResponsibility(pnodb,facnod,IsResponsible)
!          
!          !If I am responsible for the element then I Need to send it! 
!          !Not for the face, because I might be responsible for the face but not for the element
!          !But then I do not know to whom I have to send it. 
!          !So we stick to the element responsibility criteria
!          !Loop through size of RebalanceElementProcessors
!          IsResponsible = .true.
!          if (IsResponsible) then
!             !add face to the list of the ranks the element is responsible of
!             do ranki = 1,size(a%RebalanceElementProcessors(ielem)%l)
!                irank = a%RebalanceElementProcessors(ielem)%l(ranki)
!                a%RebalanceNfaceOthersNeed(irank+1) = a%RebalanceNfaceOthersNeed(irank+1)+1
!                
!                call StorageFacesOthersNeed%Add(HeaderFacesOthersNeed(irank+1),iboun)
!             enddo
!          endif
         
         !Check if I am responsible for the boundary
         !I only need to send it to the ranks of the boundary nodes
         !This is for matching EXTERNAL INFO, we do not require the hierarchy
         call a%CheckElementResponsibility(pnodb,facnod,IsResponsible)

         if (IsResponsible) then
            !add face to the list of the ranks of the boundary nodes
            do inodb = 1,pnodb
               ipoin = facnod(inodb)
               irank = a%RebalanceProcessorList(ipoin)
               if (IsSentFaceToRank(irank+1) /= iboun) then
                  IsSentFaceToRank(irank+1) = iboun
                  
                  a%RebalanceNfaceOthersNeed(irank+1) = a%RebalanceNfaceOthersNeed(irank+1)+1
                  call StorageFacesOthersNeed%Add(HeaderFacesOthersNeed(irank+1),iboun)
               endif
            enddo   
         endif
      enddo
      deallocate(IsSentFaceToRank)
      
      !Extract from linked list
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceOthersNeed(irank+1) /= 0) then
            allocate(a%RebalanceListFacesOthersNeed(irank+1)%Faces(a%RebalanceNfaceOthersNeed(irank+1)))
            call StorageFacesOthersNeed%ListToArray(HeaderFacesOthersNeed(irank+1),a%RebalanceListFacesOthersNeed(irank+1)%Faces)
            allocate(a%RebalanceListFacesOthersNeed(irank+1)%GEI(a%RebalanceNfaceOthersNeed(irank+1)))
            allocate(a%RebalanceListFacesOthersNeed(irank+1)%iface(a%RebalanceNfaceOthersNeed(irank+1)))
            do iboun = 1,a%RebalanceNfaceOthersNeed(irank+1)
               iexternalBoundary = a%RebalanceListFacesOthersNeed(irank+1)%Faces(iboun)
               ielem = a%OldExternalBoundaries(1,iexternalBoundary)
               iface = a%OldExternalBoundaries(2,iexternalBoundary)
               a%RebalanceListFacesOthersNeed(irank+1)%GEI(iboun) = a%GEI_ElementList(ielem)
               a%RebalanceListFacesOthersNeed(irank+1)%iface(iboun) = iface
            enddo
         endif
      enddo
         
      deallocate(HeaderFacesOthersNeed)
      call StorageFacesOthersNeed%Dealloc
      
      !Allocate and Initialize
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      !Communicate number of boundaries to communicate with rebalancing neighbours only
      !But I still do not know the rebalancing neighbours: AllToAll?
      call MPI_AllToAll(a%RebalanceNFaceOthersNeed, 1, MPI_INTEGER4, a%RebalanceNFaceINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
!       do irank = 0,a%MPIsize-1
!          if (a%RebalanceNelemOthersNeed(irank+1) /= 0) then
!             call MPI_ISEND(a%RebalanceNfaceOthersNeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
!             
!          endif
!          if (a%RebalanceNelemINeed(irank+1) /= 0) then
!             call MPI_IRECV(a%RebalanceNfaceINeed(irank+1),1_ip,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
!          endif
!       enddo
!       call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
!       call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Allocates
      allocate(a%RebalanceListFacesINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceINeed(irank+1) /= 0) then
            allocate(a%RebalanceListFacesINeed(irank+1)%Elements(a%RebalanceNFaceINeed(irank+1)))
            allocate(a%RebalanceListFacesINeed(irank+1)%GEI(a%RebalanceNFaceINeed(irank+1)))
            allocate(a%RebalanceListFacesINeed(irank+1)%iface(a%RebalanceNFaceINeed(irank+1)))
         endif
      enddo   
         
      !Send and receive the GEI and the face number
      !For knowing the sizes to send
      GEI_Size = sizeof(GEI_Sizer)
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL
      !Communicate number of boundaries to communicate with rebalancing neighbours only
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(a%RebalanceListFacesOthersNeed(irank+1)%GEI,a%RebalanceNfaceOthersNeed(irank+1)*GEI_Size,MPI_BYTE, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(a%RebalanceListFacesOthersNeed(irank+1)%iface,a%RebalanceNfaceOthersNeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest03(irank+1), ierr)
            
         endif
         if (a%RebalanceNfaceINeed(irank+1) /= 0) then
            call MPI_IRECV(a%RebalanceListFacesINeed(irank+1)%GEI,a%RebalanceNFaceINeed(irank+1)*GEI_Size,MPI_BYTE, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
            call MPI_IRECV(a%RebalanceListFacesINeed(irank+1)%iface,a%RebalanceNFaceINeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest04(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
      
      call a%Timers%PrepareRebalanceFaces_a%Toc
      
   end subroutine
      
   subroutine PrepareRebalanceFaces_b(a)
      !The second part of the subroutine needs to be done after the elements have been transferred
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
   
      type(LinkedListHeader), allocatable :: HeaderFacesOthersNeed(:)
      type(LinkedListStorage) :: StorageFacesOthersNeed
      
      integer(ip) :: iboun,ielem,ranki,irank,iexternalBoundary,iface
      integer(ip) :: GEI_Size
      type(GlobalElementIdentifier) :: iGEI
      integer(ip) :: ifoundElement
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      call a%Timers%PrepareRebalanceFaces_b%Tic
   
      !Find the local element numbers from the GEI
      do irank = 0,a%MPIsize-1
         do iboun = 1,a%RebalanceNfaceINeed(irank+1)
            iGEI = a%RebalanceListFacesINeed(irank+1)%GEI(iboun)
            call a%LocalElementFromGEI(iGEI,ifoundElement)
            a%RebalanceListFacesINeed(irank+1)%Elements(iboun) = ifoundElement
         enddo
      enddo
      
      call a%Timers%PrepareRebalanceFaces_b%Toc
   end subroutine
   
   subroutine GetiBoundaryMatch(a,nboun,BoundaryMatch,iBoundaryMatch)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: nboun
      integer(ip) :: BoundaryMatch(nboun), iBoundaryMatch(*)
      
      integer(ip) :: iboun
      
      iBoundaryMatch(1:size(a%OldExternalBoundaries,2)) = -1
      do iboun = 1,nboun
         if (BoundaryMatch(iboun) /= -1) then
            iBoundaryMatch(BoundaryMatch(iboun)) = iboun
         endif
      enddo
   end subroutine
   
   subroutine GetiBoundaryToElementsAndFacesMatch(a,nboun,BoundaryMatch,iBoundaryMatch,Memor)
      use Mod_Memor
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: nboun
      integer(ip) :: BoundaryMatch(2,nboun)
      type(i1p), allocatable :: iBoundaryMatch(:)
      type(MemoryMan) :: Memor
      
      integer(ip) :: iboun,ielem,iface,pface
      
      call a%Memor%alloc(a%nelem,iBoundaryMatch,'iBoundaryMatch','GetiBoundaryToElementsAndFacesMatch')
      !Allocate
      do ielem = 1,a%nelem
         pface = a%pfacs(ielem+1)-a%pfacs(ielem)
         allocate(iBoundaryMatch(ielem)%l(pface))
         iBoundaryMatch(ielem)%l = 0
      enddo
      call a%Memor%AllocObj(0,'iBoundaryMatch','GetiBoundaryToElementsAndFacesMatch',ip*(a%pfacs(a%nelem+1)-1))
      
      !Second to fill
      do iboun = 1,nboun
         if (BoundaryMatch(1,iboun) /= -1) then
            ielem = BoundaryMatch(1,iboun)
            iface = BoundaryMatch(2,iboun)
            iBoundaryMatch(ielem)%l(iface) = iboun
         endif
      enddo
   end subroutine

   subroutine RebalanceVariableIntegerBoundary1(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*)
      type(i1p)   :: inewBoundaryMatch(*)
      integer(ip), target :: OldArray(*),NewArray(*)
      
      integer(ip), pointer :: auxOldArray(:,:) => NULL(), auxNewArray(:,:) => NULL()
      
      auxOldArray(1:1,1:1) => OldArray(1:1)
      auxNewArray(1:1,1:1) => NewArray(1:1)
      
      call RebalanceVariableIntegerBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,auxOldArray,auxNewArray)
      
   end subroutine   
   
   subroutine RebalanceVariableIntegerBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*),OldArray(ndofn,*),NewArray(ndofn,*)
      type(i1p)   :: inewBoundaryMatch(*)
      
      type(i2p), allocatable :: auxCommunicatorOthersNeed(:), auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: ielem, iboun,newboun,iface,facei,iindex
      
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%l(ndofn,a%RebalanceNfaceOthersNeed(irank+1)))
            auxCommunicatorOthersNeed(irank+1)%l = 0
            do facei = 1,size(a%RebalanceListFacesOthersNeed(irank+1)%Faces)
               iface = a%RebalanceListFacesOthersNeed(irank+1)%Faces(facei)
               iindex = ioldBoundaryMatch(iface)
               if (iindex /= -1) then
                  auxCommunicatorOthersNeed(irank+1)%l(:,facei) = OldArray(:,iindex)
               endif
            enddo
         endif
         if (a%RebalanceNfaceINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%l(ndofn,a%RebalanceNfaceINeed(irank+1)))
         endif
      enddo
      
      !Communicate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNFaceOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%l,ndofn*a%RebalanceNFaceOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%RebalanceNFaceINeed(irank+1) /= 0) then
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%l,ndofn*a%RebalanceNFaceINeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      do irank = 0,a%MPIsize-1
         do iboun = 1,a%RebalanceNfaceINeed(irank+1)
            ielem = a%RebalanceListFacesINeed(irank+1)%Elements(iboun)
            iface = a%RebalanceListFacesINeed(irank+1)%iface(iboun)
            
            newboun = inewBoundaryMatch(ielem)%l(iface)
            
            !Sent boundary might not be considered a boundary by the external Mesh, check
            if (newboun /= 0) then
               NewArray(:,newboun) = auxCommunicatorINeed(irank+1)%l(:,iboun)
            endif   
         enddo
      enddo
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         if (associated(auxCommunicatorINeed(irank+1)%l)) then
            deallocate(auxCommunicatorINeed(irank+1)%l)
         endif
         if (associated(auxCommunicatorOthersNeed(irank+1)%l)) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%l)
         endif
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
   end subroutine
   
   subroutine RebalanceVariableRealBoundary1(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*)
      real(rp), target :: OldArray(*),NewArray(*)
      type(i1p)   :: inewBoundaryMatch(*)
      
      real(rp), pointer :: auxOldArray(:,:) => NULL(), auxNewArray(:,:) => NULL()
      
      auxOldArray(1:1,1:1) => OldArray(1:1)
      auxNewArray(1:1,1:1) => NewArray(1:1)
      
      call RebalanceVariableRealBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,auxOldArray,auxNewArray)
      
   end subroutine  
   
   subroutine RebalanceVariableRealBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*)
      real(rp) ::OldArray(ndofn,*),NewArray(ndofn,*)
      type(i1p)   :: inewBoundaryMatch(*)
      
      type(r2p), allocatable :: auxCommunicatorOthersNeed(:), auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: iboun,ielem,iface,newboun,iindex,facei
      
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%a(ndofn,a%RebalanceNfaceOthersNeed(irank+1)))
            auxCommunicatorOthersNeed(irank+1)%a = 0
            do facei = 1,size(a%RebalanceListFacesOthersNeed(irank+1)%Faces)
               iface = a%RebalanceListFacesOthersNeed(irank+1)%Faces(facei)
               iindex = ioldBoundaryMatch(iface)
               if (iindex /= -1) then
                  auxCommunicatorOthersNeed(irank+1)%a(:,facei) = OldArray(:,iindex)
               endif
            enddo
         endif
         if (a%RebalanceNfaceINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%a(ndofn,a%RebalanceNfaceINeed(irank+1)))
         endif   
      enddo
      
      !Communicate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNFaceOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%a,ndofn*a%RebalanceNFaceOthersNeed(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%RebalanceNFaceINeed(irank+1) /= 0) then
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%a,ndofn*a%RebalanceNFaceINeed(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      !Communications structure to new
      do irank = 0,a%MPIsize-1
         do iboun = 1,a%RebalanceNfaceINeed(irank+1)
            ielem = a%RebalanceListFacesINeed(irank+1)%Elements(iboun)
            iface = a%RebalanceListFacesINeed(irank+1)%iface(iboun)
            
            newboun = inewBoundaryMatch(ielem)%l(iface)
            
            if (newboun /= 0) then
               NewArray(:,newboun) = auxCommunicatorINeed(irank+1)%a(:,iboun)
            endif   
         enddo
      enddo
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         deallocate(auxCommunicatorINeed(irank+1)%a)
         deallocate(auxCommunicatorOthersNeed(irank+1)%a)
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
   end subroutine   
   
   subroutine RebalanceVariableTyper1pBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,OldArray,NewArray,AllocatedBytes)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*)
      type(r1p)   :: OldArray(*),NewArray(*)
      type(i1p)   :: inewBoundaryMatch(*)
      integer(ip) :: AllocatedBytes
      
      type(CSRList), allocatable :: auxCommunicatorOthersNeed(:), auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      integer(ip) :: ielem, iboun,newboun,iface,icount,boundOldi,iBoundOld,nface,isize,ispos,pnode,auxtest
      
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNfaceOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%p(a%RebalanceNfaceOthersNeed(irank+1)+1))
            !FirstToCount
            icount = 0
            do BoundOldi = 1,a%RebalanceNfaceOthersNeed(irank+1)
               iBoundOld = a%RebalanceListFacesOthersNeed(irank+1)%Faces(BoundOldi)
               iBoundOld = ioldBoundaryMatch(iBoundOld)
               if (iBoundOld /= -1 .and. associated(OldArray(iBoundOld)%a)) then
                  icount = icount + size(OldArray(iBoundOld)%a)
               endif
            enddo
            allocate(auxCommunicatorOthersNeed(irank+1)%r(icount))
            !Second to fill
            auxCommunicatorOthersNeed(irank+1)%p(1) = 1
            do BoundOldi = 1,a%RebalanceNfaceOthersNeed(irank+1)
               iBoundOld = a%RebalanceListFacesOthersNeed(irank+1)%Faces(BoundOldi)
               iBoundOld = ioldBoundaryMatch(iBoundOld)
               if (iBoundOld /= -1 .and. associated(OldArray(iBoundOld)%a)) then
                  pnode = size(OldArray(iBoundOld)%a)
               else
                  pnode = 0
               endif
               ispos = auxCommunicatorOthersNeed(irank+1)%p(BoundOldi)-1
               auxCommunicatorOthersNeed(irank+1)%p(boundOldi+1) = &
                  auxCommunicatorOthersNeed(irank+1)%p(BoundOldi) + pnode
               if (pnode /= 0) then   
                  auxCommunicatorOthersNeed(irank+1)%r(ispos+1:ispos+pnode) = OldArray(iBoundOld)%a
               endif
            enddo
         endif
      enddo
      
      !Communicate
      !First I will send the p, and then I can allocate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNFaceOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%p,a%RebalanceNFaceOthersNeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%RebalanceNFaceINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%p(a%RebalanceNFaceINeed(irank+1)+1))
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%p,a%RebalanceNFaceINeed(irank+1)+1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      !Allocate
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNFaceOthersNeed(irank+1) /= 0) then
            nface = a%RebalanceNFaceOthersNeed(irank+1)
            isize = auxCommunicatorOthersNeed(irank+1)%p(nface+1)-1
            if (isize /= 0) then
               call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%r,isize,MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            endif
         endif
         if (a%RebalanceNFaceINeed(irank+1) /= 0) then
            nface = a%RebalanceNFaceINeed(irank+1)
            isize = auxCommunicatorINeed(irank+1)%p(nface+1)-1
            if (isize /= 0) then
               allocate(auxCommunicatorINeed(irank+1)%r(isize))
               call MPI_IRECV(auxCommunicatorINeed(irank+1)%r,isize,MPI_REAL8, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
            endif
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      icount = 0
      do irank = 0,a%MPIsize-1
         do iboun = 1,a%RebalanceNfaceINeed(irank+1)
            ielem = a%RebalanceListFacesINeed(irank+1)%Elements(iboun)
            iface = a%RebalanceListFacesINeed(irank+1)%iface(iboun)
            
            ispos = auxCommunicatorINeed(irank+1)%p(iboun)-1
            pnode = auxCommunicatorINeed(irank+1)%p(iboun+1) - auxCommunicatorINeed(irank+1)%p(iboun)
            
            newboun = inewBoundaryMatch(ielem)%l(iface)
            if (newboun /= 0) then
               allocate(NewArray(newboun)%a(pnode))
               icount = icount + pnode
               if (pnode /= 0) then
                  NewArray(newboun)%a = auxCommunicatorINeed(irank+1)%r(ispos+1:ispos+pnode)
               endif
            endif
         enddo
      enddo
      AllocatedBytes = icount*rp
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         if (allocated(auxCommunicatorINeed(irank+1)%p)) then
            deallocate(auxCommunicatorINeed(irank+1)%p)
         endif
         if (allocated(auxCommunicatorINeed(irank+1)%r)) then
            deallocate(auxCommunicatorINeed(irank+1)%r)
         endif
         if (allocated(auxCommunicatorOthersNeed(irank+1)%p)) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%p)
         endif
         if (allocated(auxCommunicatorOthersNeed(irank+1)%r)) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%r)
         endif
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
      
   end subroutine
   
   subroutine Set_2_1_Balancing(a,kfl_21Balancing)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: kfl_21Balancing
      
      a%kfl_21Balancing = kfl_21Balancing
   end subroutine
   
   subroutine SetRebalancingNumbering(a,PointGlobNumber,PointProcNumber)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: PointGlobNumber(*),PointProcNumber(*)
      
      integer(ip) :: ipoin, externalIpoin
      
      if (allocated(a%RebalancePointNumbering)) then
         deallocate(a%RebalancePointNumbering)
         deallocate(a%RebalanceProcessorList)
      endif
      
      allocate(a%RebalancePointNumbering(a%npoin))
      allocate(a%RebalanceProcessorList(a%npoin))
      
      do ipoin = 1,a%npoin
         externalIpoin = a%ExportPointNumbering(ipoin)
         if (externalIpoin /= 0) then
            a%RebalancePointNumbering(ipoin) = PointGlobNumber(externalIpoin)
            a%RebalanceProcessorList(ipoin) = PointProcNumber(externalIpoin)
         endif
      enddo
      
      !GhostCommunicate new global numbering and processors for all ghost nodes
      call a%Communicator%GhostCommunicate(1_ip,a%RebalancePointNumbering)
      call a%Communicator%GhostCommunicate(1_ip,a%RebalanceProcessorList)
   end subroutine
   
   subroutine WriteTimes(a,lunout,cputotal)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: lunout
      real(rp) :: cputotal
      
      real(rp) :: timeValues(100),cpurefine
      integer(ip) :: itime
      
      call a%Timers%Initialize%GetValue(timeValues(1))
      call a%Timers%PrepareOldExternalBoundariesInfoToKeep%GetValue(timeValues(2))
      call a%Timers%Refine%GetValue(timeValues(3))
         call a%Timers%AmendMarkel%GetValue(timeValues(4))
         call a%Timers%Do21Balancing%GetValue(timeValues(5))
         call a%Timers%RefineSerial%GetValue(timeValues(6))
         call a%Timers%NumberingToParallel%GetValue(timeValues(7))
         call a%Timers%UpdateArraysToNewParallelNumbering%GetValue(timeValues(8))
         call a%Timers%HangingNodes%GetValue(timeValues(9))
         call a%Timers%SelectComponentsToExport%GetValue(timeValues(10))
         call a%Timers%MarkDeadElements%GetValue(timeValues(11))
      call a%Timers%LoadRebalance%GetValue(timeValues(12))
         call a%Timers%RebalanceRenumbering%GetValue(timeValues(13))
         call a%Timers%RebalanceElemProcs%GetValue(timeValues(14))
         call a%Timers%PrepareRebalanceFaces_a%GetValue(timeValues(15))
         call a%Timers%RebalanceSendNodesAndElements%GetValue(timeValues(16))
         call a%Timers%RebuildRefinementStructure%GetValue(timeValues(17))
         call a%Timers%PrepareRebalanceFaces_b%GetValue(timeValues(18))
         call a%Timers%PrepareRebalanceCommunicator%GetValue(timeValues(19))
         
       cpurefine = timeValues(1)+timeValues(3) + timeValues(12)  
       
       write(lunout,100) cpurefine, 100.0_rp*cpurefine/cputotal, &
         (timeValues(itime),100.0_rp*timeValues(itime)/cpurefine,itime=1,19)
      
   100 format(&            
            10x,'   REFINER                     :',F12.2,' (',F6.2,' % )',/,&
            10x,'     Initialize             :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     OldExternalBoundaries  :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Refine                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        AmendMarkel         :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        Do21Balancing       :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        RefineSerial        :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        NumberingToParallel :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        UpdateArraysParaNum :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        HangingNodes        :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        SelectComp.ToExport :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        MarkDeadElements    :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'    LoadRebalance           :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        Renumbering         :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        ElemProcs           :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        Faces_a             :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        SendNodesAndElements:   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        RebuildRefinementStr:   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        Faces_b,&           :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'        Communicator        :   ',F12.2,' (',F6.2,' % )')
   
   end subroutine
   
   subroutine Dealloc(a)
      use typre
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      
      integer(ip) :: i,irank
      logical     :: kfl_isallocated
      integer(ip), allocatable :: markel(:)
      
      !This has not been refined
      if (a%NotRefined) then
         allocate(markel(a%ExportNelem))
         markel = 0
         call a%PrepareOldExternalBoundariesInfoToKeep  
         call a%Refine(markel)
         deallocate(markel)
      endif
      a%NotRefined = .true.
      
      
      if (associated(a%OldLocalOrdering)) then
         call a%OldLocalOrdering%Deallocate(a%Memor)
         deallocate(a%OldLocalOrdering)
      endif
      
      call a%LocalOrdering%Deallocate(a%Memor)
      deallocate(a%LocalOrdering)
      call a%ElementLocalOrdering%Deallocate(a%Memor)
      call a%StorageOldElementBoundaries%Dealloc
      deallocate(a%HeaderOldElementBoundaries)
      
      
      call a%Communicator%Deallocate
      call a%ExportCommunicator%Deallocate
      
      deallocate(a%pPointerToChildren,a%lPointerToChildren,a%PointerToParent,a%pedger,a%ledger,a%OldPointerToParent)
      deallocate(a%level,a%renpo,a%irenpo,a%pnods,a%lnods,a%markel,a%renel,a%irenel,a%pfacs,a%lfacs)
      if (allocated(a%oldmarkel)) deallocate(a%oldmarkel)
      !deallocate(a%EdgerEntryPoints)
      if (allocated(a%OldGEI_ElementList)) deallocate(a%OldGEI_ElementList)
      deallocate(a%GEI_ElementList,a%pHangingList,a%lHangingList,a%rHangingList)
      deallocate(a%ElementTypeList)
      deallocate(a%HangingListForElements%p,a%HangingListForElements%l)
      deallocate(a%ExportElement,a%ExportPointNumbering,a%iRenelExport,a%OldExportPointNumbering,a%ExportElementNumbering,a%OldExportElementNumbering)
      deallocate(a%ExtraPoints,a%ExtraHangingList%p,a%ExtraHangingList%l,a%ExtraHangingList%r)
      deallocate(a%ExtraElements%p,a%ExtraElements%l)
      
      if (allocated(a%RebalancePointNumbering)) then
         deallocate(a%RebalancePointNumbering,a%RebalanceProcessorList)
         deallocate(a%RebalanceNelemINeed,a%RebalanceNelemOthersNeed)
         deallocate(a%RebalanceNfaceINeed,a%RebalanceNfaceOthersNeed)
         
         do i = 1,size(a%RebalanceElementProcessors)
            if (associated(a%RebalanceElementProcessors(i)%l)) then
               deallocate(a%RebalanceElementProcessors(i)%l)
            endif
         enddo
         deallocate(a%RebalanceElementProcessors)
         
         do i = 1,size(a%RebalanceListFacesOthersNeed)
            if (allocated(a%RebalanceListFacesOthersNeed(i)%Faces)) then
               deallocate(a%RebalanceListFacesOthersNeed(i)%Faces)
               deallocate(a%RebalanceListFacesOthersNeed(i)%GEI)
               deallocate(a%RebalanceListFacesOthersNeed(i)%iface)
            endif
         enddo
         deallocate(a%RebalanceListFacesOthersNeed)
         
         do i = 1,size(a%RebalanceListFacesINeed)
            if (allocated(a%RebalanceListFacesINeed(i)%Elements)) then
               deallocate(a%RebalanceListFacesINeed(i)%Elements)
               deallocate(a%RebalanceListFacesINeed(i)%GEI)
               deallocate(a%RebalanceListFacesINeed(i)%iface)
            endif
         enddo
         deallocate(a%RebalanceListFacesINeed)
         
         
      endif
      deallocate(a%OldExternalBoundaries)
      
      if (allocated(a%RebalanceListElementsINeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%RebalanceListElementsINeed(irank+1)%l)
         enddo
         deallocate(a%RebalanceListElementsINeed)
      endif
      
      if (allocated(a%RebalanceListElementsOthersNeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%RebalanceListElementsOthersNeed(irank+1)%l)
         enddo
         deallocate(a%RebalanceListElementsOthersNeed)
      endif
      
      if (allocated(a%ExtraElementsListOthersNeed)) then
         do irank = 0,a%MPIsize-1
            deallocate(a%ExtraElementsListOthersNeed(irank+1)%l)
         enddo
         deallocate(a%ExtraElementsListOthersNeed)
      endif
      
      if (allocated(a%ExtraNelemINeed)) then
         deallocate( a%ExtraNelemINeed,a%ExtraNelemOthersNeed)
      endif
      
      if (allocated(a%lJustRefinedElements)) deallocate(a%lJustRefinedElements)
      if (allocated(a%OldChildNumber)) deallocate(a%OldChildNumber)

      if (allocated(a%isHangingNode)) deallocate (a%isHangingNode)

      call a%RebalanceComm%IsAllocated(kfl_isallocated)
      if (kfl_isallocated .eqv. .true.) then
         call a%RebalanceComm%Deallocate
      endif
      
   end subroutine
   
   subroutine ElementalRefine(a,ARET)
      use Mod_AdaptiveElementDependentSubroutines
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerElementTransferer) :: ARET
      
      integer(ip) :: ielem,inewel,nchild,childi, ichild,ielty,ispos,ivariation,ExportIelem,ExportInewel,ExportIchild
      integer(ip) :: ExportIparent,inewparent,ioldparent

      do ielem = 1,a%oldnelem
         !Only for originally exported elements
         ExportIelem = a%OldExportElementNumbering(ielem)
         if (ExportIelem > 0) then
            inewel = a%renel(ielem)
            if (inewel > 0) then
               ExportInewel = a%ExportElementNumbering(inewel)
               nchild = a%pPointerToChildren(inewel+1)-a%pPointerToChildren(inewel)
               if (nchild == 0) then
                  !Copy
                  if (ExportInewel > 0) then
                     call ARET%CopyElement(ExportIelem,ExportInewel)
                  endif
               else
                  !Interpolate
                  ispos = a%pPointerToChildren(inewel)-1
                  ivariation = a%ElementTypeList(inewel)%ivariation
                  ielty = a%ElementTypeList(inewel)%ielty
                  do ichild = 1,nchild
                     childi = a%lPointerToChildren(ispos+ichild)
                     ExportIchild = a%ExportElementNumbering(childi)
                     if (ExportIChild > 0) then
                        a%ARET_ielty = ielty
                        a%ARET_ivariation = ivariation
                        a%ARET_ichild = ichild
                        a%ARET_InterpolationTask = 1
                        call ARET%NodallyInterpolateElement(a,ExportIelem,ExportIchild)
                     endif
                  enddo
               endif
   
            else
               !The element was unrefined, we need to contribute to the parent
               ioldparent = a%OldPointerToParent(ielem)
               inewparent = a%renel(ioldparent)
               ivariation = a%ElementTypeList(inewparent)%ivariation
               ielty = a%ElementTypeList(ivariation)%ielty
               ExportIparent = a%ExportElementNumbering(inewparent)
               !The parent might not be exported now...
               if (ExportIparent > 0) then
                  ichild = a%OldChildNumber(ielem)
                  a%ARET_ielty = ielty
                  a%ARET_ivariation = ivariation
                  a%ARET_ichild = ichild
                  a%ARET_InterpolationTask = 2
                  call ARET%NodallyInterpolateElement(a,ExportIelem,ExportIparent)
               endif
               
            endif
         endif
      enddo
      
      !Communicate the unrefine data
      call ElementalUnref(a,ARET)
      
      !Communicate the extra elemental data!!
      call ElementalExtra(a,ARET)
      
   end subroutine
   
   
   subroutine ElementalUnref(a,ARET)
      implicit none
      class(IsotropicAdaptiveRefiner), target :: a
      class(AdaptiveRefinerElementTransferer) :: ARET
      
      type(GlobalElementIdentifier) :: iGEI
      type(GlobalElementIdentifier) :: GEI_Sizer
      
      integer(ip), pointer :: lnods(:) => NULL()
      
      integer(ip) :: parentlprocs(50),parentpnode,childpnode,childlprocs(50),ielem, ielty, inode,parentnode,ioldparent,iproc,irank, ichild,ispos,ivariation,sendproc
      
      type(ListGEI), allocatable :: GEI_NelemINeed(:),GEI_NelemOthersNeed(:)
      integer(ip) :: GEI_Size
      integer(ip), allocatable :: ToReceiveElements(:),ToSendElements(:), NelemOthersNeed(:), NelemINeed(:), BufferSizeToSend(:), BufferSizeToReceive(:)
      type(r1p), allocatable :: SendBuffer(:), ReceiveBuffer(:)
      integer(ip) :: Bsize, BufferCount, ExportIparent,icount,ierr,inewparent,ilevel, exportielem
      
      
      type(LinkedListHeader), allocatable :: HeaderListToSend(:)
      type(LinkedListStorage) :: StorageListToSend
      !MPI parameters
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      
   
      GEI_Size = sizeof(GEI_Sizer)
   
      allocate(ToReceiveElements(a%MPIsize),ToSendElements(a%MPIsize),NelemOthersNeed(a%MPIsize), NelemIneed(a%MPisize),BufferSizeToSend(a%MPisize),BufferSizeToReceive(a%MPisize))
      NelemOthersNeed = 0_ip
      NelemINeed = 0_ip
      ToSendElements = -1_ip
      ToReceiveElements = -1_ip
      BufferSizeToSend = 0_ip
      BufferSizeToReceive = 0_ip
      
      
      allocate(HeaderListToSend(a%MPIsize))
      call StorageListToSend%Initialize(2*a%npoinGhost,'NonRemovable','Repeatable',a%Memor)
      
      !Unrefined elements (not dead) make sure all parents have what they need!
      !first loop to count
      do ielem = 1,a%oldnelem
         if (a%oldmarkel(ielem) == -1) then
            ioldparent = a%OldPointerToParent(ielem)
            inewparent = a%renel(ioldparent)
            ExportIparent = a%ExportElementNumbering(inewparent)
            !For each processor of the parent element, identify if I was exported or not
            parentpnode = a%oldpnods(ioldParent+1)-a%oldpnods(ioldParent)
            ispos = a%oldpnods(ioldparent)-1
            lnods => a%oldlnods(ispos+1:ispos+parentpnode)
            call a%OldLocalOrdering%GetProc(parentpnode,lnods,parentlprocs)
            
            childpnode = a%oldpnods(ielem+1)-a%oldpnods(ielem)
            ispos = a%oldpnods(ielem)-1
            lnods => a%oldlnods(ispos+1:ispos+parentpnode)
            call a%OldLocalOrdering%GetProc(childpnode,lnods,childlprocs)
            
            do inode = 1,parentpnode
               iproc = parentlprocs(inode)
               if (any(childlprocs(1:childpnode) == iproc)) then
                  !Everything is fine, no need to send info, info is already there
               else
                  !If iproc is myself, add it to the list to receive
                  if (iproc == a%MPIrank) then
                     sendproc = minval(childlprocs(1:childpnode))
                     if (ToReceiveElements(sendproc+1) /= ielem) then
                        NelemINeed(sendproc+1) = NelemINeed(sendproc+1)+1
                        ToReceiveElements(sendproc+1) = ielem
                        call ARET%GetNewElementBufferSize(exportIparent,BSize)
                        BufferSizeToReceive(sendproc+1) = BufferSizeToReceive(sendproc+1)+BSize
                     endif   
                  
                  
                  !Otherwise, if I am responsible for the element, add it to the list to send
                  elseif (minval(childlprocs(1:childpnode)) == a%MPIrank) then
                     if (ToSendElements(iproc+1) /= ielem) then
                        NelemOthersNeed(iproc+1) = NelemOthersNeed(iproc+1)+1
                        ToSendElements(iproc+1) = ielem
                        exportielem = a%OldExportElementNumbering(ielem)
                        call ARET%GetOldElementBufferSize(exportielem,BSize)
                        BufferSizeToSend(iproc+1) = BufferSizeToSend(iproc+1)+BSize
                        call StorageListToSend%Add(HeaderListToSend(iproc+1),ielem)
                     endif   
                  endif
               endif
            enddo
         endif
      enddo
      
      !Allocate fill and send the buffers
      allocate(SendBuffer(a%MPIsize),ReceiveBuffer(a%MPIsize))
      allocate(GEI_NelemINeed(a%MPIsize),GEI_NelemOthersNeed(a%MPIsize))
      
      !Allocate communications arrays
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL
      
      !Receive the buffers and process them
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1) /= 0) then
            BufferCount = 1
            allocate(SendBuffer(irank+1)%a(BufferSizeToSend(irank+1)))
            allocate(GEI_NelemOthersNeed(irank+1)%a(NelemOthersNeed(irank+1)))
            
            !Load from linked list to transfer list
            call StorageListToSend%GoToFirstOfList(HeaderListToSend(irank+1))
            call StorageListToSend%GetNext(ielem)
            icount = 1
            do while (ielem /= -1)
               GEI_NelemOthersNeed(irank+1)%a(icount) = a%OldGEI_ElementList(ielem)
               exportielem = a%OldExportElementNumbering(ielem)
               call ARET%GetOldElementBufferSize(exportielem,BSize)
               call ARET%OldElementToBuffer(exportielem,SendBuffer(irank+1)%a(BufferCount:BufferCount))
               BufferCount = BufferCount+ Bsize
               
               icount = icount +1
               call StorageListToSend%GetNext(ielem)
            enddo
            
            !Send GEI and buffer
            call MPI_ISEND(SendBuffer(irank+1)%a,BufferSizeToSend(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(GEI_NelemOthersNeed(irank+1)%a,NelemOthersNeed(irank+1)*GEI_Size,MPI_BYTE, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
            
         endif   
         
         !Receive!!
         if (NelemIneed(irank+1) /= 0) then
            allocate(ReceiveBuffer(irank+1)%a(BufferSizeToReceive(irank+1)))
            allocate(GEI_NelemINeed(irank+1)%a(NelemINeed(irank+1)))
            
            call MPI_IRECV(ReceiveBuffer(irank+1)%a,BufferSizeToReceive(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest04(irank+1), ierr)
            call MPI_IRECV(GEI_NelemINeed(irank+1)%a,NelemINeed(irank+1)*GEI_Size,MPI_byte,irank,mtag2,a%MPIcomm,irequest03(irank+1),ierr) 
         endif
         
         
      enddo
      
      !Wait for communications to finish
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
      
      
      
      do irank = 0,a%MPIsize-1
         BufferCount = 1
         do ielem = 1,NelemINeed(irank+1)
            iGEI =  GEI_NelemINeed(irank+1)%a(ielem)
            
            !which child I am
            ilevel = iGEI%level
            call GEI_GetPositionInElementForLevel(iGEI,ilevel,ichild)

            !I want to look for the parent in the current numbering
            !I modify the GEI level so that I find the parent instead of the element
            iGEI%level = iGEI%level-1
            call a%LocalElementFromGEI(iGEI,inewparent)
            iGEI%level = iGEI%level+1
            
            
            ExportIparent = a%ExportElementNumbering(inewparent)
            call ARET%GetNewElementBufferSize(ExportIparent,BSize)
            
            !The element was unrefined, we need to contribute to the parent
            ivariation = a%ElementTypeList(inewparent)%ivariation
            ielty = a%ElementTypeList(ivariation)%ielty
            ExportIparent = a%ExportElementNumbering(inewparent)
            !The parent might not be exported now...
            if (ExportIparent > 0) then
               a%ARET_ielty = ielty
               a%ARET_ivariation = ivariation
               a%ARET_ichild = ichild
               a%ARET_InterpolationTask = 2
               call ARET%NodallyInterpolateElementFromOldBuffer(a,ExportIparent,ReceiveBuffer(irank+1)%a(BufferCount:BufferCount+Bsize-1))
            endif
            BufferCount = BufferCount +Bsize
         enddo
      enddo   
         
         
      do irank = 0,a%MPIsize-1
         if (NelemOthersNeed(irank+1) /= 0) then
            deallocate(SendBuffer(irank+1)%a)
            deallocate(GEI_NelemOthersNeed(irank+1)%a)   
         endif
         if (NelemIneed(irank+1) /= 0) then
            deallocate(ReceiveBuffer(irank+1)%a)
            deallocate(GEI_NelemINeed(irank+1)%a)
         endif
      enddo
         
         
      deallocate(GEI_NelemINeed,GEI_NelemOthersNeed,SendBuffer,ReceiveBuffer)
      deallocate(ToReceiveElements,ToSendElements,NelemOthersNeed,NelemINeed,BufferSizeToSend,BufferSizeToReceive)
      deallocate(irequest01,irequest02,irequest03,irequest04)
      deallocate(HeaderListToSend)
      call StorageListToSend%Dealloc
   end subroutine
   
   
   subroutine ElementalExtra(a,ARET)
      use MPI
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerElementTransferer) :: ARET
      
      integer(ip) :: BufferSizes(a%MPIsize), BufferSizesReceive(a%MPIsize), Bsize,BufferCount
      type(r1p), allocatable :: SendBuffer(:),ReceiveBuffer(:)
      
      integer(ip), allocatable :: irequest01(:), irequest02(:)
      integer(ip) elemi,ielem,irank,ierr
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      integer(ip) :: ExtraElementCount
      
      
      !We also need to do the ExtraElements!!  :(
      !This goes as the rebalance case, we need to first know how much to communicate, then do it
      BufferSizes = 0
      BufferSizesReceive = 0
      
      !First loop to count sizes to send and receive
      do irank = 0,a%MPIsize-1
         do elemi = 1,a%ExtraNelemOthersNeed(irank+1)
            ielem = a%ExtraElementsListOthersNeed(irank+1)%l(elemi)
            ielem = a%ExportElementNumbering(ielem)
            call ARET%GetNewElementBufferSize(ielem,BSize)
            BufferSizes(irank+1) = BufferSizes(irank+1)+BSize
         enddo
      enddo
      
      !Send and receive buffer sizes
      !First I will send the p, and then I can allocate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank  =0,a%MPIsize-1
         if (a%ExtraNelemOthersNeed(irank+1)>0) then
            call MPI_ISEND(BufferSizes(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%ExtraNelemINeed(irank+1) > 0) then
            call MPI_IRECV(BufferSizesReceive(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Fill the buffers
      allocate(SendBuffer(a%MPIsize))
      do irank = 0,a%MPIsize-1
         BufferCount = 1
         allocate(SendBuffer(irank+1)%a(BufferSizes(irank+1)))
         do elemi = 1,a%ExtraNelemOthersNeed(irank+1)
            ielem = a%ExtraElementsListOthersNeed(irank+1)%l(elemi)
            ielem = a%ExportElementNumbering(ielem)
            call ARET%NewElementToBuffer(ielem,SendBuffer(irank+1)%a(BufferCount:BufferCount))
            
            call ARET%GetNewElementBufferSize(ielem,Bsize)
            BufferCount = BufferCount + Bsize
         enddo
      enddo
      
      !Send and receive the buffers 
      allocate(ReceiveBuffer(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (a%ExtraNelemOthersNeed(irank+1)>0) then
            call MPI_ISEND(SendBuffer(irank+1)%a,BufferSizes(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         
         if (a%ExtraNelemINeed(irank+1) > 0) then 
            allocate(ReceiveBuffer(irank+1)%a(BufferSizesReceive(irank+1)))
            call MPI_IRECV(ReceiveBuffer(irank+1)%a,BufferSizesReceive(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)   
         
         
      !Rebuild the structure in this side   
      ExtraElementCount = 0
      do irank = 0,a%MPIsize-1   
         BufferCount = 1
         do elemi = 1,a%ExtraNelemINeed(irank+1)
            ExtraElementCount = ExtraElementCount+1
            ielem = a%ExportNelem+ExtraElementCount
            call ARET%BufferToNewElement(ielem,ReceiveBuffer(irank+1)%a(BufferCount:BufferCount))
            call ARET%GetNewElementBufferSize(ielem,Bsize)
            BufferCount = BufferCount + Bsize
         enddo
      enddo
         
      !Deallocates   
      do irank = 0,a%MPIsize-1
         if (associated(SendBuffer(irank+1)%a)) deallocate(SendBuffer(irank+1)%a)
         if (associated(ReceiveBuffer(irank+1)%a)) deallocate(ReceiveBuffer(irank+1)%a)
      enddo
      deallocate(SendBuffer,ReceiveBuffer)
      deallocate(irequest01,irequest02)
   end subroutine
  
   subroutine ISO_InterpolateElement(a,ndofn,OldNodalArray,NewNodalArray     )
      class(IsotropicAdaptiveRefiner) :: a
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(:,:),NewNodalArray(:,:)
      
      if (a%ARET_InterpolationTask == 1) then
         call ParentToChildrenElement(ndofn,OldNodalArray, NewNodalArray,a%ARET_ielty,a%ARET_ivariation,a%ARET_ichild)
      elseif (a%ARET_InterpolationTask == 2) then
         call ChildrenToParentElement(ndofn,OldNodalArray, NewNodalArray,a%ARET_ielty,a%ARET_ivariation,a%ARET_ichild)
      endif
   end subroutine
   
   subroutine ElementalRebalance(a,ARET)
      use MPI
      implicit none
      class(IsotropicAdaptiveRefiner) :: a
      class(AdaptiveRefinerElementTransferer) :: ARET
      
      integer(ip) :: BufferSizes(a%MPIsize), BufferSizesReceive(a%MPIsize), Bsize,BufferCount
      type(r1p), allocatable :: SendBuffer(:),ReceiveBuffer(:)
      type(i1p), allocatable :: SendBufferIndexes(:), ReceiveBufferIndexes(:)
      
      integer(ip), allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer(ip) elemi,ielem,irank,ierr
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      BufferSizes = 0
      BufferSizesReceive = 0
      
      !First loop to count sizes to send and receive
      do irank = 0,a%MPIsize-1
         do elemi = 1,a%RebalanceNelemOthersNeed(irank+1)
            ielem = a%RebalanceListElementsOthersNeed(irank+1)%l(elemi)
            ielem = a%OldExportElementNumbering(ielem)
            if (ielem > 0) then
               call ARET%GetOldElementBufferSize(ielem,BSize)
               BufferSizes(irank+1) = BufferSizes(irank+1)+BSize
            endif
         enddo
      enddo
            
      !Send and receive buffer sizes
      !First I will send the p, and then I can allocate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank  =0,a%MPIsize-1
         if (a%RebalanceNelemOthersNeed(irank+1)>0) then
            call MPI_ISEND(BufferSizes(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%RebalanceNelemINeed(irank+1) > 0) then
            call MPI_IRECV(BufferSizesReceive(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Fill the buffers
      allocate(SendBuffer(a%MPIsize))
      allocate(SendBufferIndexes(a%MPIsize))
      do irank = 0,a%MPIsize-1
         BufferCount = 1
         allocate(SendBuffer(irank+1)%a(BufferSizes(irank+1)))
         allocate(SendBufferIndexes(irank+1)%l(a%RebalanceNelemOthersNeed(irank+1)))
         do elemi = 1,a%RebalanceNelemOthersNeed(irank+1)
            ielem = a%RebalanceListElementsOthersNeed(irank+1)%l(elemi)
            ielem = a%OldExportElementNumbering(ielem)
            if (ielem > 0) then
               call ARET%OldElementToBuffer(ielem,SendBuffer(irank+1)%a(BufferCount:BufferCount))
               
               !The sendbuffer indexes are required because an element might be exported in the sender, but not in the receiver
               SendBufferIndexes(irank+1)%l(elemi) = BufferCount
            
            
               call ARET%GetOldElementBufferSize(ielem,Bsize)
               BufferCount = BufferCount + Bsize
            endif
         enddo
      enddo

      !Send and receive the buffers 
      allocate(irequest03(a%MPIsize))
      allocate(irequest04(a%MPIsize))
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL
      allocate(ReceiveBuffer(a%MPIsize))
      allocate(ReceiveBufferIndexes(a%MPIsize))
      do irank = 0,a%MPIsize-1
         if (a%RebalanceNelemOthersNeed(irank+1)>0) then
            call MPI_ISEND(SendBuffer(irank+1)%a,BufferSizes(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(SendBufferIndexes(irank+1)%l,a%RebalanceNelemOthersNeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest03(irank+1), ierr)
         endif
         
         if (a%RebalanceNelemINeed(irank+1) > 0) then 
            allocate(ReceiveBuffer(irank+1)%a(BufferSizesReceive(irank+1)))
            allocate(ReceiveBufferIndexes(irank+1)%l(a%RebalanceNelemINeed(irank+1)))
            call MPI_IRECV(ReceiveBuffer(irank+1)%a,BufferSizesReceive(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
            call MPI_IRECV(ReceiveBufferIndexes(irank+1)%l,a%RebalanceNelemINeed(irank+1),MPI_INTEGER4, irank, mtag2, a%MPIcomm,irequest04(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)   
      
      deallocate(irequest01,irequest02,irequest03,irequest04)
         
         
      !Rebuild the structure in this side   
      do irank = 0,a%MPIsize-1   
         do elemi = 1,a%RebalanceNelemINeed(irank+1)
            ielem = a%RebalanceListElementsINeed(irank+1)%l(elemi)
            ielem = a%ExportElementNumbering(ielem)
            if (ielem > 0) then
               BufferCount = ReceiveBufferIndexes(irank+1)%l(elemi)
               call ARET%BufferToNewElement(ielem,ReceiveBuffer(irank+1)%a(BufferCount:BufferCount))
            endif
         enddo
      enddo
         
         
      !Deallocates   
      do irank = 0,a%MPIsize-1
         if (associated(SendBuffer(irank+1)%a)) deallocate(SendBuffer(irank+1)%a,SendBufferIndexes(irank+1)%l)
         if (associated(ReceiveBuffer(irank+1)%a)) deallocate(ReceiveBuffer(irank+1)%a,ReceiveBufferIndexes(irank+1)%l)
      enddo
      deallocate(SendBuffer,ReceiveBuffer,SendBufferIndexes,ReceiveBufferIndexes)

      !Communicate the extra elemental data!!
      call ElementalExtra(a,ARET)
      
   end subroutine
   


end module

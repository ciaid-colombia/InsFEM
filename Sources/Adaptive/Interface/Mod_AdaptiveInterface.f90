module Mod_AdaptiveInterface
   use typre
   use Mod_Memor
   use Mod_MPIObject
   implicit none
   private
   public AdaptiveRefinerInitialData, AdaptiveRefinerInterface, AdaptiveRefinerElementTransferer
   
   
   type, abstract :: AdaptiveRefinerInitialData
      integer(ip) :: ndime
      integer(ip) :: npoinlocal
      integer(ip) :: npoin
      integer(ip) :: gnpoin
      integer(ip) :: nelem
      
contains
      procedure(GetLnodsInitial), deferred :: GetLnode
      procedure(Local2Global), deferred :: Local2Global
      procedure(GetProcessor), deferred :: GetProcessor
   end type
   
   type, abstract :: AdaptiveRefinerElementTransferer
   
   contains

      procedure(CopyElement), deferred :: CopyElement
      procedure(NodallyInterpolateElement), deferred :: NodallyInterpolateElement
      procedure(NodallyInterpolateElementFromOldBuffer), deferred :: NodallyInterpolateElementFromOldBuffer
   
      procedure(GetOldElementBufferSize), deferred :: GetOldElementBufferSize
      procedure(OldElementToBuffer), deferred :: OldElementToBuffer
      
      procedure(GetNewElementBufferSize), deferred :: GetNewElementBufferSize
      procedure(BufferToNewElement), deferred :: BufferToNewElement
      procedure(NewElementToBuffer), deferred :: NewElementToBuffer
      
   
   end type
   
   abstract interface
   
      subroutine GetLnodsInitial(a,ielem,pnode,lnode)
         use typre
         import AdaptiveRefinerInitialData
         implicit none
         class(AdaptiveRefinerInitialData) :: a
         integer(ip) :: ielem,pnode
         integer(ip), pointer :: lnode(:)
      end subroutine
      
      subroutine Local2Global(a,npoin,LocalNumbering,GlobalNumbering)
         use typre
         import AdaptiveRefinerInitialData
         implicit none
         class(AdaptiveRefinerInitialData) :: a
         integer(ip) :: npoin
         integer(ip) :: LocalNumbering(npoin),GlobalNumbering(npoin)
      end subroutine
      
      subroutine GetProcessor(a,npoin,LocalNumbering,ProcessorList)
         use typre
         import AdaptiveRefinerInitialData
         implicit none
         class(AdaptiveRefinerInitialData) :: a
         integer(ip) :: npoin
         integer(ip) :: LocalNumbering(npoin),ProcessorList(npoin)
      end subroutine
   end interface

   
   
   type, abstract, extends(MpiObject) ::  AdaptiveRefinerInterface
      
   
   
contains   
      
      procedure(Initialize), deferred           :: Initialize
      procedure(SetCoordArray), deferred        :: SetCoordArray
      procedure(Refine), deferred               :: Refine
      procedure(GetHangingListDimensions), deferred :: GetHangingListDimensions
      procedure(GetHangingList), deferred       :: GetHangingList
      procedure(GetHangingFacesList), deferred       :: GetHangingFacesList
      procedure(GetPointDimensions), deferred   :: GetPointDimensions
      procedure(GetElementDimensions), deferred :: GetElementDimensions
      procedure(GetLocalOrdering), deferred     :: GetLocalOrdering
      procedure(GetLnods), deferred             :: GetLnods
      
      procedure(PrepareOldExternalBoundariesInfoToKeep), deferred ::  PrepareOldExternalBoundariesInfoToKeep
      procedure(MatchOldLboelToOldExternalBoundaries), deferred ::       MatchOldLboelToOldExternalBoundaries
      procedure(MatchNewLboelToNewElementsAndFaces), deferred ::                    MatchNewLboelToNewElementsAndFaces
      procedure(MatchOldBoundariesToNewBoundaries), deferred ::              MatchOldBoundariesToNewBoundaries
      
      procedure(UpdateVariableReal), deferred   :: UpdateVariableReal
      procedure(UpdateVariableInteger), deferred   :: UpdateVariableInteger
      procedure(UpdateVariableReal1), deferred   :: UpdateVariableReal1
      procedure(UpdateVariableInteger1), deferred   :: UpdateVariableInteger1
      procedure(RebalanceVariableReal), deferred   :: RebalanceVariableReal
      procedure(RebalanceVariableInteger), deferred   :: RebalanceVariableInteger
      procedure(RebalanceVariableReal1), deferred   :: RebalanceVariableReal1
      procedure(RebalanceVariableInteger1), deferred   :: RebalanceVariableInteger1
      procedure(UpdateVariableLogical), deferred   :: UpdateVariableLogical
      procedure(UpdateVariableLogical1), deferred   :: UpdateVariableLogical1
      procedure(RebalanceVariableLogical), deferred   :: RebalanceVariableLogical
      procedure(RebalanceVariableLogical1), deferred   :: RebalanceVariableLogical1
      procedure(GetiBoundaryMatch), deferred           :: GetiBoundaryMatch
      procedure(GetiBoundaryToElementsAndFacesMatch), deferred :: GetiBoundaryToElementsAndFacesMatch
      procedure(RebalanceVariableIntegerBoundary1), deferred :: RebalanceVariableIntegerBoundary1
      procedure(RebalanceVariableIntegerBoundary), deferred  :: RebalanceVariableIntegerBoundary
      procedure(RebalanceVariableRealBoundary1), deferred    :: RebalanceVariableRealBoundary1
      procedure(RebalanceVariableRealBoundary), deferred     :: RebalanceVariableRealBoundary
      procedure(RebalanceVariableTyper1pBoundary), deferred  :: RebalanceVariableTyper1pBoundary
      procedure(LoadRebalance), deferred :: LoadRebalance
      procedure(GetLevel), deferred :: GetLevel
      procedure(Set_2_1_Balancing), deferred    :: Set_2_1_Balancing
      procedure(SetRebalancingNumbering), deferred :: SetRebalancingNumbering
      procedure(WriteTimes), deferred :: WriteTimes

      
      generic :: UpdateVariable => UpdateVariableReal, UpdateVariableInteger,UpdateVariableReal1,UpdateVariableInteger1,UpdateVariableLogical,UpdateVariableLogical1
      generic :: RebalanceVariable => RebalanceVariableReal, RebalanceVariableInteger,RebalanceVariableReal1,RebalanceVariableInteger1,RebalanceVariableLogical,RebalanceVariableLogical1
      generic :: RebalanceVariableBoundary => RebalanceVariableRealBoundary,RebalanceVariableIntegerBoundary,RebalanceVariableRealBoundary1,RebalanceVariableIntegerBoundary1,RebalanceVariableTyper1pBoundary
      
      procedure(Dealloc), deferred :: Dealloc

      procedure(ElementalRefine), deferred :: ElementalRefine
      procedure(InterpolateElement), deferred :: InterpolateElement
      procedure(ElementalRebalance), deferred :: ElementalRebalance
   
   end type
  
  
    abstract interface
   
      subroutine Initialize(a,InitialData)
         use typre
         import AdaptiveRefinerInterface
         import AdaptiveRefinerInitialData
         implicit none
         class(AdaptiveRefinerInterface) :: a
         class(AdaptiveRefinerInitialData) :: InitialData
      end subroutine
      
      subroutine SetCoordArray(a,coord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         real(rp), target :: coord(:,:)
      end subroutine
      
      subroutine Refine(a,markelLastLevel)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: markelLastLevel(*)
      end subroutine
   
      subroutine GetElementDimensions(a,newnelem,LnodsSize)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) ::newnelem,lnodsSize
      end subroutine
        
      subroutine GetPointDimensions(a,newnpoinLocal,newnpoin,newgnpoin)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: newnelem,NewGnpoin,NewNpoin,NewNpoinLocal
      end subroutine
   
      subroutine GetLocalOrdering(a,LocalToGlobal,ProcessorList)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: NewGnpoin,NewNpoin,NewNpoinLocal
         integer(ip) :: LocalToGlobal(*),ProcessorList(*)
      end subroutine
      
      subroutine GetLnods(a,pnods,lnods)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: NewNelem
         integer(ip) :: pnods(*), lnods(*)
      end subroutine
      
      subroutine GetHangingListDimensions(a,HangingListDimensions)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: HangingListDimensions
      end subroutine
      
      subroutine GetHangingList(a,pHangingList,lHangingList,rHangingList)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: pHangingList(*),lHangingList(*)
         real(rp) :: rHangingList(*)
      end subroutine
      
      subroutine GetHangingFacesList(a,pHangingList,lHangingList)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: pHangingList(*),lHangingList(2,*)
      end subroutine
      
      subroutine PrepareOldExternalBoundariesInfoToKeep(a)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
      end subroutine
      
      subroutine LoadRebalance(a)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
      end subroutine
      
      subroutine MatchOldLboelToOldExternalBoundaries(a,nboun,pboel,lboel,BoundaryMatch,ExportOldNboun)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: nboun,pboel(nboun+1),lboel(pboel(nboun+1)-1),BoundaryMatch(nboun),ExportOldNboun
      end subroutine
      
      subroutine MatchNewLboelToNewElementsAndFaces(a,nboun,pboel,lboel,BoundaryMatch)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: nboun,pboel(nboun+1),lboel(pboel(nboun+1)-1),BoundaryMatch(2,nboun)
      end subroutine
      
      subroutine MatchOldBoundariesToNewBoundaries(a,oldnboun,oldBoundaryMatch,newnboun,newBoundaryMatch,BoundaryNewToOld)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: oldnboun, oldBoundaryMatch(oldnboun),newnboun,newBoundaryMatch(2,newnboun),BoundaryNewToOld(newnboun)
      end subroutine
      
      subroutine UpdateVariableReal(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         real(rp) :: coord(ndime,*)
         real(rp) :: newcoord(ndime,*)
      end subroutine
      
      subroutine UpdateVariableReal1(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         real(rp), target :: coord(*)
         real(rp), target :: newcoord(*)
      end subroutine
      
      subroutine UpdateVariableInteger(a,ndime,coord,newcoord,maxmin)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         integer(ip) :: coord(ndime,*)
         integer(ip) :: newcoord(ndime,*)
         character(5) :: maxmin
      end subroutine
      
      subroutine UpdateVariableInteger1(a,ndime,coord,newcoord,maxmin)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         integer(ip), target :: coord(*)
         integer(ip), target :: newcoord(*)
         character(5) :: maxmin
      end subroutine
      
      subroutine UpdateVariableLogical(a,ndime,coord,newcoord,maxmin)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         logical :: coord(ndime,*)
         logical :: newcoord(ndime,*)
         character(5) :: maxmin
      end subroutine
      
      subroutine UpdateVariableLogical1(a,ndime,coord,newcoord,maxmin)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         logical, target :: coord(*)
         logical, target :: newcoord(*)
         character(5) :: maxmin
      end subroutine
            
      subroutine RebalanceVariableReal(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         real(rp) :: coord(ndime,*)
         real(rp) :: newcoord(ndime,*)
      end subroutine
      
      subroutine RebalanceVariableReal1(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         real(rp), target :: coord(*)
         real(rp), target :: newcoord(*)
      end subroutine
      
      subroutine RebalanceVariableInteger(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         integer(ip) :: coord(ndime,*)
         integer(ip) :: newcoord(ndime,*)
      end subroutine
      
      subroutine RebalanceVariableInteger1(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         integer(ip), target :: coord(*)
         integer(ip), target :: newcoord(*)
      end subroutine

      subroutine RebalanceVariableLogical(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         logical :: coord(ndime,*)
         logical :: newcoord(ndime,*)
      end subroutine
      
      subroutine RebalanceVariableLogical1(a,ndime,coord,newcoord)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndime !array dimension
         logical, target :: coord(*)
         logical, target :: newcoord(*)
      end subroutine
  
      subroutine GetiBoundaryMatch(a,nboun,BoundaryMatch,iBoundaryMatch)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: nboun
         integer(ip) :: BoundaryMatch(nboun), iBoundaryMatch(*)
      end subroutine
      
      subroutine RebalanceVariableIntegerBoundary1(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: oldnboun,ndofn
         integer(ip) :: ioldBoundaryMatch(*)
         integer(ip), target :: OldArray(*),NewArray(*)
         type(i1p)   :: inewBoundaryMatch(*)
      end subroutine   
      
      subroutine RebalanceVariableIntegerBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: oldnboun,ndofn
         integer(ip) :: ioldBoundaryMatch(*),OldArray(ndofn,*),NewArray(ndofn,*)
         type(i1p)   :: inewBoundaryMatch(*)
      end subroutine
      
      subroutine RebalanceVariableRealBoundary1(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: oldnboun,ndofn
         integer(ip) :: ioldBoundaryMatch(*)
         real(rp), target :: OldArray(*),NewArray(*)
         type(i1p)   :: inewBoundaryMatch(*)
      end subroutine  
      
      subroutine RebalanceVariableRealBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,ndofn,OldArray,NewArray)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: oldnboun,ndofn
         integer(ip) :: ioldBoundaryMatch(*)
         real(rp) ::OldArray(ndofn,*),NewArray(ndofn,*)
         type(i1p)   :: inewBoundaryMatch(*)   
      end subroutine   
      
      subroutine GetiBoundaryToElementsAndFacesMatch(a,nboun,BoundaryMatch,iBoundaryMatch,Memor)
         use typre
         use Mod_Memor
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: nboun
         integer(ip) :: BoundaryMatch(2,nboun)
         type(i1p), allocatable   :: iBoundaryMatch(:)
         type(MemoryMan) :: Memor
      end subroutine
      
      subroutine RebalanceVariableTyper1pBoundary(a,oldnboun,ioldBoundaryMatch,inewBoundaryMatch,OldArray,NewArray,AllocatedBytes)
      use typre
      import AdaptiveRefinerInterface
      implicit none
      class(AdaptiveRefinerInterface) :: a
      integer(ip) :: oldnboun,ndofn
      integer(ip) :: ioldBoundaryMatch(*)
      type(r1p)   :: OldArray(*),NewArray(*)
      type(i1p)   :: inewBoundaryMatch(*)
      integer(ip) :: AllocatedBytes

   end subroutine
      
      subroutine GetLevel(a,level)
         use typre
         import AdaptiveRefinerInterface   
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: level(*)
      end subroutine
      
      subroutine Set_2_1_Balancing(a,kfl_21Balancing)
         use typre
         import AdaptiveRefinerInterface   
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: kfl_21Balancing
      end subroutine
      
      subroutine SetRebalancingNumbering(a,PointGlobNumber,PointProcNumber)
         use typre
         import AdaptiveRefinerInterface   
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: PointGlobNumber(*),PointProcNumber(*)
      end subroutine
      
      subroutine Dealloc(a)
         use typre
         import AdaptiveRefinerInterface   
         implicit none
         class(AdaptiveRefinerInterface) :: a
      end subroutine
      
      subroutine WriteTimes(a,lunout,cputotal)
         use typre
         import AdaptiveRefinerInterface
         implicit none
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: lunout
         real(rp)    :: cputotal
      end subroutine
      
      
      
!-------------------------------------------------------------------------
! AdaptiveRefinerElementTransferer
!-------------------------------------------------------------------------
      
      
      
      subroutine ElementalRefine(a,ARET)
         import 
         implicit none
         class(AdaptiveRefinerInterface) :: a
         class(AdaptiveRefinerElementTransferer) :: ARET
      
      
      
      end subroutine
      
      subroutine ElementalRebalance(a,ARET)
         import 
         implicit none
         class(AdaptiveRefinerInterface) :: a
         class(AdaptiveRefinerElementTransferer) :: ARET
      
      
      
      end subroutine
      
      
      !AdaptiveRefinerElementTransferer
      subroutine CopyElement(a,oldielem,newielem)
         import
         implicit none
         class(AdaptiveRefinerElementTransferer), target :: a
         integer(ip) :: oldielem,newielem
      
      end subroutine
      
      subroutine NodallyInterpolateElement(a,Refiner,oldelem,newelem)
         import
         implicit none
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: oldelem,newelem
         class(AdaptiveRefinerInterface) :: Refiner
         
         !If there is something to be nodally interpolated, then we should call
         !the refiner to do it
         
         !call Refiner%NodallyInterpolateElement(ndofn,OldNodalArray,NewNodalArray)
      end subroutine
      
      subroutine NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
         import
         implicit none
         class(AdaptiveRefinerElementTransferer) :: a
         class(AdaptiveRefinerInterface) :: Refiner
         integer(ip) :: newelem
         real(rp), target :: OldBuffer(*)
      
      end subroutine
      
      
      subroutine InterpolateElement(a,ndofn,OldNodalArray,NewNodalArray     )
         import
         class(AdaptiveRefinerInterface) :: a
         integer(ip) :: ndofn
         real(rp) :: OldNOdalArray(:,:), NewNodalArray(:,:)
      end subroutine
      
      subroutine OldElementToBuffer(a,ielem,Buffer)
         import
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: ielem
         real(rp) :: Buffer(*)
      end subroutine
      
      subroutine NewElementToBuffer(a,ielem,Buffer)
         import
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: ielem
         real(rp) :: Buffer(*)
      end subroutine
      
      subroutine BufferToNewElement(a,ielem,Buffer)
         import
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: ielem
         real(rp) :: Buffer(*)
      end subroutine
      
      subroutine GetOldElementBufferSize(a,ielem,BufferSize)
         import
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: ielem
         integer(ip) :: BufferSize
      end subroutine
      
      subroutine GetNewElementBufferSize(a,ielem,BufferSize)
         import
         class(AdaptiveRefinerElementTransferer) :: a
         integer(ip) :: ielem
         integer(ip) :: BufferSize
      end subroutine
      
      
      
   end interface

  
  
contains
   
  
  

end module

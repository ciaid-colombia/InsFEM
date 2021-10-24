module Mod_Octree
   use typre
   use Mod_Memor
   use Mod_LinkedList
   implicit none
   private
   public Octree, OctreeLinkedListHeader, OctreeLinkedListStorage, RangeGiver
   
   real(rp) :: OctreeTol = 1e-6
   integer(ip), parameter :: maxlevel = 12  !maximum number of levels
   
   !OCTREE LINKED LIST HEAD TYPE
   !This one is used only when we have a range and want all the elements in the lists
   type, extends(LinkedListHeader) :: OctreeLinkedListHeader
      real(rp) :: range(2,3)
      
      
contains

   end type
   
   !OCTREE LINKED LIST STORE TYPE
   type, extends(LinkedListStorage) :: OctreeLinkedListStorage
      real(rp) :: range(2,3)
      class(Octree), pointer :: CurrentLeaf => NULL()
      class(Octree), pointer :: TopLeaf => NULL()
      
      
      
contains
      procedure :: GoToFirstOfList => OctreeLinkedListGoToFirstOfList
      
   end type
   
   !OCTREE TYPE  
   type :: Octree
   
      !Tree structure
      integer(ip)           :: nchildren = 0    !Number of children I have
      type(Octree), pointer :: OctChildren(:) => NULL()
      class(Octree), pointer :: OctParent => NULL()
      integer(ip) :: level = 0            !In which level of the tree I am
      integer(ip) :: ichildLooper = 0     !For looping through the children (OctreeLinkedListHeader)
      
      integer(ip) :: ndime
      integer(ip) :: npoin, npoinAsParent = -1
      integer(ip) :: maxpoin
      
      !The range contains the coordinates of the box ([low,mid,up],1:ndime)
      real(rp) :: range(3,3)
      
      type(LinkedListHeader)  :: ListHead
      type(LinkedListStorage), pointer :: ListStore => NULL()
      
      !These are lists only for getting ranges
      type(OctreeLinkedListHeader) :: OctreeHead
      type(OctreeLinkedListStorage) :: OctreeStore
      
      
contains

      !Initialize procedure
      procedure :: InitializePointOctree
      procedure :: InitializeElementOctree
      procedure :: GetListForCoord
      procedure :: GetListForRange
      procedure :: Dealloc
   
   end type
   
   
   !GIVE ME RANGE TYPE
   type, abstract :: RangeGiver
      !This is an abstract type whose sole function is to give the range of elements
      !to the Octree structure
      !It is implemented in this way to respecte the hierarchy
      !Mesh is above Octree => Octree cannot know about Mesh
      !It needs to be implemented by Mesh
   
contains

      procedure(GiveMeRange), deferred :: GiveMeRange
   end type
   
   abstract interface
      subroutine GiveMeRange(a,ielem,range)
         use typre
         import RangeGiver
         implicit none
         class(RangeGiver) :: a
         integer(ip) :: ielem
         real(rp) :: range(*)
      end subroutine
   end interface
   
  
   
contains


   !Children Allocation
   subroutine AllocateChildren(a,nchildren,Memor)
      use typre
      implicit none
      class(Octree), target :: a
      integer(ip) :: nchildren
      type(MemoryMan) :: Memor
      
      integer(ip) :: istat, ichild
      

      a%nchildren = nchildren
      allocate(a%OctChildren(a%nchildren),stat=istat)
      call Memor%allocObj(istat,'children','tree',a%nchildren)

      a%OctChildren(:)%level = a%level+1
      do ichild = 1,a%nchildren
         a%OctChildren(ichild)%OctParent => a
      enddo
      
   end subroutine 
   
   !Deallocation
   recursive subroutine DeallocateDownwards(a,Memor)
      !This is a from bottom to top procedure
      use typre
      implicit none
      class(Octree) :: a
      type(MemoryMan) :: Memor
      
      integer(ip) :: ichild
      integer(ip) :: istat
      
      !First my children
      do ichild = 1,a%nchildren
         call DeallocateDownwards(a%Octchildren(ichild),Memor)
      enddo
      
      !Delete my children (if they exist)
      if (a%nchildren /= 0) then
         deallocate(a%OctChildren,stat=istat)
         call Memor%DeallocObj(istat,'children','tree',a%nchildren)
         a%nchildren = 0
      endif
   end subroutine

   !For points
   subroutine InitializePointOctree(a,ndime,npoin,coord,maxpoin,Memor)
      use typre
      implicit none
      class(Octree) :: a
      integer(ip) :: ndime,npoin,maxpoin
      real(rp) :: coord(ndime,npoin)
      type(MemoryMan), target :: Memor
      
      integer(ip) :: idime,ipoin
      integer(ip) :: istat
      
      !Initializations
      a%ndime = ndime
      a%npoin = npoin
      a%maxpoin = maxpoin
      a%OctParent => NULL()
      !Initial Box for top level leaf
      do idime = 1,a%ndime
         a%range(1,idime) = minval(coord(idime,:)) - OctreeTol             !low end
         a%range(3,idime) = maxval(coord(idime,:)) + OctreeTol             !up  end
         a%range(2,idime) = (a%range(1,idime)+a%range(3,idime))/2.0_rp  !midpoint
      enddo
      
      !Allocate LinkedListStore
      allocate(a%ListStore,stat=istat)
      call Memor%AllocObj(istat,'ListStore','InitializePointOctree',1_ip)
      call a%ListStore%Initialize(a%npoin,'Removable','NonRepeatable',Memor)
      
      !We Start by adding all the points to the list of the initial Leaf
      do ipoin = 1,a%npoin
         call a%ListStore%Add(a%ListHead,ipoin)
      enddo 
      
      !Now we build the octree from the coordinate list
      call BuildOctreeFromCoord(a,coord,Memor)
   
   end subroutine
      
   recursive subroutine BuildOctreeFromCoord(a,coord,Memor)
      use typre
      implicit none
      class(Octree) :: a
      real(rp) :: coord(a%ndime,*)
      type(MemoryMan) :: Memor
      
      integer(ip) :: ichild
      
      !This is a from top to bottom subroutine: first me, then my children
      
      !If I am 
      if (a%npoin > a%maxpoin .and. a%level < maxlevel) then
         
         !The number of children is 4 for 2D and 8 for 3D
         a%nchildren = 2**a%ndime
         call AllocateChildren(a,a%nchildren,Memor)
         
         !Computes the Range for the newly created children
         !Assigns the pointer to the linked list
         call OctreeInitializeChildren(a)
         
         call MovePointListToChildren(a,coord)
      endif 
      
      !Now my children recursively do the same
      do ichild = 1,a%nchildren
         call BuildOctreeFromCoord(a%OctChildren(ichild),coord,Memor)
      enddo   
   end subroutine
   
   subroutine OctreeInitializeChildren(a)
      !This subroutine computes the range for all the children
      !Assigns the pointer of the linked list Store 
      use typre
      implicit none
      class(Octree) :: a
      
      integer(ip) :: ichild,auxval,base,idime,ival
      integer(ip) :: coeff(3)

      
      do ichild = 1,a%nchildren
         !First we compute the range for each children
         ival = ichild
         do idime = a%ndime,1,-1
            base = 2**(idime-1)
            auxval = (ival-1)/(base)
            coeff(idime) = auxval
            ival = mod(ival-1,base)+1  !These are integers!!
         enddo
      
         do idime = 1,a%ndime
            a%OctChildren(ichild)%range(1,idime) = a%range(coeff(idime)+1,idime)
            a%OctChildren(ichild)%range(3,idime) = a%range(coeff(idime)+2,idime)
            a%OctChildren(ichild)%range(2,idime) = & 
            (a%OctChildren(ichild)%range(1,idime)+a%OctChildren(ichild)%range(3,idime))/2.0_rp
         enddo
         
         !Secondly, we assign the pointer to the parent linked list Store unit and initialize the head
         a%OctChildren(ichild)%ListStore => a%ListStore
         call a%OctChildren(ichild)%ListHead%Initialize
         
         !Other variables
         a%OctChildren(ichild)%maxpoin = a%maxpoin
         a%OctChildren(ichild)%ndime = a%ndime
         a%OctChildren(ichild)%npoin = 0
      enddo
   end subroutine 
   
   subroutine MovePointListToChildren(a,coord)
      use typre
      implicit none
      class(Octree) :: a
      real(rp) :: coord(a%ndime,*)
      
      integer(ip) :: ichild,ierr,ipoin,nextpoin
      
      !We cycle through my list of points, remove it from parent, add it to corresponding child 
      call a%ListStore%GoToFirstOfList(a%ListHead)
      
      call a%ListStore%GetNext(ipoin)
      do while (ipoin /= -1)
         call a%ListStore%GetNext(nextpoin)
         
         !Which child does the point go to
         call ChooseChildFromCoord(a,coord(:,ipoin),ichild)
         
         !Remove from my list
         call a%ListStore%Remove(a%ListHead,ipoin,ierr)
         
         !Add to corresponding child list
         call a%ListStore%Add(a%OctChildren(ichild)%ListHead,ipoin)
         
         !Continue
         ipoin = nextpoin
         
      enddo
      
      !Update the number of points in each children and myself
      do ichild = 1,a%nchildren
         call a%OctChildren(ichild)%ListHead%GetNelem(a%OctChildren(ichild)%npoin)
      enddo
      call a%ListHead%GetNelem(a%npoin)
      
   end subroutine
   
   subroutine ChooseChildFromCoord(a,coord,ichild)
      use typre
      implicit none
      class(Octree) :: a
      real(rp) :: coord(a%ndime)
      integer(ip) :: ichild
      
      integer(ip) :: idime
      
      ichild = 1
      do idime = 1,a%ndime
         if (coord(idime) > a%range(2,idime)) then
            ichild = ichild + 2**(idime-1)
         endif
      enddo
   end subroutine      
   
   subroutine GetListForCoord(a,coord,Head,Store)
      use typre
      implicit none
      class(Octree), target :: a
      real(rp) :: coord(a%ndime)
      class(LinkedListHeader), pointer  :: Head
      class(LinkedListStorage), pointer :: Store

      class(Octree), pointer :: CurrentLeaf => NULL()
      integer(ip) :: ichild
      
      !This is where the journey begins
      CurrentLeaf => a
      
      !Travel to the deepest levels in search of the magical list...
      do while (CurrentLeaf%nchildren > 0)
         call ChooseChildFromCoord(CurrentLeaf,coord,ichild)
         CurrentLeaf => CurrentLeaf%OctChildren(ichild)
      enddo
      
      !The quest was successful
      Head => CurrentLeaf%ListHead
      Store => CurrentLeaf%ListStore
   end subroutine   
   
   subroutine Dealloc(a,Memor)   
      use typre
      implicit none
      class(Octree) :: a
      type(MemoryMan) :: Memor
      
      integer(ip) :: istat
      
      !First we deallocate everything (abstract subroutine)
      call DeallocateDownwards(a,Memor)
      
      !In the level 0 I also deallocate the linked list storage
      call a%ListStore%Dealloc
      deallocate(a%ListStore,stat=istat)
      call Memor%DeAllocObj(istat,'ListStore','InitializePointOctree',1_ip)
   end subroutine
   
   
   !For elements
   subroutine InitializeElementOctree(a,ndime,npoin,MacGyver,maxpoin,Memor)
      use typre
      implicit none
      class(Octree) :: a
      integer(ip) :: ndime,npoin,maxpoin
      class(RangeGiver) :: MacGyver
      type(MemoryMan), target :: Memor
      
      real(rp) :: range(2,ndime)
      integer(ip) :: idime,ipoin,istat
      
      !Initializations
      a%ndime = ndime
      a%npoin = npoin
      a%maxpoin = maxpoin
      a%OctParent => NULL()
      
      !Initial Box for top level leaf
      call MacGyver%GiveMeRange(1,range)
      a%range(1,1:a%ndime) = range(1,1:a%ndime)
      a%range(3,1:a%ndime) = range(2,1:a%ndime)
      
      do ipoin = 2,npoin
         call MacGyver%GiveMeRange(ipoin,range)
         do idime = 1,a%ndime
            a%range(1,idime) = min(a%range(1,idime),range(1,idime))
            a%range(3,idime) = max(a%range(3,idime),range(2,idime))
         enddo
      enddo
      
      !Now apply tolerances and middpoint
      do idime = 1,a%ndime
         a%range(1,idime) = a%range(1,idime) - OctreeTol             !low end
         a%range(3,idime) = a%range(3,idime) + OctreeTol             !up  end
         a%range(2,idime) = (a%range(1,idime)+a%range(3,idime))/2.0_rp  !midpoint
      enddo
      
      !Allocate LinkedListStore
      allocate(a%ListStore,stat=istat)
      call Memor%AllocObj(istat,'ListStore','InitializeElementOctree',1_ip)
      call a%ListStore%Initialize(a%npoin,'Removable','Repeatable',Memor)
      
      !We Start by adding all the points to the list of the initial Leaf
      do ipoin = 1,a%npoin
         call a%ListStore%Add(a%ListHead,ipoin)
      enddo 
      
      !Now we build the octree from the coordinate list
      call BuildOctreeFromRangeGiver(a,MacGyver,Memor)
   
   end subroutine
   
   subroutine oct_checkOverlappingOfElements(a,MacGyver,ForceGoUp)
      implicit none
      class(Octree) :: a
      class(RangeGiver) :: MacGyver
      logical :: ForceGoUp
      
      real(rp) ::  range(2,a%ndime), auxrange(2,a%ndime)
      logical  :: DimensionalOverlap
      
      integer(ip) :: idime, ielem,auxnpoin
      
      real(rp) :: BoxSizes(a%ndime), auxBoxSizes(a%ndime)
      
      ForceGoUp = .true.

!       if (a%level > 2) then
!          if (a%npoin == a%OctParent%npoinAsParent) then
!             !If my number of elements is equal to my parent's number of elements we are not progressing
!             !Force GoUp
!             if (a%npoin > 1000) then
!                write(*,*) 'asdfasdf'
!                auxnpoin = a%npoin
!                write(*,*) auxnpoin
!             endif
!             return
!          endif
!       endif

      
      do idime = 1,a%ndime
         BoxSizes(idime) = a%range(3,idime)-a%range(1,idime)
      enddo
      
      call a%ListStore%GoToFirstOfList(a%ListHead)
      call a%ListStore%GetNext(ielem)
      DoList: do while (ielem /= -1)
         call MacGyver%GiveMeRange(ielem,auxrange)
         
         do idime = 1,a%ndime
            auxBoxSizes(idime) = auxrange(2,idime)-auxrange(1,idime)
            
            if (auxBoxSizes(idime)*5.0_rp < BoxSizes(idime)) then
               ForceGoUp = .false. 
               return
            endif
         enddo
         
         !Continue
         call a%ListStore%GetNext(ielem)
      enddo DoList
      
   !    DimensionalOverlap = .true.
   !    
   !    
   !    range(1,:) = -1e24
   !    range(2,:) = 1e24
   !    !We cycle through my list of elements, check their range and see if all of them overlap
   !    call a%ListStore%GoToFirstOfList(a%ListHead)
   !    call a%ListStore%GetNext(ielem)
   !    DoList: do while (ielem /= -1)
   !       call MacGyver%GiveMeRange(ielem,auxrange)
   !       
   !       do idime = 1,a%ndime
   !          range(1,idime) = max(range(1,idime),auxrange(1,idime))
   !          range(2,idime) = min(range(2,idime),auxrange(2,idime))
   !          
   !          !The sign seems to be important: must be a plus so that we are more restrictive
   !          !This is key for multiple elements sharing one point
   !          if (range(1,idime) > range(2,idime) + OctreeTol) then
   !             DimensionalOverlap = .false.
   !             exit DoList
   !          endif
   !       enddo
   !       
   !       !Continue
   !       call a%ListStore%GetNext(ielem)
   !    enddo DoList
   !    
   !    if (DimensionalOverlap .eqv. .true.) then
   !       write(*,*) 'aha'
   !    endif
   !    
   !    
   !    ForceGoUp = DimensionalOverlap
      
   end subroutine
   
   recursive subroutine BuildOctreeFromRangeGiver(a,MacGyver,Memor)
      use typre
      implicit none
      class(Octree) :: a
      class(RangeGiver) :: MacGyver
      type(MemoryMan) :: Memor
      
      integer(ip) :: ichild
      
      logical :: ForceGoUp
      
       
      
      
      !This is a from top to bottom subroutine: first me, then my children
      
      !If I am 
      if (a%npoin > a%maxpoin .and. a%level < maxlevel) then
      
         a%npoinAsParent = a%npoin
         
         !Check if all the boxes for all the elements overlap, if so, move up
         call oct_CheckOverlappingOfElements(a,MacGyver,ForceGoUp)
         if (ForceGoUp) return
         
         !The number of children is 4 for 2D and 8 for 3D
         a%nchildren = 2**a%ndime
         call AllocateChildren(a,a%nchildren,Memor)
         
         !Computes the Range for the newly created children
         !Assigns the pointer to the linked list
         call OctreeInitializeChildren(a)
         
         call MoveElementListToChildren(a,MacGyver)
         
         call CheckAndCorrectForNonUsefulDimensions(a)
      endif 
      
      !Now my children recursively do the same
      do ichild = 1,a%nchildren
         call BuildOctreeFromRangeGiver(a%OctChildren(ichild),MacGyver,Memor)
      enddo   
   end subroutine
   
   subroutine MoveElementListToChildren(a,MacGyver)
      use typre
      implicit none
      class(Octree) :: a
      class(RangeGiver) :: MacGyver
      
      real(rp) :: range(2,a%ndime)
      
      integer(ip) :: ichild,ierr,ipoin,nextpoin
      logical     :: IsIn
      
      !We cycle through my list of points, remove it from parent, add it to corresponding child 
      call a%ListStore%GoToFirstOfList(a%ListHead)
      
      call a%ListStore%GetNext(ipoin)
      do while (ipoin /= -1)
         call a%ListStore%GetNext(nextpoin)
         
         !Which children does the point go to
         do ichild = 1,a%nchildren
            call MacGyver%GiveMeRange(ipoin,range)
            call CheckFromRange(a%OctChildren(ichild),range,IsIn)
            !Add to corresponding child list if IsIN
            if (IsIn .eqv. .true.) call a%ListStore%Add(a%OctChildren(ichild)%ListHead,ipoin)
         enddo
         
         !Continue
         ipoin = nextpoin
      enddo
      
      !Remove All from my list, since they have been added to the children list
      call a%ListStore%RemoveAll(a%ListHead)
      
      
      !Update the number of points in each children and myself
      do ichild = 1,a%nchildren
         call a%OctChildren(ichild)%ListHead%GetNelem(a%OctChildren(ichild)%npoin)
      enddo
      call a%ListHead%GetNelem(a%npoin)
   end subroutine
   
   subroutine CheckFromRange(a,range,IsIn)
      use typre 
      implicit none
      class(Octree) :: a
      real(rp) :: range(2,a%ndime)
      logical :: IsIN
      
      integer(ip) :: idime
      
      IsIn = .true.
      !All conditions need to be satisfied
      do idime = 1,a%ndime
         if (range(2,idime) < a%range(1,idime) - OctreeTol) IsIn = .false.
         if (range(1,idime) > a%range(3,idime) + OctreeTol) IsIn = .false.
      enddo
   end subroutine
   
   
  
   
      
   
   !For Getting the list for a range
   subroutine GetListForRange(a,range,Head,Store)
      use typre
      implicit none
      class(Octree), target :: a
      real(rp) :: range(2,a%ndime)
      class(LinkedListHeader), pointer  :: Head
      class(LinkedListStorage), pointer :: Store
      
      Head => a%OctreeHead
      a%OctreeHead%range(1:2,1:a%ndime) = range(1:2,1:a%ndime)
      
      Store => a%OctreeStore
      !TopLeaf points to the entry leaf (called by the user)
      a%OctreeStore%TopLeaf => a
      a%OctreeStore%GetNext => OctreeLinkedListGetNext
   end subroutine   
   
   subroutine OctreeLinkedListGoToFirstOfList(a,Header)
      use typre
      implicit none
      class(OctreeLinkedListStorage) :: a
      class(LinkedListHeader)  :: Header
      
      logical :: found, IsIn
      
      !Copy range
      select type (Header)
      type is (LinkedListHeader)
      class is (OctreeLinkedListHeader)
         a%range = Header%range
      end select
      
      a%CurrentLeaf => a%TopLeaf
      a%CurrentLeaf%ichildLooper = 0
      call CheckFromRange(a%CurrentLeaf,a%range,IsIn)
      
      if (IsIN .eqv. .true.) then
         found = .true.
      
         call OctreeLinkedListMoveToNextList(a)
         
      else
         found = .false.
      endif
      
      
   
   end subroutine
   
   subroutine OctreeLinkedListGetNext(a,value)
      !Be careful: if this is an element list, elements can be repeated!!
      use typre
      implicit none
      class(LinkedListStorage) :: a
      integer(ip) :: value
      
      logical :: IsIn
      
      select type (a)
      type is (LinkedListStorage)
      class is (OctreeLinkedListStorage)
      
      call a%CurrentLeaf%ListStore%GetNext(value)
      !Continue with the previous process
      do while (value == -1 .and. (associated(a%CurrentLeaf%OctParent) .eqv. .true.)) 
         a%CurrentLeaf => a%CurrentLeaf%OctParent
         
         call OctreeLinkedListMoveToNextList(a)
      
         call a%CurrentLeaf%ListStore%GetNext(value)
   
      enddo

      end select
   end subroutine
   
   subroutine OctreeLinkedListMoveToNextList(a)
      class(OctreeLinkedListStorage) :: a
      
      logical :: IsIn
      
      !Travel to the deepest levels in search of the magical list...
      LeafLoop: do while (a%CurrentLeaf%nchildren > 0)
         childLoop: do while (a%CurrentLeaf%ichildLooper < a%CurrentLeaf%nchildren)
            a%CurrentLeaf%ichildLooper = a%CurrentLeaf%ichildLooper + 1 
            call CheckFromRange(a%CurrentLeaf%OctChildren(a%CurrentLeaf%ichildLooper),a%range,IsIn)
            if (IsIn .eqv. .true.) then
               a%CurrentLeaf =>a%CurrentLeaf%OctChildren(a%CurrentLeaf%ichildLooper)
               a%CurrentLeaf%ichildLooper = 0
               exit childLoop
            endif
         enddo childLoop
         if (a%CurrentLeaf%ichildLooper >= a%CurrentLeaf%nchildren) exit LeafLoop
      enddo LeafLoop
      call a%CurrentLeaf%ListStore%GoToFirstOfList(a%CurrentLeaf%ListHead)
   end subroutine
   
   
   
   
   subroutine CheckAndCorrectForNonUsefulDimensions(a)
      implicit none
      class(Octree) :: a
      
      integer(ip) :: ncheck, dimecheck
      
      
      integer(ip) :: checkarray(2,4)
      logical :: checklogicals(4)
      
      type(Octree), pointer :: OctChildren(:) => null()
      
      OctChildren => a%OctChildren
      
      
      !Check if children in both sides have exactly the same number of elements
      if (a%ndime == 2) then
         !check x
         !1 - 2
         !3 - 4
         ncheck = 2
         dimecheck = 1
         checkarray(1,1) = 1
         checkarray(2,1) = 2
         checkarray(1,2) = 3
         checkarray(2,2) = 4
         call PerformChecks
        
        
         !check y
         !1 - 3
         !2 - 4
         ncheck = 2
         dimecheck = 2
         checkarray(1,1) = 1
         checkarray(2,1) = 3
         checkarray(1,2) = 2
         checkarray(2,2) = 4
         call PerformChecks
      
      
      elseif (a%ndime == 3) then
         !check x
         !1 - 2
         !3 - 4
         !5 - 6
         !7 - 8
         
         ncheck = 4
         dimecheck = 1
         checkarray(1,1) = 1
         checkarray(2,1) = 2
         checkarray(1,2) = 3
         checkarray(2,2) = 4
         checkarray(1,3) = 5
         checkarray(2,3) = 6
         checkarray(1,4) = 7
         checkarray(2,4) = 8
         call PerformChecks
         
         !check y
         !1 - 3
         !2 - 4
         !5 - 7
         !6 - 8
         ncheck = 4
         dimecheck = 2
         checkarray(1,1) = 1
         checkarray(2,1) = 3
         checkarray(1,2) = 2
         checkarray(2,2) = 4
         checkarray(1,3) = 5
         checkarray(2,3) = 7
         checkarray(1,4) = 6
         checkarray(2,4) = 8
         call PerformChecks
         
         !check z
         !1 - 5
         !2 - 6
         !3 - 7
         !4 - 8
         ncheck = 4
         dimecheck = 3
         checkarray(1,1) = 1
         checkarray(2,1) = 5
         checkarray(1,2) = 2
         checkarray(2,2) = 6
         checkarray(1,3) = 3
         checkarray(2,3) = 7
         checkarray(1,4) = 4
         checkarray(2,4) = 8
         call PerformChecks
      endif
      
         
contains
      subroutine PerformChecks
         integer(ip) :: icheck
         
         class(LinkedListStorage), pointer :: auxStorage => NULL()
         
         checklogicals = .false.
         do icheck = 1,ncheck
            call AreEqualsChildren(a,checkarray(:,icheck),checklogicals(icheck))
         enddo   
         
         if (all(checklogicals(1:ncheck) .eqv. .true.)) then
             do icheck = 1,ncheck
                !Remove All from children 3 and 4 list
               auxStorage => a%OctChildren(checkarray(2,icheck))%ListStore
               call auxStorage%RemoveAll(a%ListHead)
               
               !1 and 2 have the same x range as the parent
               a%OctChildren(checkarray(1,icheck))%range(:,dimecheck) = a%range(:,dimecheck)
             enddo
            
            !change parent range so that we only look for elements in childs 1 and 2
            a%range(2,dimecheck) = a%range(3,dimecheck)
         endif
      end subroutine
   end subroutine
   
   subroutine AreEqualsChildren(a,checkarray,checklogical)
      implicit none
      class(Octree) :: a
      integer(ip) :: checkarray(2)
      logical :: checklogical
      
      integer(ip), allocatable :: array1(:), array2(:)
      
      checklogical = .false.
      if (a%OctChildren(checkarray(1))%npoin /= a%OctChildren(checkarray(2))%npoin) return
      !if both of them are zero also return
      if (a%OctChildren(checkarray(1))%npoin == 0) return
      
      allocate(array1(a%OctChildren(checkarray(1))%npoin))
      allocate(array2(a%OctChildren(checkarray(2))%npoin))
      
      !Take list of points for each children from the linked list storage
      call a%OctChildren(checkarray(1))%ListStore%ListToArray(a%OctChildren(checkarray(1))%ListHead,array1)
      call a%OctChildren(checkarray(2))%ListStore%ListToArray(a%OctChildren(checkarray(2))%ListHead,array2)
      
      if (all(array1.eq.array2)) checklogical = .true.
      
      deallocate(array1,array2)
   end subroutine
   
end module

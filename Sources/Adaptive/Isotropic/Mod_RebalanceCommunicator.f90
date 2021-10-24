module Mod_RebalanceCommunicator
   use typre
   use Mod_MPIObject
   implicit none
   private
   public RebalanceCommunicator

   type, extends(MPIObject) :: RebalanceCommunicator
      integer(ip) :: OldNpoinLocal, NewNpoinLocal
      integer(ip), allocatable :: NpoinOthersNeed(:),NpoinINeed(:)
      type(i1p), allocatable :: PointsINeed(:),PointsOthersNeed(:)

contains
      procedure :: IsAllocated
      procedure :: Initialize
      procedure :: CommunicateReal
      procedure :: CommunicateReal1
      procedure :: CommunicateInteger
      procedure :: CommunicateInteger1
      procedure :: CommunicateLogical
      procedure :: CommunicateLogical1
      procedure :: Deallocate => Dealloc
      
      generic :: Communicate => CommunicateReal, CommunicateReal1, CommunicateInteger, CommunicateInteger1,CommunicateLogical,CommunicateLogical1
   end type

contains
   subroutine IsAllocated(a,kfl_isAllocated)
      use typre
      class(RebalanceCommunicator) :: a
      logical :: kfl_isAllocated
      
      kfl_isAllocated = .false.
      if (allocated(a%NpoinOthersNeed)) then
         kfl_isAllocated = .true.
      endif
   end subroutine
   
   subroutine Initialize(a,OldNpoinLocal,NewNpoinLocal,RebalancePointNumbering,RebalanceProcessorList)
      use typre
      use Mod_LinkedList
      use MPI
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: OldNPoinLocal, NewNpoinLocal,RebalancePointNumbering(*),RebalanceProcessorList(*)
      
      type(LinkedListHeader), allocatable :: HeaderPointsOthersNeed(:)
      type(LinkedListStorage) :: StoragePointsOthersNeed
      
      integer(ip) :: ipoin,irank,poin0,auxpoin0
      logical :: mintaken
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      a%OldNpoinLocal = OldNpoinLocal
      a%NewNpoinLocal = NewNpoinLocal
      
      allocate(a%NpoinOthersNeed(a%MPIsize))
      allocate(a%NpoinINeed(a%MPIsize))
      a%NpoinOthersNeed = 0_ip
      allocate(HeaderPointsOthersNeed(a%MPIsize))
      call StoragePointsOthersNeed%Initialize(a%OldNpoinLocal,'NonRemove','NonRepeat',a%Memor)
      
      do ipoin = 1,a%OldNpoinLocal
         irank = RebalanceProcessorList(ipoin)
         call StoragePointsOthersNeed%Add(HeaderPointsOthersNeed(irank+1),ipoin)
         a%NpoinOthersNeed(irank+1) = a%NpoinOthersNeed(irank+1)+1
      enddo
      
      allocate(a%PointsOthersNeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(a%PointsOthersNeed(irank+1)%l(a%NpoinOthersNeed(irank+1)))
         call StoragePointsOthersNeed%ListToArray(HeaderPointsOthersNeed(irank+1),a%PointsOthersNeed(irank+1)%l)
         a%PointsOthersNeed(irank+1)%l = RebalancePointNumbering(a%PointsOthersNeed(irank+1)%l)
      enddo
      
      call MPI_AllToAll(a%NpoinOthersNeed, 1, MPI_INTEGER4, a%NpoinINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      allocate(a%PointsINeed(a%MPIsize))
      do irank = 0,a%MPIsize-1
         allocate(a%PointsINeed(irank+1)%l(a%NpoinINeed(irank+1)))
      enddo
      
      !Receive the points I Need
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(a%PointsOthersNeed(irank+1)%l,a%NpoinOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(a%PointsINeed(irank+1)%l,a%NpoinINeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !PointsOthersNeed to OldLocalNumbering
      do irank = 0,a%MPIsize-1
         !we reload it in local numbering
         call StoragePointsOthersNeed%ListToArray(HeaderPointsOthersNeed(irank+1),a%PointsOthersNeed(irank+1)%l)
      enddo
      
      call StoragePointsOthersNeed%Dealloc
      
      !PointsINeed to NewLocalNumbering
      poin0 = 0 !just in case nothing is found
      mintaken = .false.
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1)/= 0) then
            !Take the first one
            if (mintaken .eqv. .false.) then
               mintaken = .true.
               poin0 = a%PointsINeed(irank+1)%l(1)
            endif
            auxpoin0 = minval(a%PointsINeed(irank+1)%l)
            poin0 = min(poin0,auxpoin0)
         endif
      enddo
      poin0 = poin0-1
      do irank = 0,a%MPIsize-1
         a%PointsINeed(irank+1)%l = a%PointsINeed(irank+1)%l - poin0
      enddo
   end subroutine
   
   subroutine CommunicateReal1(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      real(rp), target :: OldArray(*), NewArray(*)   
      
      real(rp), pointer :: auxOldArray(:,:) => NULL(), auxNewArray(:,:) => NULL()
      
      auxOldArray(1:ndofn,1:a%OldNPoinLocal) => OldArray(1:ndofn*a%OldNpoinLocal)
      auxNewArray(1:ndofn,1:a%NewNPoinLocal) => NewArray(1:ndofn*a%NewNpoinLocal)
   
      call a%CommunicateReal(ndofn,auxOldArray,auxNewArray)
   end subroutine
   
   subroutine CommunicateReal(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      real(rp) :: OldArray(ndofn,*), NewArray(ndofn,*)
      
      type :: CommunicationData
         real(rp), allocatable :: l(:,:)
      end type
      type(CommunicationData), allocatable :: auxCommunicatorOthersNeed(:),auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      !Communications structure
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      !Old to Communications structure
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%l(ndofn,a%NpoinOthersNeed(irank+1)))
            auxCommunicatorOthersNeed(irank+1)%l = OldArray(:,a%PointsOthersNeed(irank+1)%l)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%l(ndofn,a%NpoinINeed(irank+1)))
         endif
      enddo
      
      !Communicate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%l,ndofn*a%NpoinOthersNeed(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%l,ndofn*a%NpoinINeed(irank+1),MPI_REAL8, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            NewArray(:,a%PointsINeed(irank+1)%l) = auxCommunicatorINeed(irank+1)%l
         endif
      enddo
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            deallocate(auxCommunicatorINeed(irank+1)%l)
         endif
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%l)
         endif
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
   end subroutine
   
   subroutine CommunicateInteger1(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      integer(ip), target :: OldArray(*), NewArray(*)   
      
      integer(ip), pointer :: auxOldArray(:,:) => NULL(), auxNewArray(:,:) => NULL()
      
      auxOldArray(1:ndofn,1:a%OldNPoinLocal) => OldArray(1:ndofn*a%OldNpoinLocal)
      auxNewArray(1:ndofn,1:a%NewNPoinLocal) => NewArray(1:ndofn*a%NewNpoinLocal)
   
      call a%CommunicateInteger(ndofn,auxOldArray,auxNewArray)
   end subroutine
   
   subroutine CommunicateInteger(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      integer(ip) :: OldArray(ndofn,*), NewArray(ndofn,*)
      
      type :: CommunicationData
         integer(ip), allocatable :: l(:,:)
      end type
      type(CommunicationData), allocatable :: auxCommunicatorOthersNeed(:),auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      !Communications structure
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      !Old to Communications structure
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%l(ndofn,a%NpoinOthersNeed(irank+1)))
            auxCommunicatorOthersNeed(irank+1)%l = OldArray(:,a%PointsOthersNeed(irank+1)%l)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%l(ndofn,a%NpoinINeed(irank+1)))
         endif
      enddo
      
      !Communicate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%l,ndofn*a%NpoinOthersNeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%l,ndofn*a%NpoinINeed(irank+1),MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            NewArray(:,a%PointsINeed(irank+1)%l) = auxCommunicatorINeed(irank+1)%l
         endif   
      enddo
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            deallocate(auxCommunicatorINeed(irank+1)%l)
         endif
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%l)
         endif
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
   end subroutine
   
   subroutine CommunicateLogical1(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      logical, target :: OldArray(*), NewArray(*)   
      
      logical, pointer :: auxOldArray(:,:) => NULL(), auxNewArray(:,:) => NULL()
      
      auxOldArray(1:ndofn,1:a%OldNPoinLocal) => OldArray(1:ndofn*a%OldNpoinLocal)
      auxNewArray(1:ndofn,1:a%NewNPoinLocal) => NewArray(1:ndofn*a%NewNpoinLocal)
   
      call a%CommunicateLogical(ndofn,auxOldArray,auxNewArray)
   end subroutine
   
   subroutine CommunicateLogical(a,ndofn,OldArray,NewArray)
      use typre
      use MPI 
      implicit none
      class(RebalanceCommunicator) :: a
      integer(ip) :: ndofn
      logical :: OldArray(ndofn,*), NewArray(ndofn,*)
      
      type :: CommunicationData
         logical, allocatable :: l(:,:)
      end type
      type(CommunicationData), allocatable :: auxCommunicatorOthersNeed(:),auxCommunicatorINeed(:)
      
      integer(ip) :: irank
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer, allocatable :: irequest01(:), irequest02(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr
      
      !Communications structure
      allocate(auxCommunicatorOthersNeed(a%MPIsize))
      allocate(auxCommunicatorINeed(a%MPIsize))
      
      !Old to Communications structure
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            allocate(auxCommunicatorOthersNeed(irank+1)%l(ndofn,a%NpoinOthersNeed(irank+1)))
            auxCommunicatorOthersNeed(irank+1)%l = OldArray(:,a%PointsOthersNeed(irank+1)%l)
         endif   
         if (a%NpoinINeed(irank+1) /= 0) then
            allocate(auxCommunicatorINeed(irank+1)%l(ndofn,a%NpoinINeed(irank+1)))
         endif
      enddo
      
      !Communicate
      allocate(irequest01(a%MPIsize))
      allocate(irequest02(a%MPIsize))
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(auxCommunicatorOthersNeed(irank+1)%l,ndofn*a%NpoinOthersNeed(irank+1),MPI_LOGICAL, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(auxCommunicatorINeed(irank+1)%l,ndofn*a%NpoinINeed(irank+1),MPI_LOGICAL, irank, mtag1, a%MPIcomm,irequest02(irank+1), ierr)
         endif
      enddo  
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      
      !Communications structure to new
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            NewArray(:,a%PointsINeed(irank+1)%l) = auxCommunicatorINeed(irank+1)%l
         endif
      enddo
      
      !Deallocations
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            deallocate(auxCommunicatorINeed(irank+1)%l)
         endif
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            deallocate(auxCommunicatorOthersNeed(irank+1)%l)
         endif
      enddo
      deallocate(auxCommunicatorINeed)
      deallocate(auxCommunicatorOthersNeed)
   end subroutine
   
  
   subroutine Dealloc(a)
      use typre
      implicit none
      class(RebalanceCommunicator) :: a
      
      integer(ip) :: irank
      
      do irank = 0,a%MPIsize-1
         if (associated(a%PointsINeed(irank+1)%l))       deallocate(a%PointsINeed(irank+1)%l)
         if (associated(a%PointsOthersNeed(irank+1)%l))  deallocate(a%PointsOthersNeed(irank+1)%l)
      enddo
      deallocate(a%NpoinINeed,a%NpoinOthersNeed)
      deallocate(a%PointsINeed,a%PointsOthersNeed)
      
      
      
   end subroutine

end module

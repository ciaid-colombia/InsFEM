module Mod_MyCommunicator
   use typre
   use Mod_ParallelCommunicatorInterface
   use Mod_Memor
   use Mod_LinkedList
   implicit none
   private
   public MyCommunicator
   
   
   
   type, extends(ParallelCommunicatorInterface) ::  MyCommunicator
      
      integer(ip) :: npoinLocal, npoinGhost, gnpoin
            
      integer(ip) :: MPIsize
      integer(ip) :: MPIrank
      integer(ip) :: MPIcomm
      
      type(MemoryMan), pointer :: Memor => NULL()
      
      type(LinkedListHeader) :: RanksToSendH, RanksToReceiveH
      type(LinkedListStorage) :: RanksToSendS, RanksToReceiveS
      
      integer(ip), allocatable :: NGhostsToReceive(:), NGhostsToSend(:)
      type(i1p), allocatable   :: GhostsToReceive(:), GhostsToSend(:)
      
      type(i1p), allocatable ::  iGhostsToReceive(:,:), iGhostsToSend(:,:)
      type(r1p), allocatable ::  rGhostsToReceive(:,:), rGhostsToSend(:,:)
      type(l1p), allocatable ::  lGhostsToReceive(:,:), lGhostsToSend(:,:)
      
      
      
      logical   :: ikfl_VectorInitialized(10) = .false.
      logical   :: rkfl_VectorInitialized(10) = .false.
      logical   :: lkfl_VectorInitialized(10) = .false.
   
contains   
   
      procedure :: Init => MyInitCommunicator
      procedure :: GhostCommunicateReal => MyGhostCommunicateReal
      procedure :: GhostCommunicateInteger => MyGhostCommunicateInteger
      procedure :: GhostCommunicateReal1 => MyGhostCommunicateReal1
      procedure :: GhostCommunicateInteger1 => MyGhostCommunicateInteger1
      procedure :: GhostCommunicateLogical1 => MyGhostCommunicateLogical1
      procedure :: GhostCommunicateLogical => MyGhostCommunicateLogical
      procedure :: Deallocate => MyDeallocate
      
      
   
   end type

contains

   subroutine MyInitCommunicator(a,MPIComm,MPIsize,MPIrank,npoinLocal,npoinGhost,gnpoin,ParallelOrdering,Memor)
      use typre
      use Mod_Memor
      use Mod_ParallelOrderingInterface
      use MPI 
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: MPIsize,MPIrank,MPIcomm
      class(ParallelOrderingInterface)  :: ParallelOrdering
      integer(ip) :: npoinLocal,npoinGhost,gnpoin
      type(MemoryMan),target :: Memor
      
      integer(ip) :: ipoin,jpoin,poinj
      
      integer(ip), allocatable :: iwa(:)
      logical,     allocatable :: lwa(:)
      
      integer(ip), allocatable :: ProcessorPointRange(:)
      integer(ip), allocatable :: GhostProcessor(:)
      
      integer(ip), allocatable :: GhostLocal2Global(:)
      
      !MPI
      integer :: ierr,irequest1(MPIsize), irequest2(MPIsize)
      integer status(MPI_STATUS_SIZE)
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      integer(ip) :: gipoin, irank,nrank,ranki,istat
      
      a%ikfl_VectorInitialized = .false.
      a%rkfl_VectorInitialized = .false.
      a%MPIcomm = MPIcomm
      a%MPIsize = MPIsize
      a%MPIrank = MPIrank
      
      a%Memor => Memor
      
      a%npoinLocal = npoinLocal
      a%npoinGhost = npoinGhost
      a%gnpoin     = gnpoin
      
      call a%Memor%alloc(a%MPIsize,a%NGhostsToReceive,'a%NGhostsToReceive','MyCommunicator_Init')
      call a%Memor%alloc(a%MPIsize,a%NGhostsToSend,'a%NGhostsToSend','MyCommunicator_Init')
      
      call a%Memor%alloc(a%MPIsize,a%GhostsToReceive,'a%GhostsToReceive','MyCommunicator_Init')
      call a%Memor%alloc(a%MPIsize,a%GhostsToSend,'a%GhostsToSend','MyCommunicator_Init')
      
      !Working arrays
      call a%Memor%alloc(npoinGhost,GhostLocal2Global,'GhostLocal2Global','MyCommunicator_Init')
      call a%Memor%alloc(npoinGhost,GhostProcessor,'GhostProcessor','MyCommunicator_Init')
      !call a%Memor%alloc(a%MPIsize+1,ProcessorPointRange,'ProcessorPointRange','MyCommunicator_Init')
      
      !Build the list of the global numbering of the Ghost points
      do ipoin = 1,a%npoinGhost
         GhostLocal2Global(ipoin) = a%npoinLocal + ipoin
      enddo
      call ParallelOrdering%Local2Global(a%npoinGhost,GhostLocal2Global,GhostLocal2Global)
      
      
      !We gather how many npoinLocal each process has
      !call MPI_AllGather(a%npoinLocal, 1, MPI_INTEGER4, ProcessorPointRange(2), 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
!       !We store the local points range of each processor
!       do irank = 0,a%MPIsize-1
!          ProcessorPointRange(irank+2) = ProcessorPointRange(irank+1)+ProcessorPointRange(irank+2)
!       enddo
!    
!       !Now, for each ghost point, we find to which process it belongs
!       do ipoin = 1,a%npoinGhost
!          gipoin = GhostLocal2Global(ipoin)
!          call GhostDomainSearch(a%MPIsize,ProcessorPointRange,gipoin,irank)
!          GhostProcessor(ipoin) = irank
!       enddo
      do ipoin = 1,a%npoinGhost
         GhostProcessor(ipoin) = a%npoinLocal+ipoin
      enddo
      call ParallelOrdering%GetProc(a%npoinGhost,GhostProcessor,GhostProcessor)
      
      !Now we prepare a list of the number of nodes we need to receive from each processor
      do ipoin = 1,a%npoinGhost
         irank = GhostProcessor(ipoin)
         
         if (irank > a%MPIsize-1) then
            write(*,*) 'here we are'
         endif
         
         a%NGhostsToReceive(irank+1) = a%NGhostsToReceive(irank+1)+1
      enddo
      
      !Everybody sends how many nodes they need from the other processors
      !Receives how many nodes the other processors neeed from him
      call MPI_AllToAll(a%NGhostsToReceive, 1, MPI_INTEGER4, a%NGhostsToSend, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      !Now we allocate and fill the list of nodes we need to receive from others
      !We use a list to avoid going through all ranks
      call a%RanksToReceiveH%Initialize
      call a%RanksToReceiveS%Initialize(a%MPISize,'NonRemovable','NonRepeatable',a%Memor)
      do irank = 0,a%MPIsize-1
         if (a%NGhostsToReceive(irank+1) /= 0) then
            !We add it to the list of ranks
            call a%RanksToReceiveS%Add(a%RanksToReceiveH,irank+1)
            !We allocate the space for storing which nodes
            call a%Memor%palloc(a%NGhostsToReceive(irank+1),a%GhostsToReceive(irank+1)%l,'a%GhostsToReceive','MyCommunicator_Init')
            !allocate(a%GhostsToReceive(irank+1)%l(a%NGhostsToReceive(irank+1)))
         endif
      enddo
      a%NGhostsToReceive = 0
      do ipoin = 1,a%npoinGhost
         irank = GhostProcessor(ipoin)
         a%NGhostsToReceive(irank+1) = a%NGhostsToReceive(irank+1)+1
         a%GhostsToReceive(irank+1)%l(a%NGhostsToReceive(irank+1)) = GhostLocal2Global(ipoin)
      enddo
       
      !Now we allocate the list of nodes others need us to send them
      !We use alist to avoid going through all ranks
      call a%RanksToSendH%Initialize
      call a%RanksToSendS%Initialize(a%MPISize,'NonRemovable','NonRepeatable',a%Memor)
      do irank = 0,a%MPisize-1
         if (a%NGhostsToSend(irank+1) /= 0) then
            !We add it to the list of ranks
            call a%RanksToSendS%Add(a%RanksToSendH,irank+1)
            
            !We allocate the space for storing which nodes
            call a%Memor%palloc(a%NGhostsToSend(irank+1),a%GhostsToSend(irank+1)%l,'a%GhostsToSend','MyCommunicator_Init')
            !allocate(a%GhostsToSend(irank+1)%l(a%NGhostsToSend(irank+1)))
         endif
      enddo
       
      !First we Ireceive what others need from us
      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NGhostsToSend(irank+1) /= 0) then
            call MPI_IRECV(a%GhostsToSend(irank+1)%l,a%NGhostsToSend(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
         endif
      enddo
       
      !Now we Isend what we need from others
      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NGhostsToReceive(irank+1) /= 0) then
            call MPI_ISEND(a%GhostsToReceive(irank+1)%l,a%NGhostsToReceive(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
         endif
      enddo
      
      !Now we wait for everything to be complete
      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1), status, ierr)
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo
      
      !When we have sent the toReceive nodes, we transform them to local numbering
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
         
         call ParallelOrdering%Global2Local(a%NGhostsToReceive(irank+1),a%GhostsToReceive(irank+1)%l,a%GhostsToReceive(irank+1)%l)
      enddo
      
      !We do the same with the ToSend nodes
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
         
         call ParallelOrdering%Global2Local(a%NGhostsToSend(irank+1),a%GhostsToSend(irank+1)%l,a%GhostsToSend(irank+1)%l)
      enddo
      
      
      
      !Finally, we allocate the r1p and i1p for the actual communications
      !We allocate the arrays to be used when GhostCommunicating. For integers
      allocate(a%iGhostsToReceive(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'iGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      allocate(a%iGhostsToSend(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'iGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
      !and for reals
      allocate(a%rGhostsToReceive(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'rGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      allocate(a%rGhostsToSend(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'rGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
      !and for logicals
      allocate(a%lGhostsToReceive(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'lGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      allocate(a%lGhostsToSend(a%MPIsize,10),stat=istat)
      call a%Memor%allocObj(istat,'lGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
      
      !DeallocateWorking arrays
      call a%Memor%dealloc(npoinGhost,GhostLocal2Global,'GhostLocal2Global','MyCommunicator_Init')
      call a%Memor%dealloc(npoinGhost,GhostProcessor,'GhostProcessor','MyCommunicator_Init')
      !call a%Memor%dealloc(a%MPIsize+1,ProcessorPointRange,'ProcessorPointRange','MyCommunicator_Init')      
      
      
   end subroutine
   
   subroutine GhostDomainSearch(MPIsize,ProcessorPointRange,gipoin,rank)
      use typre
      implicit none
      integer(ip) :: MPIsize,ProcessorPointRange(MPIsize+1),gipoin,rank
      
      logical :: not_found
      integer(ip) :: stepsize, modu
      
      not_found = .true.
      
      stepsize = (MPIsize-1)
      modu = mod(stepsize,2)
      stepsize = stepsize/2
      
      rank = stepsize
      do while(not_found)
         if (gipoin > ProcessorPointRange(rank+1) .and. gipoin <= ProcessorPointRange(rank+2)) then
            not_found = .false.
         else
            modu = mod(stepsize,2)
            stepsize = stepsize/2+modu
            if (gipoin <= ProcessorPointRange(rank+1)) then
               rank = rank - stepsize
            elseif (gipoin > ProcessorPointRange(rank+2)) then
               rank = rank + stepsize
            endif
         endif
      enddo
            
   end subroutine   
   
   subroutine MyGhostCommunicateReal1(a,ndime,array)
      use typre
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      real(rp), target    :: array(*)
      
      real(rp), pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call MyGhostCommunicateReal(a,ndime,aux_array)
   end subroutine
   
   subroutine MyGhostCommunicateInteger1(a,ndime,array)
      use typre
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      integer(ip),target :: array(*)
      
      integer(ip), pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call MyGhostCommunicateInteger(a,ndime,aux_array)
   end subroutine
   
   subroutine MyGhostCommunicateLogical1(a,ndime,array)
      use typre
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      logical,target :: array(*)
      
      logical, pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call MyGhostCommunicateLogical(a,ndime,aux_array)
   end subroutine
   
   subroutine MyGhostCommunicateReal(a,ndime,array)
      use typre
      use MPI
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      real(rp)    :: array(ndime,*)
      
      integer(ip) :: irank,ranki,nrank,istat,ispos,ipoin
      
      !MPI
      integer :: ierr,irequest1(a%MPIsize), irequest2(a%MPIsize)
      integer status(MPI_STATUS_SIZE)
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      !Do we need to allocate space for this type of communication
      if (a%rkfl_VectorInitialized(ndime) .eqv. .false.) then
         a%rkfl_VectorInitialized(ndime) = .true.
         
         !When we have sent the toReceive nodes, we transform them to local numbering
         call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
         call a%RanksToReceiveH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToReceiveS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToReceive(irank+1)*ndime,a%rGhostsToReceive(irank+1,ndime)%a,'rGhostsToReceive','MyGhostCommunicateReal')
         enddo
         
         call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
         call a%RanksToSendH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToSendS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToSend(irank+1)*ndime,a%rGhostsToSend(irank+1,ndime)%a,'rGhostsToSend','MyGhostCommunicateReal')
         enddo
      endif   
     
      !Now we need to send what needs to be sent, and receive what needs to be received
     
      !First we place what we need in the send arrays
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToSend(irank+1)
            ispos = (ipoin-1)*ndime
            a%rGhostsToSend(irank+1,ndime)%a(ispos+1:ispos+ndime) = array(1:ndime,a%GhostsToSend(irank+1)%l(ipoin))
         enddo
      enddo
      
      !Now we send and receive
      !First we Isend
      irequest1 = MPI_REQUEST_NULL
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         call MPI_ISEND(a%rGhostsToSend(irank+1,ndime)%a,a%NGhostsToSend(irank+1)*ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      enddo
      
      !Now we Irecv 
      irequest2 = MPI_REQUEST_NULL
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         call MPI_IRECV(a%rGhostsToReceive(irank+1,ndime)%a,a%NGhostsToReceive(irank+1)*ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      enddo
      
      
      !Now we wait for everything to be complete
      !First for Send
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest1(irank+1), status, ierr)
      enddo
      !Second for receive
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo
      
      !Finally we place what we have received in the array
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToReceive(irank+1)
            ispos = (ipoin-1)*ndime
            array(1:ndime,a%GhostsToReceive(irank+1)%l(ipoin)) = a%rGhostsToReceive(irank+1,ndime)%a(ispos+1:ispos+ndime)
         enddo
      enddo   
   
   end subroutine
   
   subroutine MyGhostCommunicateInteger(a,ndime,array)
      use typre
      use MPI
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      integer(ip) :: array(ndime,*)
      
      integer(ip) :: irank,ranki,nrank,istat,ispos,ipoin
      
      !MPI
      integer :: ierr,irequest1(a%MPIsize), irequest2(a%MPIsize)
      integer status(MPI_STATUS_SIZE)
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      !Do we need to allocate space for this type of communication
      if (a%ikfl_VectorInitialized(ndime) .eqv. .false.) then
         a%ikfl_VectorInitialized(ndime) = .true.
         
         !When we have sent the toReceive nodes, we transform them to local numbering
         call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
         call a%RanksToReceiveH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToReceiveS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToReceive(irank+1)*ndime,a%iGhostsToReceive(irank+1,ndime)%l,'iGhostsToReceive','MyGhostCommunicateInteger')
         enddo
         
         call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
         call a%RanksToSendH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToSendS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToSend(irank+1)*ndime,a%iGhostsToSend(irank+1,ndime)%l,'iGhostsToSend','MyGhostCommunicateInteger')
         enddo
      endif   
     
      !Now we need to send what needs to be sent, and receive what needs to be received
     
      !First we place what we need in the send arrays
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToSend(irank+1)
            ispos = (ipoin-1)*ndime
            a%iGhostsToSend(irank+1,ndime)%l(ispos+1:ispos+ndime) = array(1:ndime,a%GhostsToSend(irank+1)%l(ipoin))
         enddo
      enddo
      
      !Now we send and receive
      !First we Isend
      irequest1 = MPI_REQUEST_NULL
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         call MPI_ISEND(a%iGhostsToSend(irank+1,ndime)%l,a%NGhostsToSend(irank+1)*ndime,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      enddo
      
      !Now we Irecv 
      irequest2 = MPI_REQUEST_NULL
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         call MPI_IRECV(a%iGhostsToReceive(irank+1,ndime)%l,a%NGhostsToReceive(irank+1)*ndime,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      enddo
      
      
      !Now we wait for everything to be complete
      !First for Send
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest1(irank+1), status, ierr)
      enddo
      !Second for receive
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo
      
      !Finally we place what we have received in the array
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToReceive(irank+1)
            ispos = (ipoin-1)*ndime
            array(1:ndime,a%GhostsToReceive(irank+1)%l(ipoin)) = a%iGhostsToReceive(irank+1,ndime)%l(ispos+1:ispos+ndime)
         enddo
      enddo   
      
   end subroutine   
   
   subroutine MyGhostCommunicateLogical(a,ndime,array)
      use typre
      use MPI
      implicit none
      
      class(MyCommunicator) :: a
      integer(ip) :: ndime
      logical :: array(ndime,*)
      
      integer(ip) :: irank,ranki,nrank,istat,ispos,ipoin
      
      !MPI
      integer :: ierr,irequest1(a%MPIsize), irequest2(a%MPIsize)
      integer status(MPI_STATUS_SIZE)
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      
      !Do we need to allocate space for this type of communication
      if (a%lkfl_VectorInitialized(ndime) .eqv. .false.) then
         a%lkfl_VectorInitialized(ndime) = .true.
         
         !When we have sent the toReceive nodes, we transform them to local numbering
         call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
         call a%RanksToReceiveH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToReceiveS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToReceive(irank+1)*ndime,a%lGhostsToReceive(irank+1,ndime)%l,'lGhostsToReceive','MyGhostCommunicateInteger')
         enddo
         
         call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
         call a%RanksToSendH%GetNelem(nrank)
         do ranki = 1,nrank
            call a%RanksToSendS%GetNext(irank)
            irank = irank-1
         
            call a%Memor%palloc(a%NGhostsToSend(irank+1)*ndime,a%lGhostsToSend(irank+1,ndime)%l,'lGhostsToSend','MyGhostCommunicateInteger')
         enddo
      endif   
     
      !Now we need to send what needs to be sent, and receive what needs to be received
     
      !First we place what we need in the send arrays
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToSend(irank+1)
            ispos = (ipoin-1)*ndime
            a%lGhostsToSend(irank+1,ndime)%l(ispos+1:ispos+ndime) = array(1:ndime,a%GhostsToSend(irank+1)%l(ipoin))
         enddo
      enddo
      
      !Now we send and receive
      !First we Isend
      irequest1 = MPI_REQUEST_NULL
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
      
         call MPI_ISEND(a%lGhostsToSend(irank+1,ndime)%l,a%NGhostsToSend(irank+1)*ndime,MPI_LOGICAL,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      enddo
      
      !Now we Irecv 
      irequest2 = MPI_REQUEST_NULL
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         call MPI_IRECV(a%lGhostsToReceive(irank+1,ndime)%l,a%NGhostsToReceive(irank+1)*ndime,MPI_LOGICAL,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      enddo
      
      
      !Now we wait for everything to be complete
      !First for Send
      call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
      call a%RanksToSendH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToSendS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest1(irank+1), status, ierr)
      enddo
      !Second for receive
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
         
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo
      
      !Finally we place what we have received in the array
      call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
      call a%RanksToReceiveH%GetNelem(nrank)
      do ranki = 1,nrank
         call a%RanksToReceiveS%GetNext(irank)
         irank = irank-1
      
         do ipoin = 1,a%NGhostsToReceive(irank+1)
            ispos = (ipoin-1)*ndime
            array(1:ndime,a%GhostsToReceive(irank+1)%l(ipoin)) = a%lGhostsToReceive(irank+1,ndime)%l(ispos+1:ispos+ndime)
         enddo
      enddo   
      
   end subroutine   
   
   
   subroutine MyDeallocate(a)
      implicit none
      class(MyCommunicator) :: a
      integer(ip) :: ierr,ndime
      
      integer(ip) :: irank, nrank, istat,ranki
      
      !Now we allocate and fill the list of nodes we need to receive from others
      do irank = 0,a%MPIsize-1
         if (a%NGhostsToReceive(irank+1) /= 0) then
            !We allocate the space for storing which nodes
            call a%Memor%pdealloc(a%NGhostsToReceive(irank+1),a%GhostsToReceive(irank+1)%l,'a%GhostsToReceive','MyCommunicator_Init')
         endif
      enddo
      
      !Now we allocate the list of nodes others need us to send them
      do irank = 0,a%MPisize-1
         if (a%NGhostsToSend(irank+1) /= 0) then
            !We deallocate the space for storing which nodes
            call a%Memor%pdealloc(a%NGhostsToSend(irank+1),a%GhostsToSend(irank+1)%l,'a%GhostsToSend','MyCommunicator_Init')
         endif
      enddo
      
      !Do we need to deallocate
      do ndime = 1,10
         if (a%ikfl_VectorInitialized(ndime) .eqv. .true.) then
            a%ikfl_VectorInitialized(ndime) = .false.
            
            call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
            call a%RanksToReceiveH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToReceiveS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToReceive(irank+1)*ndime,a%iGhostsToReceive(irank+1,ndime)%l,'iGhostsToReceive','MyGhostCommunicateInteger')
            enddo
            
            call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
            call a%RanksToSendH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToSendS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToSend(irank+1)*ndime,a%iGhostsToSend(irank+1,ndime)%l,'iGhostsToSend','MyGhostCommunicateInteger')
            enddo
         endif   
      enddo
            
      !Do we need to allocate space for this type of communication
      do ndime = 1,10
         if (a%rkfl_VectorInitialized(ndime) .eqv. .true.) then
            a%rkfl_VectorInitialized(ndime) = .false.
            
            !When we have sent the toReceive nodes, we transform them to local numbering
            call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
            call a%RanksToReceiveH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToReceiveS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToReceive(irank+1)*ndime,a%rGhostsToReceive(irank+1,ndime)%a,'rGhostsToReceive','MyGhostCommunicateReal')
            enddo
            
            call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
            call a%RanksToSendH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToSendS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToSend(irank+1)*ndime,a%rGhostsToSend(irank+1,ndime)%a,'rGhostsToSend','MyGhostCommunicateReal')
            enddo
         endif  
      enddo
      
      !Do we need to deallocate
      do ndime = 1,10
         if (a%lkfl_VectorInitialized(ndime) .eqv. .true.) then
            a%lkfl_VectorInitialized(ndime) = .false.
            
            call a%RanksToReceiveS%GoToFirstOfList(a%RanksToReceiveH)
            call a%RanksToReceiveH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToReceiveS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToReceive(irank+1)*ndime,a%lGhostsToReceive(irank+1,ndime)%l,'lGhostsToReceive','MyGhostCommunicateInteger')
            enddo
            
            call a%RanksToSendS%GoToFirstOfList(a%RanksToSendH)
            call a%RanksToSendH%GetNelem(nrank)
            do ranki = 1,nrank
               call a%RanksToSendS%GetNext(irank)
               irank = irank-1
            
               call a%Memor%pdealloc(a%NGhostsToSend(irank+1)*ndime,a%lGhostsToSend(irank+1,ndime)%l,'lGhostsToSend','MyGhostCommunicateInteger')
            enddo
         endif   
      enddo
      
      !Deallocate lists
      call a%RanksToSendS%Dealloc
      call a%RanksToReceiveS%Dealloc
      
      !Deallocate arrays
      call a%Memor%dealloc(a%MPIsize,a%NGhostsToReceive,'a%NGhostsToReceive','MyCommunicator_Init')
      call a%Memor%dealloc(a%MPIsize,a%NGhostsToSend,'a%NGhostsToSend','MyCommunicator_Init')
      
      call a%Memor%dealloc(a%MPIsize,a%GhostsToReceive,'a%GhostsToReceive','MyCommunicator_Init')
      call a%Memor%dealloc(a%MPIsize,a%GhostsToSend,'a%GhostsToSend','MyCommunicator_Init')
      
      
      !Finally, we allocate the r1p and i1p for the actual communications
      !We allocate the arrays to be used when GhostCommunicating. For integers
      deallocate(a%iGhostsToReceive,stat=istat)
      call a%Memor%deallocObj(istat,'iGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      deallocate(a%iGhostsToSend,stat=istat)
      call a%Memor%deallocObj(istat,'iGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
      !and for reals
      deallocate(a%rGhostsToReceive,stat=istat)
      call a%Memor%deallocObj(istat,'rGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      deallocate(a%rGhostsToSend,stat=istat)
      call a%Memor%deallocObj(istat,'rGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
      !and for logicals
      deallocate(a%lGhostsToReceive,stat=istat)
      call a%Memor%deallocObj(istat,'lGhostsToReceive','MyCommunicator_Init',a%MPIsize*10)
      deallocate(a%lGhostsToSend,stat=istat)
      call a%Memor%deallocObj(istat,'lGhostsToSend','MyCommunicator_Init',a%MPIsize*10)
       
       
      a%npoinLocal = 0
      a%npoinGhost = 0
      a%gnpoin = 0
   
   end subroutine
   
   


end module

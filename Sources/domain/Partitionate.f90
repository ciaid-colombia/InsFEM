module Mod_Partitionate
   use MPI
   use typre
   use Mod_Mesh
   use Mod_NaivePartition
   implicit none

   
contains   
   
   subroutine PartitionateStoreInfoForScattering(a,NaivePartitioner,PointGlobNumber,PointProcNumber,auxPointGlobNumber,auxPointProcNumber)
      implicit none
      class(FemMesh) :: a
      type(NaivePartition) :: NaivePartitioner
      integer(ip) :: PointGlobNumber(*), PointProcNumber(*),auxPointGlobNumber(*),auxPointProcNumber(*)
      
      !Serial read data
      if (a%kfl_ReadType == 0) then
         call PartitionateStoreInfoForScatteringSerial(a,NaivePartitioner,PointGlobNumber,PointProcNumber,auxPointGlobNumber,auxPointProcNumber)
      
      !Partitionate read data
      elseif(a%kfl_ReadType == 1) then
         call PartitionateStoreInfoForScatteringPartitioned(a,NaivePartitioner,PointGlobNumber,PointProcNumber,auxPointGlobNumber,auxPointProcNumber)
         
      endif
   end subroutine
   
   subroutine PartitionateStoreInfoForScatteringSerial(a,NaivePartitioner,PointGlobNumber,PointProcNumber,auxPointGlobNumber,auxPointProcNumber)
      class(FemMesh) :: a
      type(NaivePartition) :: NaivePartitioner
      integer(ip) :: PointGlobNumber(*),PointProcNumber(*),auxPointGlobNumber(*),auxPointProcNumber(*)
      integer(ip) :: inipoin,endpoin,npoinLocal
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer status(MPI_STATUS_SIZE),irequest1,irequest2
      integer :: ierr
      integer(ip) :: irank
      
      irequest1 = MPI_REQUEST_NULL
      irequest2 = MPI_REQUEST_NULL
      
      !Send the partitioning result from all to root (for reading nodes and boundary conditions)
      if (a%npoinLocal > 0) then
         call MPI_ISEND(PointProcNumber,a%npoinLocal,MPI_INTEGER4,a%MPIroot,mtag1,a%MPIcomm,irequest1,ierr) 
         call MPI_ISEND(PointGlobNumber,a%npoinLocal,MPI_INTEGER4,a%MPIroot,mtag2,a%MPIcomm,irequest2,ierr)
      endif

      !Receive the data in root
      call a%Memor%alloc(a%gnpoin,a%gPointGlobNumber,'gPointGlobNumber','Partitionate')
      call a%Memor%alloc(a%gnpoin,a%gPointProcNumber,'gPointProcNumber','Partitionate') 
      
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call NaivePartitioner%GetIniEndPoint(irank,inipoin,endpoin)
            npoinLocal = endpoin-inipoin+1
            if (npoinLocal > 0) then
               call MPI_RECV(a%gPointProcNumber(inipoin),endpoin-inipoin+1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr)
               call MPI_RECV(a%gPointGlobNumber(inipoin),endpoin-inipoin+1,MPI_INTEGER4,irank,mtag2,a%MPIcomm,status,ierr)
            endif
         enddo
      endif  
      
      call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
      a%RP_poin0 = inipoin-1
      a%RP_npoinLocal = a%npoinLocal
      a%RP_npoinGhost = a%npoinGhost
      a%RP_npoin = a%npoin
      
      call MPI_WAIT(irequest1, status, ierr)
      call MPI_WAIT(irequest2, status, ierr)
   
      CALL MPI_BCAST(a%gPointGlobNumber,a%gnpoin,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%gPointProcNumber,a%gnpoin,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
   end subroutine   
   
   subroutine PartitionateStoreInfoForScatteringPartitioned(a,NaivePartitioner,PointGlobNumber,PointProcNumber,auxPointGlobNumber,auxPointProcNumber)
      use MPI
      class(FemMesh) :: a
      type(NaivePartition) :: NaivePartitioner
      integer(ip) :: PointGlobNumber(*),PointProcNumber(*),auxPointGlobNumber(*),auxPointProcNumber(*)
      
      integer(ip) :: ipoin,irank,inipoin,endpoin
      integer :: ierr
      
      a%RP_npoinLocal = a%npoinLocal
      a%RP_npoinGhost = a%npoinGhost
      a%RP_npoin = a%npoin
      
      call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
      call NaivePartitioner%getiniendpoint(a%MPIrank,inipoin,endpoin)
      a%RP_poin0 = inipoin-1
      
      !Create and Initialize the LocalOrdering
      call a%ParallelLibrary%CreateOrdering(a%RP_LocalOrdering,a%Memor)  !FACTORY
      call a%RP_LocalOrdering%Init(a%MPIcomm,a%npoinLocal+a%npoinGhost,PointGlobNumber,PointProcNumber,a%Memor)
      
      call a%ParallelLibrary%CreateOrdering(a%RP_InitialOrdering,a%Memor)  !FACTORY
      call a%RP_InitialOrdering%Init(a%MPIcomm,a%npoinGhost,auxPointGlobNumber,auxPointProcNumber,a%Memor)
      
      !Allocate auxs
      call a%Memor%alloc(a%MPIsize,a%RP_NpoinINeed,'RP_NpoinINeed','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,a%RP_NpoinOthersNeed,'RP_NpoinOthersNeed','reaCoordPartitioned')
      
      !Count the number of points others need
      do ipoin = 1,a%RP_npoinLocal
         call a%RP_LocalOrdering%GetProc(ipoin,irank)
         a%RP_NpoinOthersNeed(irank+1) = a%RP_NpoinOthersNeed(irank+1)+1
      enddo
      
      !Amount of required communications
      call MPI_AllToAll(a%RP_NpoinOthersNeed, 1, MPI_INTEGER4, a%RP_NpoinINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr)
      
   end subroutine
   
   
   subroutine StatisticsAndElementPostprocessing(a)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      
      integer(ip), allocatable :: Elem0Array(:)
      integer(ip) :: inelem
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize),irequest4(a%MPIsize),irequest5(a%MPIsize),irequest6(a%MPIsize)
      integer :: ierr
      integer(ip) :: irank

      !For Statistics and elements postprocessing
      call MPI_ISEND(a%nelem,1,MPI_INTEGER4,a%MPIroot,mtag1,a%MPIcomm,irequest1(a%MPIroot+1),ierr) 
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%alloc(a%MPIsize,Elem0Array,'Elem0Array','Partitionate')
         Elem0Array(1) = 0
         a%MaxNelemPerProc = 0
         do irank = 0,a%MPIsize-1
            call MPI_RECV(inelem,1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr) 
            a%MaxNelemPerProc = max(a%MaxNelemPerProc,inelem)
            if (irank < a%MPIsize-1) then
               Elem0Array(irank+2) = Elem0Array(irank+1) + inelem
            endif
            call MPI_ISEND(Elem0Array(irank+1),1,MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest2(irank+1),ierr) 
         enddo
      endif
      call MPI_WAIT(irequest1(a%MPIroot+1), status, ierr)
      call MPI_RECV(a%elem0,1_ip, MPI_INTEGER4, a%MPIroot, mtag2, a%MPIcomm,status, ierr) 
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest2(a%MPIroot+1), status, ierr)
         enddo
         call a%Memor%dealloc(a%MPIsize,Elem0Array,'Elem0Array','Partitionate')
      endif
   end subroutine
   
   
   
   
   !For scattering elements
   subroutine PartitionateStoreInfoForScatteringElements(a)
      implicit none
      class(FemMesh) :: a
      type(NaivePartition) :: NaivePartitioner
      
      !Serial read data
      if (a%kfl_ReadType == 0) then
         call PartitionateStoreInfoForScatteringSerialElements(a)
      
      !Partitionate read data
      elseif(a%kfl_ReadType == 1) then
         call PartitionateStoreInfoForScatteringPartitionedElements(a)
         
      endif
   end subroutine
   
   subroutine PartitionateStoreInfoForScatteringSerialElements(a)
      use MPI
      implicit none
      class(FemMesh) :: a
      integer(ip), allocatable :: Local2Initial(:)
      type(i1p), allocatable :: InitialElements(:)
      
      integer(ip), allocatable :: ElementSendCount(:)
      integer(ip) :: elemi,ielem,ierr,irequest1,irank
      integer :: status(MPI_STATUS_SIZE)
      integer(ip), parameter :: mtag1=1, mtag2=2,mtag3=3, mtag20 = 20
      integer(ip), allocatable :: nelems(:)
      integer(ip) :: auxcount
      
      !ElementNaiveSendToProcessor is no longer needed in this case
      auxcount = 0
      do ielem = 1,size(a%ElementNaiveSendToProcessor)
         if (associated(a%ElementNaiveSendToProcessor(ielem)%l)) then
            auxcount = auxcount + size(a%ElementNaiveSendToProcessor(ielem)%l)
            deallocate(a%ElementNaiveSendToProcessor(ielem)%l)
         endif
      enddo
      call a%Memor%deallocObj(0,'gElementInitialSendToProcessor','Partitionate',auxcount*ip)
      call a%Memor%dealloc(size(a%ElementNaiveSendToProcessor),a%ElementNaiveSendToProcessor,'gElementInitialSendToProcessor','Partitionate')
      
      call a%Memor%alloc(a%nelem,Local2Initial,'Local2Initial','Partitionate')
      call a%Memor%alloc(a%MPIsize,nelems,'nelems','Partitionate')
            
      do ielem = 1,a%nelem
         Local2Initial(ielem) = ielem
      enddo   
      
      call a%ElementLocal2Initial%Local2Global(a%nelem,Local2Initial,Local2Initial)
      
      call MPI_Gather( a%nelem, 1, MPI_INTEGER4, nelems,1, MPI_INTEGER4, a%MPIroot, a%MPIcomm,ierr); 
      
      irequest1 =  MPI_REQUEST_NULL
      call MPI_ISEND(Local2Initial, a%nelem, MPI_INTEGER4, a%MPIroot, mtag20, a%MPIcomm,irequest1, ierr)
      
      call a%Memor%alloc(a%MPIsize,InitialElements,'nelems','Partitionate')
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            allocate(InitialElements(irank+1)%l(nelems(irank+1)))
         
            call MPI_RECV(InitialElements(irank+1)%l,nelems(irank+1),MPI_INTEGER4,irank,mtag20,a%MPIcomm,status,ierr) 
         enddo
      endif
      call MPI_WAIT(irequest1, status, ierr)   
      call a%Memor%dealloc(a%nelem,Local2Initial,'Local2Initial','Partitionate')
      
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%alloc(a%gnelem,ElementSendCount,'ElementSendCount','Partitionate')
         ElementSendCount = 0
         
         allocate(a%gElementInitialSendToProcessor(a%gnelem))
         call a%Memor%allocObj(0,'gElementInitialSendToProcessor','Partitionate',a%gnelem*ip)
         
         !First to count
         do irank = 0,a%MPIsize-1
            do ielem = 1,size(InitialElements(irank+1)%l)
               elemi = InitialElements(irank+1)%l(ielem)
               ElementSendCount(elemi) = ElementSendCount(elemi)+1
            enddo
         enddo
         auxcount = 0
         do ielem = 1,a%gnelem
            allocate(a%gElementInitialSendToProcessor(ielem)%l(ElementSendCount(ielem)))
            auxcount = auxcount + ElementSendCount(ielem)
         enddo
         call a%Memor%allocObj(0,'gElementInitialSendToProcessor','Partitionate',auxcount*ip)
         
         !Second to fill
         ElementSendCount = 0
         do irank = 0,a%MPIsize-1
            do ielem = 1,size(InitialElements(irank+1)%l)
               elemi = InitialElements(irank+1)%l(ielem)
               ElementSendCount(elemi) = ElementSendCount(elemi)+1
               a%gElementInitialSendToProcessor(elemi)%l(ElementSendCount(elemi)) = irank
            enddo
            deallocate(InitialElements(irank+1)%l)
         enddo
         call a%Memor%dealloc(a%gnelem,ElementSendCount,'ElementSendCount','Partitionate')
      endif  
      call a%Memor%dealloc(a%MPIsize,InitialElements,'nelems','Partitionate')
      
      call a%Memor%dealloc(a%MPIsize,nelems,'nelems','Partitionate')

      
      
   end subroutine
   
   subroutine PartitionateStoreInfoForScatteringPartitionedElements(a)
      class(FemMesh), target :: a
      
      !We just need to copy 
      a%gElementInitialSendToProcessor => a%ElementNaiveSendToProcessor
      
   end subroutine


end module


subroutine Partitionate(a)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_BuildGraph
   use Mod_NaivePartition
   use Mod_ParallelPartitionerInterface
   use MPI
   use Mod_Partitionate
   implicit none
   class(FemMesh) :: a
   
   class(ParallelPartitionerInterface), pointer :: Partitioner => NULL()
   
   integer(ip) :: inipoin,endpoin,ielem,inode,ipoin,ispos,jv,pnode
   integer(ip), allocatable :: ia(:),ja(:),lelpo(:),pelpo(:),PointGlobNumber(:),PointProcNumber(:)
   integer(ip), allocatable :: Renumbering(:),iRenumbering(:),RankTouch(:)
   integer(ip), allocatable :: auxPointProcNumber(:), auxPointGlobNumber(:)
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3, mtag4 = 4
   integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize),irequest4(a%MPIsize),irequest5(a%MPIsize),irequest6(a%MPIsize),irequest7(a%MPIsize),irequest8(a%MPIsize)
   integer :: ierr
   integer(ip) :: irank
   
   integer(ip) :: alloc_count,alloc_size,jelem,jspos,lsize,nelem,elemi, istatus
   
   integer(ip), allocatable :: SizesINeed(:), SizesOthersNeed(:)
   integer(ip), allocatable :: iNodesINeed(:),jNodesINeed(:),jNodesINeed2(:),iNodesOthersNeed(:),jNodesOthersNeed(:),jNodesOthersNeed2(:)  !1: Glob number, 2: Processor number
   
   integer(ip), allocatable :: SizesElementsINeed(:,:), SizesElementsOthersNeed(:,:)
   type(i1p), allocatable :: iElementsINeed(:), jElementsINeed(:),iElementsOthersNeed(:), jElementsOthersNeed(:)
   
   integer(ip), allocatable :: auxnelty(:),auxlnods(:),auxpnods(:),lstart(:),pstart(:)
   integer(ip) :: ielty,maxnelem
   
   integer(ip), allocatable :: RankList(:)
   integer(ip) :: icount_rank,ranki,auxcount
   
   real(rp) :: time2,time1
   
   integer(ip) :: inelem
   integer(ip), allocatable :: Elem0Array(:), auxElementRenumbering(:),auxLocal2InitialElementNumbering(:), ProcList(:)

   type(NaivePartition) :: NaivePartitioner
   integer(ip) :: initialElem,auxnelem
   type(i1p), allocatable :: gInitialNumberingElementsINeed(:),gInitialNumberingElementsOthersNeed(:)

   call a%Timer%Partitionate%Tic
   
   a%npoinGhost = 0
   
   !Initialize NaivePartition 
   call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
   
   !Renumbering for Building the graph 
   !This is so that there are no empty node numbers when we send to Buildgraph
   !(npoin = npoinLocal + npoinGhost, not gnpoin)
   !We undo the renumbering after the graph has been built
   call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
   a%npoinLocal = endpoin - inipoin + 1
   
   call a%Memor%alloc(a%gnpoin,Renumbering,'Renumbering','Partitionate')
   call a%Memor%alloc(a%gnpoin,iRenumbering,'iRenumbering','Partitionate')
   
   do ielem = 1,a%nelem
      call getlnods(ielem,pnode,ispos,a%gpnods,a%glnods,a%nelty,a%nnode)
      do inode = 1,pnode
         ipoin = a%glnods(ispos+inode)
         if (Renumbering(ipoin) == 0) then
            if (ipoin >= inipoin .and. ipoin <= endpoin) then
               Renumbering(ipoin) = ipoin-inipoin+1
               iRenumbering(ipoin-inipoin+1) = ipoin
            else
               a%npoinGhost = a%npoinGhost + 1
               Renumbering(ipoin) = a%npoinGhost + a%npoinLocal
               iRenumbering(a%npoinLocal+a%npoinGhost) = ipoin
            endif
         endif
         a%glnods(ispos+inode) = Renumbering(ipoin)
      enddo
   enddo
   
   a%npoin = a%npoinLocal + a%npoinGhost
   
   
   !Build the Graph
   call BuildElpo(pelpo,lelpo,a%npoin,a%nelem,a%gpnods,a%glnods,a%nelty,a%nnode,a%Memor)     !Compute glelpo and gpelpo
   
   call BuildGraph(ia,ja,pelpo,lelpo,a%npoin,a%nelem,a%gpnods,a%glnods,a%nelty,a%nnode,a%Memor)
   
   !Undo the Renumbering
   do ielem = 1,a%nelem
      call getlnods(ielem,pnode,ispos,a%gpnods,a%glnods,a%nelty,a%nnode)
      do inode = 1,pnode
         ipoin = a%glnods(ispos+inode)
         a%glnods(ispos+inode) = iRenumbering(ipoin)
      enddo
   enddo
   
   !Undo the renumbering for ja
   do jv = 1,ia(a%npoin+1)-1
      ipoin = ja(jv)
      ja(jv) = iRenumbering(ipoin)
   enddo
   
   !Memory allocations
   call a%Memor%alloc(a%npoin,PointProcNumber,'PointProcNumber','Partitionate')
   call a%Memor%alloc(a%npoin,PointGlobNumber,'PointGlobNumber','Partitionate')
   
   !Copy iRenumbering to Ghost points positions in PointGlobNumber
   PointGlobNumber(a%npoinLocal+1:a%npoin) = iRenumbering(a%npoinLocal+1:a%npoinLocal+a%npoinGhost)
   call a%Memor%dealloc(a%gnpoin,iRenumbering,'iRenumbering','Partitionate')
   
   !-----------------------------------------------------------------------------------------------------------
   !Partitionate
   call a%ParallelLibrary%CreatePartitioner(Partitioner,a%Memor)
   call Partitioner%Initialize(a%Memor)
   call Partitioner%GraphPart(a%MPIcomm,a%npoinLocal,a%gnpoin,ia,ja,PointProcNumber,PointGlobNumber)
   call a%ParallelLibrary%DeallocatePartitioner(Partitioner,a%Memor)
   !-----------------------------------------------------------------------------------------------------------
   
   call a%Memor%dealloc(ia(a%npoin+1)-1,ja,'ja','partitionate')
   call a%Memor%dealloc(a%npoin+1,ia,'ia','partitionate')
   
   call a%Memor%dealloc(pelpo(a%npoin+1)-1,lelpo,'lelpo','Partitionate')
   call a%Memor%dealloc(a%npoin+1,pelpo,'pelpo','Partitionate')
   
   !Decide how many data I need to receive and send from/to each processor
   call a%Memor%alloc(a%MPIsize,SizesINeed,'SizesINeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,SizesOthersNeed,'SizesOthersNeed','Partitionate')
   
   !For ghost points I still do not know the global numbering and processor
   !How many nodes I do not know the numbering and the processor of
   !In the ghost nodes, PointGlobNumber still stores the original numbering!
   SizesINeed = 0
   do inode = a%npoinLocal+1,a%npoin
      ipoin = PointGlobNumber(inode)
      call NaivePartitioner%GetPointRank(ipoin,irank)
      SizesINeed(irank+1) = SizesINeed(irank+1)+1
   enddo   
   
   !Send how many nodes I need from others, receive how many nodes others need from me
   call MPI_AllToAll(SizesINeed, 1, MPI_INTEGER4, SizesOthersNeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
   
   !Some memory allocations
   !This needed to be done this way in order to avoid crash with MPIGNU release
   !call a%Memor%alloc(a%MPIsize+1,iNodesINeed,'iNodesINeed','Partitionate')
   !call a%Memor%alloc(a%npoinGhost,jNodesINeed,'jNodesINeed','Partitionate')
   !call a%Memor%alloc(a%MPIsize+1,iNodesOthersNeed,'iNodesOthersNeed','Partitionate')
   !call a%Memor%alloc(sum(SizesOthersNeed),jNodesOthersNeed,'jNodesOthersNeed','Partitionate')
   allocate(iNodesINeed(a%MPIsize+1),STAT=istatus)
   call a%Memor%allocObj(istatus,'iNodesINeed','Partitionate',(a%MPIsize+1)*ip)
   allocate(iNodesOthersNeed(a%MPIsize+1),STAT=istatus)
   call a%Memor%allocObj(istatus,'iNodesOthersNeed','Partitionate',(a%MPIsize+1)*ip)
   allocate(jNodesINeed(a%npoinGhost),STAT=istatus)
   call a%Memor%allocObj(istatus,'jNodesINeed2','Partitionate',a%npoinGhost*ip)
   allocate(jNodesINeed2(a%npoinGhost),STAT=istatus)
   call a%Memor%allocObj(istatus,'jNodesINeed','Partitionate',a%npoinGhost*ip)

   allocate(jNodesOthersNeed(sum(SizesOthersNeed)),STAT=istatus)
   call a%Memor%allocObj(istatus,'jNodesOthersNeed','Partitionate',sum(SizesOthersNeed)*ip)
   allocate(jNodesOthersNeed2(sum(SizesOthersNeed)),STAT=istatus)
   call a%Memor%allocObj(istatus,'jNodesOthersNeed2','Partitionate',sum(SizesOthersNeed)*ip)
   
   !--------------------------------------------------------------------------------------------------------------------------
   !Send which ghost nodes I need from others 
   !Starting point
   iNodesINeed(1) = 1
   do irank = 0,a%MPIsize -1
      iNodesINeed(irank+2) = iNodesINeed(irank+1) + SizesINeed(irank+1)
   enddo
   
   !List of nodes
   do inode = a%npoinLocal+1,a%npoin
      ipoin = PointGlobNumber(inode)
      call NaivePartitioner%GetPointRank(ipoin,irank)
      jNodesINeed(iNodesINeed(irank+1)) = ipoin
      iNodesINeed(irank+1) = iNodesINeed(irank+1) + 1
   enddo
   
   !Reset starting point
   iNodesINeed(1) = 1
   do irank = 0,a%MPIsize -1
      iNodesINeed(irank+2) = iNodesINeed(irank+1) + SizesINeed(irank+1)
   enddo
   
   !Send  
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesINeed(irank+1) /= 0) then
         ispos = iNodesINeed(irank+1)
         call MPI_ISEND(jNodesINeed(ispos),SizesINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
     
   
   !----------------------------------------------------------------------------------------------------------------------------
   !Receive which ghost nodes others need from me
   !Starting point
   iNodesOthersNeed(1) = 1
   do irank = 0,a%MPIsize -1
      iNodesOthersNeed(irank+2) = iNodesOthersNeed(irank+1) + SizesOthersNeed(irank+1)
   enddo
   
   !Receive the nodes
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesOthersNeed(irank+1) /= 0) then
         ispos = iNodesOthersNeed(irank+1)
         call MPI_IRECV(jNodesOthersNeed(ispos),SizesOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !-----------------------------------------------------------------------------------------------------------------------------
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   !-----------------------------------------------------------------------------------------------------------------------------
   !Send what others need from me
   call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
   do inode = 1,iNodesOthersNeed(a%MPIsize+1)-1
      ipoin = jNodesOthersNeed(inode)
      ipoin = ipoin - inipoin + 1
      jNodesOthersNeed(inode) = PointGlobNumber(ipoin)  !1 Global POint
      jNodesOthersNeed2(inode) = PointProcNumber(ipoin)  !2 Processor number
   enddo
   
   !Send the nodes
   irequest1 = MPI_REQUEST_NULL
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesOthersNeed(irank+1) /= 0) then
         ispos = iNodesOthersNeed(irank+1)
         call MPI_ISEND(jNodesOthersNeed(ispos),SizesOthersNeed(irank+1),MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest1(irank+1),ierr) 
         call MPI_ISEND(jNodesOthersNeed2(ispos),SizesOthersNeed(irank+1),MPI_INTEGER4,irank,mtag3,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !Receive what I need from others
   irequest3 = MPI_REQUEST_NULL
   irequest4 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesINeed(irank+1) /= 0) then
         ispos = iNodesINeed(irank+1)
         call MPI_IRECV(jNodesINeed(ispos),SizesINeed(irank+1),MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest3(irank+1),ierr) 
         call MPI_IRECV(jNodesINeed2(ispos),SizesINeed(irank+1),MPI_INTEGER4,irank,mtag3,a%MPIcomm,irequest4(irank+1),ierr) 
      endif
   enddo

   !-----------------------------------------------------------------------------------------------------------------------------
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
      call MPI_WAIT(irequest3(irank+1), status, ierr)
      call MPI_WAIT(irequest4(irank+1), status, ierr)
   enddo
   
   !----------------------------------------------------------------------------------------------------------------------------
   !Use the received data to transform lnods to new Global Numbering
   
   !Modify PointGlobNumber(ghostpoints) so that it is in new Global Numbering
   !SizesINeed is just a counter so that I keep updating properly
   
   !Aux info which needs to be stored
   call a%Memor%alloc(a%npoinghost,auxPointGlobNumber,'auxPointGlobNumber','Partitionate')
   call a%Memor%alloc(a%npoinghost,auxPointProcNumber,'auxPointProcNumber','Partitionate')
   auxPointProcNumber = PointProcNumber(a%npoinLocal+1:a%npoin)
   auxPointGlobNumber = PointGlobNumber(a%npoinLocal+1:a%npoin)
   
   SizesINeed = 0
   do inode = a%npoinLocal+1,a%npoin
      ipoin = PointGlobNumber(inode)
      call NaivePartitioner%GetPointRank(ipoin,irank)
      SizesINeed(irank+1) = SizesINeed(irank+1) +1
      PointGlobNumber(inode) = jNodesINeed(iNodesINeed(irank+1)+SizesINeed(irank+1)-1)
      PointProcNumber(inode) = jNodesINeed2(iNodesINeed(irank+1)+SizesINeed(irank+1)-1)
   enddo  
   
   !------------------------------------------------------------------------------------------------!
   !                                                                                                !   
   ! Info on the original numbering to parallel numbering needs to be stored                        !
       call PartitionateStoreInfoForScattering(a,NaivePartitioner,PointGlobNumber,PointProcNumber, &!
                                                 auxPointGlobNumber, auxPointProcNumber)            !
   !                                                                                                ! 
   !------------------------------------------------------------------------------------------------!
   
   
   !Deallocate some working arrays
   call a%Memor%dealloc(a%npoinghost,auxPointGlobNumber,'auxPointGlobNumber','Partitionate')
   call a%Memor%dealloc(a%npoinghost,auxPointProcNumber,'auxPointProcNumber','Partitionate')
   call a%Memor%dealloc(a%MPIsize+1,iNodesINeed,'iNodesINeed','Partitionate')
   call a%Memor%dealloc(a%npoinGhost,jNodesINeed,'jNodesINeed','Partitionate')
   call a%Memor%dealloc(a%npoinGhost,jNodesINeed2,'jNodesINeed2','Partitionate')
   
   call a%Memor%dealloc(a%MPIsize+1,iNodesOthersNeed,'iNodesOthersNeed','Partitionate')
   call a%Memor%dealloc(size(jNodesOthersNeed),jNodesOthersNeed,'jNodesOthersNeed','Partitionate')
   call a%Memor%dealloc(size(jNodesOthersNeed2),jNodesOthersNeed2,'jNodesOthersNeed2','Partitionate')
   
   call a%Memor%dealloc(a%MPIsize,SizesINeed,'SizesINeed','Partitionate')
   call a%Memor%dealloc(a%MPIsize,SizesOthersNeed,'SizesOthersNeed','Partitionate')
   
   
   
   !-------------------------------------------------------------------------------------------------------------------------
   !Decide which processors I have elements of   
   call a%Memor%alloc(2,a%MPIsize,SizesElementsINeed,'SizesElementsINeed','Partitionate') !(1): number of elements (2) size of glnods
   call a%Memor%alloc(2,a%MPIsize,SizesElementsOthersNeed,'SizesElementsOthersNeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,iElementsINeed,'iElementsINeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,iElementsOthersNeed,'iElementsOthersNeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,jElementsINeed,'jElementsINeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,jElementsOthersNeed,'jElementsOthersNeed','Partitionate')

   
   !First to count
   call a%Memor%alloc(a%MPIsize,RankTouch,'RankTouch','Partitionate')
   RankTouch = -1
   call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
   do ielem = 1,a%nelem
      call getlnods(ielem,pnode,ispos,a%gpnods,a%glnods,a%nelty,a%nnode)
      !Only one of the processors who has the element sends it to everybody
      if (minval(a%glnods(ispos+1:ispos+pnode))>=inipoin) then
         do inode = 1,pnode
            ipoin = a%glnods(ispos+inode)
            irank = PointProcNumber(Renumbering(ipoin))
      
            if (RankTouch(irank+1) /= ielem) then
               RankTouch(irank+1) = ielem
               SizesElementsOthersNeed(1,irank+1) = SizesElementsOthersNeed(1,irank+1) + 1
               SizesElementsOthersNeed(2,irank+1) = SizesElementsOthersNeed(2,irank+1) + pnode
            endif
         enddo
      endif
   enddo
   
   !Now send their sizes to all, receive the sizes of  what I need from all
   call MPI_AllToAll(SizesElementsOthersNeed, 2, MPI_INTEGER4, SizesElementsINeed, 2, MPI_INTEGER4, a%MPIcomm,ierr);
   
   !Initial numbering of elements needs also to be sent
   call a%Memor%alloc(a%MPIsize,gInitialNumberingElementsOthersNeed,'gInitialNumberingElementsOthersNeed','Partitionate')
   call a%Memor%alloc(a%MPIsize,gInitialNumberingElementsINeed,'gInitialNumberingElementsOthersNeed','Partitionate')
   do irank = 0,a%MPIsize-1
      call a%Memor%palloc(SizesElementsOthersNeed(1,irank+1),gInitialNumberingElementsOthersNeed(irank+1)%l,'gInitialNumberingElementsOthersNeed','Partitionate')
      call a%Memor%palloc(SizesElementsINeed(1,irank+1),gInitialNumberingElementsINeed(irank+1)%l,'gInitialNumberingElementsINeed','Partitionate')
   enddo
   
   
   !-------------------------------------------------------------------------------------------------------------
   !In the second Loop I do 2 things:
   !  1) Fill the list of elements which need to be sent to other processors
   !  2) Update lnods to the new Global Numbering
   !-------------------------------------------------------------------------------------------------------------
   
   !Allocate as counted
   alloc_count = 0
   do irank = 0,a%MPIsize-1
      alloc_size = SizesElementsOthersNeed(1,irank+1)
      allocate(iElementsOthersNeed(irank+1)%l(alloc_size+1))
      alloc_count = alloc_count + alloc_size + 1
      
      alloc_size = SizesElementsOthersNeed(2,irank+1)
      allocate(jElementsOthersNeed(irank+1)%l(alloc_size))
      alloc_count = alloc_count + alloc_size
   enddo
   call a%Memor%allocObj(0,'ElementsOthersNeed','Partitionate',alloc_count*ip)
      
   !Initialize Number of elements to sent
   SizesElementsOthersNeed(1,:) = 0
   do irank = 0,a%MPIsize-1
      iElementsOthersNeed(irank+1)%l(1) = 1
   enddo
  
   !Start the loop to: 1) Fill the list of elements to be sent
   !                   2) Update lnods to the new Global Numbering
   
   !Auxiliary list
   call a%Memor%alloc(a%MPIsize,RankList,'RankList','Partitionate') 
   
   !we count gElementInitialSendToProcessor, because later it is pointed here
   call a%Memor%alloc(a%nelem,a%ElementNaiveSendToProcessor,'gElementInitialSendToProcessor','Partitionate')
   
   RankTouch = -1
   auxcount = 0
   do ielem = 1,a%nelem
      icount_rank = 0
      call getlnods(ielem,pnode,ispos,a%gpnods,a%glnods,a%nelty,a%nnode)
      !Only one of the processors who has the element sends it to everybody
      if (minval(a%glnods(ispos+1:ispos+pnode))>=inipoin) then
         do inode = 1,pnode
            ipoin = a%glnods(ispos+inode)

            if (ipoin >= inipoin .and. ipoin <= endpoin) then
               !For filling the list of elements
               irank = PointProcNumber(ipoin-inipoin+1)
               
               !For updating lnods
               a%glnods(ispos+inode) = pointGlobNumber(ipoin-inipoin+1)
            else
               ipoin = Renumbering(ipoin) 
               
               !For Filling the list of elements
               irank = PointProcNumber(ipoin)
               
               !For updating lnods
               a%glnods(ispos+inode) = PointGlobNumber(ipoin)
            endif
            
            !For updating lnods
            if (RankTouch(irank+1) /= ielem) then
               RankTouch(irank+1) = ielem
               icount_rank = icount_rank + 1 
               RankList(icount_rank) = irank
            endif
            
         enddo
         
         !it is not convenient to use Memor here
         !we count gElementInitialSendToProcessor, because later it is pointed here
         allocate(a%ElementNaiveSendToProcessor(ielem)%l(icount_rank))
         auxcount = auxcount + icount_rank
         !call a%Memor%palloc(icount_rank,a%ElementNaiveSendToProcessor(ielem)%l,'ElementNaiveSendToProcessor','Partitionate')
         
         
         a%ElementNaiveSendToProcessor(ielem)%l = RankList(1:icount_rank)
      
         !For all the processors I need to send  the element
         do ranki = 1,icount_rank
            irank = RankList(ranki)
            SizesElementsOthersNeed(1,irank+1) = SizesElementsOthersNeed(1,irank+1) + 1
            jelem = SizesElementsOthersNeed(1,irank+1)
            iElementsOthersNeed(irank+1)%l(jelem+1) = iElementsOthersNeed(irank+1)%l(jelem) + pnode
            jspos = iElementsOthersNeed(irank+1)%l(jelem)-1
            jElementsOthersNeed(irank+1)%l(jspos+1:jspos+pnode) = a%glnods(ispos+1:ispos+pnode)
            
            call a%ElementNaive2Initial%Local2Global(ielem,initialElem)
            
            !Initial Element numbering
            gInitialNumberingElementsOthersNeed(irank+1)%l(jelem) = initialElem
         
         enddo
      endif
   enddo
   !we count gElementInitialSendToProcessor, because later it is pointed here
   call a%Memor%AllocObj(0,'gElementInitialSendToProcessor','Partitionate',auxcount*ip)
   
   !Deallocate Auxiliary list
   call a%Memor%dealloc(a%MPIsize,RankList,'RankList','Partitionate') 
   
   !Deallocate some working arrays
   
   call a%Memor%dealloc(a%gnpoin,Renumbering,'Renumbering','Partitionate')

   call a%Memor%dealloc(a%npoin,PointProcNumber,'PointProcNumber','Partitionate')
   call a%Memor%dealloc(a%npoin,PointGlobNumber,'PointGlobNumber','Partitionate')
   
   !Deallocate some working arrays
   call a%Memor%dealloc(a%MPIsize,RankTouch,'RankTouch','Partitionate')
   
   !I can now deallocate glnods and gpnods
   call a%Memor%dealloc(a%gpnods(a%nelem+1)-1,a%glnods,'glnods','Partitionate')
   call a%Memor%dealloc(a%nelem+1,a%gpnods,'gpnods','Partitionate')

   
   !Allocate what I need to receive
   alloc_count = 0
   do irank = 0,a%MPIsize-1
      alloc_size = SizesElementsINeed(1,irank+1)
      allocate(iElementsINeed(irank+1)%l(alloc_size+1))
      alloc_count = alloc_count + alloc_size + 1
      
      alloc_size = SizesElementsINeed(2,irank+1)
      allocate(jElementsINeed(irank+1)%l(alloc_size))
      alloc_count = alloc_count + alloc_size
   enddo
   call a%Memor%allocObj(0,'ElementsINeed','Partitionate',alloc_count*ip)
    
   !Send everything
   irequest1 = MPI_REQUEST_NULL
   irequest2 = MPI_REQUEST_NULL
   irequest3 = MPI_REQUEST_NULL
   irequest7 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesElementsOthersNeed(1,irank+1) /= 0) then
         nelem = SizesElementsOthersNeed(1,irank+1)
         lsize = SizesElementsOthersNeed(2,irank+1)
         
         
         call MPI_ISEND(iElementsOthersNeed(irank+1)%l,nelem+1,MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest2(irank+1),ierr) 
         call MPI_ISEND(jElementsOthersNeed(irank+1)%l,lsize,MPI_INTEGER4,irank,mtag3,a%MPIcomm,irequest3(irank+1),ierr) 
         
         call MPI_ISEND(gInitialNumberingElementsOthersNeed(irank+1)%l,nelem,MPI_INTEGER4,irank,mtag4,a%MPIcomm,irequest7(irank+1),ierr) 
      endif
   enddo
   
   !Receive Everything
   irequest4 = MPI_REQUEST_NULL
   irequest5 = MPI_REQUEST_NULL
   irequest6 = MPI_REQUEST_NULL
   irequest8 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesElementsINeed(1,irank+1) /= 0) then
         nelem = SizesElementsINeed(1,irank+1)
         lsize = SizesElementsINeed(2,irank+1)

         call MPI_IRECV(iElementsINeed(irank+1)%l,nelem+1,MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest5(irank+1),ierr) 
         call MPI_IRECV(jElementsINeed(irank+1)%l,lsize,MPI_INTEGER4,irank,mtag3,a%MPIcomm,irequest6(irank+1),ierr) 
         
         call MPI_IRECV(gInitialNumberingElementsINeed(irank+1)%l,nelem,MPI_INTEGER4,irank,mtag4,a%MPIcomm,irequest8(irank+1),ierr) 
      endif
   enddo
   
   !Wait
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
      call MPI_WAIT(irequest3(irank+1), status, ierr)
      call MPI_WAIT(irequest4(irank+1), status, ierr)
      call MPI_WAIT(irequest5(irank+1), status, ierr)
      call MPI_WAIT(irequest6(irank+1), status, ierr)
      call MPI_WAIT(irequest7(irank+1), status, ierr)
      call MPI_WAIT(irequest8(irank+1), status, ierr)
   enddo
   
   !Now I can Deallocate the Elements the Others Need, as they have already been sent
   alloc_count = 0
   do irank = 0,a%MPIsize-1
      alloc_size = SizesElementsOthersNeed(1,irank+1)
      deallocate(iElementsOthersNeed(irank+1)%l)
      alloc_count = alloc_count + alloc_size + 1
      
      alloc_size = SizesElementsOthersNeed(2,irank+1)
      deallocate(jElementsOthersNeed(irank+1)%l)
      alloc_count = alloc_count + alloc_size
   enddo
   call a%Memor%deallocObj(0,'ElementsOthersNeed','Partitionate',alloc_count*ip)
   call a%Memor%dealloc(2,a%MPIsize,SizesElementsOthersNeed,'SizesElementsOthersNeed','Partitionate')
   call a%Memor%dealloc(a%MPIsize,iElementsOthersNeed,'iElementsOthersNeed','Partitionate')
   call a%Memor%dealloc(a%MPIsize,jElementsOthersNeed,'jElementsOthersNeed','Partitionate')
   
  
   !-------------------------------------------------------------------------------------------------------

   !Now build pnods and lnods from what I have received. 
   !I know that no element will be repeated because only one of the naive subdomains sends each element
   
   !First to count
   a%nelem = 0
   lsize = 0
   do irank = 0,a%MPIsize-1
      if (SizesElementsINeed(1,irank+1) /= 0) then
         !Loop through received elements
         nelem = SizesElementsINeed(1,irank+1)
         do ielem = 1,nelem
            ispos = iElementsINeed(irank+1)%l(ielem)-1
            pnode = iElementsINeed(irank+1)%l(ielem+1)-iElementsINeed(irank+1)%l(ielem)
            
            a%nelem = a%nelem+1
            lsize = lsize + pnode
         enddo
      endif
   enddo
   
   !Allocate
   call a%Memor%alloc(a%nelem+1,a%pnods,'pnods','Partitionate')
   call a%Memor%alloc(lsize,a%lnods,'lnods','Partitionate')
   
   
   !Second to fill
   a%nelem = 0
   a%pnods(1) = 1
   do irank = 0,a%MPIsize-1
      if (SizesElementsINeed(1,irank+1) /= 0) then
         !Loop through received elements
         nelem = SizesElementsINeed(1,irank+1)
         do ielem = 1,nelem
            ispos = iElementsINeed(irank+1)%l(ielem)-1
            pnode = iElementsINeed(irank+1)%l(ielem+1)-iElementsINeed(irank+1)%l(ielem)
            
            a%nelem = a%nelem+1
            
            jspos = a%pnods(a%nelem)-1
            a%lnods(jspos+1:jspos+pnode) = jElementsINeed(irank+1)%l(ispos+1:ispos+pnode)
            a%pnods(a%nelem+1) = a%pnods(a%nelem)+pnode
         enddo
      endif
   enddo
   
   call a%Memor%alloc(a%nelem,auxElementRenumbering,'auxElementRenumbering','ScatterData')
   
   if (a%nelty > 1) then
           
      !---------------------------------------------------------------------------------------------
      !Reorder pnods and lnods so that all elements of the same type are together
      call a%Memor%alloc(a%nelem+1,auxpnods,'auxpnods','ScatterData')
      call a%Memor%alloc(size(a%lnods),auxlnods,'auxlnods','ScatterData')
      call a%Memor%alloc(a%nelty,auxnelty,'auxnelty','ScatterData')
      call a%Memor%alloc(a%nelty+1,pstart,'pstart','ScatterData')
      call a%Memor%alloc(a%nelty+1,lstart,'lstart','ScatterData')
      
      
      
      auxnelty = 0
      !count the number of elements of each type
      do ielem = 1,a%nelem
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
         auxnelty(ielty) = auxnelty(ielty)+1
      enddo
      
      pstart(1) = 1
      lstart(1) = 1
      auxpnods(1) = 1
      do ielty = 1,a%nelty
         pstart(ielty+1) = pstart(ielty) + auxnelty(ielty)
         pnode = a%nnode(ielty)
         lstart(ielty+1) = lstart(ielty) + auxnelty(ielty)*pnode
         auxpnods(pstart(ielty+1)) = lstart(ielty+1)
      enddo   
      
      do ielem = 1,a%nelem
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
         elemi = pstart(ielty)
         auxElementRenumbering(ielem) = elemi
         auxpnods(elemi+1) = auxpnods(elemi) + a%nnode(ielty)
         pstart(ielty) = pstart(ielty)+1
         
         auxlnods(lstart(ielty):lstart(ielty)+a%nnode(ielty)-1) = a%lnods(a%pnods(ielem):a%pnods(ielem+1)-1)
         lstart(ielty) = lstart(ielty)+a%nnode(ielty)
      enddo   
      
      do ielem = 1,a%nelem
         a%pnods(ielem) = auxpnods(ielem)
         a%lnods(auxpnods(ielem):auxpnods(ielem+1)-1) = auxlnods(auxpnods(ielem):auxpnods(ielem+1)-1)
      enddo   
      
      call a%Memor%dealloc(a%nelem+1,auxpnods,'auxpnods','ScatterData')
      call a%Memor%dealloc(size(a%lnods),auxlnods,'auxlnods','ScatterData')
      call a%Memor%dealloc(a%nelty,auxnelty,'auxnelty','ScatterData')
      call a%Memor%dealloc(a%nelty+1,pstart,'pstart','ScatterData')
      call a%Memor%dealloc(a%nelty+1,lstart,'lstart','ScatterData')
      
      !---------------------------------------------------------------------------------------------
      
   else
      call a%Memor%dealloc(a%nelem+1,a%pnods,'pnods','Partitionate')
      call a%Memor%alloc(1,a%pnods,'pnods','Partitionate')
      a%pnods(1) = a%nnode(1)
      
      do ielem = 1,a%nelem
         auxElementRenumbering(ielem) = ielem
      enddo
   endif
   
   
   
   !Statistics and Element Postprocessing
   call StatisticsAndElementPostprocessing(a)
   
   
   
   !Before we are done here, I need to build a data structure for storing the initial element numberings and their correspondance with multiple processors
   
   !First let me store the processors an element needs to be sent to
   

   !Secondly let me store the correspondance of an initial element number with the new local element numbering
   call a%Memor%alloc(a%nelem,auxLocal2InitialElementNumbering,'auxLocal2InitialElementNumbering','Partitionate')
   
   !Second to fill
   auxnelem = 0
   do irank = 0,a%MPIsize-1
      if (SizesElementsINeed(1,irank+1) /= 0) then
         !Loop through received elements
         nelem = SizesElementsINeed(1,irank+1)
         do ielem = 1,nelem
            
            auxnelem = auxnelem+1
            
            auxLocal2InitialElementNumbering(auxElementRenumbering(auxnelem)) = gInitialNumberingElementsINeed(irank+1)%l(ielem)
         enddo
      endif
   enddo
   
   call a%Memor%dealloc(a%nelem,auxElementRenumbering,'auxElementRenumbering','ScatterData')
   
   !Deallocate working arrays
   !DeAllocate what I need to receive
   alloc_count = 0
   do irank = 0,a%MPIsize-1
      alloc_size = SizesElementsINeed(1,irank+1)

      deallocate(iElementsINeed(irank+1)%l)
      alloc_count = alloc_count + alloc_size + 1
      
      alloc_size = SizesElementsINeed(2,irank+1)
      deallocate(jElementsINeed(irank+1)%l)
      alloc_count = alloc_count + alloc_size
   enddo
   call a%Memor%dealloc(2,a%MPIsize,SizesElementsINeed,'SizesElementsINeed','Partitionate') !(1): number of elements (2) size of glnods
   call a%Memor%deallocObj(0,'ElementsINeed','Partitionate',alloc_count*ip)
   call a%Memor%dealloc(a%MPIsize,iElementsINeed,'iElementsINeed','Partitionate')
   call a%Memor%dealloc(a%MPIsize,jElementsINeed,'jElementsINeed','Partitionate')
   
   
   
   call a%Memor%alloc(auxnelem,ProcList,'ProcList','Partitionate')
   !Create and Initialize the LocalOrdering
   call a%ParallelLibrary%CreateOrdering(a%ElementLocal2Initial,a%Memor)  !FACTORY
   call a%ElementLocal2Initial%Init(a%MPIcomm,auxnelem,auxLocal2InitialElementNumbering,ProcList,a%Memor)
   call a%Memor%dealloc(auxnelem,ProcList,'ProcList','Partitionate')
   call a%Memor%dealloc(a%nelem,auxLocal2InitialElementNumbering,'auxLocal2InitialElementNumbering','Partitionate')
   
   !Initial numbering of elements needs also to be sent
   do irank = 0,a%MPIsize-1
      call a%Memor%pdealloc(size(gInitialNumberingElementsOthersNeed(irank+1)%l),gInitialNumberingElementsOthersNeed(irank+1)%l,'gInitialNumberingElementsOthersNeed','Partitionate')
      call a%Memor%pdealloc(size(gInitialNumberingElementsINeed(irank+1)%l),gInitialNumberingElementsINeed(irank+1)%l,'gInitialNumberingElementsINeed','Partitionate')
   enddo
   call a%Memor%dealloc(a%MPIsize,gInitialNumberingElementsOthersNeed,'gInitialNumberingElementsOthersNeed','Partitionate')
   call a%Memor%dealloc(a%MPIsize,gInitialNumberingElementsINeed,'gInitialNumberingElementsOthersNeed','Partitionate')
   
   
   call PartitionateStoreInfoForScatteringElements(a)
   
   call a%Timer%Partitionate%Toc
   
   
end subroutine


module Mod_PeriodicBoundaryConditionsData
   use typre
   implicit none
   
   real(rp), parameter :: periodicTol = 1e-5
   
   type :: PeriodicBoundaryConditionsData
      real(rp) :: range(2,3)
      real(rp), allocatable :: SubdomainRange(:,:,:)
      
      !Data initially in my subdomain
      integer(ip), allocatable :: MasterSlave(:,:)
      !Component 1: -1 means master, 0 normal, positive means slave and points to master
      !Component 2: Which is the subdomain of the master
      !list for looping through slaves
      integer(ip) :: nslave
      integer(ip), allocatable :: SlaveList(:)
      
      !Data transferred from SlaveDomains
      !Nodes transferred from Slave Domains  1: gslavepoin, 2: gmasterpoin
      integer(ip) :: SlaveNpoin
      integer(ip), allocatable :: SlaveMasterSlave(:,:)
      !The domain information for SlamveMasterSlave: 1: rank of the gslavepoin, 2: rank of the gmasterpoin
      integer(ip), allocatable :: SlaveMasterSlaveProc(:,:)
      !Elements transferred from Slave Domains
      integer(ip) :: SlaveNelem
      integer(ip), allocatable :: PnodsSlave(:), LnodsSlave(:)
   end type
end module

subroutine PeriodicBCModify(a)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_PeriodicBoundaryConditionsData
   use MPI
   implicit none
   class(FemMesh) :: a
   
   real(rp) :: cputim1, cputim2
   
   interface
      subroutine FindRange(a,PBCD)
         use Mod_Mesh
         use Mod_PeriodicBoundaryConditionsData
         implicit none
         class(FemMesh) :: a
         type(PeriodicBoundaryConditionsData) :: PBCD
      end subroutine
      
      subroutine PeriodicizeMesh(a,PBCD)
         use Mod_Mesh
         use Mod_PeriodicBoundaryConditionsData
         implicit none
         class(FemMesh) :: a
         type(PeriodicBoundaryConditionsData) :: PBCD
      end subroutine

      subroutine FindMasterSlave(a,PBCD)
         use Mod_Mesh
         use Mod_PeriodicBoundaryConditionsData
         implicit none
         class(FemMesh) :: a
         type(PeriodicBoundaryConditionsData) :: PBCD
      end subroutine
      
      subroutine SendElementsFromSlaveToMaster(a,PBCD)
         use Mod_Mesh
         use Mod_PeriodicBoundaryConditionsData
         implicit none
         class(FemMesh) :: a
         type(PeriodicBoundaryConditionsData) :: PBCD
      end subroutine
      
      subroutine PBC_RebuildLocalOrdering(a,PBCD)
         use Mod_Mesh
         use Mod_PeriodicBoundaryConditionsData
         implicit none
         class(FemMesh) :: a
         type(PeriodicBoundaryConditionsData) :: PBCD
      end subroutine
      
      subroutine GlobalMasterSlaveList(a)
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
   end interface 
   
   type(PeriodicBoundaryConditionsData) :: PBCD
   
   integer :: ierr
   
   
   call a%Timer%PeriodicBoundaryConditions%Tic
   
   call cpu_time(cputim1)
   call FindRange(a,PBCD)
   call cpu_time(cputim2)
   
   if (a%kfl_ForcePeriodicMesh == 1) then
      call PeriodicizeMesh(a,PBCD)
   endif
   
   call cpu_time(cputim1)
   call FindMasterSlave(a,PBCD)
   call cpu_time(cputim2)
   
   call cpu_time(cputim1)
   call SendElementsFromSlaveToMaster(a,PBCD)
   call cpu_time(cputim2)
   
   call cpu_time(cputim1)
   call PBC_RebuildLocalOrdering(a,PBCD)
   call cpu_time(cputim2)
   
   call cpu_time(cputim1)
   call GlobalMasterSlaveList(a)
   call cpu_time(cputim2)
   
   call a%Timer%PeriodicBoundaryConditions%Toc
   
end subroutine
   



subroutine FindRange(a,PBCD)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use MPI
   use Mod_PeriodicBoundaryConditionsData
   implicit none
   class(FemMesh) :: a
   type(PeriodicBoundaryConditionsData) :: PBCD
   
   
   real(rp) :: myrange(2,a%ndime)
   integer(ip) :: idime
   
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer :: ierr
   integer(ip) :: irank
   
   real(rp) :: cputim1, cputim2
   
   real(rp), allocatable :: sendSubdomainRange(:,:,:)
   
   !-------------------------------------------------------------------
   !First compute my range
   call a%GetRange(myrange)
   
   !Now send my range to everybody, and receive everybody's range
   !First allocate SubdomainRange 
   call a%Memor%alloc(2,a%ndime,a%MPIsize,PBCD%SubdomainRange,'SubdomainRange','FindRange')
   
   call MPI_AllGather(myrange, 2*a%ndime, MPI_REAL8, PBCD%SubdomainRange, 2*a%ndime, MPI_REAL8, a%MPIcomm,ierr);
   
   !Now get the Global Range
   call a%GetGrange(PBCD%range)
   

end subroutine

module Mod_FindMasterSlaveGiver
   use typre
   use Mod_LinkedList
   use Mod_Octree
   implicit none
   

   type, extends(RangeGiver) :: SubdomainRangeGiver
      real(rp), pointer :: range(:,:,:) => NULL()
      integer(ip) :: ndime
contains
      procedure :: GiveMeRange => GiveRange
      procedure :: Initialize
   
   end type
   
contains
   
   subroutine GiveRange(a,ielem,range)
      implicit none
      class(SubdomainRangeGiver) :: a
      integer(ip) :: ielem
      real(rp) :: range(*)
      
      integer(ip) :: icount, idime, j
      
      icount = 0
      do idime = 1,a%ndime
         do j = 1,2
            icount = icount +1
            range(icount) = a%range(j,idime,ielem)
         enddo
      enddo
   end subroutine
   
   subroutine Initialize(a,nelem,ndime,range)
      implicit none
      class(SubdomainRangeGiver) :: a
      integer(ip) :: ndime
      integer(ip) :: nelem
      real(rp), target :: range(2,ndime,nelem)
      
      a%ndime = ndime
      a%range => range
   end subroutine
   
   
end module

module mod_AuxPeriodicize
   use typre
   use Mod_MeshInterpolator
   implicit none
   
   type, extends(Interpolator)   :: RoundInterpolator

contains
      procedure :: ModifyShafun => RoundModifyShafun
   end type
   
contains
   subroutine RoundModifyShafun(a,newshape)
      use Mod_Element      
class(RoundInterpolator) :: a
      real(rp) :: newshape(:)
      
      !This is put in a subroutine so that it can be extended differently later
      !(used in periodic BC)
      newshape(:) = nint(newshape(:))

   end subroutine
end module
      
subroutine PeriodicizeMesh(a,PBCD)
   use typre
   use Mod_Mesh
   use Mod_PeriodicBoundaryConditionsData
   use mod_AuxPeriodicize
   implicit none
   class(FemMesh) :: a
   type(PeriodicBoundaryConditionsData) :: PBCD
   
   type(RoundInterpolator) :: RIP
  
   real(rp) , allocatable :: coord(:,:)
   real(rp) :: grange(2,3)
   integer(ip) :: idime, ipoin
   
   !first we build an auxiliary coord array
   !in which all nodes are moved on top of the periodic surface (minimum one)
   call a%Memor%alloc(a%ndime,a%npoin,coord,'coord','pm')
   coord = a%coord
   call a%GetGrange(grange)
   do idime = 1,a%ndime
      if (a%PerDime(idime)) then
         coord(idime,:) = grange(1,idime)
      endif
   enddo
   
   !We use the interpoaltor to look for which are the generator nodes
   !The interpolator is rounded so that it does not interpolate
   !but it looks for the closest node (Round Interpolator)
   call RIP%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call RIP%SetTol(a%ForcePeriodicMeshTol)
   call RIP%Initialize_Interp(a,coord)

   
   !We look for the coordinates
   call RIP%Interpolate(a%ndime,a%coord,coord)
   
   

   do idime = 1,a%ndime
      if (a%PerDime(idime) .eqv. .false.) then		
         a%coord(idime,:) = coord(idime,:)
      endif
   enddo
   
   call RIP%Finalize
   
   call a%Memor%dealloc(a%ndime,a%npoin,coord,'coord','pm')
   
end subroutine
   
      
subroutine AreCoordPeriodic(coord,PBCD,ndime,PerDime,IsPeriodic,MasterCoord)
   use typre
   use Mod_Mesh
   use Mod_PeriodicBoundaryConditionsData
   implicit none
   integer(ip) :: ndime
   real(rp)    :: coord(ndime)
   type(PeriodicBoundaryConditionsData) :: PBCD
   logical :: PerDime(ndime)
   logical :: IsPeriodic
   real(rp) :: MasterCoord(ndime)
   
   integer(ip) :: idime
   
   IsPeriodic = .false.
   
   do idime = 1,ndime
      if ((PerDime(idime) .eqv. .true.) .and. abs(coord(idime) - PBCD%range(2,idime)) < PeriodicTol ) then
         IsPeriodic = .true.
         MasterCoord(idime) = PBCD%range(1,idime)
      else
         MasterCoord(idime) = coord(idime)
      endif
      
   enddo
end subroutine
      
subroutine FindMasterSlave(a,PBCD)
   use MPI
   use typre
   use Mod_LinkedList
   use Mod_Octree
   use Mod_Mesh
   use Mod_FindMasterSlaveGiver
   use Mod_PeriodicBoundaryConditionsData
   implicit none
   class(FemMesh) :: a
   type(PeriodicBoundaryConditionsData) :: PBCD
   !At the end of this subroutine I want to know, which of my nodes are slaves,
   !who is their master, and which is the subdomain of the master
   
   type(LinkedListHeader), allocatable :: TestNodesForSubdomainH(:)
   type(LinkedListStorage)             :: TestNodesForSubdomainS 
   
   class(LinkedListHeader), pointer  :: PerHead => NULL()
   class(LinkedListStorage), pointer :: PerStore => NULL()
   
   type(Octree) :: SubdomainOctree
   type(SubdomainRangeGiver) :: MacGyver
   
   type(Octree) :: PointOctree
   
   integer(ip), allocatable :: TestNpoin(:)
   type(i1p), allocatable   :: TestNodes(:)
   type(r2p), allocatable   :: TestCoord(:)
   
   integer(ip), allocatable :: TestNpoinForOthers(:)
   type(r2p), allocatable   :: TestCoordForOthers(:)
   
   type(i1p), allocatable :: TestResultForOthers(:)
   type(i1p), allocatable :: TestResult(:)
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize)
   integer :: ierr
   integer(ip) :: irank
   
   integer(ip) :: icount, idime,ipoin,poini
   logical :: IsFound, IsIn, IsPeriodic
   integer(ip) :: istatus
   
   real(rp) :: MasterCoord(a%ndime)
   
   integer(ip) :: iSlaveCount
   
   real(rp) :: range(2,3)
   
   !Identify which of my nodes are slaves
   call a%Memor%alloc(2,a%npoin,PBCD%MasterSlave,'MasterSlave','FindMasterSlave')
   
   !Allocate the List of Nodes to be Tested by each subdomain
   allocate(TestNodesForSubdomainH(a%MPIsize),stat=istatus)
   call a%Memor%allocObj(istatus,'TestNodesForSubdomainH','FindMasterSlave',a%MPIsize)
   do irank = 0,a%MPIsize-1
      call TestNodesForSubdomainH(irank+1)%Initialize
   enddo
   call TestNodesForSubdomainS%Initialize(a%npoin,'NoRemove','Repeatable',a%Memor)
   
   !Prepare an Octree of Subdomains
   call MacGyver%Initialize(a%MPIsize,a%ndime,PBCD%SubdomainRange)
   call SubdomainOctree%InitializeElementOctree(a%ndime,a%MPIsize,MacGyver,10,a%Memor)
   
   PBCD%nslave = 0
   do ipoin = 1,a%npoin
      call AreCoordPeriodic(a%coord(:,ipoin),PBCD,a%ndime,a%PerDime,IsPeriodic,MasterCoord)
      if (isPeriodic .eqv. .true.) then
         PBCD%nslave = PBCD%nslave+1
         !Get the list of subdomains I need to test
         call SubdomainOctree%GetListForCoord(MasterCoord,PerHead,PerStore)
         call PerStore%GoToFirstOfList(PerHead)
         call PerStore%GetNext(irank)
         do while (irank /= -1)
            irank = irank-1
            
            !Check if I really need to test it
            IsIn = .true.
            do idime = 1,a%ndime
               if (MasterCoord(idime) >= PBCD%SubdomainRange(1,idime,irank+1)-PeriodicTol .and. MasterCoord(idime) <= PBCD%SubdomainRange(2,idime,irank+1)+PeriodicTol ) then  
                  !The node is not out
               else
                  IsIN = .false.
               endif
            enddo
            
            !If it is in, then add the subdomain to the list of check subdomains
            if (IsIn .eqv. .true.) then
               call TestNodesForSubdomainS%Add(TestNodesForSubdomainH(irank+1),ipoin)
            endif
            
            call PerStore%GetNext(irank)
         enddo
         
      endif
   enddo
   
   !SubdomainOctree is no longer needed
   call SubdomainOctree%Dealloc(a%Memor)
   
   !Prepare the list of nodes to be tested, so that we can send it to other subdomains
   call a%Memor%alloc(a%MPIsize,TestNodes,'TestNodes','FindMasterSlave')
   call a%Memor%alloc(a%MPIsize,TestCoord,'TestCoord','FindMasterSlave')
   call a%Memor%alloc(a%MPIsize,TestNpoin,'TestNpoin','FindMasterSlave')
   
   do irank = 0,a%MPIsize-1
      call TestNodesForSubdomainH(irank+1)%GetNelem(TestNpoin(irank+1))
      call a%Memor%palloc(TestNpoin(irank+1),TestNodes(irank+1)%l,'TestNodes','FindMasterSlave')
      call a%Memor%palloc(a%ndime,TestNpoin(irank+1),TestCoord(irank+1)%a,'TestNodes','FindMasterSlave')
   
      call TestNodesForSubdomainS%GoToFirstOfList(TestNodesForSubdomainH(irank+1))
      
      do icount = 1,TestNpoin(irank+1)
         call TestNodesForSubdomainS%GetNext(ipoin)
         call a%Local2Global(ipoin,poini)
         TestNodes(irank+1)%l(icount) = poini
         TestCoord(irank+1)%a(:,icount) = a%coord(:,ipoin)
      enddo
   enddo
   
   !TestNodesForSubdomainList is no longer needed
   call TestNodesForSubdomainS%Dealloc
   deallocate(TestNodesForSubdomainH,stat=istatus)
   call a%Memor%deallocObj(istatus,'TestNodesForSubdomainS','FindMasterSlave',a%MPIsize)
   
   !Now I have what I want to test in:
   !TestNpoinOthersNeed
   !TestNodesOthersNeed
   !TestCoordOthersNeed
   
   call a%Memor%alloc(a%MPIsize,TestNpoinForOthers,'TestNpoinForOthers','FindMasterSlave')
   call a%Memor%alloc(a%MPIsize,TestCoordForOthers,'TestCoordForOthers','FindMasterSlave')
   
   call MPI_AllToAll(TestNpoin, 1, MPI_INTEGER4, TestNpoinForOthers, 1, MPI_INTEGER4, a%MPIcomm,ierr);

   
   !First I receive
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (TestNpoinForOthers(irank+1) /= 0) then
         !call a%Memor%palloc(TestNpoinForOthers(irank+1),TestNodesForOthers(irank+1)%a,'TestNodesForOthers','FindMasterSlave')
         call a%Memor%palloc(a%ndime,TestNpoinForOthers(irank+1),TestCoordForOthers(irank+1)%a,'TestCoordForOthers','FindMasterSlave')
         call MPI_IRECV(TestCoordForOthers(irank+1)%a,TestNpoinForOthers(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !Second I send
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (TestNpoin(irank+1) /= 0) then
         call MPI_ISEND(TestCoord(irank+1)%a,TestNpoin(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   !Now I have what I need to test in TestCoordForOthers, TestNpoinForOthers
   !Allocate the array for the solution
   call a%Memor%alloc(a%MPIsize,TestResultForOthers,'TestResultForOthers','FindMasterSlave')
   do irank = 0,a%MPIsize-1
      call a%Memor%palloc(TestNpoinForOthers(irank+1),TestResultForOthers(irank+1)%l,'TestResultForOthers','FindMasterSlave')
   enddo
   
   !I need to build an octree of points for fast searching
   call PointOctree%InitializePointOctree(a%ndime,a%npoinLocal,a%coord,10_ip,a%Memor)
   
   !Now search each point, for each processor 
   do irank = 0,a%MPIsize-1
      do icount = 1,TestNpoinForOthers(irank+1)
         call AreCoordPeriodic(TestCoordForOthers(irank+1)%a(:,icount),PBCD,a%ndime,a%PerDime,IsPeriodic,MasterCoord)
         do idime = 1,a%ndime
            range(1,idime) = MasterCoord(idime)-PeriodicTol
            range(2,idime) = MasterCoord(idime)+PeriodicTol
         enddo
         call PointOctree%GetListForRange(range,PerHead,PerStore)
         call PerStore%GoToFirstOfList(PerHead)
         call PerStore%GetNext(ipoin)
         listloop: do while (ipoin /= -1)
            IsFound = .true.
            do idime = 1,a%ndime
               if (abs(MasterCoord(idime) - a%coord(idime,ipoin)) > PeriodicTol) IsFound = .false.
            enddo
            if (IsFound .eqv. .true.) then
               call a%Local2Global(ipoin,poini)
               TestResultForOthers(irank+1)%l(icount) = poini
               exit listloop
            endif
            
            call PerStore%GetNext(ipoin)
         enddo listloop
      enddo
   enddo
   
   !Now I do not need the Octree anymore
   call PointOctree%Dealloc(a%Memor)
   
   !Deallocate also many other things
   do irank = 0,a%MPIsize-1
      if (TestNpoinForOthers(irank+1) /= 0) then
         !call a%Memor%palloc(TestNpoinForOthers(irank+1),TestNodesForOthers(irank+1)%a,'TestNodesForOthers','FindMasterSlave')
         call a%Memor%pdealloc(a%ndime,TestNpoinForOthers(irank+1),TestCoordForOthers(irank+1)%a,'TestCoordForOthers','FindMasterSlave')
      endif
   enddo
   
   !Now I receive my results
   irequest2 = MPI_REQUEST_NULL
   call a%Memor%alloc(a%MPIsize,TestResult,'TestResult','FindMasterSlave')
   do irank = 0,a%MPIsize-1
      if (TestNpoin(irank+1) /= 0) then
         call a%Memor%palloc(TestNpoin(irank+1),TestResult(irank+1)%l,'TestResult','FindMasterSlave')
         call MPI_IRECV(TestResult(irank+1)%l,TestNpoin(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !Now I send
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (TestNpoinForOthers(irank+1) /= 0) then
         call MPI_ISEND(TestResultForOthers(irank+1)%l,TestNpoinForOthers(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   !Now fill the MasterSlave from this list
   !Also fill the SlaveList (for looping only through slaves)
   call a%Memor%alloc(PBCD%nslave,PBCD%SlaveList,'SlaveList','FindMasterSlave')
   
   iSlaveCount = 0
   do irank = 0,a%MPIsize-1
      do icount = 1,TestNpoin(irank+1)
         if (TestResult(irank+1)%l(icount) /= 0) then
            ipoin = TestNodes(irank+1)%l(icount)
            call a%Global2Local(ipoin,poini)
            
            !MasterSlave
            PBCD%MasterSlave(1,poini) = TestResult(irank+1)%l(icount)
            PBCD%MasterSlave(2,poini) = irank
            
            if (islaveCount == size(PBCD%SlaveList)) then
               call runend('at least one has been found more than once')
            endif
            
            !SlaveList
            iSlaveCount = iSlaveCount+1
            PBCD%SlaveList(iSlaveCount) = poini
         endif
      enddo
   enddo
   
   !Deallocate lots of things here
   do irank = 0,a%MPIsize-1
      call a%Memor%pdealloc(TestNpoin(irank+1),TestNodes(irank+1)%l,'TestNodes','FindMasterSlave')
      call a%Memor%pdealloc(a%ndime,TestNpoin(irank+1),TestCoord(irank+1)%a,'TestNodes','FindMasterSlave')
   enddo   
   call a%Memor%dealloc(a%MPIsize,TestNodes,'TestNodes','FindMasterSlave')
   call a%Memor%dealloc(a%MPIsize,TestCoord,'TestCoord','FindMasterSlave')
   do irank = 0,a%MPIsize-1
      if (TestNpoin(irank+1) /= 0) then
         call a%Memor%pdealloc(TestNpoin(irank+1),TestResult(irank+1)%l,'TestResult','FindMasterSlave')
      endif
   enddo
   call a%Memor%dealloc(a%MPIsize,TestNpoin,'TestNpoin','FindMasterSlave')
   call a%Memor%dealloc(a%MPIsize,TestResult,'TestResult','FindMasterSlave')
   do irank = 0,a%MPIsize-1
      call a%Memor%pdealloc(TestNpoinForOthers(irank+1),TestResultForOthers(irank+1)%l,'TestResultForOthers','FindMasterSlave')
   enddo
   call a%Memor%dealloc(a%MPIsize,TestResultForOthers,'TestResultForOthers','FindMasterSlave')
   call a%Memor%dealloc(a%MPIsize,TestNpoinForOthers,'TestNpoinForOthers','FindMasterSlave')
   call a%Memor%dealloc(a%MPIsize,TestCoordForOthers,'TestCoordForOthers','FindMasterSlave')
   
   call a%Memor%dealloc(2,a%ndime,a%MPIsize,PBCD%SubdomainRange,'SubdomainRange','FindRange')
   
   !And now we are done
   
end subroutine

subroutine SendElementsFromSlaveToMaster(a,PBCD)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   use Mod_LinkedList
   use Mod_PeriodicBoundaryConditionsData
   use MPI
   implicit none
   class(FemMesh) :: a
   type(PeriodicBoundaryConditionsData) :: PBCD
   
   type(LinkedListHeader), allocatable :: ToSendListH(:)
   type(LinkedListStorage)             :: ToSendListS
   
   integer(ip) :: istatus
   class(FiniteElement), pointer :: e => NULL()
   
   logical, allocatable :: ElementIsDone(:)
   integer(ip), allocatable :: ElementToSubdomain(:), SizesToSend(:),SizesToReceive(:)
   type(i1p), allocatable :: PnodsToSend(:),PnodsToReceive(:),LnodsToSend(:),LnodsToReceive(:)
   
   integer(ip), allocatable :: NelemToReceive(:), NelemToSend(:),NpoinToReceive(:),NpoinToSend(:)
   type(i2p), allocatable :: MasterSlaveToReceive(:), MasterSlaveToSend(:)
   
   integer(ip), allocatable :: glnods(:),IsSentMasterSlave(:),iwa(:)
   logical, allocatable :: lwa(:)
   
   logical :: ElementAlreadyHave, AllNodesAreGhost
   integer(ip) :: elemi,gipoin,gjpoin,icount,ielem,inode,ipoin,ispos,isubdomain,jnode,jpoin,giproc
   integer(ip) :: MaxSlaveNpoin,maxSlaveLnodsSize,maxSlaveNelem,minGlobNumber,minslave,pnode
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize),irequest4(a%MPIsize)
   integer :: ierr
   integer(ip) :: irank
   
   integer(ip) :: jspos,islave
   integer(ip) :: aux_ipoin, aux_lnods(5)
   
   real(rp) :: cputim2, cputim1
   
   allocate(ToSendListH(a%MPIsize),stat = istatus)
   call a%Memor%allocObj(istatus,'ToSendListH','SendElementsFromSlaveToMaster',a%MPIsize)
   do irank = 0,a%MPIsize-1
      call ToSendListH(irank+1)%Initialize
   enddo
   call ToSendListS%Initialize(a%nelem,'Removable','Repeatable',a%Memor)
   
   !For Element Data
   call a%ElementAlloc(e,a%Memor,'DefaultRule','nsm_elmope')
   call a%Memor%alloc(a%nelem,ElementIsDone,'ElementIsDone','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,ElementToSubdomain,'ElementToSubdomain','SendElementsFromSlaveToMaster')
   
   call a%Memor%alloc(a%MPIsize,SizesToSend,'SizesToSend','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,SizesToReceive,'SizesToReceive','SendElementsFromSlaveToMaster')
   
   !Loop through all my slaves
   SlaveLoop: do icount = 1,PBCD%nslave
      ipoin = PBCD%SlaveList(icount) 
      
      if (ipoin == 0) call runend('FindMasterSlave: Master not found. Check your mesh and make sure it is periodic')
      
      call a%Local2Global(ipoin,aux_ipoin)
      
      !Now loop through my pelpo and see if I need to add the element to the master domain
      pelpoLoop: do ispos = a%pelpo(ipoin),a%pelpo(ipoin+1)-1
         ielem = a%lelpo(ispos)
         !Do each element only once
         if (ElementIsDone(ielem) .eqv. .false.) then
            ElementIsDone(ielem) = .true.
            
            call a%ElementLoad(ielem,e)  
            
            !Check if the element needs to be sent
            !(it only needs to be sent if the slave with the lowest global number is mine)
            !This ensures that every element is sent only from one subdomain 
            
            !Now I want to find the slave with the minum Global number and check if it is mine
            minslave = 0
            minGlobNumber = a%gnpoin+1
            
            do jnode = 1,e%pnode
               jpoin = e%lnods(jnode)
               if (PBCD%MasterSlave(1,jpoin) > 0) then
                  call a%Local2Global(jpoin,gjpoin)
                  if (gjpoin < MinGlobNumber) then
                     minslave = jpoin
                     minGlobNumber = gjpoin
                  endif
               endif
            enddo
            
            !Check if the lowest gobally numbered slave is mine, and if so, send
            if (minslave <= a%npoinLocal) then
               !The element needs to be sent
               do jnode = 1,e%pnode
                  jpoin = e%lnods(jnode)
                  if (PBCD%MasterSlave(1,jpoin) > 0) then
                     isubdomain = PBCD%MasterSlave(2,jpoin)
                     !Ensure that the element is only sent to a subdomain once
                     !(Several slaves belonging to the same element can have master's in the same subdomain)
                     if (ElementToSubdomain(isubdomain+1) /= ielem) then
                        ElementToSubdomain(isubdomain+1) = ielem
                        
                        !Also store the size of the list of elements to be sent
                        SizesToSend(isubdomain+1) = SizesToSend (isubdomain+1)+ e%pnode
                        
                        !Add Element to the list of elements to be sent to isubdomain
                        call ToSendListS%Add(ToSendListH(isubdomain+1),ielem)
                     endif
                  endif
               enddo
            endif
         endif
         
      enddo pelpoLoop
   enddo SlaveLoop
      
   !Dealloc unnecessary data
   call a%ElementDealloc(e,a%Memor,'DefaultRule','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%nelem,ElementIsDone,'ElementIsDone','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,ElementToSubdomain,'ElementToSubdomain','SendElementsFromSlaveToMaster')
      
      
   !Prepare the data to be sent to the other subdomains
   call a%Memor%alloc(a%MPIsize,NelemToSend,'NelemToSend','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,NelemToReceive,'NelemToReceive','SendElementsFromSlaveToMaster')
   
   call a%Memor%alloc(a%MPIsize,LnodsToSend,'LnodsToSend','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,PnodsToSend,'PnodsToSend','SendElementsFromSlaveToMaster')
   
   call a%ElementAlloc(e,a%Memor,'DefaultRule','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(e%mnode,gLnods,'gLnods','SendElementsFromSlaveToMaster')
   
   do irank = 0,a%MPIsize-1
      call ToSendListH(irank+1)%GetNelem(NelemToSend(irank+1))
      if (NelemToSend(irank+1) /= 0) then
         call a%Memor%palloc(NelemToSend(irank+1)+1,PnodsToSend(irank+1)%l,'PnodsToSend','SendElementsFromSlaveToMaster')
         call a%Memor%palloc(SizesToSend(irank+1),LnodsToSend(irank+1)%l,'LnodsToSend','SendElementsFromSlaveToMaster')
         
         call ToSendListS%GoToFirstOfList(ToSendListH(irank+1))
         PnodsToSend(irank+1)%l(1) = 1
         do ielem = 1,NelemToSend(irank+1)
            call ToSendListS%GetNext(elemi)
            call a%ElementLoad(elemi,e)
            
            call a%Local2Global(e%pnode,e%lnods,gLnods)
            ispos = PnodsToSend(irank+1)%l(ielem)-1
            LnodsToSend(irank+1)%l(ispos+1:ispos+e%pnode) = gLnods(1:e%pnode)
            PnodsToSend(irank+1)%l(ielem+1) = PnodsToSend(irank+1)%l(ielem)+e%pnode
         enddo
      endif
   enddo
   
   call a%Memor%dealloc(e%mnode,gLnods,'gLnods','SendElementsFromSlaveToMaster')
   call a%ElementDealloc(e,a%Memor,'DefaultRule','nsm_elmope')
   
   call MPI_AllToAll(NelemToSend, 1, MPI_INTEGER4, NelemToReceive, 1, MPI_INTEGER4, a%MPIcomm,ierr);
   call MPI_AllToAll(SizesToSend, 1, MPI_INTEGER4, SizesToReceive, 1, MPI_INTEGER4, a%MPIcomm,ierr);
   
   !Allocate what I need to receive
   call a%Memor%alloc(a%MPIsize,LnodsToReceive,'LnodsToReceive','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,PnodsToReceive,'PnodsToSend','SendElementsFromSlaveToMaster')
   
   do irank = 0,a%MPIsize-1
      if (NelemToReceive(irank+1) /= 0) then
         call a%Memor%palloc(NelemToReceive(irank+1)+1,PnodsToReceive(irank+1)%l,'PnodsToReceive','SendElementsFromSlaveToMaster')
         call a%Memor%palloc(SizesToReceive(irank+1),LnodsToReceive(irank+1)%l,'LnodsToReceive','SendElementsFromSlaveToMaster')
      endif
   enddo
   
   !-----------------------------------------
   !Send Pnods 
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         call MPI_ISEND(PnodsToSend(irank+1)%l,NelemToSend(irank+1)+1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Receive Pnods
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToReceive(irank+1) /= 0) then
         call MPI_IRECV(PnodsToReceive(irank+1)%l,NelemToReceive(irank+1)+1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !-----------------------------------------
   !Send Lnods
   irequest3 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         call MPI_ISEND(LnodsToSend(irank+1)%l,SizesToSend(irank+1),MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest3(irank+1),ierr) 
      endif
   enddo
   
   !Receive Lnods
   irequest4 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToReceive(irank+1) /= 0) then
         call MPI_IRECV(LnodsToReceive(irank+1)%l,SizesToReceive(irank+1)+1,MPI_INTEGER4,irank,mtag2,a%MPIcomm,irequest4(irank+1),ierr) 
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
   
   !I do also need to send a list of the masters of nodes involved in LnodsToSend
   call a%Memor%alloc(a%npoin,IsSentMasterSlave,'IsSentMasterSlave','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%npoin,iwa,'iwa','SendElementsFromSlaveToMaster')
   
   call a%Memor%alloc(a%MPIsize,MasterSlaveToSend,'MasterSlaveToSend','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(a%MPIsize,NpoinToSend,'NpoinToSend','SendElementsFromSlaveToMaster')
      
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         icount = 0
         do ielem = 1,NelemToSend(irank+1)
            ispos = PnodsToSend(irank+1)%l(ielem)-1
            pnode = PnodsToSend(irank+1)%l(ielem+1)-PnodsToSend(irank+1)%l(ielem)
            do inode = 1,pnode
               ipoin = LnodsToSend(irank+1)%l(ispos+inode)
               call a%Global2Local(ipoin,ipoin)
               if (IsSentMasterSlave(ipoin) /= irank+1) then
                  IsSentMasterSlave(ipoin) = irank+1
                  icount = icount + 1
                  iwa(icount) = ipoin
               endif
            enddo
         enddo
         
         npoinToSend(irank+1) = icount
         
         !MasterSlaveToSend: pos1: global node number. pos2: global master number. pos3: master domain
         call a%Memor%palloc(3,NpoinToSend(irank+1),MasterSlaveToSend(irank+1)%l,'MasterSlaveToSend','SendElementsFromSlaveToMaster')
         do icount = 1,npoinToSend(irank+1)
            ipoin = iwa(icount)
            call a%Local2Global(ipoin,gipoin)
            MasterSlaveToSend(irank+1)%l(1,icount) = gipoin
            MasterSlaveToSend(irank+1)%l(2:3,icount) = PBCD%MasterSlave(:,ipoin)
         enddo
      endif
   enddo
   
   call a%Memor%dealloc(a%npoin,iwa,'iwa','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%npoin,IsSentMasterSlave,'IsSentMasterSlave','SendElementsFromSlaveToMaster')
   
   call a%Memor%alloc(a%MPIsize,NpoinToReceive,'NpoinToReceive','SendElementsFromSlaveToMaster')
   
   !Send NpoinToSend
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         call MPI_ISEND(NpoinToSend(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Receive NpoinToReceive
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NelemToReceive(irank+1) /= 0) then
         call MPI_IRECV(NpoinToReceive(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo

   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   !Now allocate and communicate MasterSlaveToReceive
   call a%Memor%alloc(a%MPIsize,MasterSlaveToReceive,'MasterSlaveToReceive','SendElementsFromSlaveToMaster')
   do irank = 0,a%MPIsize-1
      if (NpoinToReceive(irank+1) /= 0) then
         call a%Memor%palloc(3,NpoinToReceive(irank+1),MasterSlaveToReceive(irank+1)%l,'MasterSlaveToReceive','SendElementsFromSlaveToMaster')
      endif   
   enddo
   
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NpoinToSend(irank+1) /= 0) then
         call MPI_ISEND(MasterSlaveToSend(irank+1)%l,3*NpoinToSend(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (NpoinToReceive(irank+1) /= 0) then
         call MPI_IRECV(MasterSlaveToReceive(irank+1)%l,3*NpoinToReceive(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   !Deallocate
   
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         call a%Memor%pdealloc(3,NpoinToSend(irank+1),MasterSlaveToSend(irank+1)%l,'MasterSlaveToSend','SendElementsFromSlaveToMaster')
      endif   
   enddo      
   call a%Memor%dealloc(a%MPIsize,MasterSlaveToSend,'MasterSlaveToSend','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,NpoinToSend,'NpoinToSend','SendElementsFromSlaveToMaster')
   do irank = 0,a%MPIsize-1
      if (NelemToSend(irank+1) /= 0) then
         call a%Memor%pdealloc(NelemToSend(irank+1)+1,PnodsToSend(irank+1)%l,'PnodsToSend','SendElementsFromSlaveToMaster')
         call a%Memor%pdealloc(SizesToSend(irank+1),LnodsToSend(irank+1)%l,'LnodsToSend','SendElementsFromSlaveToMaster')
      endif
   enddo
   call a%Memor%dealloc(a%MPIsize,NelemToSend,'NelemToSend','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,LnodsToSend,'LnodsToSend','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,PnodsToSend,'PnodsToSend','SendElementsFromSlaveToMaster')
   
   deallocate(ToSendListH,stat = istatus)
   call a%Memor%deallocObj(istatus,'ToSendListH','SendElementsFromSlaveToMaster',a%MPIsize)
   call ToSendListS%Dealloc
   
   
   !Now I have all the UNREPEATED, BUT MAYBE ALREADY PRESENT elements and REPEATED points
   !Let us go for the new element list (in global coordinates)
   
   maxSlaveNelem = sum(NelemToReceive)
   maxSlaveLnodsSize = 0
   do irank = 0,a%MPIsize-1
      maxSlaveLnodsSize = maxSlaveLnodsSize + SizesToReceive(irank+1)
   enddo
   
   call a%Memor%alloc(maxSlaveNelem+1,PBCD%PnodsSlave,'PBCD%PnodsSlave','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(maxSlaveLnodsSize,PBCD%LnodsSlave,'PBCD%LnodsSlave','SendElementsFromSlaveToMaster')
   
   PBCD%PnodsSlave(1) = 1
   PBCD%SlaveNelem = 0
   do irank = 0,a%MPIsize-1
      do ielem = 1,NelemToReceive(irank+1)
         ispos = PnodsToReceive(irank+1)%l(ielem)-1
         pnode = PnodsToReceive(irank+1)%l(ielem+1) - PnodsToReceive(irank+1)%l(ielem)
      
         !now check the element and see if I already have
         !Already have the element means all points have a local numbering for me
         ElementAlreadyHave = .true.
         AllNodesAreGhost = .true.
         nodeloop : do inode = 1,pnode
            ipoin = LnodsToReceive(irank+1)%l(ispos+inode)
           
            call a%Global2Local(ipoin,ipoin)
            if (ipoin == 0) then
               ElementAlreadyHave = .false.
               exit nodeloop
            endif
            if (ipoin <= a%npoinLocal) then
               AllNodesAreGhost = .false.
            endif
            
         enddo nodeloop
         
         !Add it if necessary
         if ((ElementAlreadyHave .eqv. .false.) .or. (AllnodesAreGhost .eqv. .true.)) then
            PBCD%SlaveNelem = PBCD%SlaveNelem + 1
            PBCD%PnodsSlave(PBCD%SlaveNelem+1) = PBCD%PnodsSlave(PBCD%SlaveNelem) + pnode
            jspos = PBCD%PnodsSlave(PBCD%SlaveNelem)-1
            PBCD%LnodsSlave(jspos+1:jspos+pnode) = LnodsToReceive(irank+1)%l(ispos+1:ispos+pnode)
         endif
      enddo
   enddo
   
   !Deallocate
   do irank = 0,a%MPIsize-1
      if (NelemToReceive(irank+1) /= 0) then
         call a%Memor%pdealloc(NelemToReceive(irank+1)+1,PnodsToReceive(irank+1)%l,'PnodsToReceive','SendElementsFromSlaveToMaster')
         call a%Memor%pdealloc(SizesToReceive(irank+1),LnodsToReceive(irank+1)%l,'LnodsToReceive','SendElementsFromSlaveToMaster')
      endif
   enddo
   call a%Memor%dealloc(a%MPIsize,NelemToReceive,'NelemToReceive','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,LnodsToReceive,'LnodsToReceive','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,PnodsToReceive,'PnodsToSend','SendElementsFromSlaveToMaster')
   
   
   !Reallocate arrays
   call a%Memor%realloc(PBCD%SlaveNelem+1,PBCD%PnodsSlave,'PBCD%PnodsSlave','SendElementsFromSlaveToMaster')
   call a%Memor%realloc(PBCD%PnodsSlave(PBCD%SlaveNelem+1)-1,PBCD%LnodsSlave,'PBCD%LnodsSlave','SendElementsFromSlaveToMaster')

   !Now I need to prepare a list of MasterSlave for the new nodes
   !new nodes (SlaveNpoin) are those which:
   !  1. Are masters of one of my nodes and did not belong to the domain previously
   !  2. Are connected to a slave in another domain, and the corresponding master is 
   !     in my domain
   
   call a%Memor%alloc(a%gnpoin,lwa,'lwa','SendElementsFromSlaveToMaster')
   do ipoin = 1,a%npoin
      call a%Local2Global(ipoin,gipoin)
      lwa(gipoin) = .true.
   enddo
   
   MaxSlaveNpoin = PBCD%nslave !All the masters of my slave could be new nodes
   do irank = 0,a%MPIsize-1
      MaxSlaveNpoin = MaxSlaveNpoin + 2*NpoinToReceive(irank+1)  !All the points I receive could be new points, and have masters out of my domain
   enddo
   call a%Memor%alloc(2,MaxSlaveNpoin,PBCD%SlaveMasterSlave,'SlaveMasterSlave','SendElementsFromSlaveToMaster')
   call a%Memor%alloc(2,MaxSlaveNpoin,PBCD%SlaveMasterSlaveProc,'SlaveMasterSlaveProc','SendElementsFromSlaveToMaster')
   
   !Initialize
   PBCD%SlaveNpoin = 0
   
   !First I loop my slave nodes and add their masters if necessary
   do ipoin = 1,PBCD%nslave
      islave = PBCD%SlaveList(ipoin)
      gipoin = PBCD%MasterSlave(1,islave)
      giproc = PBCD%MasterSlave(2,islave)
      if (lwa(gipoin) .eqv. .false.) then
         lwa(gipoin) = .true.
         
         PBCD%SlaveNpoin = PBCD%SlaveNpoin+1
         PBCD%SlaveMasterSlave(1,PBCD%SlaveNpoin) = gipoin
         PBCD%SlaveMasterSlaveProc(1,PBCD%SlaveNpoin) = giproc
      endif
   enddo   
      
   !Now I loop the received nodes from other subdomains
   do irank  = 0,a%MPIsize-1
      do ipoin = 1,NpoinToReceive(irank+1)
         !First the slave
         gipoin = MasterSlaveToReceive(irank+1)%l(1,ipoin)
         if (lwa(gipoin) .eqv. .false.) then
            lwa(gipoin) = .true.
            
            PBCD%SlaveNpoin = PBCD%SlaveNpoin+1
            PBCD%SlaveMasterSlave(1:2,PBCD%SlaveNpoin) = MasterSlaveToReceive(irank+1)%l(1:2,ipoin)
            PBCD%SlaveMasterSlaveProc(1,PBCD%SlaveNpoin) = irank
            PBCD%SlaveMasterSlaveProc(2,PBCD%SlaveNpoin) = MasterSlaveToReceive(irank+1)%l(3,ipoin)
         endif
         
         !If it has one master, it needs to be added too!
         gipoin = MasterSlaveToReceive(irank+1)%l(2,ipoin)
         if (gipoin /= 0) then
            if (lwa(gipoin) .eqv. .false.) then
               lwa(gipoin) = .true.
               
               PBCD%SlaveNpoin = PBCD%SlaveNpoin+1
               PBCD%SlaveMasterSlave(1,PBCD%SlaveNpoin) = gipoin
               
               giproc = MasterSlaveToReceive(irank+1)%l(3,ipoin)
               PBCD%SlaveMasterSlaveProc(1,PBCD%SlaveNpoin) = giproc
               
            endif
         endif   
         
      enddo
   enddo
   call a%Memor%realloc(2_ip,PBCD%SlaveNpoin,PBCD%SlaveMasterSlave,'SlaveMasterSlave','SendElementsFromSlaveToMaster')
   call a%Memor%realloc(2_ip,PBCD%SlaveNpoin,PBCD%SlaveMasterSlaveProc,'SlaveMasterSlaveProc','SendElementsFromSlaveToMaster')
   
   !Deallocate
   call a%Memor%dealloc(PBCD%nslave,PBCD%SlaveList,'SlaveList','FindMasterSlave')
   call a%Memor%dealloc(a%gnpoin,lwa,'lwa','SendElementsFromSlaveToMaster')
   do irank = 0,a%MPIsize-1
      if (NpoinToReceive(irank+1) /= 0) then
         call a%Memor%pdealloc(3,NpoinToReceive(irank+1),MasterSlaveToReceive(irank+1)%l,'MasterSlaveToReceive','SendElementsFromSlaveToMaster')
      endif   
   enddo      
   call a%Memor%dealloc(a%MPIsize,MasterSlaveToReceive,'MasterSlaveToReceive','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,NpoinToReceive,'NpoinToReceive','SendElementsFromSlaveToMaster')
   
   call a%Memor%dealloc(a%MPIsize,SizesToSend,'SizesToSend','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(a%MPIsize,SizesToReceive,'SizesToReceive','SendElementsFromSlaveToMaster')
   
end subroutine

subroutine PBC_RebuildLocalOrdering(a,PBCD)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_PeriodicBoundaryConditionsData
   use Mod_BuildGraph
   implicit none
   class(FemMesh) :: a
   type(PeriodicBoundaryConditionsData) :: PBCD
   
   integer(ip) :: new_npoin, new_npoinGhost, new_nelem,new_nslave,new_pnodsSize
   integer(ip), allocatable :: iwa(:),iwa2(:)
   
   integer(ip), allocatable :: aux_pnods(:),aux_lnods(:), aux_pelpo(:),aux_lelpo(:)
   integer(ip) :: aux_nelem
   
   integer(ip) :: icount,sizeLnods,ielem,imaster,ipoin,iSlaveCount
   integer(ip) :: old_lnodssize,sizeLnodsSlave
   
   new_npoin = a%npoin + PBCD%SlaveNpoin
   new_npoinGhost= a%npoinGhost + PBCD%SlaveNpoin
   new_nelem = a%nelem + PBCD%SlaveNelem
   
   
   call a%Memor%alloc(new_npoin,iwa,'iwa','PBC_RebuildLocalOrdering')
   call a%Memor%alloc(new_npoin,iwa2,'iwa2','PBC_RebuildLocalOrdering')
  
   !First the already existing points
   do ipoin = 1,a%npoin
      iwa(ipoin) = ipoin
   enddo
   call a%LocalOrdering%GetProc(a%npoin,iwa,iwa2)
   call a%Local2Global(a%npoin,iwa,iwa)
   
   !now for the new points
   do ipoin = 1,PBCD%SlaveNpoin
      iwa(a%npoin+ipoin) = PBCD%SlaveMasterSlave(1,ipoin)
      iwa2(a%npoin+ipoin) = PBCD%SlaveMasterSlaveProc(1,ipoin)
   enddo
   
   !Now rebuild local ordering
   call a%LocalOrdering%Deallocate(a%Memor)
   
   call a%LocalOrdering%Init(a%MPIcomm,new_npoin,iwa,iwa2,a%Memor)
   
   !Deallcoate
   call a%Memor%dealloc(new_npoin,iwa,'iwa','PBC_RebuildLocalOrdering')
   call a%Memor%dealloc(new_npoin,iwa2,'iwa2','PBC_RebuildLocalOrdering')
   
   !Fom Here On, we start to overwrite existing structure
   
   !Rebuild Lnods accounting for the new elements wich arrived from other subdomains
   !LnodsSlave from global to local
   sizeLnodsSlave = size(PBCD%LnodsSlave)
   call a%Global2Local(sizeLnodsSlave,PBCD%LnodsSlave,PBCD%LnodsSlave)
   
   old_lnodssize = size(a%lnods)
   call a%Memor%realloc(size(a%lnods)+sizeLnodsSlave,a%lnods,'lnods','PBC_RebuildLocalOrdering')
   
   a%lnods(old_lnodssize+1:old_lnodssize+sizeLnodsSlave) = PBCD%LnodsSlave
   if (a%nelty == 1) then
      !pnods is ok
   else
      new_pnodsSize = a%nelem+PBCD%SlaveNelem+1
      call a%Memor%realloc(new_pnodsSize,a%pnods,'pnods','PBC_RebuildLocalOrdering')
      a%pnods(a%nelem+2:new_pnodsSize) = a%pnods(a%nelem+1)+PBCD%PnodsSlave(2:PBCD%SlaveNelem+1)
   endif
   
   !Reallocate Coord (more ghost positions available for the new nodes, coord still
   !                  needs to be ghostCommunicated)
   call a%Memor%realloc(a%ndime,new_npoin,a%coord,'coord','PBC_RebuildLocalOrdering')
      
   !Build local MasterSlave (all masters are now ghost or were already local)
   call a%Memor%alloc(new_npoin,a%MasterSlave,'MasterSlave','PBC_RebuildLocalOrdering')
   !First already existing points
   new_nslave = 0
   do ipoin = 1,a%npoin
      if (PBCD%MasterSlave(1,ipoin) == 0) then
         a%MasterSlave(ipoin) = ipoin
      else
         call a%Global2Local(PBCD%MasterSlave(1,ipoin),imaster)
         a%MasterSlave(ipoin) = imaster
         new_nslave = new_nslave+1
      endif
   enddo
   !Second new points (and count slaves)
   !new_nslave = PBCD%nslave  (I recompute it)
   do ipoin = 1,PBCD%SlaveNpoin
      if (PBCD%SlaveMasterSlave(2,ipoin) == 0) then
         a%MasterSlave(ipoin+a%npoin) = ipoin+a%npoin
      else
         call a%Global2Local(PBCD%SlaveMasterSlave(2,ipoin),imaster)
         a%MasterSlave(ipoin+a%npoin) = imaster
         new_nslave = new_nslave + 1
      endif
   enddo
   
   !Build SlaveList (local numbering)
   a%nslave = new_nslave
   call a%Memor%alloc(a%nslave,a%SlaveList,'SlaveList','PBC_RebuildLocalOrdering')
   iSlaveCount = 0
   do ipoin = 1,new_npoin
      if (a%MasterSlave(ipoin) /= ipoin) then
         iSlaveCount = iSlaveCount + 1
         a%SlaveList(iSlaveCount) = ipoin
      endif
   enddo
   
   !Build New Local Graph which:
      !1. Includes new elements
      !2. Slaves are replaced by their masters
      !3. Includes connections Master-Slave and Slave-Master
   call a%Memor%alloc(new_nelem+new_nslave+1,aux_pnods,'aux_pnods','PBC_RebuildLocalOrdering')
   
   !First Pnods
   if (a%nelty == 1) then
      aux_pnods(1) = 1
      do ielem = 1,new_nelem
         aux_pnods(ielem+1) = aux_pnods(ielem)+a%nnode(1)
      enddo
   else
      aux_pnods(1:new_nelem+1) = a%pnods(1:new_nelem+1)
   endif
   aux_nelem = new_nelem+new_nslave
   do ielem = new_nelem+2,aux_nelem+1
      aux_pnods(ielem) = aux_pnods(ielem-1)+2
   enddo   
  
   !Second Lnods
   call a%Memor%alloc(aux_pnods(new_nelem+1)-1+2*new_nslave,aux_lnods,'aux_lnods','PBC_RebuildLocalOrdering')
   sizeLnods = size(a%lnods)
   aux_lnods(1:sizeLnods) = a%lnods
   !Replace slaves by their master
   do icount = 1,sizeLnods
      aux_lnods(icount) = a%MasterSlave(aux_lnods(icount))
   enddo
   
   !Connect slaves to their master
   do iSlaveCount = 1,a%nslave
      aux_Lnods(sizeLnods+(iSlaveCount-1)*2+1) = a%SlaveList(iSlaveCount)
      aux_Lnods(sizeLnods+(iSlaveCount-1)*2+2) = a%MasterSlave(a%SlaveList(iSlaveCount))
   enddo
   
   !Pelpo and Lelpo
   !call a%Memor%dealloc(size(a%pelpo),a%pelpo,'pelpo','PBC_RebuildLocalOrdering')
   !call a%Memor%dealloc(size(a%lelpo),a%lelpo,'lelpo','PBC_RebuildLocalOrdering')
   call BuildElpo(aux_pelpo,aux_lelpo,new_npoin,aux_nelem,aux_pnods,aux_lnods,2_ip,a%nnode,a%Memor)     !Compute glelpo and gpelpo
   
   !Local graph
   call a%Memor%dealloc(size(a%ia),a%ia,'ia','PBC_RebuildLocalOrdering')
   call a%Memor%dealloc(size(a%ja),a%ja,'ja','PBC_RebuildLocalOrdering')
   call BuildGraph(a%ia,a%ja,aux_pelpo,aux_lelpo,new_npoin,aux_nelem,aux_pnods,aux_lnods,2_ip,a%nnode,a%Memor)
   
   call a%Memor%dealloc(aux_pnods(new_nelem+1)-1+2*new_nslave,aux_lnods,'aux_lnods','PBC_RebuildLocalOrdering')
   call a%Memor%dealloc(new_nelem+new_nslave+1,aux_pnods,'aux_pnods','PBC_RebuildLocalOrdering')
   
   !Dealloc auxiliary arrays
   call a%Memor%dealloc(size(aux_pelpo),aux_pelpo,'pelpo','PBC_RebuildLocalOrdering')
   call a%Memor%dealloc(size(aux_lelpo),aux_lelpo,'lelpo','PBC_RebuildLocalOrdering')
   
   call a%Memor%dealloc(size(PBCD%PnodsSlave),PBCD%PnodsSlave,'PBCD%PnodsSlave','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(size(PBCD%LnodsSlave),PBCD%LnodsSlave,'PBCD%LnodsSlave','SendElementsFromSlaveToMaster')
   
   
   !Build New Array Communicator
   call a%ArrayCommunicator%Deallocate
   
   call a%ArrayCommunicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,new_npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)

   !Deallocate lots of things here
   call a%Memor%dealloc(2,a%npoin,PBCD%MasterSlave,'MasterSlave','FindMasterSlave')
   call a%Memor%dealloc(2_ip,PBCD%SlaveNpoin,PBCD%SlaveMasterSlave,'SlaveMasterSlave','SendElementsFromSlaveToMaster')
   call a%Memor%dealloc(2_ip,PBCD%SlaveNpoin,PBCD%SlaveMasterSlaveProc,'SlaveMasterSlaveProc','SendElementsFromSlaveToMaster')
   
   !Rebuild Pelpo and Lelpo so that they include the new elements
   call a%Memor%dealloc(size(a%pelpo),a%pelpo,'pelpo','PBC_RebuildLocalOrdering')
   call a%Memor%dealloc(size(a%lelpo),a%lelpo,'lelpo','PBC_RebuildLocalOrdering')
   call BuildElpo(a%pelpo,a%lelpo,new_npoin,new_nelem,a%pnods,a%lnods,a%nelty,a%nnode,a%Memor)   
   
   call a%Memor%dealloc(a%npoin,a%IsExnor,'isExnor','PeriodicBC')
   call a%Memor%dealloc(a%npoin,a%ExternalNormal,'ExternalNormal','PeriodicBC')
   
   call a%Memor%alloc(new_npoin,a%IsExnor,'isExnor','PeriodicBC')
   call a%Memor%alloc(new_npoin,a%ExternalNormal,'ExternalNormal','PeriodicBC')
   
   a%npoin = new_npoin
   a%npoinGhost = new_npoinGhost
   a%nelem = new_nelem
   
end subroutine   
   
subroutine GlobalMasterSlaveList(a)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use MPI
   implicit none
   class(FemMesh) :: a
   
   integer(ip), allocatable :: aux_gMasterSlave(:),npoinLocals(:)
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer status(MPI_STATUS_SIZE),irequest1,irequest2
   integer :: ierr
   integer(ip) :: irank
   
   integer(ip) :: ispos
   
   
   if (a%MPIrank == a%MPIroot) then
      call a%Memor%alloc(a%gnpoin,a%gMasterSlave,'gMasterSlave','GlobalMasterSlaveList')
   endif
   call a%Memor%alloc(a%npoinLocal,aux_gMasterSlave,'aux_gMasterSlave','GlobalMasterSlaveList')
   
   call a%Local2Global(a%npoinLocal,a%MasterSlave,aux_gMasterSlave)
   
   !First send sizes (npoinLocal)
   irequest1 = MPI_REQUEST_NULL
   irequest2 = MPI_REQUEST_NULL
   call MPI_ISEND(a%npoinLocal,1,MPI_INTEGER4,a%MPIroot,mtag1,a%MPIcomm,irequest1,ierr) 
   call MPI_ISEND(aux_gMasterSlave,a%npoinLocal,MPI_INTEGER4,a%MPIroot,mtag2,a%MPIcomm,irequest2,ierr) 
   
   if (a%MPIrank == a%MPIroot) then
      call a%Memor%alloc(a%MPIsize,npoinLocals,'npoinLocals','GlobalMasterSlaveList')
      ispos = 1
      do irank = 0,a%MPIsize-1
         call MPI_RECV(npoinLocals(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr)
         call MPI_RECV(a%gMasterSlave(ispos),npoinLocals(irank+1),MPI_INTEGER4,irank,mtag2,a%MPIcomm,status,ierr) 
         ispos = ispos + npoinLocals(irank+1)
      enddo
      call a%Memor%dealloc(a%MPIsize,npoinLocals,'npoinLocals','GlobalMasterSlaveList')
   endif   
      
   call MPI_WAIT(irequest1, status, ierr)   
   call MPI_WAIT(irequest2, status, ierr)
   
   !Deallocate
   call a%Memor%dealloc(a%npoinLocal,aux_gMasterSlave,'aux_gMasterSlave','GlobalMasterSlaveList')
   
   

   
   
   
   
   
end subroutine




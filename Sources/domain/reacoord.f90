module Mod_Reacoord
   use typre
   use Mod_MPIObject
   use Mod_ParallelOrderingInterface
   use Mod_Mesh
   implicit none
   
contains


   subroutine reaCoordSerial(a)
      use MPI
      use typre
      use Mod_Int2Str
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2, mtag3 = 3, mtag4 = 4
      integer            :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize),irequest4(a%MPIsize)
      integer            :: status(MPI_STATUS_SIZE)
      integer(ip) :: irank
      
      integer(ip) :: kfl_ContinueReading,myIostat1
      
      integer(ip) :: ipoin,gipoin,auxipoin,poini,auxpoin2,npoinBlock
      real(rp)    :: coord(a%ndime)
      
      integer(ip), allocatable :: SBuffer_Npoin(:), SBuffer_Lpoin(:,:),SBuffer_Block(:,:)
      real(rp), allocatable    :: SBuffer_Coord(:,:,:)
      
      integer(ip) :: RBuffer_Npoin
      integer(ip), allocatable :: RBuffer_Lpoin(:)
      integer(ip), allocatable :: RBuffer_Block(:)
      real(rp), allocatable :: RBuffer_Coord(:,:)
      
      integer(ip) :: SBBytesPoint,SBBytesProc,SBNpoinProc,SBipoin
      integer(ip) :: kfl_doreceiver
      
      
      
      call a%Memor%alloc(a%ndime,a%npoin,a%coord,'coord','reaCoord')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','reaCoord')
      a%geoblock=0
      
      !Allocate a 5 MBytes buffer for sending and receiving coordinates
      SBBytesPoint = ip+a%ndime*rp
      SBBytesProc = SBBytesPoint*a%MPIsize
      SBNpoinProc = min(5e6/SBBytesProc,1.0_rp*a%gnpoin)
      
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%alloc(a%MPIsize,SBuffer_Npoin,'SBuffer_Npoin','reageo')
         call a%Memor%alloc(SBNpoinProc,a%MPIsize,SBuffer_Lpoin,'SBuffer_Lpoin','reageo')
         call a%Memor%alloc(SBNpoinProc,a%MPIsize,SBuffer_Block,'SBuffer_Block','reageo')
         call a%Memor%alloc(a%ndime,SBNpoinProc,a%MPIsize,SBuffer_Coord,'SBuffer_Coord','reageo')
      endif
      
      call a%Memor%alloc(SBNpoinProc,RBuffer_Lpoin,'RBuffer_Lpoin','reageo')
      call a%Memor%alloc(SBNpoinProc,RBuffer_Block,'RBuffer_Block','reageo')
      call a%Memor%alloc(a%ndime,SBNpoinProc,RBuffer_Coord,'RBuffer_Coord','reageo')
      
      
      if (a%MPIrank == a%MPIroot) then
         call ReachCoordinatesSection(a)
      endif
      
      !Still not done flag
      kfl_ContinueReading = 0
      
      !Number of read elements
      if (a%MPIrank == a%MPIroot) gipoin = 0
      
      irequest1 = MPI_REQUEST_NULL
      irequest2 = MPI_REQUEST_NULL
      irequest3 = MPI_REQUEST_NULL
      irequest4 = MPI_REQUEST_NULL
      
      do while (kfl_ContinueReading /= -1)
         
         !Default is do the receive part
         kfl_doreceiver = 1
      
         !To be done by root
         if (a%MPIrank == a%MPIroot) then
            gipoin = gipoin +1
            
            !If I am root, I do not do the receive part by default
            kfl_doreceiver = 0
            
            !If we are done, first pass, send what remains in the buffers
            if (gipoin == a%gnpoin+1) then
               do irank = 0,a%MPIsize-1
                  if (SBuffer_Npoin(irank+1) > 0 .and. SBuffer_Npoin(irank+1) /= SBNpoinProc) then
                     call reacoord_sendbuffers(a%MPIcomm,irank,a%ndime,a%kfl_ReadBlock,SBuffer_Npoin(irank+1),SBuffer_Lpoin(1,irank+1),SBuffer_Coord(1,1,irank+1),SBuffer_Block(1,irank+1),irequest1(irank+1),irequest2(irank+1),irequest3(irank+1),irequest4(irank+1))
                     !If I sent it to root, then root must receive
                     if (irank == a%MPIroot) kfl_doreceiver = 1
                  endif
               enddo
            
            !If we are done, second pass, we tell everybody we are done
            elseif (gipoin == a%gnpoin+2) then
               do irank = 0,a%MPIsize-1
                  call MPI_WAIT(irequest1(irank+1), status, ierr)
                  call MPI_ISEND(-1_ip,1, MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest1(irank+1), ierr) 
               enddo
               !Root must also receive
               kfl_doreceiver = 1
            
            !Continue sending
            else   
               !First read
               npoinBlock = 0
               call ReadPointCoordinate(a,ipoin,coord,npoinBlock)
               !if(a%kfl_ReadBlock) call ReadPointBlock(a,npoinBlock)
               
               irank = a%gPointProcNumber(ipoin)
               if (irank /= a%gPointProcNumber(ipoin)) then
                  write(*,*) 'ja esta liat altre cop: ','irank: ',irank,'ipoin: ', ipoin,'agPointProcNumber(ipoin): ',a%gPointProcNumber(ipoin)
               endif
               poini = a%gPointGlobNumber(ipoin)
               
               !If we just sent the buffer, we wait for the send to be complete
               if (SBuffer_Npoin(irank+1) == SBNpoinProc) then
                  call MPI_WAIT(irequest1(irank+1), status, ierr)
                  call MPI_WAIT(irequest2(irank+1), status, ierr)
                  call MPI_WAIT(irequest3(irank+1), status, ierr)
                  if(a%kfl_ReadBlock)call MPI_WAIT(irequest4(irank+1), status, ierr)
                  
                  !We start the buffer again
                  SBuffer_Npoin(irank+1) = 0
               endif
               
               !Then place in the buffer
               SBuffer_Npoin(irank+1) = SBuffer_Npoin(irank+1)+1
               
               SBipoin = SBuffer_Npoin(irank+1)
               SBuffer_Lpoin(SBipoin,irank+1) = poini
               SBuffer_Block(SBipoin,irank+1) = npoinBlock
               SBuffer_Coord(:,SBipoin,irank+1) = coord
               
               !If the Buffer is full, we send it
               if (SBuffer_npoin(irank+1) == SBNpoinProc) then
                  call reacoord_sendbuffers(a%MPIcomm,irank,a%ndime,a%kfl_ReadBlock,SBuffer_Npoin(irank+1),SBuffer_Lpoin(1,irank+1),SBuffer_Coord(1,1,irank+1),SBuffer_Block(1,irank+1),irequest1(irank+1),irequest2(irank+1),irequest3(irank+1),irequest4(irank+1))
                  !If I sent it to root, then root must receive
                  if (irank == a%MPIroot) kfl_doreceiver = 1
               endif
            endif
         endif   
         
         !The receive part, to be done by everybody
         if (kfl_doreceiver == 1) then
            !To be done by everybody
            call MPI_RECV(kfl_ContinueReading,1_ip, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,status, ierr) 
            
            if (kfl_continueReading > 0_ip) then
               RBuffer_Npoin = kfl_ContinueReading
               
               call MPI_RECV(RBuffer_Lpoin,RBuffer_Npoin, MPI_INTEGER4, a%MPIroot, mtag2, a%MPIcomm,status, ierr) 
               if(a%kfl_ReadBlock)call MPI_RECV(RBuffer_Block,RBuffer_Npoin, MPI_INTEGER4, a%MPIroot, mtag4, a%MPIcomm,status, ierr) 
               call MPI_RECV(RBuffer_Coord,RBuffer_Npoin*a%ndime, MPI_REAL8, a%MPIroot, mtag3, a%MPIcomm,status, ierr) 
               
               do ipoin = 1,RBuffer_Npoin
                  poini = RBuffer_Lpoin(ipoin)
                  call a%Global2Local(poini,poini)
                  a%coord(:,poini) = RBuffer_Coord(:,ipoin)
                  if(a%kfl_ReadBlock)a%geoblock(poini) = RBuffer_Block(ipoin)
               enddo
            endif
         
         endif
         
      enddo
      
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest1(irank+1), status, ierr)
            call MPI_WAIT(irequest2(irank+1), status, ierr)
            call MPI_WAIT(irequest3(irank+1), status, ierr)
            if(a%kfl_ReadBlock)call MPI_WAIT(irequest4(irank+1), status, ierr)
         enddo
      endif   
      
      !Deallocate the 5 MBytes buffer
      !Allocate a 5 MBytes buffer for sending and receiving coordinates
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%dealloc(a%MPIsize,SBuffer_Npoin,'SBuffer_Npoin','reageo')
         call a%Memor%dealloc(SBNpoinProc,a%MPIsize,SBuffer_Lpoin,'SBuffer_Lpoin','reageo')
         call a%Memor%dealloc(SBNpoinProc,a%MPIsize,SBuffer_Block,'SBuffer_Block','reageo')
         call a%Memor%dealloc(a%ndime,SBNpoinProc,a%MPIsize,SBuffer_Coord,'SBuffer_Coord','reageo')
      endif
      
      call a%Memor%dealloc(SBNpoinProc,RBuffer_Lpoin,'RBuffer_Lpoin','reageo')
      call a%Memor%dealloc(SBNpoinProc,RBuffer_Block,'RBuffer_Block','reageo')
      call a%Memor%dealloc(a%ndime,SBNpoinProc,RBuffer_Coord,'RBuffer_Coord','reageo')
      
   end subroutine

   subroutine reacoord_sendbuffers(MPIcomm,irank,ndime,kfl_ReadBlock,SBuffer_Npoin,SBuffer_Lpoin,SBuffer_Coord,SBuffer_Block,irequest1,irequest2,irequest3,irequest4)
      use MPI
      use typre
      implicit none
      integer(ip)             :: ndime,MPIcomm
      integer(ip), intent(in) :: SBuffer_Npoin,irank
      integer(ip), intent(in) :: SBuffer_Lpoin(*),SBuffer_Block(*)
      real(rp), intent(in)    :: SBuffer_Coord(ndime,*)
      integer(ip), intent(in) :: irequest1,irequest2,irequest3,irequest4
      logical:: kfl_ReadBlock
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3,mtag4 = 4
      integer     :: ierr   
      
      integer(ip) :: lnods_size
      
      call MPI_ISEND(SBuffer_Npoin,1, MPI_INTEGER4, irank, mtag1, MPIcomm,irequest1, ierr)

      if(kfl_ReadBlock)call MPI_ISEND(SBuffer_Block(1),SBuffer_Npoin, MPI_INTEGER4, irank, mtag4, MPIcomm,irequest4, ierr) 
      !Lpoin
      call MPI_ISEND(SBuffer_Lpoin(1),SBuffer_Npoin, MPI_INTEGER4, irank, mtag2, MPIcomm,irequest2, ierr) 
      !Coord
      call MPI_ISEND(SBuffer_Coord(1,1),ndime*SBuffer_Npoin, MPI_REAL8, irank, mtag3, MPIcomm,irequest3, ierr) 
   end subroutine
      
   subroutine reaCoordPartitioned(a)
      use MPI
      use typre
      use Mod_NaivePartition
      use Mod_MPIObject
      use Mod_ParallelOrderingInterface
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      
      type(r2p), allocatable   :: CoordINeed(:), CoordOthersNeed(:)
      type(i1p), allocatable   :: PointsINeed(:), PointsOthersNeed(:)
      
      type(NaivePartition)     :: NaivePartitioner
      integer(ip), allocatable :: countPoint(:)
      integer(ip) :: inipoin,endpoin,gipoin,ipoin,LocalIpoin
      integer(ip) :: OriginalGlobalIpoin,OriginalPoin0,OriginalLocalIpoin
      real(rp) :: coord(a%ndime)
      
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3, mtag4 = 4
      integer, allocatable :: irequest01(:), irequest02(:),irequest03(:),irequest04(:)
      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr,aux
      integer(ip) :: irank
      
      
      !Allocate communications structure
      call a%Memor%alloc(a%MPIsize,PointsOthersNeed,'PointsOthersNeed','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,CoordOthersNeed,'CoordOthersNeed','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,PointsINeed,'PointsINeed','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,CoordINeed,'CoordINeed','reaCoordPartitioned')
      do irank = 0,a%MPIsize-1
         if (a%RP_NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%RP_NpoinOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','reaCoordPartitioned')
            call a%Memor%palloc(a%ndime,a%RP_NpoinOthersNeed(irank+1),CoordOthersNeed(irank+1)%a,'CoordOthersNeed','reaCoordPartitioned')
         endif
         
         if (a%RP_NpoinINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%RP_NpoinINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','reaCoordPartitioned')
            call a%Memor%palloc(a%ndime,a%RP_NpoinINeed(irank+1),CoordINeed(irank+1)%a,'CoordINeed','reaCoordPartitioned')
         endif
      enddo
      
      !Reach the coordinates section
      call ReachCoordinatesSection(a)
      
      !Initialize NaivePartition and get poin0
      call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
      call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
      OriginalPoin0 = inipoin-1
      
      call a%Memor%alloc(a%MPIsize,countPoint,'countPoint','reaCoordPartitioned')
      countPoint = 0
      !Read coordinates and fill the communications structure
      do ipoin = 1,a%RP_npoinLocal
         call ReadPointCoordinate(a,OriginalGlobalIpoin,coord,aux)
         
         !To Local Numbering (NaivePartition)
         OriginalLocalIpoin = OriginalGlobalIpoin - OriginalPoin0
         !To Global Numbering (Final)
         call a%RP_LocalOrdering%Local2Global(OriginalLocalIpoin,gipoin)
         call a%RP_LocalOrdering%GetProc(OriginalLocalIpoin,irank)
         
         countPoint(irank+1) = countPoint(irank+1)+1
         
         CoordOthersNeed(irank+1)%a(:,countPoint(irank+1)) = coord(1:a%ndime)
         PointsOthersNeed(irank+1)%l(countPoint(irank+1)) = gipoin
      enddo
      call a%Memor%dealloc(a%MPIsize,countPoint,'countPoint','reaCoordPartitioned')
      
      
      call a%Memor%alloc(a%MPIsize,irequest01,'irequest01','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,irequest02,'irequest01','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,irequest03,'irequest01','reaCoordPartitioned')
      call a%Memor%alloc(a%MPIsize,irequest04,'irequest01','reaCoordPartitioned')
      
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      irequest03 = MPI_REQUEST_NULL
      irequest04 = MPI_REQUEST_NULL
      !Communicate the coordinates
      do irank = 0,a%MPIsize-1
         if (a%RP_NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(PointsOthersNeed(irank+1)%l, a%RP_NpoinOthersNeed(irank+1), MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
            call MPI_ISEND(CoordOthersNeed(irank+1)%a, a%ndime*a%RP_NpoinOthersNeed(irank+1), MPI_REAL8, irank, mtag2, a%MPIcomm,irequest02(irank+1), ierr)
         endif
         if (a%RP_NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(PointsINeed(irank+1)%l, a%RP_NpoinINeed(irank+1), MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest03(irank+1), ierr)
            call MPI_IRECV(CoordINeed(irank+1)%a, a%ndime*a%RP_NpoinINeed(irank+1), MPI_REAL8, irank, mtag2, a%MPIcomm,irequest04(irank+1), ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest02, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest03, MPI_STATUSES_IGNORE, ierr)
      call MPI_WAITALL(a%MPIsize,irequest04, MPI_STATUSES_IGNORE, ierr)
         
      !Copy from buffers to coordinates
      call a%Memor%alloc(a%ndime,a%npoin,a%coord,'coord','reaCoord')
      call a%Memor%alloc(a%npoin,a%geoblock,'geoblock','reaCoord')
      a%geoblock=0
      
      
      do irank = 0,a%MPIsize-1
         if (a%RP_NpoinINeed(irank+1) /= 0) then
            do ipoin = 1,a%RP_NpoinINeed(irank+1)
               gipoin = PointsINeed(irank+1)%l(ipoin)
               call a%LocalOrdering%Global2Local(gipoin,LocalIpoin)
               a%coord(:,LocalIpoin) = CoordINeed(irank+1)%a(:,ipoin)
            enddo
         endif
      enddo
               
      !Deallocate buffers
      !Allocate communications structure
      do irank = 0,a%MPIsize-1
         if (a%RP_NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%RP_NpoinOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','reaCoordPartitioned')
            call a%Memor%pdealloc(a%ndime,a%RP_NpoinOthersNeed(irank+1),CoordOthersNeed(irank+1)%a,'CoordOthersNeed','reaCoordPartitioned')
         endif
         
         if (a%RP_NpoinINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%RP_NpoinINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','reaCoordPartitioned')
            call a%Memor%pdealloc(a%ndime,a%RP_NpoinINeed(irank+1),CoordINeed(irank+1)%a,'CoordINeed','reaCoordPartitioned')
         endif
      enddo
      call a%Memor%dealloc(a%MPIsize,PointsOthersNeed,'PointsOthersNeed','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,CoordOthersNeed,'CoordOthersNeed','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,PointsINeed,'PointsINeed','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,CoordINeed,'CoordINeed','reaCoordPartitioned')
      
      call a%Memor%dealloc(a%MPIsize,irequest01,'irequest01','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,irequest02,'irequest01','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,irequest03,'irequest01','reaCoordPartitioned')
      call a%Memor%dealloc(a%MPIsize,irequest04,'irequest01','reaCoordPartitioned')
      
   end subroutine
   
   subroutine ReachCoordinatesSection(a)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      
      call a%Listener%listen('REAGEO')
      do while(a%Listener%words(1)/='COORD')
         call a%Listener%listen('REAGEO')
      end do
   end subroutine
   
   subroutine ReadPointCoordinate(a,ipoin,coord,npoinBlock)
      use Mod_Int2Str
      use Mod_MpiObject
      use Mod_Mesh
      
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ipoin,npoinBlock
      real(rp) :: coord(a%ndime)
      
      integer(ip) :: myIOStat1
      
      if(a%kfl_ReadBlock) then
          read(a%Listener%nunit,*,IOStat=myIOStat1) ipoin, coord(1:a%ndime),npoinBlock
      else
          read(a%Listener%nunit,*,IOStat=myIOStat1) ipoin, coord(1:a%ndime)
      endif
      if (myIOStat1/=0) call runend('Reageo: Error reading coordinates in point '//adjustl(trim(int2str(ipoin))))
      
      if(ipoin<0.or.ipoin>a%gnpoin) &
                     call runend('REAGEO: WRONG NUMBER OF NODES') 
   end subroutine   

   
end module

subroutine reaCoord(a)
   use typre
   use Mod_Reacoord
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   call a%Timer%reaCoord%Tic
   
   if (a%kfl_ReadType == 0) then
      call reaCoordSerial(a)
      
   elseif (a%kfl_ReadType == 1) then
      call reaCoordPartitioned(a)
      
   endif
   
   call a%Timer%reaCoord%Toc
   
end subroutine


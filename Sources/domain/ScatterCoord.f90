subroutine ScatterCoord(a)
   use MPI
   use typre
   use Mod_Mesh
   use Mod_NaivePartition
   implicit none
   class(FemMesh) :: a
   
   type(NaivePartition) :: NaivePartitioner
   integer(ip) :: inipoin,endpoin,ipoin,TempNpoinLocal,poini
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer            :: ierr
   integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize)
   integer(ip) :: irank
   
   integer(ip), allocatable :: SizesOthersNeed(:)
   type(i1p), allocatable :: PointsOthersNeed(:)
   type(r2p), allocatable    :: CoordOthersNeed(:)
   
   integer(ip), allocatable :: SizesINeed(:)
   type(i1p), allocatable :: PointsINeed(:)
   type(r2p), allocatable    :: CoordINeed(:)
   
   integer(ip) :: istat
   
   
   !Initialize NaivePartitioner
   !THIS IS ONLY USED TO KNOW HOW MANY NODES I HAVE (TempNpoinLocal)
   call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
   call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
   TempNpoinLocal = endpoin-inipoin+1
   
   call a%Memor%alloc(a%MPIsize,SizesOthersNeed,'SizesOthersNeed','reaCoord')
   call a%Memor%alloc(a%MPIsize,SizesINeed,'SizesOthersNeed','reaCoord')
   
   SizesOthersNeed = 0
   !Count
   do ipoin = 1,TempNpoinLocal
      irank = a%pointProcNumber(ipoin)
      SizesOthersNeed(irank+1) = SizesOthersNeed(irank+1)+1
   enddo    
   
   !Send how many nodes I need from others, receive how many nodes others need from me
   call MPI_AllToAll(SizesOthersNeed, 1, MPI_INTEGER4, SizesINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr);
   
   !Others need
   allocate(PointsOthersNeed(a%MPIsize),stat=istat)
   call a%Memor%allocObj(istat,'PointsOthersNeed','ScatterCoord',a%MPIsize*ip)
   allocate(CoordOthersNeed(a%MPIsize),stat=istat)
   call a%Memor%allocObj(istat,'CoordOthersNeed','ScatterCoord',a%MPIsize*ip)
   do irank = 0,a%MPIsize-1
      call a%Memor%palloc(SizesOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','ScatterCoord')
      call a%Memor%palloc(a%ndime,SizesOthersNeed(irank+1),CoordOthersNeed(irank+1)%a,'CoordOthersNeed','ScatterCoord')
   enddo   
   
   !I need
   allocate(PointsINeed(a%MPIsize),stat=istat)
   call a%Memor%allocObj(istat,'PointsINeed','ScatterCoord',a%MPIsize*ip)
   allocate(CoordINeed(a%MPIsize),stat=istat)
   call a%Memor%allocObj(istat,'CoordINeed','ScatterCoord',a%MPIsize*ip)
   do irank = 0,a%MPIsize-1
      call a%Memor%palloc(SizesINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','ScatterCoord')
      call a%Memor%palloc(a%ndime,SizesINeed(irank+1),CoordINeed(irank+1)%a,'CoordINeed','ScatterCoord')
   enddo   
   
   SizesOthersNeed = 0
   !Fill others need
   do ipoin = 1,TempNpoinLocal
      irank = a%pointProcNumber(ipoin)
      SizesOthersNeed(irank+1) = SizesOthersNeed(irank+1)+1
      PointsOthersNeed(irank+1)%l(SizesOthersNeed(irank+1)) = a%pointGlobNumber(ipoin)
      CoordOthersNeed(irank+1)%a(1:a%ndime,SizesOthersNeed(irank+1)) = a%gcoord(1:a%ndime,ipoin)
   enddo
   
   !Send PointsOthersNeed
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesOthersNeed(irank+1) /= 0) then
         call MPI_ISEND(PointsOthersNeed(irank+1)%l,SizesOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Receive the Points others need
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesINeed(irank+1) /= 0) then
         call MPI_IRECV(PointsINeed(irank+1)%l,SizesINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !-----------------------------------------------------------------------------------------------------------------------------
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   
   !Now the coordinates
   !Send CoordOthersNeed
   irequest1 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesOthersNeed(irank+1) /= 0) then
         call MPI_ISEND(CoordOthersNeed(irank+1)%a,SizesOthersNeed(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr) 
      endif
   enddo
   
   !Receive the Coord I need
   irequest2 = MPI_REQUEST_NULL
   do irank = 0,a%MPIsize-1
      if (SizesINeed(irank+1) /= 0) then
         call MPI_IRECV(CoordINeed(irank+1)%a,SizesINeed(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
      endif
   enddo
   
   !-----------------------------------------------------------------------------------------------------------------------------
   !Make sure everything is sent and received
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   
   !Allocate coord
   call a%Memor%alloc(a%ndime,a%npoin,a%coord,'coord','reaCoord')
   
   !Now move the data to coord
   do irank = 0,a%MPIsize-1
      do ipoin = 1,SizesINeed(irank+1)
         poini = PointsINeed(irank+1)%l(ipoin)
         call a%Global2Local(poini,poini)
         a%coord(:,poini) = CoordINeed(irank+1)%a(:,ipoin)
      enddo
   enddo
   
   !Now I can deallocate
   
   !Others need
   do irank = 0,a%MPIsize-1
      call a%Memor%pdealloc(SizesOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','ScatterCoord')
      call a%Memor%pdealloc(a%ndime,SizesOthersNeed(irank+1),CoordOthersNeed(irank+1)%a,'CoordOthersNeed','ScatterCoord')
   enddo   
   deallocate(PointsOthersNeed,stat=istat)
   call a%Memor%deallocObj(istat,'PointsOthersNeed','ScatterCoord',a%MPIsize*ip)
   deallocate(CoordOthersNeed,stat=istat)
   call a%Memor%deallocObj(istat,'CoordOthersNeed','ScatterCoord',a%MPIsize*ip)
   
   
   !I need
   do irank = 0,a%MPIsize-1
      call a%Memor%pdealloc(SizesINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','ScatterCoord')
      call a%Memor%pdealloc(a%ndime,SizesINeed(irank+1),CoordINeed(irank+1)%a,'CoordINeed','ScatterCoord')
   enddo  
   deallocate(PointsINeed,stat=istat)
   call a%Memor%deallocObj(istat,'PointsINeed','ScatterCoord',a%MPIsize*ip)
   deallocate(CoordINeed,stat=istat)
   call a%Memor%deallocObj(istat,'CoordINeed','ScatterCoord',a%MPIsize*ip)
    
   
   call a%Memor%dealloc(a%MPIsize,SizesOthersNeed,'SizesOthersNeed','reaCoord')
   call a%Memor%dealloc(a%MPIsize,SizesINeed,'SizesOthersNeed','reaCoord')
   
   
   !Deallocate gcoord (no longer needed)
   call a%Memor%dealloc(a%ndime,TempNpoinLocal,a%gcoord,'gcoord','ScatterCoord')
      
   call a%Memor%dealloc(size(a%pointProcNumber),a%pointProcNumber,'pointProcNumber','ScatterCoord')
   call a%Memor%dealloc(size(a%pointGlobNumber),a%pointGlobNumber,'pointGlobNumber','ScatterCoord')
      
end subroutine

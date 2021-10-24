subroutine TestCommunications(a)
   use MPI
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip) :: i, ierr
   
   real(rp), allocatable :: AlltoAll1(:), AllToAll2(:)
   real(rp) :: rwa1, rwa2
   integer(ip) :: npoinGhost, npoinGhostGlobal
   
   integer(ip), allocatable :: Gbyte1(:), GByte2(:)
   
   integer, parameter :: mtag1 = 1, mtag2 = 2, mtag3 = 3
   integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize),irequest2(a%MPIsize)
   
   integer(ip) :: irank
   
   real(rp) :: time1, time2
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
   
   !--------------------------------------------------------------------
   !Test the communication of ghost points, and time it
   call a%Timer%TestGhostCommunicate%Tic
   do i = 1,100
      call a%ArrayCommunicator%GhostCommunicate(a%ndime,a%coord)
   enddo
   call MPI_BARRIER(a%MPIcomm,ierr)  
   call a%Timer%TestGhostCommunicate%Toc
   !--------------------------------------------------------------------
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
      
   !Test how long does it take to communicate 1 Mbit between all processors
   call a%Timer%TestBandwidth%Tic
   
   a%BandwidthTestSize = 1e7_rp/4_rp
   call a%Memor%alloc(a%BandwidthTestSize,Gbyte1,'gbyte1','TestCommunications')
   call a%Memor%alloc(a%BandwidthTestSize,Gbyte2,'gbyte2','TestCommunications')
   
   do irank = 0,a%MPIsize-1
      call MPI_ISEND(Gbyte1, a%BandwidthTestSize, MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest1(irank+1), ierr)
   enddo
      
   do irank = 0,a%MPIsize-1
      call MPI_IRECV(Gbyte2,a%BandwidthTestSize,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
   enddo
   
   do irank = 0,a%MPIsize-1
      call MPI_WAIT(irequest1(irank+1), status, ierr)
      call MPI_WAIT(irequest2(irank+1), status, ierr)
   enddo
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
   call a%Timer%TestBandwidth%Toc
   
   call a%Timer%TestBandwidth%GetValue(time1)
   a%Bandwidth = (4.0_rp*8.0_rp*a%BandwidthTestSize) * a%MPIsize * a%MPIsize /time1 /1e9
   
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
   
   !--------------------------------------------------------------------
   !Test the communications of an AllToAll array, and time it
   call a%Memor%alloc(a%MPIsize,AllToAll1,'AllToAll1','TestCommunications')
   call a%Memor%alloc(a%MPIsize,AllToAll2,'AllToAll2','TestCommunications')
   
   call a%Timer%TestAllToAll%Tic
   do i = 1,100
      call MPI_AllToAll(AllToAll1, 1, MPI_REAL8, AllToAll2, 1, MPI_REAL8, a%MPIcomm,ierr);
   enddo
   call MPI_BARRIER(a%MPIcomm,ierr)   
   call a%Timer%TestAllToAll%Toc
   
   call a%Memor%dealloc(a%MPIsize,AllToAll1,'AllToAll1','TestCommunications')
   call a%Memor%dealloc(a%MPIsize,AllToAll2,'AllToAll2','TestCommunications')
   !--------------------------------------------------------------------
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
   
   !--------------------------------------------------------------------
   !Test the communications of MPI_REDUCE, and time it
   rwa1 = a%MPIrank
   call a%Timer%TestReduce%Tic
    do i = 1,100
      call  MPI_REDUCE( rwa1, rwa2, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   enddo
   call MPI_BARRIER(a%MPIcomm,ierr)   
   call a%Timer%TestReduce%Toc
   !--------------------------------------------------------------------
   
   call MPI_BARRIER(a%MPIcomm,ierr)   
   
   !--------------------------------------------------------------------
   !Test the communications of MPI_BCAST, and time it

   call a%Timer%TestBcast%Tic
    do i = 1,100
      CALL MPI_BCAST(rwa1, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   enddo
   call MPI_BARRIER(a%MPIcomm,ierr)  
   call a%Timer%TestBcast%Toc
   
end subroutine

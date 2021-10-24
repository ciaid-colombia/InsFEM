module Mod_GatherScatterDtcri
contains
subroutine php_GatherScatterDtcri(a)
   !This subroutine gathers dtcri from all processes, finds the minimum dtcri, and scatters 
   !the minimum value to all processes
   use MPI
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: irank
   real(rp)    :: aux_dtinv(a%MPIsize)
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer status(MPI_STATUS_SIZE)
   integer :: ierr,irequest,irequest01(a%MPIsize)
   
   irequest = MPI_REQUEST_NULL
   irequest01 = MPI_REQUEST_NULL
   
   call MPI_ISEND(a%dtcri, 1, MPI_REAL8, a%MPIroot, mtag1, a%MPIcomm,irequest, ierr)
   if (a%MPIrank == a%MPIroot) then
      !Blocking
      if (a%kfl_MPIComType == 0) then
         do irank = 0,a%MPIsize-1
            call MPI_RECV(aux_dtinv(irank+1), 1, MPI_REAL8, irank, mtag1, a%MPIcomm, status, ierr)
         enddo
      !Non-Blocking
      else
         do irank = 0,a%MPIsize-1
            call MPI_IRECV(aux_dtinv(irank+1), 1, MPI_REAL8, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr)
         enddo
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest01(irank+1), status,ierr)
         enddo
      endif
      
      a%dtcri = minval(aux_dtinv)
   endif  
   call MPI_WAIT(irequest, status, ierr)
   CALL MPI_BCAST(a%dtcri, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
end subroutine
end module
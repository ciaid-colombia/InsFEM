subroutine InitializeMPI
   use def_master
   use MPI
   implicit none
   
   integer(ip) :: ierr
   real(rp) :: start, endt
   
   call cpu_start(4)%Tic2
   call MPI_INIT( ierr )
   MPIcomm = MPI_COMM_WORLD
   if (ierr /= 0) call runend('Error in MPI_INIT, InitializeMPI')
   call MPI_COMM_RANK( MPIcomm, MPIrank, ierr )
   call MPI_COMM_SIZE( MPIcomm, MPIsize, ierr )

   call MPI_BARRIER(MPIcomm,ierr)   
   call cpu_start(4)%Toc2
   
end subroutine

subroutine FinalizeMPI
   use def_master
   use MPI

   implicit none
   
   integer(ip) :: ierr
   call MPI_FINALIZE(ierr)
   
end subroutine   
   
   

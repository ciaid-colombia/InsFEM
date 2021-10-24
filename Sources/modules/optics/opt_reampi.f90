subroutine opt_reampi(a)
   use typre
   use MPI
   use Mod_Optics
   implicit none
   
   class(OpticsProblem) :: a
   integer(ip) :: ierr,istat,ibeam
   
   !Communicate opt_reaphy
   CALL MPI_BCAST(a%kfl_velocity, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_pressure, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_temperature, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)   
   
   CALL MPI_BCAST(a%densi, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   CALL MPI_BCAST(a%sphea, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   CALL MPI_BCAST(a%tcond, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   CALL MPI_BCAST(a%visco, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   CALL MPI_BCAST(a%prtur, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   
   !Communicate opt_reanut
   CALL MPI_BCAST(a%kfl_dissi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%turbu,size(a%turbu),MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   CALL MPI_BCAST(a%a_obukhov, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   
   !Communicate opt_reaous
   call MPI_BCAST(a%nbeams,1,MPI_INTEGER4,a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%kfl_Avg1DCn2,1,MPI_INTEGER4,a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%Avg1DIdime,1,MPI_INTEGER4,a%MPIroot, a%MPIcomm, ierr)
   
   call MPI_BCAST(a%lambda,1,MPI_REAL8,a%MPIroot, a%MPIcomm, ierr)
   
   if (a%MPIrank /= a%MPIroot) then
      allocate(a%beams(a%nbeams),STAT=istat)
      call a%Memor%allocObj(istat,'beams','opt_reampi',6_ip*rp*a%nbeams)
   endif
   
   do ibeam = 1,a%nbeams
      call MPI_BCAST(a%beams(ibeam)%origin,size(a%beams(ibeam)%origin),MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      call MPI_BCAST(a%beams(ibeam)%direct,size(a%beams(ibeam)%direct),MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      call MPI_BCAST(a%beams(ibeam)%length,1,MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   enddo
   
   

end subroutine
subroutine php_Initializations(a)
   use MPI
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: ndime, nptra,ierr

   call a%Mesh%GetNdime(ndime)
   if(a%nptra>0) then
       nptra = 0
       CALL MPI_BCAST(nptra, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
       if (a%MPIrank == a%MPIroot) nptra = a%nptra

       if (a%MPIrank .ne. a%MPIroot) then
           call a%Memor%alloc(ndime,a%nptra,a%cptra,'cptra','php_reaous')
       endif

       call a%TrackingInterpolator%Initialize_Interp(a%Mesh,a%cptra(:,1:nptra))

   endif
   
end subroutine 

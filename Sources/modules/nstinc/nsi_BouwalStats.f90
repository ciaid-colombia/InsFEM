module Mod_nsi_BouwalStats
   use MPI
   use typre
   use Mod_NavierStokes
   implicit none

contains   
   subroutine nsi_BouwalStats(a,itask)
      class(NavierStokesProblem) :: a
      integer(ip) :: itask
      
      real(rp) :: ypmean, ypmax, ypmin
      real(rp) :: tamean, tamax, tamin
      integer :: ierr
      
      !Begste
      if (itask == 0) then
         a%ypmean = 0
         a%nmean_walllaw = 0
         a%ypmax = -1e12
         a%ypmin = 1e12
      !Endste   
      elseif (itask == 1) then
      
         if (a%nmean_walllaw /= 0) a%ypmean = a%ypmean/a%nmean_walllaw
         call MPI_REDUCE( a%ypmean, ypmean, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
         call MPI_REDUCE( a%ypmin, ypmin, 1, MPI_REAL8, MPI_MIN, a%MPIroot,a%MPIcomm, ierr )
         call MPI_REDUCE( a%ypmax, ypmax, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
         
         call MPI_REDUCE( a%tamea, tamean, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
         call MPI_REDUCE( a%tamin, tamin, 1, MPI_REAL8, MPI_MIN, a%MPIroot,a%MPIcomm, ierr )
         call MPI_REDUCE( a%tamax, tamax, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

         if (a%MPIrank == a%MPIroot) then
            a%ypmean = ypmean/a%MPIsize
            a%ypmax = ypmax
            a%ypmin = ypmin
            a%tamea = tamean/a%MPIsize
            a%tamax = tamax
            a%tamin = tamin
         endif
      endif
   end subroutine
 
end module

subroutine ComputeRanges(a)
   use MPI
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   real(rp) :: maxrange(a%ndime),minrange(a%ndime),maxgrange(a%ndime),mingrange(a%ndime)
   
   integer(ip) :: idime
   integer     :: ierr

   !First the local range
   do idime = 1,a%ndime
      a%range(1,idime) = minval(a%coord(idime,:))
      a%range(2,idime) = maxval(a%coord(idime,:))
   enddo

   minrange = a%range(1,1:a%ndime)
   maxrange = a%range(2,1:a%ndime)
   
   !Now the global one
   call  MPI_AllREDUCE( minrange, mingrange, a%ndime, MPI_REAL8, MPI_MIN,a%MPIcomm, ierr )
   call  MPI_AllREDUCE( maxrange, maxgrange, a%ndime, MPI_REAL8, MPI_MAX,a%MPIcomm, ierr )
   
   a%grange(1,1:a%ndime) = mingrange
   a%grange(2,1:a%ndime) = maxgrange

end subroutine

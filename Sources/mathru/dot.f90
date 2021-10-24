subroutine dot(v1,v2,n,vdot)
   use typre
   implicit none
   integer(ip) :: n
   real(rp) :: v1(n),v2(n),vdot

   integer(ip) :: idime
   
   vdot = 0.0_rp
   do idime = 1,n
      vdot = vdot + v1(idime)*v2(idime)
   enddo
end subroutine

subroutine MPIdot(v1,v2,n,vdot,mpicomm,mpiroot)
   use typre
   use MPI  
   implicit none
   integer(ip) :: n, mpicomm,mpiroot
   real(rp) :: v1(n), v2(n), vdot

   real(rp) :: rwa, rwa2
   
   integer                  :: ierr
   
   rwa =  dot_product(v1(1:n),v2(1:n))
   
   call  MPI_REDUCE( rwa, vdot, 1, MPI_REAL8, MPI_SUM, mpiroot,MPIcomm, ierr )
   
end subroutine
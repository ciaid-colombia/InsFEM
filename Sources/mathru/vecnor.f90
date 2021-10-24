subroutine vecnor(v,n,vnor,inor)

!-----------------------------------------------------------------------
!
! Compute the L1, L2 and L-inf norms of a vector V of length N
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n,inor
  real(rp),    intent(in)  :: v(n)
  real(rp),    intent(out) :: vnor
  integer(ip)              :: i

  vnor=0.0_rp
  if(inor==0.or.inor==3) then
     do i=1,n
        if(abs(v(i))>vnor) vnor=abs(v(i))
     end do
  else if(inor==1) then
     do i=1,n
        vnor=vnor+abs(v(i))
     end do
  else if(inor==2) then
     do i=1,n
        vnor=vnor+v(i)*v(i)
     end do
     vnor=sqrt(vnor)
  end if
  
end subroutine vecnor

subroutine vecnorALLMPI(v,ndime,n,vnor,inor,MPIComm,MPIrank,MPIroot,MPIsize) 

!-----------------------------------------------------------------------
!
! Compute the L1, L2 and L-inf norms of a vector V of length N
!
!-----------------------------------------------------------------------
!Returns correct result only for rank 0!!
  use      typre
  use MPI
  implicit none
  integer(ip), intent(in)  :: n,inor,ndime
  real(rp),    intent(in)  :: v(ndime,n)
  real(rp),    intent(out) :: vnor
  integer(ip) :: i
  integer(ip) :: MPIComm,MPIrank,MPIroot,MPIsize
  integer(ip) :: ierr,idime
  real(rp)    :: vnor0,va

  vnor0=0.0_rp
  if(inor==0.or.inor==3) then
     do i = 1,n
        do idime = 1,ndime
         va = v(idime,i)
         if(abs(va)>vnor0) vnor0=abs(va)
        enddo
     end do

     call  MPI_ALLREDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_MAX, MPIroot,MPIcomm, ierr )

  else if(inor==1) then
     do i = 1,n
        do idime = 1,ndime
         va = v(idime,i)
         vnor0=vnor0+abs(va)
        enddo
     end do

     call  MPI_ALLREDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )

  else if(inor==2) then

     do i = 1,n
        do idime = 1,ndime
         va = v(idime,i)
         vnor0 = vnor0 + va*va
        enddo
     end do

     call  MPI_ALLREDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )

     vnor=sqrt(vnor)
  end if
  
end subroutine vecnorALLMPI

subroutine vecnorMPI(v,n,vnor,inor,MPIComm,MPIrank,MPIroot,MPIsize) 

!-----------------------------------------------------------------------
!
! Compute the L1, L2 and L-inf norms of a vector V of length N
!
!-----------------------------------------------------------------------
!Returns correct result only for rank 0!!
  use      typre
  use MPI
  implicit none
  integer(ip), intent(in)  :: n,inor
  real(rp),    intent(in)  :: v(n)
  real(rp),    intent(out) :: vnor
  
  integer(ip)              :: i
  integer(ip) :: MPIComm,MPIrank,MPIroot,MPIsize
  integer(ip) :: ierr
  real(rp) :: vnor0

  vnor0=0.0_rp
  if(inor==0.or.inor==3) then
     do i=1,n
        if(abs(v(i))>vnor0) vnor0=abs(v(i))
     end do
     
     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_MAX, MPIroot,MPIcomm, ierr )
  else if(inor==1) then
     do i=1,n
        vnor0=vnor0+abs(v(i))
     end do
     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
  else if(inor==2) then
     do i=1,n
        vnor0=vnor0+v(i)*v(i)
     end do
     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
     vnor=sqrt(vnor)
  end if
  
end subroutine vecnorMPI

subroutine vecnorGenericMPI(v,n,ndime,vnor,inor,MPIcomm,MPIrank,MPIroot,MPIsize) 

!-----------------------------------------------------------------------
!
! Compute the L1, L2 and L-inf norms of a vector V of length N
!
!-----------------------------------------------------------------------
!Returns correct result only for rank 0!!
  use      typre
  use MPI
  implicit none
  integer(ip), intent(in)  :: n,inor,ndime
  real(rp),    intent(in)  :: v(ndime,n)
  real(rp),    intent(out) :: vnor
  
  integer(ip)              :: i,idime
  integer(ip) :: MPIcomm,MPIrank,MPIroot,MPIsize
  integer(ip) :: ierr
  real(rp) :: vnor0,val

  val  =0.0_rp
  vnor0=0.0_rp
  if(inor==0.or.inor==3) then

      do i=1,n
          do idime = 1,ndime
              if(abs(v(idime,i))>vnor0) vnor0=abs(v(idime,i))
          enddo
      end do
     
     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_MAX, MPIroot,MPIcomm, ierr )

  else if(inor==1) then

     do i=1,n
         do idime = 1,ndime
             vnor0=vnor0+abs(v(idime,i))
         enddo
     end do
     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )

  else if(inor==2) then

     do i = 1,n
        do idime = 1,ndime
         val = v(idime,i)
         vnor0 = vnor0 + val*val
        enddo
     end do

     call  MPI_REDUCE( vnor0, vnor, 1, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
     vnor=sqrt(vnor)
  end if
  
end subroutine vecnorGenericMPI

subroutine matNormF(ndime,matrix,matnor)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the norm of a matrix
   !    norm = [sum_j sum_i |A_ij|^2]^0.5
   !
   !-----------------------------------------------------------------------
   use      typre
   implicit none
   integer(ip), intent(in)    :: ndime
   real(rp),    intent(in)    :: matrix(ndime,ndime)
   real(rp),    intent(out)   :: matnor

   real(rp)                   :: aux
   integer(ip)                :: idime,jdime

   aux = 0.0_rp      
   do idime = 1,ndime
      do jdime = 1,ndime
         aux = aux + (matrix(idime,jdime)*matrix(idime,jdime))
      enddo
   enddo

   matnor = sqrt(aux)

end subroutine matNormF

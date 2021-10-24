subroutine vecres(norm,n,v1,v2,redif,zero)

!-----------------------------------------------------------------------
!
! Compute the relative difference between two vectors:
! 
! redif = 1*||v1 - v2|| / ||v1||      
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n,norm
  real(rp),    intent(in)  :: v1(n),v2(n),zero
  real(rp),    intent(out) :: redif
  integer(ip)              :: i
  real(rp)                 :: numer,denom,va,vo

  redif = 0.0_rp
  numer = 0.0_rp
  denom = 0.0_rp

  select case(norm)

  case(0)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = max(numer,abs(va-vo))
        denom = max(denom,abs(va))
     end do
     if(denom.gt.zero) redif = 1.0_rp*numer/denom

  case(1)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + abs(va-vo)
        denom = denom + abs(va)
     end do
     if(denom>zero) redif = 1.0_rp*numer/denom

  case(2)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + (va-vo)*(va-vo)
        denom = denom + va*va
     end do
     if(denom>zero) redif = 1.0_rp*sqrt(numer/denom)

  end select

end subroutine vecres


subroutine vecresMPI(norm,n,v1,v2,redif,zero,kfl_MPIComType,MPIcomm,MPIroot,MPIrank,MPIsize)

!-----------------------------------------------------------------------
!
! Compute the relative difference between two vectors:
! 
! redif = 100*||v1 - v2|| / ||v1||      
!
!-----------------------------------------------------------------------
  use typre
  use MPI
  implicit none
  integer(ip), intent(in)  :: n,norm
  real(rp),    intent(in)  :: v1(n),v2(n),zero
  real(rp),    intent(out) :: redif
  integer(ip), intent(in)  :: kfl_MPIComType,MPIcomm,MPIroot,MPIrank,MPIsize
  integer(ip)              :: i
  real(rp)                 :: numer,denom,va,vo,rwa(2),rwa2(2)
  integer                  :: ierr

  redif = 0.0_rp
  numer = 0.0_rp
  denom = 0.0_rp

  select case(norm)

  case(0)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = max(numer,abs(va-vo))
        denom = max(denom,abs(va))
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call MPI_REDUCE_Implementation(2_ip,rwa,rwa2,kfl_MPIComType,MPIcomm,MPIrank,MPIroot,MPIsize)
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom.gt.zero) redif = 1.0_rp*numer/denom
  case(1)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + abs(va-vo)
        denom = denom + abs(va)
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call MPI_REDUCE_Implementation(2_ip,rwa,rwa2,kfl_MPIComType,MPIcomm,MPIrank,MPIroot,MPIsize)
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom>zero) redif = 1.0_rp*numer/denom
  case(2)
     do i = 1,n
        va = v1(i)
        vo = v2(i)
        numer = numer + (va-vo)*(va-vo)
        denom = denom + va*va
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call MPI_REDUCE_Implementation(2_ip,rwa,rwa2,kfl_MPIComType,MPIcomm,MPIrank,MPIroot,MPIsize)
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom>zero) redif = 1.0_rp*sqrt(numer/denom)
  end select

contains
   
  
end subroutine vecresMPI

subroutine MPI_REDUCE_Implementation(rwasize,rwa,rwa2,kfl_MPIComType,MPIcomm,MPIrank,MPIroot,MPIsize)
   use typre
   use MPI
   implicit none
   integer(ip) :: rwasize, MPIcomm,MPIrank,MPIroot,MPIsize,kfl_MPIComType
   real(rp) :: rwa(rwasize),rwa2(rwasize)
   
   real(rp) :: auxrwa2(2,MPIsize)
   integer(ip) :: irank,isize
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer :: irequest, irequest01(MPIsize)
   integer status(MPI_STATUS_SIZE)
   integer :: ierr
   
   !Blocking
   if (kfl_MPIComType == 0) then
      call  MPI_REDUCE( rwa, rwa2, rwasize, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
   !Non-Blocking
   else
      !Everybody sends to root
      call MPI_ISEND(rwa,rwasize,MPI_REAL8, MPIroot, mtag1, MPIcomm,irequest, ierr)
   
      !Root receives
      if (MPIrank == MPIroot) then
         do irank = 0,MPIsize-1
            call MPI_IRECV(auxrwa2(:,irank+1), rwasize, MPI_REAL8, irank, mtag1, MPIcomm,irequest01(irank+1), ierr)
         enddo
         do irank = 0,MPIsize-1
            call MPI_WAIT(irequest01(irank+1), status,ierr)
         enddo
         do isize = 1,rwasize
            rwa2(isize) = sum(auxrwa2(isize,:))
         enddo
      endif  
      call MPI_WAIT(irequest, status, ierr)
   endif
end subroutine



subroutine vecresMPIHeterogeneous(norm,ndime1,ndime2,npoinLocal,startcomp1,startcomp2,ncomp,v1,v2,redif,zero,MPIcomm,MPIroot)

!-----------------------------------------------------------------------
!
! Compute the relative difference between two vectors:
! 
! redif = 100*||v1 - v2|| / ||v1||      
!
!-----------------------------------------------------------------------
  use typre
  use MPI
  implicit none
  integer(ip), intent(in)  :: ndime1,ndime2,ncomp,startcomp1,startcomp2,norm,npoinLocal
  real(rp),    intent(in)  :: v1(ndime1,npoinLocal),v2(ndime2,npoinLocal),zero
  real(rp),    intent(out) :: redif
  integer(ip)              :: i,MPIcomm,MPIroot,idime
  real(rp)                 :: numer,denom,va,vo,rwa(2),rwa2(2)
  integer                  :: ierr

  redif = 0.0_rp
  numer = 0.0_rp
  denom = 0.0_rp

  select case(norm)

  case(0)
     do i = 1,npoinLocal
        do idime = 0,ncomp-1
         va = v1(startcomp1+idime,i)
         vo = v2(startcomp2+idime,i)
         numer = max(numer,abs(va-vo))
         denom = max(denom,abs(va))
        enddo
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call  MPI_REDUCE( rwa, rwa2, 2, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom.gt.zero) redif = 1.0_rp*numer/denom

  case(1)
      do i = 1,npoinLocal
        do idime = 0,ncomp-1
         va = v1(startcomp1+idime,i)
         vo = v2(startcomp2+idime,i)
         numer = numer + abs(va-vo)
         denom = denom + abs(va)
        enddo
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call  MPI_REDUCE( rwa, rwa2, 2, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom>zero) redif = 1.0_rp*numer/denom

  case(2)
     do i = 1,npoinLocal
        do idime = 0,ncomp-1
         va = v1(startcomp1+idime,i)
         vo = v2(startcomp2+idime,i)
         numer = numer + (va-vo)*(va-vo)
         denom = denom + va*va
        enddo
     end do
     
     rwa(1) = numer
     rwa(2) = denom
     call  MPI_REDUCE( rwa, rwa2, 2, MPI_REAL8, MPI_SUM, MPIroot,MPIcomm, ierr )
     numer = rwa2(1)
     denom = rwa2(2)
     
     if(denom>zero) redif = 1.0_rp*sqrt(numer/denom)

  end select

end subroutine vecresMPIHeterogeneous

subroutine mbvab1(a,b,c,n1,n2,en1,en2)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix-vector product A = B C, where
! A -> R(n1), B -> Mat(n1,n2), C -> R(n2), and also computes  the
! Euclidian norms of A and C
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n1,n2
  real(rp),    intent(in)  :: b(n1,n2), c(n2)
  real(rp),    intent(out) :: a(n1),en1,en2
  integer(ip)              :: i,k

  do i=1,n1
     a(i)=0.0_rp
     do k=1,n2
        a(i)=a(i)+b(i,k)*c(k)
     end do
  end do
  
  en1=0.0_rp
  en2=0.0_rp
  do i=1,n1
     en1=en1+a(i)*a(i)
  end do
  do i=1,n2
     en2=en2+c(i)*c(i)
  end do
  en1=sqrt(en1)
  en2=sqrt(en2)

end subroutine mbvab1



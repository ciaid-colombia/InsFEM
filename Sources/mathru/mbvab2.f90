subroutine mbvab2(a,b,c,n1,n2)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix-vector product A = B^t C, where
! A -> R(n1), B -> Mat(n2,n1), C -> R(n2)
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n1,n2
  real(rp),    intent(in)  :: b(n2,n1),c(n2)
  real(rp),    intent(out) :: a(n1)
  integer(ip)              :: i,k

  do i=1,n1
     a(i)=0.0_rp
     do k=1,n2
        a(i)=a(i)+b(k,i)*c(k)
     end do
  end do
  
end subroutine mbvab2

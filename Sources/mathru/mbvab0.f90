subroutine mbvab0(a,b,c,n1,n2)

!------------------------------------------------------------------------
!
! This routine evaluates the matrix-vector product A = B C, where
! A -> R(n1), B -> Mat(n1,n2), C -> R(n2)
!
!------------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n1,n2
  real(rp),    intent(in)  :: b(n1,n2),c(n2)
  real(rp),    intent(out) :: a(n1)
  integer(ip)              :: i,j

  do i=1,n1
     a(i)=0.0_rp
     do j=1,n2
        a(i)=a(i)+b(i,j)*c(j)
     end do
  end do
  
end subroutine mbvab0

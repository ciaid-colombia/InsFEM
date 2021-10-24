subroutine mbmabt(a,b,c,n1,n2,n3)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix product A = B Ct, where
! A -> Mat(n1,n2), B -> Mat(n1,n3), C -> Mat(n2,n3)
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip) :: n1,n2,n3,i,j,k
  real(rp)    :: a(n1,n2), b(n1,n3), c(n2,n3)
    
  do i=1,n1
     do j=1,n2
        a(i,j)=0.0_rp
        do k=1,n3
           a(i,j)=a(i,j)+b(i,k)*c(j,k)
        end do
     end do
  end do

end subroutine mbmabt

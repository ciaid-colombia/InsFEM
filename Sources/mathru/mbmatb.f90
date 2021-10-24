subroutine mbmatb(a,b,c,n1,n2,n3)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix product A = Bt C, where
! A -> Mat(n1,n2), B -> Mat(n3,n1), C -> Mat(n3,n2)
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip) :: n1,n2,n3,i,j,k
  real(rp)    :: a(n1,n2), b(n3,n1), c(n3,n2)    

  do i=1,n1
     do j=1,n2
        a(i,j)=0.0_rp
        do k=1,n3
           a(i,j)=a(i,j)+b(k,i)*c(k,j)
        end do
     end do
  end do

end subroutine mbmatb


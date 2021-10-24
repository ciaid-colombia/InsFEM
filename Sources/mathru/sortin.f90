subroutine sortin(n,a)

!-----------------------------------------------------------------------
!
! Sort a vector
!
!-----------------------------------------------------------------------
  use       typre
  implicit none
  integer(ip) :: n,i,j,t
  integer(ip) :: a(n)
  
  do i=1,n-1
     do j=i+1,n
        if (a(i)>a(j)) then
           t   =a(i)
           a(i)=a(j) 
           a(j)=t
        end if
     end do
  end do
  
end subroutine sortin

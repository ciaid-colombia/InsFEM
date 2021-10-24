subroutine ordena(n,v)

!-----------------------------------------------------------------------
!****f* Mathru/ordena
! NAME 
! DESCRIPTION
!    This routine orders array v cording to its first column
! USES
! USED BY
!***
!-----------------------------------------------------------------------
  use       typre
  implicit   none
  integer(ip), intent(in)    :: n
  real(rp),    intent(inout) :: v(2,n)
  integer(ip)                :: i,j
  real(rp)                   :: v1,v2 

  do j=2,n
     v1=v(1,j)
     v2=v(2,j)
     do i=j-1,1,-1
        if(v(1,i)<=v1) go to 30
        v(1,i+1)=v(1,i)
        v(2,i+1)=v(2,i)
     end do
     i=0
30   v(1,i+1)=v1
     v(2,i+1)=v2
  end do

end subroutine ordena

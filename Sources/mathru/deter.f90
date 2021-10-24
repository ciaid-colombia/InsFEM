subroutine deter(a,det,nsize)

!-----------------------------------------------------------------------
!
! This routine calculates the determinant DETER of a matrix a
!
!    
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: nsize
  real(rp),    intent(in)  :: a(nsize,nsize)
  real(rp),    intent(out) :: det
  integer(ip)              :: isize,jsize
  real(rp)                 :: denom,t1,t2,t3,t4

  select case (nsize)
     
  case(1)
     det=a(1,1)

  case(2)
     det=a(1,1)*a(2,2)-a(2,1)*a(1,2)

  case(3)
     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
     det= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3

  case(4)
     t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
     t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
     t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
     t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
     det= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4

  end select

end subroutine deter

subroutine vecuni(n,v,mod)

!-----------------------------------------------------------------------
!
! This routine computes the length of vector V and converts it to  
! a unit one 
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in) :: n
  real(rp), intent(out)   :: v(n),mod
  integer(ip)             :: i

  mod=0.0_rp
  do i=1,n
     mod=mod + v(i)*v(i)
  end do
  mod=sqrt(mod)
  if(mod>1.0e-8) v=v/mod
 
end subroutine vecuni
      

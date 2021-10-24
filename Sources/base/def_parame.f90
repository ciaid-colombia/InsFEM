module def_parame

!-----------------------------------------------------------------------
!
!     Parameters
!
!-----------------------------------------------------------------------  
  use typre

! Frequently used mathematical constants (with precision to spare):
  real(rp),    parameter :: pi=3.141592653589793238462643383279502884197_rp
  real(rp),    parameter :: pio2=1.57079632679489661923132169163975144209858_rp
  real(rp),    parameter :: twopi=6.283185307179586476925286766559005768394_rp
  real(rp),    parameter :: sqrt2=1.41421356237309504880168872420969807856967_rp
  real(rp),    parameter :: euler=0.5772156649015328606065120900824024310422_rp
  real(rp),    parameter :: zero_rp=epsilon(0.0_rp)

! Yes/No, Load/Unload and figures
  integer(ip), parameter :: yes=1,no=0
  integer(ip), parameter :: mone=-1,zero=0,one=1,two=2,three=3,four=4
  integer(ip), parameter :: five=5,six=6,seven=7,eight=8,nine=9


! ! Initial allocation parameter arrays to avoid unassociated pointers
!   real(rp),    target :: raux1(1)=0.0_rp
!   real(rp),    target :: raux2(1,1)=0.0_rp
!   real(rp),    target :: raux3(1,1,1)=0.0_rp
!   real(rp),    target :: raux4(1,1,1,1)=0.0_rp
!   integer(ip), target :: iaux1(1)=0

end module def_parame

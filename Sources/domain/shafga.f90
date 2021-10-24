subroutine shafga(xieta,ndime,ltopo,ngaus,shaga)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions associated to gauss points
! for linear and quadratic isoparametric elements
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: ndime,ltopo,ngaus
  real(rp),    intent(in)  :: xieta(ndime)
  real(rp),    intent(out) :: shaga(ngaus) 
  
  select case(ndime)
  case(1)
     call shaga1(xieta(1),ngaus,shaga)
  case(2)
     call shaga2(xieta(1),xieta(ndime),ltopo,ngaus,shaga)
  case(3)
     call shaga3(xieta(1),xieta(2),xieta(ndime),ltopo,ngaus,shaga)
  end select
  
end subroutine shafga

subroutine rulepw(ndime,ngaus,nrule,posgp,weigp)

!-----------------------------------------------------------------------
! 
! This routine sets up the integration constants
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  ::  ndime,ngaus,nrule
  real(rp),    intent(out) ::  posgp(ndime,ngaus), weigp(ngaus)

  posgp=0.0_rp
  weigp=0.0_rp

  select case (nrule)
    
  case(1) 
     call ruqope(ndime,ngaus,posgp,weigp)
  case(2) 
     call ruqclo(ndime,ngaus,posgp,weigp)
  case(3) 
     call rutope(ndime,ngaus,posgp,weigp)
  case(4) 
     call rutclo(ndime,ngaus,posgp,weigp)
  case(5) 
     call rupope(ndime,ngaus,posgp,weigp)
  case(6) 
     call rupclo(ndime,ngaus,posgp,weigp)

  end select

end subroutine rulepw

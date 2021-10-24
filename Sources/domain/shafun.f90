subroutine shafun(&
     posgp,ndime,nnode,ntens,shape,deriv,heslo)

!-----------------------------------------------------------------------
!
!    This routine evaluates shape functions and their derivatives
!    for linear and quadratic isoparametric elements
!
!-----------------------------------------------------------------------

  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,nnode,ntens
  real(rp),    intent(in)  :: posgp(ndime)
  real(rp),    intent(out) :: shape(nnode),deriv(ndime,nnode)
  real(rp),    intent(out) :: heslo(ntens,nnode)

!
! Initializations
!
  shape=0.0_rp
  deriv=0.0_rp
  heslo=0.0_rp
!
! Evaluation of the shape functions
!
  if(ndime==1) then
     call shape1(posgp(1),nnode,shape,deriv)
  else if(ndime==2) then
     call shape2(posgp(1),posgp(2),nnode,shape,deriv,heslo)
  else if(ndime==3) then
     call shape3(posgp(1),posgp(2),posgp(3),nnode,shape,deriv,heslo)
  end if
  
end subroutine shafun

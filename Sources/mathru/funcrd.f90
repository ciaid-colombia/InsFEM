function funcrd(param,npara,ifuge,timev)
!------------------------------------------------------------------------
!
! This function yields the derivative of a parabolic or periodic
! evolution
!
!------------------------------------------------------------------------
  use      typre
  implicit none
  real(rp)                 :: funcrd
  integer(ip), intent(in)  :: npara,ifuge
  real(rp),    intent(in)  :: param(npara),timev
  real(rp)                 :: timea,zerom

  funcrd=0.0_rp 
  zerom=epsilon(1.0_rp)

  if(ifuge==0) then
!
! No time dependence 
!
     funcrd=0.0_rp
   
  else if(ifuge==1) then
!
! Parabolic evolution
!
     if(param(1)-zerom<=timev.and.timev<=param(2)+zerom) then 
        funcrd=2.0_rp*param(3)*timev+param(4)
     else if (timev>param(2)+zerom) then
        timea=param(2)
        funcrd=2.0_rp*param(3)*timea+param(4)
     else if (timev<param(1)-zerom) then
        timea=param(1)
        funcrd=2.0_rp*param(3)*timea+param(4)
     end if
  else if(ifuge==2) then
!
! Periodic evolution
!
     if(param(1)-zerom<=timev.and.timev>=param(2)+zerom) then 
        funcrd=param(3)*param(4)*cos(param(4)*timev+param(5))
     else if (timev>param(2)+zerom) then
        timea=param(2)
        funcrd=param(3)*param(4)*cos(param(4)*timea+param(5))
     else if (timev<param(1)-zerom) then
        timea=param(1)
        funcrd=param(3)*param(4)*cos(param(4)*timea+param(5))
     end if

  end if

end function funcrd

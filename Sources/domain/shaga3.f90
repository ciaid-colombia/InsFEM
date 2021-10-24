subroutine shaga3(s,t,z,ltopo,ngaus,shaga)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions associated to gauss points
! for 3D. The cases available so far are:
!
! BRICKS    -->   NGAUS =   1   8  27  
! SIMPLICES -->   NGAUS =   1   4
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: ngaus,ltopo
  real(rp),    intent(in)  :: s,t,z
  real(rp),    intent(out) :: shaga(ngaus)
  real(rp)                 :: r3,r6,a,b,c,d
!
! Non-available cases
!
  if(ltopo==0) then
     if((ngaus/=1).and.(ngaus/=8).and.(ngaus/=27))&
          call runend('SHAGA3: INTERPOLATION NOT AVAILABLE')
  else if(ltopo==1) then
     if((ngaus/=1).and.(ngaus/=4))&
          call runend('SHAGA3: INTERPOLATION NOT AVAILABLE')
  end if
!
! Available cases
!
  if (ngaus==1) then
     shaga(1)=1.0_rp
  else if (ltopo==0 .and. ngaus==8) then            ! Brick
     r3=1.0_rp/sqrt(3.0_rp)
     shaga(1)=-1.125_rp*r3*(s-r3)*(t-r3)*(z-r3)
     shaga(2)= 1.125_rp*r3*(s-r3)*(t-r3)*(z+r3)
     shaga(3)= 1.125_rp*r3*(s-r3)*(t+r3)*(z-r3)
     shaga(4)=-1.125_rp*r3*(s-r3)*(t+r3)*(z+r3)
     shaga(5)= 1.125_rp*r3*(s+r3)*(t-r3)*(z-r3)
     shaga(6)=-1.125_rp*r3*(s+r3)*(t-r3)*(z+r3)
     shaga(7)=-1.125_rp*r3*(s+r3)*(t+r3)*(z-r3)
     shaga(8)= 1.125_rp*r3*(s+r3)*(t+r3)*(z+r3)
  else if (ltopo==0 .and. ngaus==27) then       ! Brick
     r6=sqrt(.6)
     shaga( 1)= 125.0_rp/216.0_rp*(s-r6)*s*(t-r6)*t*(z-r6)*z
     shaga( 2)=-125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*t*(z-r6)*(z+r6)
     shaga( 3)= 125.0_rp/216.0_rp*(s-r6)*s*(t-r6)*t*(z+r6)*z
     shaga( 4)= 125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z-r6)*z
     shaga( 5)=-125.0_rp/54.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z-r6)*(z+r6)
     shaga( 6)= 125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z+r6)*z
     shaga( 7)= 125.0_rp/216.0_rp*(s-r6)*s*(t+r6)*t*(z-r6)*z
     shaga( 8)=-125.0_rp/108.0_rp*(s-r6)*s*(t+r6)*t*(z-r6)*(z+r6)
     shaga( 9)= 125.0_rp/216.0_rp*(s-r6)*s*(t+r6)*t*(z+r6)*z
     shaga(10)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z-r6)*z
     shaga(11)= 125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z-r6)*(z+r6)
     shaga(12)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z+r6)*z
     shaga(13)=-125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z-r6)*z
     shaga(14)= 125.0_rp/27.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z-r6)*(z+r6)
     shaga(15)=-125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z+r6)*z
     shaga(16)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z-r6)*z
     shaga(17)= 125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z-r6)*(z+r6)
     shaga(18)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z+r6)*z
     shaga(19)= 125.0_rp/216.0_rp*(s+r6)*s*(t-r6)*t*(z-r6)*z
     shaga(20)=-125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*t*(z-r6)*(z+r6)
     shaga(21)= 125.0_rp/216.0_rp*(s+r6)*s*(t-r6)*t*(z+r6)*z
     shaga(22)= 125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z-r6)*z
     shaga(23)=-125.0_rp/54.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z-r6)*(z+r6)
     shaga(24)= 125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z+r6)*z
     shaga(25)= 125.0_rp/216.0_rp*(s+r6)*s*(t+r6)*t*(z-r6)*z
     shaga(26)=-125.0_rp/108.0_rp*(s+r6)*s*(t+r6)*t*(z-r6)*(z+r6)
     shaga(27)= 125.0_rp/216.0_rp*(s+r6)*s*(t+r6)*t*(z+r6)*z
  else if (ltopo==1 .and. ngaus==4) then
     a= 1.927050983124842_rp
     b=-2.236067977499790_rp
     c=-2.236067977499789_rp   
     d=-2.236067977499789_rp
     shaga(1)=a+b*s+c*t+d*z
     a=-0.3090169943749473_rp
     b= 2.236067977499790_rp
     c= 0.0_rp
     d= 0.0_rp
     shaga(2)=a+b*s+c*t+d*z
     a=-0.3090169943749473_rp
     b= 0.0_rp
     c= 2.236067977499790_rp    
     d= 0.0_rp
     shaga(3)=a+b*s+c*t+d*z
     a=-0.3090169943749473_rp
     b= 0.0_rp
     c= 0.0_rp
     d= 2.236067977499790_rp
     shaga(4)=a+b*s+c*t+d*z
  end if

end subroutine shaga3

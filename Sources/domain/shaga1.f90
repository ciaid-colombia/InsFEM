
subroutine shaga1(    s,ngaus,shaga)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions associated to gauss points
! for 1D : NGAUS = 1, 2, 3, 4, 5
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: ngaus
  real(rp),    intent(in)  :: s
  real(rp),    intent(out) :: shaga(ngaus)
  real(rp)                 :: x1,x2,x3,x4,x5

  if(ngaus==1) then
     shaga(1)=1.0_rp
  else if(ngaus==2) then
     shaga(1)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)-s)
     shaga(2)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)+s)
  else if(ngaus==3) then
     shaga(1)= 5.0_rp*(s-sqrt(0.6_rp))*s/6.0_rp
     shaga(2)=-5.0_rp*(s-sqrt(0.6_rp))*(s+sqrt(0.6_rp))/3.0_rp
     shaga(3)= 5.0_rp*(s+sqrt(0.6_rp))*s/6.0_rp

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!*                                                           * 
!* Shaga con lagrange, en puntos de Gauss tomados de ruqope  *
!*                                                           *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  else if(ngaus==4) then
     x1=-0.861136311594053_rp
     x2=-0.339981043584856_rp
     x3= 0.339981043584856_rp
     x4= 0.861136311594053_rp     
     shaga(1)= ((s-x2)*(s-x3)*(s-x4))/((x1-x2)*(x1-x3)*(x1-x4))
     shaga(2)= ((s-x1)*(s-x3)*(s-x4))/((x2-x1)*(x2-x3)*(x2-x4))
     shaga(3)= ((s-x1)*(s-x2)*(s-x4))/((x3-x1)*(x3-x2)*(x3-x4))     
     shaga(4)= ((s-x1)*(s-x2)*(s-x3))/((x4-x1)*(x4-x2)*(x4-x3))
  else if(ngaus==5) then
     x1 = -0.906179845938664_rp
     x2 = -0.538469310105683_rp
     x3 =  0.0_rp
     x4 =  0.538469310105683_rp
     x5 =  0.906179845938664_rp     
     shaga(1)= ((s-x2)*(s-x3)*(s-x4)*(s-x5))/((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5))
     shaga(2)= ((s-x1)*(s-x3)*(s-x4)*(s-x5))/((x2-x1)*(x2-x3)*(x2-x4)*(x2-x5))
     shaga(3)= ((s-x1)*(s-x2)*(s-x4)*(s-x5))/((x3-x1)*(x3-x2)*(x3-x4)*(x3-x5))     
     shaga(4)= ((s-x1)*(s-x2)*(s-x3)*(s-x5))/((x4-x1)*(x4-x2)*(x4-x3)*(x4-x5))
     shaga(5)= ((s-x1)*(s-x2)*(s-x3)*(s-x4))/((x5-x1)*(x5-x2)*(x5-x3)*(x5-x4))
  end if
  
end subroutine shaga1

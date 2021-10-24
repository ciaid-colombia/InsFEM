subroutine shape3(s,t,z,nnode,shape,deriv,heslo)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions and their first and
! second derivatives 3-d standar continuous interpolation
! elements.
! 
! TETRAHEDRA:  4  10  &  20  nodes
! HEXAHEDRA:   8  27  &  64  nodes
! PRISM:       6             nodes
!
!-----------------------------------------------------------------------

  use      typre
  implicit none
  integer(ip), intent(in)  :: nnode
  real(rp),    intent(in)  :: s,t,z
  real(rp),    intent(out) :: deriv(3,nnode),shape(nnode),heslo(6,nnode)
  integer(ip)              :: istop,i
  real(rp)                 :: a1,a2,a3,a4,a,p1,p2,p3,z1,z2,z3,z4,s1,s2,s3,s4
  real(rp)                 :: t1,t2,t3,t4,sm,tm,zm,sq,tp,zp,s11,s21,s31,s41
  real(rp)                 :: t11,t21,t31,t41,z11,z21,z31,s12,s22,s32,s42
  real(rp)                 :: t12,t22,t32,t42,z41,z12,z22,z32,z42,sl,tl,zl

  istop=0
!
! Linear tetrahedron 
!
  if(nnode==4) then
     shape(   1) = 1.0_rp-s-t-z
     deriv(1, 1) =-1.0_rp
     deriv(2, 1) =-1.0_rp
     deriv(3, 1) =-1.0_rp
     shape(   2) = s
     deriv(1, 2) = 1.0_rp
     shape(   3) = t
     deriv(2, 3) = 1.0_rp
     shape(   4) = z
     deriv(3, 4) = 1.0_rp
!
! Quadratic tetrahedron 
!
  else if(nnode==10) then
     a1= 1.0_rp-s-t-z
     a2=s
     a3=t
     a4=z
     shape(   1) = (2.0_rp*a1-1.0_rp)*a1
     deriv(1, 1) = 1.0_rp-4.0_rp*a1
     deriv(2, 1) = 1.0_rp-4.0_rp*a1
     deriv(3, 1) = 1.0_rp-4.0_rp*a1
     do i = 1,6
        heslo(i, 1) = 4.0_rp
     enddo
     shape(   2) = (2.0_rp*a2-1.0_rp)*a2
     deriv(1, 2) = 4.0_rp*a2-1.0_rp
     heslo(1, 2) = 4.0_rp
     shape(   3) = (2.0_rp*a3-1.0_rp)*a3
     deriv(2, 3) = 4.0_rp*a3-1.0_rp
     heslo(2, 3) = 4.0_rp
     shape(   4) = (2.0_rp*a4-1.0_rp)*a4
     deriv(3, 4) = 4.0_rp*a4-1.0_rp
     heslo(3, 4) = 4.0_rp
     shape(   5) = 4.0_rp*a1*a2
     deriv(1, 5) = 4.0_rp*(a1-a2)
     deriv(2, 5) =-4.0_rp*a2
     deriv(3, 5) =-4.0_rp*a2
     heslo(1, 5) =-8.0_rp
     heslo(4, 5) =-4.0_rp
     heslo(5, 5) =-4.0_rp
     shape(   6) = 4.0_rp*a2*a3
     deriv(1, 6) = 4.0_rp*a3
     deriv(2, 6) = 4.0_rp*a2
     heslo(4, 6) = 4.0_rp
     shape(   7) = 4.0_rp*a1*a3
     deriv(1, 7) =-4.0_rp*a3
     deriv(2, 7) = 4.0_rp*(a1-a3)
     deriv(3, 7) =-4.0_rp*a3
     heslo(2, 7) =-8.0_rp
     heslo(4, 7) =-4.0_rp
     heslo(6, 7) =-4.0_rp
     shape(   8) = 4.0_rp*a1*a4
     deriv(1, 8) =-4.0_rp*a4
     deriv(2, 8) =-4.0_rp*a4
     deriv(3, 8) = 4.0_rp*(a1-a4)
     heslo(3, 8) =-8.0_rp
     heslo(5, 8) =-4.0_rp
     heslo(6, 8) =-4.0_rp
     shape(   9) = 4.0_rp*a2*a4
     deriv(1, 9) = 4.0_rp*a4
     deriv(3, 9) = 4.0_rp*a2
     heslo(5, 9) = 4.0_rp
     shape(  10) = 4.0_rp*a3*a4
     deriv(2,10) = 4.0_rp*a4
     deriv(3,10) = 4.0_rp*a3
     heslo(6,10) = 4.0_rp
!
! Cubic tetrahedron
!        
  else if(nnode==20) then
     a=4.5_rp
     p1= 1-s-t-z
     p2= 2.0_rp/3.0_rp-s-t-z
     p3= 1.0_rp/3.0_rp-s-t-z
     z1= 2.0_rp/3.0_rp-z
     z2= 1.0_rp/3.0_rp-z
     s1= 2.0_rp/3.0_rp-s
     s2= 1.0_rp/3.0_rp-s
     t1= 2.0_rp/3.0_rp-t
     t2= 1.0_rp/3.0_rp-t
     
     shape(   1) = a*p1*p2*p3
     deriv(1, 1) =-a*(p1*p2+p1*p3+p2*p3)
     deriv(2, 1) =-a*(p1*p2+p1*p3+p2*p3)
     deriv(3, 1) =-a*(p1*p2+p1*p3+p2*p3)
     shape(   2) = a*z*z1*z2
     deriv(1, 2) = 0.0_rp
     deriv(2, 2) = 0.0_rp
     deriv(3, 2) =-a*(z*z1+z*z2-z1*z2)
     shape(   3) = a*s*s1*s2
     deriv(1, 3) =-a*(s*s1+s*s2-s1*s2)
     deriv(2, 3) = 0.0_rp
     deriv(3, 3) = 0.0_rp
     shape(   4) = a*t*t1*t2
     deriv(1, 4) = 0.0_rp
     deriv(2, 4) =-a*(t*t1+t*t2-t1*t2)
     deriv(3, 4) = 0.0_rp
     shape(   5) = 3.0_rp*a*p1*p2*t
     deriv(1, 5) = 3.0_rp*a*(-p2*t-p1*t)
     deriv(2, 5) = 3.0_rp*a*(-p2*t-p1*t+p1*p2)
     deriv(3, 5) = 3.0_rp*a*(-p2*t-p1*t)      
     shape(   6) = 3.0_rp*a*p1*p2*z
     deriv(1, 6) = 3.0_rp*a*(-p2*z-p1*z)      
     deriv(2, 6) = 3.0_rp*a*(-p2*z-p1*z)
     deriv(3, 6) = 3.0_rp*a*(-p2*z-p1*z+p1*p2)      
     shape(   7) = 3.0_rp*a*p1*p2*s
     deriv(1, 7) = 3.0_rp*a*(-p2*s-p1*s+p1*p2)
     deriv(2, 7) = 3.0_rp*a*(-p2*s-p1*s)
     deriv(3, 7) = 3.0_rp*a*(-p2*s-p1*s)
     shape(   8) =-3.0_rp*a*p1*t2*t
     deriv(1, 8) = 3.0_rp*a*(t2*t)
     deriv(2, 8) = 3.0_rp*a*(t2*t+p1*t-p1*t2)
     deriv(3, 8) = 3.0_rp*a*(t2*t)
     shape(   9) =-3.0_rp*a*p1*s2*s
     deriv(1, 9) = 3.0_rp*a*(s2*s+p1*s-p1*s2)
     deriv(2, 9) = 3.0_rp*a*(s2*s)
     deriv(3, 9) = 3.0_rp*a*(s2*s)
     shape(  10) =-3.0_rp*a*p1*z2*z
     deriv(1,10) = 3.0_rp*a*(z2*z)
     deriv(2,10) = 3.0_rp*a*(z2*z)
     deriv(3,10) = 3.0_rp*a*(z2*z+p1*z-p1*z2)
     shape(  11) =-3.0_rp*a*t2*t*z
     deriv(1,11) = 0.0_rp
     deriv(2,11) =-3.0_rp*a*(t2*z-t*z)
     deriv(3,11) =-3.0_rp*a*t2*t
     shape(  12) =-3.0_rp*a*z2*t*z
     deriv(1,12) = 0.0_rp
     deriv(2,12) =-3.0_rp*a*z2*z
     deriv(3,12) =-3.0_rp*a*(t*z2-t*z)
     shape(  13) =-3.0_rp*a*z2*s*z
     deriv(1,13) =-3.0_rp*a*z2*z
     deriv(2,13) = 0.0_rp
     deriv(3,13) =-3.0_rp*a*(s*z2-s*z)
     shape(  14) =-3.0_rp*a*s2*s*z
     deriv(1,14) =-3.0_rp*a*(s2*z-s*z)
     deriv(2,14) = 0.0_rp
     deriv(3,14) =-3.0_rp*a*s2*s
     shape(  15) =-3.0_rp*a*s2*s*t
     deriv(1,15) =-3.0_rp*a*(s2*t-s*t)
     deriv(2,15) =-3.0_rp*a*s2*s
     deriv(3,15) = 0.0_rp
     shape(  16) =-3.0_rp*a*t2*t*s
     deriv(1,16) =-3.0_rp*a*t2*t
     deriv(2,16) =-3.0_rp*a*(t2*s-s*t)
     deriv(3,16) = 0.0_rp
     shape(  17) = 27.0_rp*p1*t*z
     deriv(1,17) =-27.0_rp*t*z
     deriv(2,17) = 27.0_rp*(p1*z-t*z)
     deriv(3,17) = 27.0_rp*(p1*t-t*z)
     shape(  18) = 27.0_rp*p1*s*z
     deriv(1,18) = 27.0_rp*(p1*z-s*z)
     deriv(2,18) =-27.0_rp*s*z
     deriv(3,18) = 27.0_rp*(p1*s-s*z)
     shape(  19) = 27.0_rp*p1*s*t
     deriv(1,19) = 27.0_rp*(p1*t-s*t)
     deriv(2,19) = 27.0_rp*(p1*s-s*t)
     deriv(3,19) =-27.0_rp*s*t
     shape(  20) = 27.0_rp*s*t*z
     deriv(1,20) = 27.0_rp*z*t
     deriv(2,20) = 27.0_rp*s*z
     deriv(3,20) = 27.0_rp*s*t
       
     heslo(1, 1) = 2.0_rp*a*(p1+p2+p3)
     heslo(2, 1) = 2.0_rp*a*(p1+p2+p3)
     heslo(3, 1) = 2.0_rp*a*(p1+p2+p3)
     heslo(4, 1) = 2.0_rp*a*(p1+p2+p3)
     heslo(5, 1) = 2.0_rp*a*(p1+p2+p3)
     heslo(6, 1) = 2.0_rp*a*(p1+p2+p3)
     
     heslo(3, 2) = 2.0_rp*a*(z-z1-z2)
     
     heslo(1, 3) = 2.0_rp*a*(s-s1-s2)
                  
     heslo(2, 4) = 2.0_rp*a*(t-t1-t2)
     
     heslo(1, 5) = 6.0_rp*a*t
     heslo(2, 5) = 6.0_rp*a*(-p1-p2+t)
     heslo(3, 5) = 6.0_rp*a*t
     heslo(4, 5) = 3.0_rp*a*(-p2-p1+2*t)
     heslo(5, 5) = 6.0_rp*a*t
     heslo(6, 5) = 3.0_rp*a*(-p2-p1+2*t)
                  
     heslo(1, 6) = 6.0_rp*a*z
     heslo(2, 6) = 6.0_rp*a*z
     heslo(3, 6) = 6.0_rp*a*(-p1-p2+z)
     heslo(4, 6) = 6.0_rp*a*z
     heslo(5, 6) = 3.0_rp*a*(-p2-p1+2*z)
     heslo(6, 6) = 3.0_rp*a*(-p2-p1+2*z)
                  
     heslo(1, 7) = 6.0_rp*a*(-p1-p2+s)
     heslo(2, 7) = 6.0_rp*a*s
     heslo(3, 7) = 6.0_rp*a*s
     heslo(4, 7) = 3.0_rp*a*(-p2-p1+2*s)
     heslo(5, 7) = 3.0_rp*a*(-p2-p1+2*s)
     heslo(6, 7) = 6.0_rp*a*s
                  
     heslo(1, 8) = 0.0_rp
     heslo(2, 8) = 6.0_rp*a*(t2-t+p1)
     heslo(3, 8) = 0.0_rp
     heslo(4, 8) = 3.0_rp*a*(t2-t)
     heslo(5, 8) = 0.0_rp
     heslo(6, 8) = 3.0_rp*a*(t2-t)
                  
     heslo(1, 9) = 6.0_rp*a*(s2-s+p1)
     heslo(2, 9) = 0.0_rp
     heslo(3, 9) = 0.0_rp
     heslo(4, 9) = 3.0_rp*a*(s2-s)
     heslo(5, 9) = 3.0_rp*a*(s2-s)
     heslo(6, 9) = 0.0_rp
                  
     heslo(1,10) = 0.0_rp
     heslo(2,10) = 0.0_rp
     heslo(3,10) = 6.0_rp*a*(z2-z+p1)
     heslo(4,10) = 0.0_rp
     heslo(5,10) = 3.0_rp*a*(z2-z)           
     heslo(6,10) = 3.0_rp*a*(z2-z)
                  
     heslo(1,11) = 0.0_rp
     heslo(2,11) = 6.0_rp*a*z
     heslo(3,11) = 0.0_rp
     heslo(4,11) = 0.0_rp
     heslo(5,11) = 0.0_rp
     heslo(6,11) = 3.0_rp*a*(t-t2)
                  
     heslo(1,12) = 0.0_rp
     heslo(2,12) = 0.0_rp
     heslo(3,12) = 6.0_rp*a*t
     heslo(4,12) = 0.0_rp
     heslo(5,12) = 0.0_rp
     heslo(6,12) = 3.0_rp*a*(z-z2)
                  
     heslo(1,13) = 0.0_rp
     heslo(2,13) = 0.0_rp
     heslo(3,13) = 6.0_rp*a*s
     heslo(4,13) = 0.0_rp
     heslo(5,13) = 3.0_rp*a*(z-z2)
     heslo(6,13) = 0.0_rp
                  
     heslo(1,14) = 6.0_rp*a*z
     heslo(2,14) = 0.0_rp
     heslo(3,14) = 0.0_rp
     heslo(4,14) = 0.0_rp
     heslo(5,14) = 3.0_rp*a*(s-s2)
     heslo(6,14) = 0.0_rp
                  
     heslo(1,15) = 6.0_rp*a*t
     heslo(2,15) = 0.0_rp
     heslo(3,15) = 0.0_rp
     heslo(4,15) = 3.0_rp*a*(s-s2)
     heslo(5,15) = 0.0_rp
     heslo(6,15) = 0.0_rp
                  
     heslo(1,16) = 0.0_rp
     heslo(2,16) = 6.0_rp*a*s
     heslo(3,16) = 0.0_rp
     heslo(4,16) = 3.0_rp*a*(t-t2)
     heslo(5,16) = 0.0_rp
     heslo(6,16) = 0.0_rp
                  
     heslo(1,17) = 0.0_rp
     heslo(2,17) =-54.0_rp*z
     heslo(3,17) =-54.0_rp*t
     heslo(4,17) =-27.0_rp*z
     heslo(5,17) =-27.0_rp*t
     heslo(6,17) = 27.0_rp*(p1-z-t)
                  
     heslo(1,18) =-54.0_rp*z
     heslo(2,18) = 0.0_rp
     heslo(3,18) =-54.0_rp*s
     heslo(4,18) =-27.0_rp*z
     heslo(5,18) =-27.0_rp*s
     heslo(6,18) = 27.0_rp*(p1-z-s)
                  
     heslo(1,19) =-54.0_rp*t
     heslo(2,19) =-54.0_rp*s
     heslo(3,19) = 0.0_rp
     heslo(4,19) = 27.0_rp*(p1-t-s)
     heslo(5,19) =-27.0_rp*t
     heslo(6,19) =-27.0_rp*s
                  
     heslo(1,20) = 0.0_rp
     heslo(2,20) = 0.0_rp
     heslo(3,20) = 0.0_rp
     heslo(4,20) = 27.0_rp*z
     heslo(5,20) = 27.0_rp*t
     heslo(6,20) = 27.0_rp*s
!
! Trilinear brick 
!        
  else if(nnode==8) then
     sm = 0.5_rp*(1.0_rp-s)
     tm = 0.5_rp*(1.0_rp-t)
     zm = 0.5_rp*(1.0_rp-z)
     sq = 0.5_rp*(1.0_rp+s)
     tp = 0.5_rp*(1.0_rp+t)
     zp = 0.5_rp*(1.0_rp+z)
     shape(   1) = sm*tm*zm
     deriv(1, 1) =-0.5_rp*tm*zm
     deriv(2, 1) =-0.5_rp*sm*zm
     deriv(3, 1) =-0.5_rp*sm*tm
     heslo(4, 1) = 0.25_rp*zm
     heslo(5, 1) = 0.25_rp*tm
     heslo(6, 1) = 0.25_rp*sm
     shape(   2) = sq*tm*zm
     deriv(1, 2) = 0.5_rp*tm*zm
     deriv(2, 2) =-0.5_rp*sq*zm
     deriv(3, 2) =-0.5_rp*sq*tm
     heslo(4, 2) =-0.25_rp*zm
     heslo(5, 2) =-0.25_rp*tm
     heslo(6, 2) = 0.25_rp*sq
     shape(   3) = sq*tp*zm
     deriv(1, 3) = 0.5_rp*tp*zm
     deriv(2, 3) = 0.5_rp*sq*zm
     deriv(3, 3) =-0.5_rp*sq*tp
     heslo(4, 3) = 0.25_rp*zm
     heslo(5, 3) =-0.25_rp*tp
     heslo(6, 3) =-0.25_rp*sq
     shape(   4) = sm*tp*zm
     deriv(1, 4) =-0.5_rp*tp*zm
     deriv(2, 4) = 0.5_rp*sm*zm
     deriv(3, 4) =-0.5_rp*sm*tp
     heslo(4, 4) =-0.25_rp*zm
     heslo(5, 4) = 0.25_rp*tp
     heslo(6, 4) =-0.25_rp*sm
     shape(   5) = sm*tm*zp
     deriv(1, 5) =-0.5_rp*tm*zp
     deriv(2, 5) =-0.5_rp*sm*zp
     deriv(3, 5) = 0.5_rp*sm*tm
     heslo(4, 5) = 0.25_rp*zp
     heslo(5, 5) =-0.25_rp*tm
     heslo(6, 5) =-0.25_rp*sm
     shape(   6) = sq*tm*zp 
     deriv(1, 6) = 0.5_rp*tm*zp
     deriv(2, 6) =-0.5_rp*sq*zp
     deriv(3, 6) = 0.5_rp*sq*tm
     heslo(4, 6) =-0.25_rp*zp
     heslo(5, 6) = 0.25_rp*tm
     heslo(6, 6) =-0.25_rp*sq
     shape(   7) = sq*tp*zp
     deriv(1, 7) = 0.5_rp*tp*zp
     deriv(2, 7) = 0.5_rp*sq*zp
     deriv(3, 7) = 0.5_rp*sq*tp
     heslo(4, 7) = 0.25_rp*zp
     heslo(5, 7) = 0.25_rp*tp
     heslo(6, 7) = 0.25_rp*sq
     shape(   8) = sm*tp*zp
     deriv(1, 8) =-0.5_rp*tp*zp
     deriv(2, 8) = 0.5_rp*sm*zp
     deriv(3, 8) = 0.5_rp*sm*tp
     heslo(4, 8) =-0.25_rp*zp
     heslo(5, 8) =-0.25_rp*tp
     heslo(6, 8) = 0.25_rp*sm
!
! Triquadratic brick
!        
  else if(nnode==27) then
     sl=s*(s-1.0_rp)
     tl=t*(t-1.0_rp)
     zl=z*(z-1.0_rp)
     sq=s*(s+1.0_rp)
     tp=t*(t+1.0_rp)
     zp=z*(z+1.0_rp)
     s1= 2.0_rp*s-1.0_rp
     t1= 2.0_rp*t-1.0_rp
     z1= 2.0_rp*z-1.0_rp
     s2= 1.0_rp-s*s
     t2= 1.0_rp-t*t
     z2= 1.0_rp-z*z
     s3= 1.0_rp+2.0_rp*s
     t3= 1.0_rp+2.0_rp*t
     z3= 1.0_rp+2.0_rp*z
     s4=-2.0_rp*s
     t4=-2.0_rp*t
     z4=-2.0_rp*z
     shape(   1) = 0.125_rp*sl*tl*zl
     deriv(1, 1) = 0.125_rp*s1*tl*zl
     deriv(2, 1) = 0.125_rp*sl*t1*zl
     deriv(3, 1) = 0.125_rp*sl*tl*z1
     heslo(1, 1) = 0.25_rp*tl*zl
     heslo(2, 1) = 0.25_rp*sl*zl
     heslo(3, 1) = 0.25_rp*sl*tl
     heslo(4, 1) = 0.125_rp*s1*t1*zl
     heslo(5, 1) = 0.125_rp*s1*tl*z1
     heslo(6, 1) = 0.125_rp*sl*t1*z1
     shape(   2) = 0.125_rp*sq*tl*zl
     deriv(1, 2) = 0.125_rp*s3*tl*zl
     deriv(2, 2) = 0.125_rp*sq*t1*zl
     deriv(3, 2) = 0.125_rp*sq*tl*z1
     heslo(1, 2) = 0.25_rp*tl*zl
     heslo(2, 2) = 0.25_rp*sq*zl
     heslo(3, 2) = 0.25_rp*sq*tl
     heslo(4, 2) = 0.125_rp*s3*t1*zl
     heslo(5, 2) = 0.125_rp*s3*tl*z1
     heslo(6, 2) = 0.125_rp*sq*t1*z1
     shape(   3) = 0.125_rp*sq*tp*zl
     deriv(1, 3) = 0.125_rp*s3*tp*zl
     deriv(2, 3) = 0.125_rp*sq*t3*zl
     deriv(3, 3) = 0.125_rp*sq*tp*z1
     heslo(1, 3) = 0.25_rp*tp*zl
     heslo(2, 3) = 0.25_rp*sq*zl
     heslo(3, 3) = 0.25_rp*sq*tp
     heslo(4, 3) = 0.125_rp*s3*t3*zl
     heslo(5, 3) = 0.125_rp*s3*tp*z1
     heslo(6, 3) = 0.125_rp*sq*t3*z1
     shape(   4) = 0.125_rp*sl*tp*zl
     deriv(1, 4) = 0.125_rp*s1*tp*zl
     deriv(2, 4) = 0.125_rp*sl*t3*zl
     deriv(3, 4) = 0.125_rp*sl*tp*z1
     heslo(1, 4) = 0.25_rp*tp*zl
     heslo(2, 4) = 0.25_rp*sl*zl
     heslo(3, 4) = 0.25_rp*sl*tp
     heslo(4, 4) = 0.125_rp*s1*t3*zl
     heslo(5, 4) = 0.125_rp*s1*tp*z1
     heslo(6, 4) = 0.125_rp*sl*t3*z1
     shape(   5) = 0.125_rp*sl*tl*zp
     deriv(1, 5) = 0.125_rp*s1*tl*zp
     deriv(2, 5) = 0.125_rp*sl*t1*zp
     deriv(3, 5) = 0.125_rp*sl*tl*z3
     heslo(1, 5) = 0.25_rp*tl*zp
     heslo(2, 5) = 0.25_rp*sl*zp
     heslo(3, 5) = 0.25_rp*sl*tl
     heslo(4, 5) = 0.125_rp*s1*t1*zp
     heslo(5, 5) = 0.125_rp*s1*tl*z3
     heslo(6, 5) = 0.125_rp*sl*t1*z3
     shape(   6) = 0.125_rp*sq*tl*zp
     deriv(1, 6) = 0.125_rp*s3*tl*zp
     deriv(2, 6) = 0.125_rp*sq*t1*zp
     deriv(3, 6) = 0.125_rp*sq*tl*z3
     heslo(1, 6) = 0.25_rp*tl*zp
     heslo(2, 6) = 0.25_rp*sq*zp
     heslo(3, 6) = 0.25_rp*sq*tl
     heslo(4, 6) = 0.125_rp*s3*t1*zp
     heslo(5, 6) = 0.125_rp*s3*tl*z3
     heslo(6, 6) = 0.125_rp*sq*t1*z3
     shape(   7) = 0.125_rp*sq*tp*zp
     deriv(1, 7) = 0.125_rp*s3*tp*zp
     deriv(2, 7) = 0.125_rp*sq*t3*zp
     deriv(3, 7) = 0.125_rp*sq*tp*z3
     heslo(1, 7) = 0.25_rp*tp*zp
     heslo(2, 7) = 0.25_rp*sq*zp
     heslo(3, 7) = 0.25_rp*sq*tp
     heslo(4, 7) = 0.125_rp*s3*t3*zp
     heslo(5, 7) = 0.125_rp*s3*tp*z3
     heslo(6, 7) = 0.125_rp*sq*t3*z3
     shape(   8) = 0.125_rp*sl*tp*zp
     deriv(1, 8) = 0.125_rp*s1*tp*zp
     deriv(2, 8) = 0.125_rp*sl*t3*zp
     deriv(3, 8) = 0.125_rp*sl*tp*z3
     heslo(1, 8) = 0.25_rp*tp*zp
     heslo(2, 8) = 0.25_rp*sl*zp
     heslo(3, 8) = 0.25_rp*sl*tp
     heslo(4, 8) = 0.125_rp*s1*t3*zp
     heslo(5, 8) = 0.125_rp*s1*tp*z3
     heslo(6, 8) = 0.125_rp*sl*t3*z3
     shape(   9) = 0.25_rp*s2*tl*zl
     deriv(1, 9) = 0.25_rp*s4*tl*zl
     deriv(2, 9) = 0.25_rp*s2*t1*zl
     deriv(3, 9) = 0.25_rp*s2*tl*z1
     heslo(1, 9) =-0.5_rp*tl*zl
     heslo(2, 9) = 0.5_rp*s2*zl
     heslo(3, 9) = 0.5_rp*s2*tl
     heslo(4, 9) = 0.25_rp*s4*t1*zl
     heslo(5, 9) = 0.25_rp*s4*tl*z1
     heslo(6, 9) = 0.25_rp*s2*t1*z1
     shape(  10) = 0.25_rp*sq*t2*zl
     deriv(1,10) = 0.25_rp*s3*t2*zl
     deriv(2,10) = 0.25_rp*sq*t4*zl
     deriv(3,10) = 0.25_rp*sq*t2*z1
     heslo(1,10) = 0.5_rp*t2*zl
     heslo(2,10) =-0.5_rp*sq*zl
     heslo(3,10) = 0.5_rp*sq*t2
     heslo(4,10) = 0.25_rp*s3*t4*zl
     heslo(5,10) = 0.25_rp*s3*t2*z1
     heslo(6,10) = 0.25_rp*sq*t4*z1
     shape(  11) = 0.25_rp*s2*tp*zl
     deriv(1,11) = 0.25_rp*s4*tp*zl
     deriv(2,11) = 0.25_rp*s2*t3*zl
     deriv(3,11) = 0.25_rp*s2*tp*z1
     heslo(1,11) =-0.5_rp*tp*zl
     heslo(2,11) = 0.5_rp*s2*zl
     heslo(3,11) = 0.5_rp*s2*tp
     heslo(4,11) = 0.25_rp*s4*t3*zl
     heslo(5,11) = 0.25_rp*s4*tp*z1
     heslo(6,11) = 0.25_rp*s2*t3*z1
     shape(  12) = 0.25_rp*sl*t2*zl
     deriv(1,12) = 0.25_rp*s1*t2*zl
     deriv(2,12) = 0.25_rp*sl*t4*zl
     deriv(3,12) = 0.25_rp*sl*t2*z1
     heslo(1,12) = 0.5_rp*t2*zl
     heslo(2,12) =-0.5_rp*sl*zl
     heslo(3,12) = 0.5_rp*sl*t2
     heslo(4,12) = 0.25_rp*s1*t4*zl
     heslo(5,12) = 0.25_rp*s1*t2*z1
     heslo(6,12) = 0.25_rp*sl*t4*z1
     shape(  13) = 0.25_rp*sl*tl*z2
     deriv(1,13) = 0.25_rp*s1*tl*z2
     deriv(2,13) = 0.25_rp*sl*t1*z2
     deriv(3,13) = 0.25_rp*sl*tl*z4
     heslo(1,13) = 0.5_rp*tl*z2
     heslo(2,13) = 0.5_rp*sl*z2
     heslo(3,13) =-0.5_rp*sl*tl
     heslo(4,13) = 0.25_rp*s1*t1*z2
     heslo(5,13) = 0.25_rp*s1*tl*z4
     heslo(6,13) = 0.25_rp*sl*t1*z4
     shape(  14) = 0.25_rp*sq*tl*z2
     deriv(1,14) = 0.25_rp*s3*tl*z2
     deriv(2,14) = 0.25_rp*sq*t1*z2
     deriv(3,14) = 0.25_rp*sq*tl*z4
     heslo(1,14) = 0.5_rp*tl*z2
     heslo(2,14) = 0.5_rp*sq*z2
     heslo(3,14) =-0.5_rp*sq*tl
     heslo(4,14) = 0.25_rp*s3*t1*z2
     heslo(5,14) = 0.25_rp*s3*tl*z4
     heslo(6,14) = 0.25_rp*sq*t1*z4
     shape(  15) = 0.25_rp*sq*tp*z2
     deriv(1,15) = 0.25_rp*s3*tp*z2
     deriv(2,15) = 0.25_rp*sq*t3*z2
     deriv(3,15) = 0.25_rp*sq*tp*z4
     heslo(1,15) = 0.5_rp*tp*z2
     heslo(2,15) = 0.5_rp*sq*z2
     heslo(3,15) =-0.5_rp*sq*tp
     heslo(4,15) = 0.25_rp*s3*t3*z2
     heslo(5,15) = 0.25_rp*s3*tp*z4
     heslo(6,15) = 0.25_rp*sq*t3*z4
     shape(  16) = 0.25_rp*sl*tp*z2
     deriv(1,16) = 0.25_rp*s1*tp*z2
     deriv(2,16) = 0.25_rp*sl*t3*z2
     deriv(3,16) = 0.25_rp*sl*tp*z4
     heslo(1,16) = 0.5_rp*tp*z2
     heslo(2,16) = 0.5_rp*sl*z2
     heslo(3,16) =-0.5_rp*sl*tp
     heslo(4,16) = 0.25_rp*s1*t3*z2
     heslo(5,16) = 0.25_rp*s1*tp*z4
     heslo(6,16) = 0.25_rp*sl*t3*z4
     shape(  17) = 0.25_rp*s2*tl*zp
     deriv(1,17) = 0.25_rp*s4*tl*zp
     deriv(2,17) = 0.25_rp*s2*t1*zp
     deriv(3,17) = 0.25_rp*s2*tl*z3
     heslo(1,17) =-0.5_rp*tl*zp
     heslo(2,17) = 0.5_rp*s2*zp
     heslo(3,17) = 0.5_rp*s2*tl
     heslo(4,17) = 0.25_rp*s4*t1*zp
     heslo(5,17) = 0.25_rp*s4*tl*z3
     heslo(6,17) = 0.25_rp*s2*t1*z3
     shape(  18) = 0.25_rp*sq*t2*zp
     deriv(1,18) = 0.25_rp*s3*t2*zp
     deriv(2,18) = 0.25_rp*sq*t4*zp
     deriv(3,18) = 0.25_rp*sq*t2*z3
     heslo(1,18) = 0.5_rp*t2*zp
     heslo(2,18) =-0.5_rp*sq*zp
     heslo(3,18) = 0.5_rp*sq*t2
     heslo(4,18) = 0.25_rp*s3*t4*zp
     heslo(5,18) = 0.25_rp*s3*t2*z3
     heslo(6,18) = 0.25_rp*sq*t4*z3
     shape(  19) = 0.25_rp*s2*tp*zp
     deriv(1,19) = 0.25_rp*s4*tp*zp
     deriv(2,19) = 0.25_rp*s2*t3*zp
     deriv(3,19) = 0.25_rp*s2*tp*z3
     heslo(1,19) =-0.5_rp*tp*zp
     heslo(2,19) = 0.5_rp*s2*zp
     heslo(3,19) = 0.5_rp*s2*tp
     heslo(4,19) = 0.25_rp*s4*t3*zp
     heslo(5,19) = 0.25_rp*s4*tp*z3
     heslo(6,19) = 0.25_rp*s2*t3*z3
     shape(  20) = 0.25_rp*sl*t2*zp
     deriv(1,20) = 0.25_rp*s1*t2*zp
     deriv(2,20) = 0.25_rp*sl*t4*zp
     deriv(3,20) = 0.25_rp*sl*t2*z3
     heslo(1,20) = 0.5_rp*t2*zp
     heslo(2,20) =-0.5_rp*sl*zp
     heslo(3,20) = 0.5_rp*sl*t2
     heslo(4,20) = 0.25_rp*s1*t4*zp
     heslo(5,20) = 0.25_rp*s1*t2*z3
     heslo(6,20) = 0.25_rp*sl*t4*z3
     shape(  21) = 0.5_rp*s2*t2*zl
     deriv(1,21) = 0.5_rp*s4*t2*zl
     deriv(2,21) = 0.5_rp*s2*t4*zl
     deriv(3,21) = 0.5_rp*s2*t2*z1
     heslo(1,21) =-t2*zl
     heslo(2,21) =-s2*zl
     heslo(3,21) = s2*t2
     heslo(4,21) = 0.5_rp*s4*t4*zl
     heslo(5,21) = 0.5_rp*s4*t2*z1
     heslo(6,21) = 0.5_rp*s2*t4*z1
     shape(  22) = 0.5_rp*s2*tl*z2
     deriv(1,22) = 0.5_rp*s4*tl*z2
     deriv(2,22) = 0.5_rp*s2*t1*z2
     deriv(3,22) = 0.5_rp*s2*tl*z4
     heslo(1,22) =-tl*z2
     heslo(2,22) = s2*z2
     heslo(3,22) =-s2*tl
     heslo(4,22) = 0.5_rp*s4*t1*z2
     heslo(5,22) = 0.5_rp*s4*tl*z4
     heslo(6,22) = 0.5_rp*s2*t1*z4
     shape(  23) = 0.5_rp*sq*t2*z2
     deriv(1,23) = 0.5_rp*s3*t2*z2
     deriv(2,23) = 0.5_rp*sq*t4*z2
     deriv(3,23) = 0.5_rp*sq*t2*z4
     heslo(1,23) = t2*z2
     heslo(2,23) =-sq*z2
     heslo(3,23) =-sq*t2
     heslo(4,23) = 0.5_rp*s3*t4*z2
     heslo(5,23) = 0.5_rp*s3*t2*z4
     heslo(6,23) = 0.5_rp*sq*t4*z4
     shape(  24) = 0.5_rp*s2*tp*z2
     deriv(1,24) = 0.5_rp*s4*tp*z2
     deriv(2,24) = 0.5_rp*s2*t3*z2
     deriv(3,24) = 0.5_rp*s2*tp*z4
     heslo(1,24) =-tp*z2
     heslo(2,24) = s2*z2
     heslo(3,24) =-s2*tp
     heslo(4,24) = 0.5_rp*s4*t3*z2
     heslo(5,24) = 0.5_rp*s4*tp*z4
     heslo(6,24) = 0.5_rp*s2*t3*z4
     shape(  25) = 0.5_rp*sl*t2*z2
     deriv(1,25) = 0.5_rp*s1*t2*z2
     deriv(2,25) = 0.5_rp*sl*t4*z2
     deriv(3,25) = 0.5_rp*sl*t2*z4
     heslo(1,25) = t2*z2
     heslo(2,25) =-sl*z2
     heslo(3,25) =-sl*t2
     heslo(4,25) = 0.5_rp*s1*t4*z2
     heslo(5,25) = 0.5_rp*s1*t2*z4
     heslo(6,25) = 0.5_rp*sl*t4*z4
     shape(  26) = 0.5_rp*s2*t2*zp
     deriv(1,26) = 0.5_rp*s4*t2*zp
     deriv(2,26) = 0.5_rp*s2*t4*zp
     deriv(3,26) = 0.5_rp*s2*t2*z3
     heslo(1,26) =-t2*zp
     heslo(2,26) =-s2*zp
     heslo(3,26) = s2*t2
     heslo(4,26) = 0.5_rp*s4*t4*zp
     heslo(5,26) = 0.5_rp*s4*t2*z3
     heslo(6,26) = 0.5_rp*s2*t4*z3
     shape(  27) = s2*t2*z2
     deriv(1,27) = s4*t2*z2
     deriv(2,27) = s2*t4*z2
     deriv(3,27) = s2*t2*z4
     heslo(1,27) =-2.0_rp*t2*z2
     heslo(2,27) =-2.0_rp*s2*z2
     heslo(3,27) =-2.0_rp*s2*t2
     heslo(4,27) = s4*t4*z2
     heslo(5,27) = s4*t2*z4
     heslo(6,27) = s2*t4*z4
!
! Tricubic brick
!        
    else if (nnode==64) then
       a=729.0_rp/4096.0_rp
       s1=(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)*(1.0_rp-s)                           
       s2=(1.0_rp+s)*(1.0_rp/3.0_rp-s)*(1.0_rp-s)                              
       s3=(1.0_rp+s)*(1.0_rp/3.0_rp+s)*(1.0_rp-s)
       s4=(1.0_rp+s)*(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)                                
       t1=(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)*(1.0_rp-t)         
       t2=(1.0_rp+s)*(1.0_rp/3.0_rp-t)*(1.0_rp-t)
       t3=(1.0_rp+t)*(1.0_rp/3.0_rp+t)*(1.0_rp-t)                         
       t4=(1.0_rp+t)*(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)
       z1=(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)*(1.0_rp-z)
       z2=(1.0_rp+z)*(1.0_rp/3.0_rp-z)*(1.0_rp-z)
       z3=(1.0_rp+z)*(1.0_rp/3.0_rp+z)*(1.0_rp-z)
       z4=(1.0_rp+z)*(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)
       s11=(1.0_rp/3.0_rp-s)*(1.0_rp-s)-(1.0_rp/3.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)
       s21=(1.0_rp/3.0_rp-s)*(1.0_rp-s)-(1.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp-s)
       s31=(1.0_rp/3.0_rp+s)*(1.0_rp-s)+(1.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp+s)
       s41=(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)+(1.0_rp+s)*(1.0_rp/3.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp+s)
       t11=(1.0_rp/3.0_rp-t)*(1.0_rp-t)-(1.0_rp/3.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)
       t21=(1.0_rp/3.0_rp-t)*(1.0_rp-t)-(1.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp-t)    
       t31=(1.0_rp/3.0_rp+t)*(1.0_rp-t)+(1.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp+t)    
       t41=(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)+(1.0_rp+t)*(1.0_rp/3.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp+t)
       z11=(1.0_rp/3.0_rp-z)*(1.0_rp-z)-(1.0_rp/3.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)
       z21=(1.0_rp/3.0_rp-z)*(1.0_rp-z)-(1.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp-z)    
       z31=(1.0_rp/3.0_rp+z)*(1.0_rp-z)+(1.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp+z)    
       z41=(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)+(1.0_rp+z)*(1.0_rp/3.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp+z)
       shape(   1) =   a*s1*t1*z1
       deriv(1, 1) =   a*s11*t1*z1
       deriv(2, 1) =   a*s1*t11*z1
       deriv(3, 1) =   a*s1*t1*z11
       shape(   2) =   a*s4*t1*z1
       deriv(1, 2) =   a*s41*t1*z1
       deriv(2, 2) =   a*s4*t11*z1
       deriv(3, 2) =   a*s4*t1*z11
       shape(   3) =   a*s4*t4*z1
       deriv(1, 3) =   a*s41*t4*z1
       deriv(2, 3) =   a*s4*t41*z1
       deriv(3, 3) =   a*s4*t4*z11
       shape(   4) =   a*s1*t4*z1
       deriv(1, 4) =   a*s11*t4*z1
       deriv(2, 4) =   a*s1*t41*z1
       deriv(3, 4) =   a*s1*t4*z11
       shape(   5) =   a*s1*t1*z4
       deriv(1, 5) =   a*s11*t1*z4
       deriv(2, 5) =   a*s1*t11*z4
       deriv(3, 5) =   a*s1*t1*z41
       shape(   6) =   a*s4*t1*z4
       deriv(1, 6) =   a*s41*t1*z4
       deriv(2, 6) =   a*s4*t11*z4
       deriv(3, 6) =   a*s4*t1*z41
       shape(   7) =   a*s4*t4*z4
       deriv(1, 7) =   a*s41*t4*z4
       deriv(2, 7) =   a*s4*t41*z4
       deriv(3, 7) =   a*s4*t4*z41
       shape(   8) =   a*s1*t4*z4
       deriv(1, 8) =   a*s11*t4*z4
       deriv(2, 8) =   a*s1*t41*z4
       deriv(3, 8) =   a*s1*t4*z41
       shape(   9) = 3.0_rp*a*s2*t1*z1
       deriv(1, 9) = 3.0_rp*a*s21*t1*z1
       deriv(2, 9) = 3.0_rp*a*s2*t11*z1
       deriv(3, 9) = 3.0_rp*a*s2*t1*z11
       shape(  10) = 3.0_rp*a*s3*t1*z1
       deriv(1,10) = 3.0_rp*a*s31*t1*z1
       deriv(2,10) = 3.0_rp*a*s3*t11*z1
       deriv(3,10) = 3.0_rp*a*s3*t1*z11        
       shape(  11) = 3.0_rp*a*s4*t2*z1
       deriv(1,11) = 3.0_rp*a*s41*t2*z1
       deriv(2,11) = 3.0_rp*a*s4*t21*z1
       deriv(3,11) = 3.0_rp*a*s4*t2*z11        
       shape(  12) = 3.0_rp*a*s4*t3*z1
       deriv(1,12) = 3.0_rp*a*s41*t3*z1        
       deriv(2,12) = 3.0_rp*a*s4*t31*z1        
       deriv(3,12) = 3.0_rp*a*s4*t3*z11                
       shape(  13) = 3.0_rp*a*s3*t4*z1
       deriv(1,13) = 3.0_rp*a*s31*t4*z1
       deriv(2,13) = 3.0_rp*a*s3*t41*z1
       deriv(3,13) = 3.0_rp*a*s3*t4*z11        
       shape(  14) = 3.0_rp*a*s2*t4*z1
       deriv(1,14) = 3.0_rp*a*s21*t4*z1
       deriv(2,14) = 3.0_rp*a*s2*t41*z1
       deriv(3,14) = 3.0_rp*a*s2*t4*z11        
       shape(  15) = 3.0_rp*a*s1*t3*z1
       deriv(1,15) = 3.0_rp*a*s11*t3*z1
       deriv(2,15) = 3.0_rp*a*s1*t31*z1
       deriv(3,15) = 3.0_rp*a*s1*t3*z11        
       shape(  16) = 3.0_rp*a*s1*t2*z1
       deriv(1,16) = 3.0_rp*a*s11*t2*z1
       deriv(2,16) = 3.0_rp*a*s1*t21*z1
       deriv(3,16) = 3.0_rp*a*s1*t2*z11        
       shape(  17) = 3.0_rp*a*s1*t1*z2
       deriv(1,17) = 3.0_rp*a*s11*t1*z2
       deriv(2,17) = 3.0_rp*a*s1*t11*z2
       deriv(3,17) = 3.0_rp*a*s1*t1*z21        
       shape(  18) = 3.0_rp*a*s4*t1*z2
       deriv(1,18) = 3.0_rp*a*s41*t1*z2
       deriv(2,18) = 3.0_rp*a*s4*t11*z2
       deriv(3,18) = 3.0_rp*a*s4*t1*z21        
       shape(  19) = 3.0_rp*a*s4*t4*z2
       deriv(1,19) = 3.0_rp*a*s41*t4*z2
       deriv(2,19) = 3.0_rp*a*s4*t41*z2
       deriv(3,19) = 3.0_rp*a*s4*t4*z21        
       shape(  20) = 3.0_rp*a*s1*t4*z2
       deriv(1,20) = 3.0_rp*a*s11*t4*z2
       deriv(2,20) = 3.0_rp*a*s1*t41*z2
       deriv(3,20) = 3.0_rp*a*s1*t4*z21        
       shape(  21) = 3.0_rp*a*s1*t1*z3
       deriv(1,21) = 3.0_rp*a*s11*t1*z3
       deriv(2,21) = 3.0_rp*a*s1*t11*z3
       deriv(3,21) = 3.0_rp*a*s1*t1*z31        
       shape(  22) = 3.0_rp*a*s4*t1*z3
       deriv(1,22) = 3.0_rp*a*s41*t1*z3
       deriv(2,22) = 3.0_rp*a*s4*t11*z3
       deriv(3,22) = 3.0_rp*a*s4*t1*z31
       shape(  23) = 3.0_rp*a*s4*t4*z3
       deriv(1,23) = 3.0_rp*a*s41*t4*z3
       deriv(2,23) = 3.0_rp*a*s4*t41*z3
       deriv(3,23) = 3.0_rp*a*s4*t4*z31        
       shape(  24) = 3.0_rp*a*s1*t4*z3
       deriv(1,24) = 3.0_rp*a*s11*t4*z3
       deriv(2,24) = 3.0_rp*a*s1*t41*z3
       deriv(3,24) = 3.0_rp*a*s1*t4*z31        
       shape(  25) = 3.0_rp*a*s2*t1*z4
       deriv(1,25) = 3.0_rp*a*s21*t1*z4
       deriv(2,25) = 3.0_rp*a*s2*t11*z4
       deriv(3,25) = 3.0_rp*a*s2*t1*z41        
       shape(  26) = 3.0_rp*a*s3*t1*z4
       deriv(1,21) = 3.0_rp*a*s31*t1*z4
       deriv(2,21) = 3.0_rp*a*s3*t11*z4
       deriv(3,21) = 3.0_rp*a*s3*t1*z41        
       shape(  27) = 3.0_rp*a*s4*t2*z4
       deriv(1,27) = 3.0_rp*a*s41*t2*z4
       deriv(2,27) = 3.0_rp*a*s4*t21*z4
       deriv(3,27) = 3.0_rp*a*s4*t2*z41        
       shape(  28) = 3.0_rp*a*s4*t3*z4
       deriv(1,28) = 3.0_rp*a*s41*t3*z4
       deriv(2,28) = 3.0_rp*a*s4*t31*z4
       deriv(3,28) = 3.0_rp*a*s4*t3*z41        
       shape(  29) = 3.0_rp*a*s3*t4*z4
       deriv(1,29) = 3.0_rp*a*s31*t4*z4
       deriv(2,29) = 3.0_rp*a*s3*t41*z4
       deriv(3,29) = 3.0_rp*a*s3*t4*z41        
       shape(  30) = 3.0_rp*a*s2*t4*z4
       deriv(1,30) = 3.0_rp*a*s21*t4*z4  
       deriv(2,30) = 3.0_rp*a*s2*t41*z4  
       deriv(3,30) = 3.0_rp*a*s2*t4*z41          
       shape(  31) = 3.0_rp*a*s1*t3*z4
       deriv(1,31) = 3.0_rp*a*s11*t3*z4
       deriv(2,31) = 3.0_rp*a*s1*t31*z4
       deriv(3,31) = 3.0_rp*a*s1*t3*z41        
       shape(  32) = 3.0_rp*a*s1*t2*z4
       deriv(1,32) = 3.0_rp*a*s11*t2*z4
       deriv(2,32) = 3.0_rp*a*s1*t21*z4
       deriv(3,32) = 3.0_rp*a*s1*t2*z41        
       shape(  33) =-9.0_rp*a*s2*t2*z1
       deriv(1,33) =-9.0_rp*a*s21*t2*z1
       deriv(2,33) =-9.0_rp*a*s2*t21*z1
       deriv(3,33) =-9.0_rp*a*s2*t2*z11        
       shape(  34) =-9.0_rp*a*s3*t2*z1
       deriv(1,34) =-9.0_rp*a*s31*t2*z1
       deriv(2,34) =-9.0_rp*a*s3*t21*z1
       deriv(3,34) =-9.0_rp*a*s3*t2*z11        
       shape(  35) =-9.0_rp*a*s3*t3*z1
       deriv(1,35) =-9.0_rp*a*s31*t3*z1
       deriv(2,35) =-9.0_rp*a*s3*t31*z1
       deriv(3,35) =-9.0_rp*a*s3*t3*z11        
       shape(  36) =-9.0_rp*a*s2*t3*z1
       deriv(1,36) =-9.0_rp*a*s21*t3*z1
       deriv(2,36) =-9.0_rp*a*s2*t31*z1
       deriv(3,36) =-9.0_rp*a*s2*t3*z11        
       shape(  37) =-9.0_rp*a*s2*t1*z2
       deriv(1,37) =-9.0_rp*a*s21*t1*z2
       deriv(2,37) =-9.0_rp*a*s2*t11*z2
       deriv(3,37) =-9.0_rp*a*s2*t1*z21        
       shape(  38) =-9.0_rp*a*s3*t1*z2
       deriv(1,38) =-9.0_rp*a*s31*t1*z2        
       deriv(2,38) =-9.0_rp*a*s3*t11*z2        
       deriv(3,38) =-9.0_rp*a*s3*t1*z21                
       shape(  39) =-9.0_rp*a*s4*t2*z2
       deriv(1,39) =-9.0_rp*a*s41*t2*z2        
       deriv(2,39) =-9.0_rp*a*s4*t21*z2        
       deriv(3,39) =-9.0_rp*a*s4*t2*z21                
       shape(  40) =-9.0_rp*a*s4*t3*z2
       deriv(1,40) =-9.0_rp*a*s41*t3*z2        
       deriv(2,40) =-9.0_rp*a*s4*t31*z2        
       deriv(3,40) =-9.0_rp*a*s4*t3*z21                
       shape(  41) =-9.0_rp*a*s3*t4*z2
       deriv(1,41) =-9.0_rp*a*s31*t4*z2        
       deriv(2,41) =-9.0_rp*a*s3*t41*z2        
       deriv(3,41) =-9.0_rp*a*s3*t4*z21                
       shape(  42) =-9.0_rp*a*s2*t4*z2
       deriv(1,42) =-9.0_rp*a*s21*t4*z2             
       deriv(2,42) =-9.0_rp*a*s2*t41*z2             
       deriv(3,42) =-9.0_rp*a*s2*t4*z21                     
       shape(  43) =-9.0_rp*a*s1*t3*z2
       deriv(1,43) =-9.0_rp*a*s11*t3*z2             
       deriv(2,43) =-9.0_rp*a*s1*t31*z2             
       deriv(3,43) =-9.0_rp*a*s1*t3*z21                     
       shape(  44) =-9.0_rp*a*s1*t2*z2
       deriv(1,44) =-9.0_rp*a*s11*t2*z2             
       deriv(2,44) =-9.0_rp*a*s1*t21*z2             
       deriv(3,44) =-9.0_rp*a*s1*t2*z21                     
       shape(  45) =-9.0_rp*a*s2*t1*z3
       deriv(1,45) =-9.0_rp*a*s21*t1*z3             
       deriv(2,45) =-9.0_rp*a*s2*t11*z3             
       deriv(3,45) =-9.0_rp*a*s2*t1*z31                     
       shape(  46) =-9.0_rp*a*s3*t1*z3
       deriv(1,46) =-9.0_rp*a*s31*t1*z3
       deriv(2,46) =-9.0_rp*a*s3*t11*z3
       deriv(3,46) =-9.0_rp*a*s3*t1*z31        
       shape(  47) =-9.0_rp*a*s4*t2*z3
       deriv(1,47) =-9.0_rp*a*s41*t2*z3
       deriv(2,47) =-9.0_rp*a*s4*t21*z3
       deriv(3,47) =-9.0_rp*a*s4*t2*z31
       shape(  48) =-9.0_rp*a*s4*t3*z3
       deriv(1,48) =-9.0_rp*a*s41*t3*z3
       deriv(2,48) =-9.0_rp*a*s4*t31*z3
       deriv(3,48) =-9.0_rp*a*s4*t3*z31        
       shape(  49) =-9.0_rp*a*s3*t4*z3
       deriv(1,49) =-9.0_rp*a*s31*t4*z3
       deriv(2,49) =-9.0_rp*a*s3*t41*z3
       deriv(3,49) =-9.0_rp*a*s3*t4*z31        
       shape(  50) =-9.0_rp*a*s2*t4*z3
       deriv(1,50) =-9.0_rp*a*s21*t4*z3
       deriv(2,50) =-9.0_rp*a*s2*t41*z3
       deriv(3,50) =-9.0_rp*a*s2*t4*z31        
       shape(  51) =-9.0_rp*a*s1*t3*z3
       deriv(1,51) =-9.0_rp*a*s11*t3*z3
       deriv(2,51) =-9.0_rp*a*s1*t31*z3
       deriv(3,51) =-9.0_rp*a*s1*t3*z31        
       shape(  52) =-9.0_rp*a*s1*t2*z3
       deriv(1,52) =-9.0_rp*a*s11*t2*z3
       deriv(2,52) =-9.0_rp*a*s1*t21*z3
       deriv(3,52) =-9.0_rp*a*s1*t2*z31        
       shape(  53) =-9.0_rp*a*s2*t2*z4
       deriv(1,53) =-9.0_rp*a*s21*t2*z4
       deriv(2,53) =-9.0_rp*a*s2*t21*z4
       deriv(3,53) =-9.0_rp*a*s2*t2*z41        
       shape(  54) =-9.0_rp*a*s3*t2*z4
       deriv(1,54) =-9.0_rp*a*s31*t2*z4
       deriv(2,54) =-9.0_rp*a*s3*t21*z4
       deriv(3,54) =-9.0_rp*a*s3*t2*z41        
       shape(  55) =-9.0_rp*a*s3*t3*z4
       deriv(1,55) =-9.0_rp*a*s31*t3*z4
       deriv(2,55) =-9.0_rp*a*s3*t31*z4
       deriv(3,55) =-9.0_rp*a*s3*t3*z41        
       shape(  56) =-9.0_rp*a*s2*t3*z4
       deriv(1,56) =-9.0_rp*a*s21*t3*z4
       deriv(2,56) =-9.0_rp*a*s2*t31*z4
       deriv(3,56) =-9.0_rp*a*s2*t3*z41        
       shape(  57) = 27.0_rp*a*s2*t2*z2
       deriv(1,57) = 27.0_rp*a*s21*t2*z2
       deriv(2,57) = 27.0_rp*a*s2*t21*z2
       deriv(3,57) = 27.0_rp*a*s2*t2*z21        
       shape(  58) = 27.0_rp*a*s3*t2*z2
       deriv(1,58) = 27.0_rp*a*s31*t2*z2
       deriv(2,58) = 27.0_rp*a*s3*t21*z2
       deriv(3,58) = 27.0_rp*a*s3*t2*z21        
       shape(  59) = 27.0_rp*a*s3*t3*z2
       deriv(1,59) = 27.0_rp*a*s31*t3*z2
       deriv(2,59) = 27.0_rp*a*s3*t31*z2
       deriv(3,59) = 27.0_rp*a*s3*t3*z21        
       shape(  60) = 27.0_rp*a*s2*t3*z2
       deriv(1,60) = 27.0_rp*a*s21*t3*z2
       deriv(2,60) = 27.0_rp*a*s2*t31*z2
       deriv(3,60) = 27.0_rp*a*s2*t3*z21        
       shape(  61) = 27.0_rp*a*s2*t2*z3
       deriv(1,61) = 27.0_rp*a*s21*t2*z3
       deriv(2,61) = 27.0_rp*a*s2*t21*z3
       deriv(3,61) = 27.0_rp*a*s2*t2*z31        
       shape(  62) = 27.0_rp*a*s3*t2*z3
       deriv(1,62) = 27.0_rp*a*s31*t2*z3
       deriv(2,62) = 27.0_rp*a*s3*t21*z3
       deriv(3,62) = 27.0_rp*a*s3*t2*z31        
       shape(  63) = 27.0_rp*a*s3*t3*z3
       deriv(1,63) = 27.0_rp*a*s31*t3*z3
       deriv(2,63) = 27.0_rp*a*s3*t31*z3
       deriv(3,63) = 27.0_rp*a*s3*t3*z31        
       shape(  64) = 27.0_rp*a*s2*t3*z3
       deriv(1,64) = 27.0_rp*a*s21*t3*z3
       deriv(2,64) = 27.0_rp*a*s2*t31*z3
       deriv(3,64) = 27.0_rp*a*s2*t3*z31
       s12=-2.0_rp*((1.0_rp-s)+(1.0_rp/3.0_rp-s)-(1.0_rp/3.0_rp+s))
       s22=-2.0_rp*((1.0_rp-s)+(1.0_rp/3.0_rp-s)-(1.0_rp+s))
       s32= 2.0_rp*((1.0_rp-s)-(1.0_rp/3.0_rp+s)-(1.0_rp+s))
       s42= 2.0_rp*((1.0_rp/3.0_rp-s)-(1.0_rp/3.0_rp+s)-(1.0_rp+s))
       t12=-2.0_rp*((1.0_rp-t)+(1.0_rp/3.0_rp-t)-(1.0_rp/3.0_rp+t))
       t22=-2.0_rp*((1.0_rp-t)+(1.0_rp/3.0_rp-t)-(1.0_rp+t)) 
       t32= 2.0_rp*((1.0_rp-t)-(1.0_rp/3.0_rp+t)-(1.0_rp+t))  
       t42= 2.0_rp*((1.0_rp/3.0_rp-t)-(1.0_rp/3.0_rp+t)-(1.0_rp+t))
       z12=-2.0_rp*((1.0_rp-z)+(1.0_rp/3.0_rp-z)-(1.0_rp/3.0_rp+z))
       z22=-2.0_rp*((1.0_rp-z)+(1.0_rp/3.0_rp-z)-(1.0_rp+z)) 
       z32= 2.0_rp*((1.0_rp-z)-(1.0_rp/3.0_rp+z)-(1.0_rp+z))  
       z42= 2.0_rp*((1.0_rp/3.0_rp-z)-(1.0_rp/3.0_rp+z)-(1.0_rp+z))

       heslo(1, 1) =   a*s12*t1*z1
       heslo(2, 1) =   a*s1*t12*z1
       heslo(3, 1) =   a*s1*t1*z12
       heslo(4, 1) =   a*s11*t11*z1
       heslo(5, 1) =   a*s11*t1*z11
       heslo(6, 1) =   a*s1*t11*z11

       heslo(1, 2) =   a*s42*t1*z1  
       heslo(2, 2) =   a*s4*t12*z1 
       heslo(3, 2) =   a*s4*t1*z12 
       heslo(4, 2) =   a*s41*t11*z1
       heslo(5, 2) =   a*s41*t1*z11
       heslo(6, 2) =   a*s4*t11*z11          

       heslo(1, 3) =   a*s42*t4*z1           
       heslo(2, 3) =   a*s4*t42*z1  
       heslo(3, 3) =   a*s4*t4*z12      
       heslo(4, 3) =   a*s41*t41*z1
       heslo(5, 3) =   a*s41*t4*z11
       heslo(6, 3) =   a*s4*t41*z11      

       heslo(1, 4) =   a*s12*t4*z1 
       heslo(2, 4) =   a*s1*t42*z1 
       heslo(3, 4) =   a*s1*t4*z12 
       heslo(4, 4) =   a*s11*t41*z1
       heslo(5, 4) =   a*s11*t4*z11
       heslo(6, 4) =   a*s1*t41*z11

       heslo(1, 5) =   a*s12*t1*z4 
       heslo(2, 5) =   a*s1*t12*z4 
       heslo(3, 5) =   a*s1*t1*z42 
       heslo(4, 5) =   a*s11*t11*z4
       heslo(5, 5) =   a*s11*t1*z41
       heslo(6, 5) =   a*s1*t11*z41

       heslo(1, 6) =   a*s42*t4*z4 
       heslo(2, 6) =   a*s4*t42*z4 
       heslo(3, 6) =   a*s4*t4*z42 
       heslo(4, 6) =   a*s41*t41*z4
       heslo(5, 6) =   a*s41*t4*z41
       heslo(6, 6) =   a*s4*t41*z41

       heslo(1, 7) =   a*s42*t4*z4 
       heslo(2, 7) =   a*s4*t42*z4 
       heslo(3, 7) =   a*s4*t4*z42 
       heslo(4, 7) =   a*s41*t41*z4
       heslo(5, 7) =   a*s41*t4*z41
       heslo(6, 7) =   a*s4*t41*z41

       heslo(1, 8) =   a*s12*t4*z4 
       heslo(2, 8) =   a*s1*t42*z4 
       heslo(3, 8) =   a*s1*t4*z42 
       heslo(4, 8) =   a*s11*t41*z4
       heslo(5, 8) =   a*s11*t4*z41
       heslo(6, 8) =   a*s1*t41*z41

       heslo(1, 9) = 3.0_rp*a*s22*t1*z1 
       heslo(2, 9) = 3.0_rp*a*s2*t12*z1 
       heslo(3, 9) = 3.0_rp*a*s2*t1*z12 
       heslo(4, 9) = 3.0_rp*a*s21*t11*z1
       heslo(5, 9) = 3.0_rp*a*s21*t1*z11
       heslo(6, 9) = 3.0_rp*a*s2*t11*z11

       heslo(1,10) = 3.0_rp*a*s32*t1*z1 
       heslo(2,10) = 3.0_rp*a*s3*t12*z1 
       heslo(3,10) = 3.0_rp*a*s3*t1*z12 
       heslo(4,10) = 3.0_rp*a*s31*t11*z1
       heslo(5,10) = 3.0_rp*a*s31*t1*z11
       heslo(6,10) = 3.0_rp*a*s3*t11*z11

       heslo(1,11) = 3.0_rp*a*s42*t2*z1 
       heslo(2,11) = 3.0_rp*a*s4*t22*z1 
       heslo(3,11) = 3.0_rp*a*s4*t2*z12 
       heslo(4,11) = 3.0_rp*a*s41*t21*z1
       heslo(5,11) = 3.0_rp*a*s41*t2*z11
       heslo(6,11) = 3.0_rp*a*s4*t21*z11

       heslo(1,12) = 3.0_rp*a*s42*t3*z1 
       heslo(2,12) = 3.0_rp*a*s4*t32*z1 
       heslo(3,12) = 3.0_rp*a*s4*t3*z12 
       heslo(4,12) = 3.0_rp*a*s41*t31*z1
       heslo(5,12) = 3.0_rp*a*s41*t3*z11
       heslo(6,12) = 3.0_rp*a*s4*t31*z11

       heslo(1,13) = 3.0_rp*a*s32*t4*z1 
       heslo(2,13) = 3.0_rp*a*s3*t42*z1 
       heslo(3,13) = 3.0_rp*a*s3*t4*z12 
       heslo(4,13) = 3.0_rp*a*s31*t41*z1
       heslo(5,13) = 3.0_rp*a*s31*t4*z11
       heslo(6,13) = 3.0_rp*a*s3*t41*z11

       heslo(1,14) = 3.0_rp*a*s22*t4*z1 
       heslo(2,14) = 3.0_rp*a*s2*t42*z1 
       heslo(3,14) = 3.0_rp*a*s2*t4*z12 
       heslo(4,14) = 3.0_rp*a*s21*t41*z1
       heslo(5,14) = 3.0_rp*a*s21*t4*z11
       heslo(6,14) = 3.0_rp*a*s2*t41*z11

       heslo(1,15) = 3.0_rp*a*s12*t3*z1 
       heslo(2,15) = 3.0_rp*a*s1*t32*z1 
       heslo(3,15) = 3.0_rp*a*s1*t3*z12 
       heslo(4,15) = 3.0_rp*a*s11*t31*z1
       heslo(5,15) = 3.0_rp*a*s11*t3*z11
       heslo(6,15) = 3.0_rp*a*s1*t31*z11

       heslo(1,16) = 3.0_rp*a*s12*t2*z1 
       heslo(2,16) = 3.0_rp*a*s1*t22*z1 
       heslo(3,16) = 3.0_rp*a*s1*t2*z12 
       heslo(4,16) = 3.0_rp*a*s11*t21*z1
       heslo(5,16) = 3.0_rp*a*s11*t2*z11
       heslo(6,16) = 3.0_rp*a*s1*t21*z11

       heslo(1,17) = 3.0_rp*a*s12*t1*z2 
       heslo(2,17) = 3.0_rp*a*s1*t12*z2 
       heslo(3,17) = 3.0_rp*a*s1*t1*z22 
       heslo(4,17) = 3.0_rp*a*s11*t11*z2
       heslo(5,17) = 3.0_rp*a*s11*t1*z21
       heslo(6,17) = 3.0_rp*a*s1*t11*z21

       heslo(1,18) = 3.0_rp*a*s42*t1*z2 
       heslo(2,18) = 3.0_rp*a*s4*t12*z2 
       heslo(3,18) = 3.0_rp*a*s4*t1*z22 
       heslo(4,18) = 3.0_rp*a*s41*t11*z2
       heslo(5,18) = 3.0_rp*a*s41*t1*z21
       heslo(6,18) = 3.0_rp*a*s4*t11*z21

       heslo(1,19) = 3.0_rp*a*s42*t4*z2 
       heslo(2,19) = 3.0_rp*a*s4*t42*z2 
       heslo(3,19) = 3.0_rp*a*s4*t4*z22 
       heslo(4,19) = 3.0_rp*a*s41*t41*z2
       heslo(5,19) = 3.0_rp*a*s41*t4*z21
       heslo(6,19) = 3.0_rp*a*s4*t41*z21

       heslo(1,20) = 3.0_rp*a*s12*t4*z2 
       heslo(2,20) = 3.0_rp*a*s1*t42*z2 
       heslo(3,20) = 3.0_rp*a*s1*t4*z22 
       heslo(4,20) = 3.0_rp*a*s11*t41*z2
       heslo(5,20) = 3.0_rp*a*s11*t4*z21
       heslo(6,20) = 3.0_rp*a*s1*t41*z21

       heslo(1,21) = 3.0_rp*a*s12*t1*z3 
       heslo(2,21) = 3.0_rp*a*s1*t12*z3 
       heslo(3,21) = 3.0_rp*a*s1*t1*z32 
       heslo(4,21) = 3.0_rp*a*s11*t11*z3
       heslo(5,21) = 3.0_rp*a*s11*t1*z31
       heslo(6,21) = 3.0_rp*a*s1*t11*z31

       heslo(1,22) = 3.0_rp*a*s42*t1*z3 
       heslo(2,22) = 3.0_rp*a*s4*t12*z3 
       heslo(3,22) = 3.0_rp*a*s4*t1*z32 
       heslo(4,22) = 3.0_rp*a*s41*t11*z3
       heslo(5,22) = 3.0_rp*a*s41*t1*z31
       heslo(6,22) = 3.0_rp*a*s4*t11*z31

       heslo(1,23) = 3.0_rp*a*s42*t4*z3 
       heslo(2,23) = 3.0_rp*a*s4*t42*z3 
       heslo(3,23) = 3.0_rp*a*s4*t4*z32 
       heslo(4,23) = 3.0_rp*a*s41*t41*z3
       heslo(5,23) = 3.0_rp*a*s41*t4*z31
       heslo(6,23) = 3.0_rp*a*s4*t41*z31

       heslo(1,24) = 3.0_rp*a*s12*t4*z3 
       heslo(2,24) = 3.0_rp*a*s1*t42*z3 
       heslo(3,24) = 3.0_rp*a*s1*t4*z32 
       heslo(4,24) = 3.0_rp*a*s11*t41*z3
       heslo(5,24) = 3.0_rp*a*s11*t4*z31
       heslo(6,24) = 3.0_rp*a*s1*t41*z31

       heslo(1,25) = 3.0_rp*a*s22*t4*z4 
       heslo(2,25) = 3.0_rp*a*s2*t42*z4 
       heslo(3,25) = 3.0_rp*a*s2*t4*z42 
       heslo(4,25) = 3.0_rp*a*s21*t41*z4
       heslo(5,25) = 3.0_rp*a*s21*t4*z41
       heslo(6,25) = 3.0_rp*a*s2*t41*z41

       heslo(1,26) = 3.0_rp*a*s32*t1*z4 
       heslo(2,26) = 3.0_rp*a*s3*t12*z4 
       heslo(3,26) = 3.0_rp*a*s3*t1*z42 
       heslo(4,26) = 3.0_rp*a*s31*t11*z4
       heslo(5,26) = 3.0_rp*a*s31*t1*z41
       heslo(6,26) = 3.0_rp*a*s3*t11*z41

       heslo(1,27) = 3.0_rp*a*s42*t2*z4 
       heslo(2,27) = 3.0_rp*a*s4*t22*z4 
       heslo(3,27) = 3.0_rp*a*s4*t2*z42 
       heslo(4,27) = 3.0_rp*a*s41*t21*z4
       heslo(5,27) = 3.0_rp*a*s41*t2*z41
       heslo(6,27) = 3.0_rp*a*s4*t21*z41

       heslo(1,28) = 3.0_rp*a*s42*t3*z4 
       heslo(2,28) = 3.0_rp*a*s4*t32*z4 
       heslo(3,28) = 3.0_rp*a*s4*t3*z42 
       heslo(4,28) = 3.0_rp*a*s41*t31*z4
       heslo(5,28) = 3.0_rp*a*s41*t3*z41
       heslo(6,28) = 3.0_rp*a*s4*t31*z41

       heslo(1,29) = 3.0_rp*a*s32*t4*z4 
       heslo(2,29) = 3.0_rp*a*s3*t42*z4 
       heslo(3,29) = 3.0_rp*a*s3*t4*z42 
       heslo(4,29) = 3.0_rp*a*s31*t41*z4
       heslo(5,29) = 3.0_rp*a*s31*t4*z41
       heslo(6,29) = 3.0_rp*a*s3*t41*z41

       heslo(1,30) = 3.0_rp*a*s22*t4*z4 
       heslo(2,30) = 3.0_rp*a*s2*t42*z4 
       heslo(3,30) = 3.0_rp*a*s2*t4*z42 
       heslo(4,30) = 3.0_rp*a*s21*t41*z4
       heslo(5,30) = 3.0_rp*a*s21*t4*z41
       heslo(6,30) = 3.0_rp*a*s2*t41*z31

       heslo(1,31) = 3.0_rp*a*s12*t3*z4 
       heslo(2,31) = 3.0_rp*a*s1*t32*z4 
       heslo(3,31) = 3.0_rp*a*s1*t3*z42 
       heslo(4,31) = 3.0_rp*a*s11*t31*z4
       heslo(5,31) = 3.0_rp*a*s11*t3*z41
       heslo(6,31) = 3.0_rp*a*s1*t31*z41

       heslo(1,32) = 3.0_rp*a*s12*t2*z4 
       heslo(2,32) = 3.0_rp*a*s1*t22*z4 
       heslo(3,32) = 3.0_rp*a*s1*t2*z42 
       heslo(4,32) = 3.0_rp*a*s11*t21*z4
       heslo(5,32) = 3.0_rp*a*s11*t2*z41
       heslo(6,32) = 3.0_rp*a*s1*t21*z41

       heslo(1,33) =-9.0_rp*a*s22*t2*z1 
       heslo(2,33) =-9.0_rp*a*s2*t22*z1 
       heslo(3,33) =-9.0_rp*a*s2*t2*z12 
       heslo(4,33) =-9.0_rp*a*s21*t21*z1
       heslo(5,33) =-9.0_rp*a*s21*t2*z11
       heslo(6,33) =-9.0_rp*a*s2*t21*z11

       heslo(1,34) =-9.0_rp*a*s32*t2*z1 
       heslo(2,34) =-9.0_rp*a*s3*t22*z1 
       heslo(3,34) =-9.0_rp*a*s3*t2*z12 
       heslo(4,34) =-9.0_rp*a*s31*t21*z1
       heslo(5,34) =-9.0_rp*a*s31*t2*z11
       heslo(6,34) =-9.0_rp*a*s3*t21*z11

       heslo(1,35) =-9.0_rp*a*s32*t3*z1 
       heslo(2,35) =-9.0_rp*a*s3*t32*z1 
       heslo(3,35) =-9.0_rp*a*s3*t3*z12 
       heslo(4,35) =-9.0_rp*a*s31*t31*z1
       heslo(5,35) =-9.0_rp*a*s31*t3*z11
       heslo(6,35) =-9.0_rp*a*s3*t31*z11

       heslo(1,36) =-9.0_rp*a*s22*t3*z1 
       heslo(2,36) =-9.0_rp*a*s2*t32*z1 
       heslo(3,36) =-9.0_rp*a*s2*t3*z12 
       heslo(4,36) =-9.0_rp*a*s21*t31*z1
       heslo(5,36) =-9.0_rp*a*s21*t3*z11
       heslo(6,36) =-9.0_rp*a*s2*t31*z11

       heslo(1,37) =-9.0_rp*a*s22*t1*z2 
       heslo(2,37) =-9.0_rp*a*s2*t12*z2 
       heslo(3,37) =-9.0_rp*a*s2*t1*z22 
       heslo(4,37) =-9.0_rp*a*s21*t11*z2
       heslo(5,37) =-9.0_rp*a*s21*t1*z21
       heslo(6,37) =-9.0_rp*a*s2*t11*z21

       heslo(1,38) =-9.0_rp*a*s32*t1*z2 
       heslo(2,38) =-9.0_rp*a*s3*t12*z2 
       heslo(3,38) =-9.0_rp*a*s3*t1*z22 
       heslo(4,38) =-9.0_rp*a*s31*t11*z2
       heslo(5,38) =-9.0_rp*a*s31*t1*z21
       heslo(6,38) =-9.0_rp*a*s3*t11*z21

       heslo(1,39) =-9.0_rp*a*s42*t2*z2 
       heslo(2,39) =-9.0_rp*a*s4*t22*z2 
       heslo(3,39) =-9.0_rp*a*s4*t2*z22 
       heslo(4,39) =-9.0_rp*a*s41*t21*z2
       heslo(5,39) =-9.0_rp*a*s41*t2*z21
       heslo(6,39) =-9.0_rp*a*s4*t21*z21

       heslo(1,40) =-9.0_rp*a*s42*t3*z2 
       heslo(2,40) =-9.0_rp*a*s4*t32*z2 
       heslo(3,40) =-9.0_rp*a*s4*t3*z22 
       heslo(4,40) =-9.0_rp*a*s41*t31*z2
       heslo(5,40) =-9.0_rp*a*s41*t3*z21
       heslo(6,40) =-9.0_rp*a*s4*t31*z21

       heslo(1,41) =-9.0_rp*a*s32*t4*z2 
       heslo(2,41) =-9.0_rp*a*s3*t42*z2 
       heslo(3,41) =-9.0_rp*a*s3*t4*z22 
       heslo(4,41) =-9.0_rp*a*s31*t41*z2
       heslo(5,41) =-9.0_rp*a*s31*t4*z21
       heslo(6,41) =-9.0_rp*a*s3*t41*z21

       heslo(1,42) =-9.0_rp*a*s22*t4*z2 
       heslo(2,42) =-9.0_rp*a*s2*t42*z2 
       heslo(3,42) =-9.0_rp*a*s2*t4*z22 
       heslo(4,42) =-9.0_rp*a*s21*t41*z2
       heslo(5,42) =-9.0_rp*a*s21*t4*z21
       heslo(6,42) =-9.0_rp*a*s2*t41*z21

       heslo(1,43) =-9.0_rp*a*s12*t3*z2 
       heslo(2,43) =-9.0_rp*a*s1*t32*z2 
       heslo(3,43) =-9.0_rp*a*s1*t3*z22 
       heslo(4,43) =-9.0_rp*a*s11*t31*z2
       heslo(5,43) =-9.0_rp*a*s11*t3*z21
       heslo(6,43) =-9.0_rp*a*s1*t31*z21

       heslo(1,44) =-9.0_rp*a*s12*t2*z2 
       heslo(2,44) =-9.0_rp*a*s1*t22*z2 
       heslo(3,44) =-9.0_rp*a*s1*t2*z22 
       heslo(4,44) =-9.0_rp*a*s11*t21*z2
       heslo(5,44) =-9.0_rp*a*s11*t2*z21
       heslo(6,44) =-9.0_rp*a*s1*t21*z21

       heslo(1,45) =-9.0_rp*a*s22*t1*z3 
       heslo(2,45) =-9.0_rp*a*s2*t12*z3 
       heslo(3,45) =-9.0_rp*a*s2*t1*z32 
       heslo(4,45) =-9.0_rp*a*s21*t11*z3
       heslo(5,45) =-9.0_rp*a*s21*t1*z31
       heslo(6,45) =-9.0_rp*a*s2*t11*z31

       heslo(1,46) =-9.0_rp*a*s32*t1*z3 
       heslo(2,46) =-9.0_rp*a*s3*t12*z3 
       heslo(3,46) =-9.0_rp*a*s3*t1*z32 
       heslo(4,46) =-9.0_rp*a*s31*t11*z3
       heslo(5,46) =-9.0_rp*a*s31*t1*z31
       heslo(6,46) =-9.0_rp*a*s3*t11*z31

       heslo(1,47) =-9.0_rp*a*s42*t2*z3 
       heslo(2,47) =-9.0_rp*a*s4*t22*z3 
       heslo(3,47) =-9.0_rp*a*s4*t2*z32 
       heslo(4,47) =-9.0_rp*a*s41*t21*z3
       heslo(5,47) =-9.0_rp*a*s41*t2*z31
       heslo(6,47) =-9.0_rp*a*s4*t21*z31

       heslo(1,48) =-9.0_rp*a*s42*t3*z3 
       heslo(2,48) =-9.0_rp*a*s4*t32*z3 
       heslo(3,48) =-9.0_rp*a*s4*t3*z32 
       heslo(4,48) =-9.0_rp*a*s41*t31*z3
       heslo(5,48) =-9.0_rp*a*s41*t3*z31
       heslo(6,48) =-9.0_rp*a*s4*t31*z31

       heslo(1,49) =-9.0_rp*a*s32*t4*z3 
       heslo(2,49) =-9.0_rp*a*s3*t42*z3 
       heslo(3,49) =-9.0_rp*a*s3*t4*z32 
       heslo(4,49) =-9.0_rp*a*s31*t41*z3
       heslo(5,49) =-9.0_rp*a*s31*t4*z31
       heslo(6,49) =-9.0_rp*a*s3*t41*z31

       heslo(1,50) =-9.0_rp*a*s22*t4*z3 
       heslo(2,50) =-9.0_rp*a*s2*t42*z3 
       heslo(3,50) =-9.0_rp*a*s2*t4*z32 
       heslo(4,50) =-9.0_rp*a*s21*t41*z3
       heslo(5,50) =-9.0_rp*a*s21*t4*z31
       heslo(6,50) =-9.0_rp*a*s2*t41*z31

       heslo(1,51) =-9.0_rp*a*s12*t3*z3 
       heslo(2,51) =-9.0_rp*a*s1*t32*z3 
       heslo(3,51) =-9.0_rp*a*s1*t3*z32 
       heslo(4,51) =-9.0_rp*a*s11*t31*z3
       heslo(5,51) =-9.0_rp*a*s11*t3*z31
       heslo(6,51) =-9.0_rp*a*s1*t31*z31

       heslo(1,52) =-9.0_rp*a*s12*t2*z3 
       heslo(2,52) =-9.0_rp*a*s1*t22*z3 
       heslo(3,52) =-9.0_rp*a*s1*t2*z32 
       heslo(4,52) =-9.0_rp*a*s11*t21*z3
       heslo(5,52) =-9.0_rp*a*s11*t2*z31
       heslo(6,52) =-9.0_rp*a*s1*t21*z31

       heslo(1,53) =-9.0_rp*a*s22*t2*z4 
       heslo(2,53) =-9.0_rp*a*s2*t22*z4 
       heslo(3,53) =-9.0_rp*a*s2*t2*z42 
       heslo(4,53) =-9.0_rp*a*s21*t21*z4
       heslo(5,53) =-9.0_rp*a*s21*t2*z41
       heslo(6,53) =-9.0_rp*a*s2*t21*z41

       heslo(1,54) =-9.0_rp*a*s32*t2*z4 
       heslo(2,54) =-9.0_rp*a*s3*t22*z4 
       heslo(3,54) =-9.0_rp*a*s3*t2*z42 
       heslo(4,54) =-9.0_rp*a*s31*t21*z4
       heslo(5,54) =-9.0_rp*a*s31*t2*z41
       heslo(6,54) =-9.0_rp*a*s3*t21*z41

       heslo(1,55) =-9.0_rp*a*s32*t3*z4 
       heslo(2,55) =-9.0_rp*a*s3*t32*z4 
       heslo(3,55) =-9.0_rp*a*s3*t3*z42 
       heslo(4,55) =-9.0_rp*a*s31*t31*z4
       heslo(5,55) =-9.0_rp*a*s31*t3*z41
       heslo(6,55) =-9.0_rp*a*s3*t31*z41

       heslo(1,56) =-9.0_rp*a*s22*t3*z4 
       heslo(2,56) =-9.0_rp*a*s2*t32*z4 
       heslo(3,56) =-9.0_rp*a*s2*t3*z42 
       heslo(4,56) =-9.0_rp*a*s21*t31*z4
       heslo(5,56) =-9.0_rp*a*s21*t3*z41
       heslo(6,56) =-9.0_rp*a*s2*t31*z41

       heslo(1,57) = 27.0_rp*a*s22*t2*z2 
       heslo(2,57) = 27.0_rp*a*s2*t22*z2 
       heslo(3,57) = 27.0_rp*a*s2*t2*z22 
       heslo(4,57) = 27.0_rp*a*s21*t21*z2
       heslo(5,57) = 27.0_rp*a*s21*t2*z21
       heslo(6,57) = 27.0_rp*a*s2*t21*z21

       heslo(1,58) = 27.0_rp*a*s32*t2*z2 
       heslo(2,58) = 27.0_rp*a*s3*t22*z2 
       heslo(3,58) = 27.0_rp*a*s3*t2*z22 
       heslo(4,58) = 27.0_rp*a*s31*t21*z2
       heslo(5,58) = 27.0_rp*a*s31*t2*z21
       heslo(6,58) = 27.0_rp*a*s3*t21*z21

       heslo(1,59) = 27.0_rp*a*s32*t3*z2 
       heslo(2,59) = 27.0_rp*a*s3*t32*z2 
       heslo(3,59) = 27.0_rp*a*s3*t3*z22 
       heslo(4,59) = 27.0_rp*a*s31*t31*z2
       heslo(5,59) = 27.0_rp*a*s31*t3*z21
       heslo(6,59) = 27.0_rp*a*s3*t31*z21

       heslo(1,60) = 27.0_rp*a*s22*t3*z2 
       heslo(2,60) = 27.0_rp*a*s2*t32*z2 
       heslo(3,60) = 27.0_rp*a*s2*t3*z22 
       heslo(4,60) = 27.0_rp*a*s21*t31*z2
       heslo(5,60) = 27.0_rp*a*s21*t3*z21
       heslo(6,60) = 27.0_rp*a*s2*t31*z21

       heslo(1,61) = 27.0_rp*a*s22*t2*z3 
       heslo(2,61) = 27.0_rp*a*s2*t22*z3 
       heslo(3,61) = 27.0_rp*a*s2*t2*z32 
       heslo(4,61) = 27.0_rp*a*s21*t21*z3
       heslo(5,61) = 27.0_rp*a*s21*t2*z31
       heslo(6,61) = 27.0_rp*a*s2*t21*z31

       heslo(1,62) = 27.0_rp*a*s32*t2*z3 
       heslo(2,62) = 27.0_rp*a*s3*t22*z3 
       heslo(3,62) = 27.0_rp*a*s3*t2*z32 
       heslo(4,62) = 27.0_rp*a*s31*t21*z3
       heslo(5,62) = 27.0_rp*a*s31*t2*z31
       heslo(6,62) = 27.0_rp*a*s3*t21*z31

       heslo(1,63) = 27.0_rp*a*s32*t3*z3 
       heslo(2,63) = 27.0_rp*a*s3*t32*z3 
       heslo(3,63) = 27.0_rp*a*s3*t3*z32 
       heslo(4,63) = 27.0_rp*a*s31*t31*z3
       heslo(5,63) = 27.0_rp*a*s31*t3*z31
       heslo(6,63) = 27.0_rp*a*s3*t31*z31

       heslo(1,64) = 27.0_rp*a*s22*t3*z3 
       heslo(2,64) = 27.0_rp*a*s2*t32*z3 
       heslo(3,64) = 27.0_rp*a*s2*t3*z32 
       heslo(4,64) = 27.0_rp*a*s21*t31*z3
       heslo(5,64) = 27.0_rp*a*s21*t3*z31
       heslo(6,64) = 27.0_rp*a*s2*t31*z31
!
! Linear Prism
!
  else if(nnode==6) then
     shape(   1) = (1.0_rp-s-t)*(1.0_rp-z)
     deriv(1, 1) = z-1.0_rp
     deriv(2, 1) = z-1.0_rp
     deriv(3, 1) = s+t-1.0_rp
     shape(   2) = s*(1.0_rp-z)
     deriv(1, 2) = 1.0_rp-z
     deriv(3, 2) = -s
     shape(   3) = t*(1.0_rp-z)
     deriv(2, 3) = 1.0_rp-z
     deriv(3, 3) = -t
     shape(   4) = (1.0_rp-s-t)*z
     deriv(1, 4) = -z
     deriv(2, 4) = -z
     deriv(3, 4) = 1.0_rp-s-t
     shape(   5) = s*z
     deriv(1, 5) = z
     deriv(3, 5) = s
     shape(   6) = t*z
     deriv(2, 6) = z
     deriv(3, 6) = t
     
  else
     istop = 1
  end if
!
! Errors
!      
      if(istop==1) call runend('SHAPE3: NOT AVAILABLE ELEMENT INERPOLATION')
      
    end subroutine shape3

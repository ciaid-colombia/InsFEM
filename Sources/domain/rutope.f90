
subroutine rutope(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
! 
!     This routine sets up the integration constants of open rules for
!     triangles and tetrahedra
! 
!             NDIME = 2             NDIME = 3
! 
!          NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!          -----  ----------     -----  ----------
!            1       p1            1       p1
!            3       p2            4       p2
!            4       p3            5       p3
!            6       p4           11       p4
!            7       p5           14       p5
!           13       p9
! 
!-----------------------------------------------------------------------
  
  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: a,b,c,d,e,f,g,h,p,q,r,w1,w2,w3,w4
  real(rp)                 :: ex1,et1,ez1,ex2,et2,ez2
  integer(ip)              :: istop
  
  istop=0
!
! Line integral (the same as for brick elements)
!
  if(ndime==1) then
     if(ngaus==1) then
        posgp(1,1)=0.0_rp
        weigp(  1)=2.0_rp
     else if(ngaus==2) then
        posgp(1,1)=-0.577350269189626_rp
        posgp(1,2)= 0.577350269189626_rp
        weigp(  1)= 1.0_rp
        weigp(  2)= 1.0_rp
     else if(ngaus==3) then
        posgp(1,1)=-0.774596669241483_rp
        posgp(1,2)= 0.0_rp
        posgp(1,3)= 0.774596669241483_rp
        weigp(  1)= 0.555555555555556_rp
        weigp(  2)= 0.888888888888889_rp
        weigp(  3)= 0.555555555555556_rp
     else if(ngaus==4) then
        posgp(1,1)=-0.861136311594053_rp
        posgp(1,2)=-0.339981043584856_rp
        posgp(1,3)= 0.339981043584856_rp
        posgp(1,4)= 0.861136311594053_rp
        weigp(1)= 0.347854845137454_rp
        weigp(2)= 0.652145154862546_rp
        weigp(3)= 0.652145154862546_rp
        weigp(4)= 0.347854845137454_rp  
!* * * * * * * * * * * * * * * * * * * * * * * * * *
!* 
!* Valores añadidos para ngaus==5 tomados de ruqope
!*
!* * * * * * * * * * * * * * * * * * * * * * * * * *
     else if(ngaus==5) then
            posgp(1,1) = -0.906179845938664_rp
            posgp(1,2) = -0.538469310105683_rp
            posgp(1,3) =  0.0_rp
            posgp(1,4) =  0.538469310105683_rp
            posgp(1,5) =  0.906179845938664_rp
            weigp(1) =  0.236926885056189_rp
            weigp(2) =  0.478628670499366_rp
            weigp(3) =  0.568888888888889_rp
            weigp(4) =  0.478628670499366_rp
            weigp(5) =  0.236926885056189_rp   
     else
        istop=1
     end if
     
!
! Area integral (triangles)
!
  else if(ndime==2) then
     if(ngaus==1) then
        posgp(1,1)= 1.0_rp/3.0_rp
        posgp(2,1)= 1.0_rp/3.0_rp
        weigp(  1)= 1.0_rp/2.0_rp
       else if(ngaus==3) then   ! abierta para p2
        posgp(1,1)= 2.0_rp/3.0_rp
        posgp(2,1)= 1.0_rp/6.0_rp
        posgp(1,2)= 1.0_rp/6.0_rp
        posgp(2,2)= 2.0_rp/3.0_rp
        posgp(1,3)= 1.0_rp/6.0_rp
        posgp(2,3)= 1.0_rp/6.0_rp
        weigp(  1)= 1.0_rp/6.0_rp
        weigp(  2)= 1.0_rp/6.0_rp
        weigp(  3)= 1.0_rp/6.0_rp  
      else if(ngaus==4) then   !abierta para p3
        posgp(1,1)= 1.0_rp/3.0_rp
        posgp(2,1)= 1.0_rp/3.0_rp
        posgp(1,2)= 1.0_rp/5.0_rp
        posgp(2,2)= 1.0_rp/5.0_rp
        posgp(1,3)= 3.0_rp/5.0_rp
        posgp(2,3)= 1.0_rp/5.0_rp
        posgp(1,4)= 1.0_rp/5.0_rp
        posgp(2,4)= 3.0_rp/5.0_rp
        weigp(  1)=-27.0_rp/96.0_rp
        weigp(  2)= 25.0_rp/96.0_rp
        weigp(  3)= 25.0_rp/96.0_rp
        weigp(  4)= 25.0_rp/96.0_rp 
     else if(ngaus==6) then   !abierta para p4
        ex1 = 0.816847572980459_rp
        et1 = 0.091576213509771_rp
        ez1 = 0.091576213509771_rp
        ex2 = 0.108103018168070_rp
        et2 = 0.445948490915965_rp
        ez2 = 0.445948490915965_rp
        posgp(1,1)= ex1
        posgp(2,1)= et1
        posgp(1,2)= et1
        posgp(2,2)= ez1
        posgp(1,3)= ez1
        posgp(2,3)= ex1
        posgp(1,4)= ex2
        posgp(2,4)= et2
        posgp(1,5)= et2
        posgp(2,5)= ez2
        posgp(1,6)= ez2
        posgp(2,6)= ex2
        a = 0.054975870996713638_rp
        b = 0.1116907969117165_rp    
        weigp(1)  = a
        weigp(2)  = a
        weigp(3)  = a
        weigp(4)  = b
        weigp(5)  = b
        weigp(6)  = b
     else if(ngaus==7) then !abierta para p5
        e = 1.0_rp / 3.0_rp
        a = ( 6.0_rp +          sqrt ( 15.0_rp ) ) / 21.0_rp
        b = ( 6.0_rp -          sqrt ( 15.0_rp ) ) / 21.0_rp
        c = ( 9.0_rp + 2.0_rp * sqrt ( 15.0_rp ) ) / 21.0_rp
        d = ( 9.0_rp - 2.0_rp * sqrt ( 15.0_rp ) ) / 21.0_rp
        w1 = 0.1125_rp
        w2 = ( 155.0_rp + sqrt ( 15.0_rp ) ) / 2400.0_rp
        w3 = ( 155.0_rp - sqrt ( 15.0_rp ) ) / 2400.0_rp
        posgp(1,1)= e
        posgp(2,1)= e
        posgp(1,2)= a
        posgp(2,2)= a
        posgp(1,3)= d
        posgp(2,3)= a
        posgp(1,4)= a
        posgp(2,4)= d
        posgp(1,5)= b
        posgp(2,5)= b
        posgp(1,6)= c
        posgp(2,6)= b
        posgp(1,7)= b
        posgp(2,7)= c
        weigp(  1)= w1
        weigp(  2)= w2
        weigp(  3)= w2
        weigp(  4)= w2
        weigp(  5)= w3
        weigp(  6)= w3
        weigp(  7)= w3
     else if(ngaus==12) then   !abierta para p6
        posgp(1,1) =  -0.501426509658180_rp               
        posgp(1,2) =  -0.501426509658180_rp               
        posgp(1,3) =   0.002853019316358_rp               
        posgp(1,4) =  -0.873821971016996_rp               
        posgp(1,5) =  -0.873821971016996_rp               
        posgp(1,6) =   0.747643942033992_rp               
        posgp(1,7) =  -0.379295097932432_rp               
        posgp(1,8) =   0.273004998242798_rp               
        posgp(1,9) =  -0.893709900310366_rp               
        posgp(1,10)=  -0.379295097932432_rp               
        posgp(1,11)=   0.273004998242798_rp               
        posgp(1,12)=  -0.893709900310366_rp               
        posgp(2,1) =  -0.501426509658180_rp               
        posgp(2,2) =   0.002853019316358_rp               
        posgp(2,3) =  -0.501426509658180_rp               
        posgp(2,4) =  -0.873821971016996_rp               
        posgp(2,5) =   0.747643942033992_rp               
        posgp(2,6) =  -0.873821971016996_rp               
        posgp(2,7) =   0.273004998242798_rp               
        posgp(2,8) =  -0.893709900310366_rp               
        posgp(2,9) =  -0.379295097932432_rp               
        posgp(2,10)=  -0.893709900310366_rp               
        posgp(2,11)=  -0.379295097932432_rp               
        posgp(2,12)=   0.273004998242798_rp               
        weigp( 1) = 0.233572551452758*0.25_rp
        weigp( 2) = 0.233572551452758*0.25_rp  
        weigp( 3) = 0.233572551452758*0.25_rp
        weigp( 4) = 0.101689812740414*0.25_rp  
        weigp( 5) = 0.101689812740414*0.25_rp  
        weigp( 6) = 0.101689812740414*0.25_rp  
        weigp( 7) = 0.165702151236748*0.25_rp  
        weigp( 8) = 0.165702151236748*0.25_rp  
        weigp( 9) = 0.165702151236748*0.25_rp  
        weigp(10) = 0.165702151236748*0.25_rp  
        weigp(11) = 0.165702151236748*0.25_rp  
        weigp(12) = 0.165702151236748*0.25_rp  
!        do igaus =1, ngaus
        posgp = 0.5_rp*(posgp + 1.0_rp)
!        end do
     else if(ngaus==13) then  ! abierta para p7
        a = 0.333333333333333_rp
        b = 0.479308067841920_rp
        c = 0.869739794195568_rp
        d = 0.638444188569810_rp
        e = 0.260345966079040_rp
        f = 0.065130102902216_rp
        g = 0.312865496004874_rp
        h = 0.048690315425316_rp
        w1=-0.149570044467670_rp/2.0_rp
        w2= 0.175615257433204_rp/2.0_rp
        w3= 0.053347235608839_rp/2.0_rp
        w4= 0.077113760890257_rp/2.0_rp
        posgp(1, 1)= a
        posgp(2, 1)= a         
        posgp(1, 2)= e
        posgp(2, 2)= e
        posgp(1, 3)= b
        posgp(2, 3)= e        
        posgp(1, 4)= e
        posgp(2, 4)= b        
        posgp(1, 5)= f
        posgp(2, 5)= f        
        posgp(1, 6)= c
        posgp(2, 6)= f        
        posgp(1, 7)= f
        posgp(2, 7)= c        
        posgp(1, 8)= d
        posgp(2, 8)= g        
        posgp(1, 9)= d
        posgp(2, 9)= h        
        posgp(1,10)= g
        posgp(2,10)= d        
        posgp(1,11)= g
        posgp(2,11)= h        
        posgp(1,12)= h
        posgp(2,12)= d        
        posgp(1,13)= h
        posgp(2,13)= g
        weigp( 1) = w1
        weigp( 2) = w2
        weigp( 3) = w2
        weigp( 4) = w2
        weigp( 5) = w3
        weigp( 6) = w3
        weigp( 7) = w3
        weigp( 8) = w4
        weigp( 9) = w4
        weigp(10) = w4
        weigp(11) = w4
        weigp(12) = w4
        weigp(13) = w4
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *
! *  Regla cerrada para ngaus==15, p5 exacto  
! *   
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
     else if(ngaus==15) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(1, 4)= 1.0_rp/4.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(1, 5)= 1.0_rp/2.0_rp
        posgp(2, 5)= 0.0_rp
        posgp(1, 6)= 3.0_rp/4.0_rp 
        posgp(2, 6)= 0.0_rp
        posgp(1, 7)= 3.0_rp/4.0_rp
        posgp(2, 7)= 1.0_rp/4.0_rp
        posgp(1, 8)= 1.0_rp/2.0_rp
        posgp(2, 8)= 1.0_rp/2.0_rp
        posgp(1, 9)= 1.0_rp/4.0_rp
        posgp(2, 9)= 3.0_rp/4.0_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 3.0_rp/4.0_rp
        posgp(1,11)= 0.0_rp
        posgp(2,11)= 1.0_rp/2.0_rp
        posgp(1,12)= 0.0_rp
        posgp(2,12)= 1.0_rp/4.0_rp
        posgp(1,13)= (7.0_rp-sqrt(7.0_rp))/21.0_rp
        posgp(2,13)= (7.0_rp-sqrt(7.0_rp))/21.0_rp
        posgp(1,14)= (7.0_rp+2.0_rp*sqrt(7.0_rp))/21.0_rp
        posgp(2,14)= (7.0_rp-sqrt(7.0_rp))/21.0_rp
        posgp(1,15)= (7.0_rp-sqrt(7.0_rp))/21.0_rp
        posgp(2,15)= (7.0_rp+2.0_rp*sqrt(7.0_rp))/21.0_rp
        weigp(   1)= 11.0_rp*sqrt(7.0_rp)/15120.0_rp+1.0_rp/216.0_rp
        weigp(   2)= 11.0_rp*sqrt(7.0_rp)/15120.0_rp+1.0_rp/216.0_rp
        weigp(   3)= 11.0_rp*sqrt(7.0_rp)/15120.0_rp+1.0_rp/216.0_rp
        weigp(   4)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(   5)= 11.0_rp*sqrt(7.0_rp)/630.0_rp-1.0_rp/30.0_rp
        weigp(   6)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(   7)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(   8)= 11.0_rp*sqrt(7.0_rp)/630.0_rp-1.0_rp/30.0_rp
        weigp(   9)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(  10)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(  11)= 11.0_rp*sqrt(7.0_rp)/630.0_rp-1.0_rp/30.0_rp
        weigp(  12)= 4.0_rp/135.0_rp-4.0_rp*sqrt(7.0_rp)/945.0_rp
        weigp(  13)= 49.0_rp/360.0_rp-7.0_rp*sqrt(7.0_rp)/720.0_rp
        weigp(  14)= 49.0_rp/360.0_rp-7.0_rp*sqrt(7.0_rp)/720.0_rp
        weigp(  15)= 49.0_rp/360.0_rp-7.0_rp*sqrt(7.0_rp)/720.0_rp 
     else
        istop=1
     end if
!
! Volume integral ( tetrahedra )
!
  else if(ndime==3) then
     if(ngaus==1) then
        posgp(1,1)= 1.0_rp/4.0_rp
        posgp(2,1)= 1.0_rp/4.0_rp
        posgp(3,1)= 1.0_rp/4.0_rp
        weigp(1)  = 1.0_rp/6.0_rp
     else if(ngaus==4) then
        a=0.5854101966249685_rp
        b=0.1381966011250105_rp
        posgp(1,1)= b
        posgp(2,1)= b
        posgp(3,1)= b
        posgp(1,2)= a
        posgp(2,2)= b
        posgp(3,2)= b
        posgp(1,3)= b
        posgp(2,3)= a
        posgp(3,3)= b
        posgp(1,4)= b
        posgp(2,4)= b
        posgp(3,4)= a
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 1.0_rp/24.0_rp
     else if(ngaus==5) then
        posgp(1,1)= 1.0_rp/4.0_rp
        posgp(2,1)= 1.0_rp/4.0_rp
        posgp(3,1)= 1.0_rp/4.0_rp
        posgp(1,2)= 1.0_rp/6.0_rp
        posgp(2,2)= 1.0_rp/6.0_rp
        posgp(3,2)= 1.0_rp/6.0_rp
        posgp(1,3)= 1.0_rp/2.0_rp
        posgp(2,3)= 1.0_rp/6.0_rp
        posgp(3,3)= 1.0_rp/6.0_rp
        posgp(1,4)= 1.0_rp/6.0_rp
        posgp(2,4)= 1.0_rp/2.0_rp
        posgp(3,4)= 1.0_rp/6.0_rp
        posgp(1,5)= 1.0_rp/6.0_rp
        posgp(2,5)= 1.0_rp/6.0_rp
        posgp(3,5)= 1.0_rp/2.0_rp
        weigp(  1)=-2.0_rp/15.0_rp
        weigp(  2)= 1.5_rp/20.0_rp
        weigp(  3)= 1.5_rp/20.0_rp
        weigp(  4)= 1.5_rp/20.0_rp
        weigp(  5)= 1.5_rp/20.0_rp
     else if(ngaus==11) then
        a=0.3994035761667992_rp
        b=0.1005964238332008_rp
        c=343.0_rp/7500.0_rp/6.0_rp
        d=56.0_rp/375.0_rp/6.0_rp
        posgp(1,1) = 1.0_rp/4.0_rp
        posgp(2,1) = 1.0_rp/4.0_rp
        posgp(3,1) = 1.0_rp/4.0_rp
        posgp(1,2) = 11.0_rp/14.0_rp
        posgp(2,2) = 1.0_rp/14.0_rp
        posgp(3,2) = 1.0_rp/14.0_rp
        posgp(1,3) = 1.0_rp/14.0_rp
        posgp(2,3) = 11.0_rp/14.0_rp
        posgp(3,3) = 1.0_rp/14.0_rp
        posgp(1,4) = 1.0_rp/14.0_rp
        posgp(2,4) = 1.0_rp/14.0_rp
        posgp(3,4) = 11.0_rp/14.0_rp
        posgp(1,5) = 1.0_rp/14.0_rp
        posgp(2,5) = 1.0_rp/14.0_rp
        posgp(3,5) = 1.0_rp/14.0_rp
        posgp(1,6) = a
        posgp(2,6) = a
        posgp(3,6) = b
        posgp(1,7) = a
        posgp(2,7) = b
        posgp(3,7) = a
        posgp(1,8) = a
        posgp(2,8) = b
        posgp(3,8) = b
        posgp(1,9) = b
        posgp(2,9) = a
        posgp(3,9) = a
        posgp(1,10)= b
        posgp(2,10)= a
        posgp(3,10)= b
        posgp(1,11)= b
        posgp(2,11)= b
        posgp(3,11)= a
        weigp(1)   =-148.0_rp/1875.0_rp/6.0_rp
        weigp(2)   = c
        weigp(3)   = c
        weigp(4)   = c
        weigp(5)   = c
        weigp(6)   = d
        weigp(7)   = d
        weigp(8)   = d
        weigp(9)   = d
        weigp(10)  = d
        weigp(11)  = d
     else if(ngaus==14) then
        a=0.0673422422100983_rp
        b=0.3108859192633005_rp
        c=0.7217942490673264_rp
        d=0.0927352503108912_rp
        e=0.4544962958743506_rp
        f=0.0455037041256494_rp
        p=0.1126879257180162_rp/6.0_rp
        q=0.0734930431163619_rp/6.0_rp
        r=0.0425460207770812_rp/6.0_rp
        posgp(1,1) = a
        posgp(2,1) = b
        posgp(3,1) = b
        posgp(1,2) = b
        posgp(2,2) = a
        posgp(3,2) = b
        posgp(1,3) = b
        posgp(2,3) = b
        posgp(3,3) = a
        posgp(1,4) = b
        posgp(2,4) = b
        posgp(3,4) = b
        posgp(1,5) = c
        posgp(2,5) = d
        posgp(3,5) = d
        posgp(1,6) = d
        posgp(2,6) = c
        posgp(3,6) = d
        posgp(1,7) = d
        posgp(2,7) = d
        posgp(3,7) = c
        posgp(1,8) = d
        posgp(2,8) = d
        posgp(3,8) = d
        posgp(1,9) = e
        posgp(2,9) = e
        posgp(3,9) = f
        posgp(1,10)= e
        posgp(2,10)= f
        posgp(3,10)= e
        posgp(1,11)= e
        posgp(2,11)= f
        posgp(3,11)= f
        posgp(1,12)= f
        posgp(2,12)= e
        posgp(3,12)= e
        posgp(1,13)= f
        posgp(2,13)= e
        posgp(3,13)= f
        posgp(1,14)= f
        posgp(2,14)= f
        posgp(3,14)= e
        weigp(1)   = p
        weigp(2)   = p
        weigp(3)   = p
        weigp(4)   = p
        weigp(5)   = q
        weigp(6)   = q
        weigp(7)   = q
        weigp(8)   = q
        weigp(9)   = r
        weigp(10)  = r
        weigp(11)  = r
        weigp(12)  = r
        weigp(13)  = r
        weigp(14)  = r
     else
        istop=1
     end if
  end if
!      
! Errors
!
  if(istop==1) call runend('RUTOPE: NOT AVAILABLE QUADRATURE')

end subroutine rutope

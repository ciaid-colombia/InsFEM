
subroutine rutclo(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
! 
!     This routine sets up the integration constants of closed rules for
!     triangles and tetrahedra 
! 
!             NDIME = 2             NDIME = 3
! 
!          NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!          -----  ----------     -----  ----------
!            3       p1            4       p1
!            4       p2            5       p2
!            6       p1           10       p2
!            7       p3           11       p3
!           10       p3           15       p3
!           15       p5           20       p2 & x^4,...      
! 
!-----------------------------------------------------------------------
  
  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: w1,w2,w3
  integer(ip)              :: istop
  
!
! Line integral
!
  istop=0
  if(ndime==1) then
     if(ngaus==2) then
        posgp(1,1)=-1.0_rp
        posgp(1,2)= 1.0_rp
        weigp(1)= 1.0_rp
        weigp(2)= 1.0_rp
     else if(ngaus==3) then
        posgp(1,1)=-1.0_rp
        posgp(1,2)= 0.0_rp
        posgp(1,3)= 1.0_rp
        weigp(1)= 1.0_rp/3.0_rp
        weigp(2)= 4.0_rp/3.0_rp
        weigp(3)= 1.0_rp/3.0_rp
     else if(ngaus==4) then
        posgp(1,1)=-1.0_rp
        posgp(1,2)=-1.0_rp/3.0_rp
        posgp(1,3)= 1.0_rp/3.0_rp
        posgp(1,4)= 1.0_rp
        weigp(1)= 1.0_rp/4.0_rp
        weigp(2)= 3.0_rp/4.0_rp
        weigp(3)= 3.0_rp/4.0_rp
        weigp(4)= 1.0_rp/4.0_rp
!* * * * * * * * * * * * * * * * * * 
!* 
!* Valores calculados para ngaus==5
!*
!* * * * * * * * * * * * * * * * * * 
     else if(ngaus==5) then
     posgp(1,1)=-1.0_rp
     posgp(1,2)=-1.0_rp/2.0_rp
     posgp(1,3)= 0.0_rp
     posgp(1,4)= 1.0_rp/2.0_rp
     posgp(1,5)= 1.0_rp
     weigp(1)=  7.0_rp/45.0_rp
     weigp(2)= 32.0_rp/45.0_rp
     weigp(3)=  4.0_rp/15.0_rp
     weigp(4)= 32.0_rp/45.0_rp
     weigp(5)=  7.0_rp/45.0_rp
     else
        istop=1
     end if
!
! Area integral (triangles)
!
  else if(ndime==2) then
     if(ngaus==3) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        weigp(  1)= 1.0_rp/6.0_rp
        weigp(  2)= 1.0_rp/6.0_rp
        weigp(  3)= 1.0_rp/6.0_rp
     else if(ngaus==4) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 1.0_rp/3.0_rp
        posgp(2,4)= 1.0_rp/3.0_rp
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 3.0_rp/8.0_rp
     else if(ngaus==6) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 0.5_rp
        posgp(2,4)= 0.0_rp
        posgp(1,5)= 0.5_rp
        posgp(2,5)= 0.5_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 0.5_rp
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 1.0_rp/8.0_rp
        weigp(  5)= 1.0_rp/8.0_rp
        weigp(  6)= 1.0_rp/8.0_rp
     else if(ngaus==7) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 0.5_rp
        posgp(2,4)= 0.0_rp
        posgp(1,5)= 0.5_rp
        posgp(2,5)= 0.5_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 0.5_rp
        posgp(1,7)= 1.0_rp/3.0_rp
        posgp(2,7)= 1.0_rp/3.0_rp
        weigp(  1)= 1.0_rp/40.0_rp
        weigp(  2)= 1.0_rp/40.0_rp
        weigp(  3)= 1.0_rp/40.0_rp
        weigp(  4)= 1.0_rp/15.0_rp
        weigp(  5)= 1.0_rp/15.0_rp
        weigp(  6)= 1.0_rp/15.0_rp
        weigp(  7)= 9.0_rp/40.0_rp
!* * * * * * * * * * * * * * * * * * * * * * 
!* 
!* Regla calculada para ngaus=10, p3 exacto
!*
!* * * * * * * * * * * * * * * * * * * * * * 
     else if(ngaus==10) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(1, 4)= 1.0_rp/3.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(1, 5)= 2.0_rp/3.0_rp
        posgp(2, 5)= 0.0_rp
        posgp(1, 6)= 2.0_rp/3.0_rp 
        posgp(2, 6)= 1.0_rp/3.0_rp
        posgp(1, 7)= 1.0_rp/3.0_rp
        posgp(2, 7)= 2.0_rp/3.0_rp
        posgp(1, 8)= 0.0_rp
        posgp(2, 8)= 2.0_rp/3.0_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 1.0_rp/3.0_rp
        posgp(1,10)= 1.0_rp/3.0_rp
        posgp(2,10)= 1.0_rp/3.0_rp
        weigp(   1)= 1.0_rp/60.0_rp
        weigp(   2)= 1.0_rp/60.0_rp
        weigp(   3)= 1.0_rp/60.0_rp
        weigp(   4)= 3.0_rp/80.0_rp
        weigp(   5)= 3.0_rp/80.0_rp
        weigp(   6)= 3.0_rp/80.0_rp
        weigp(   7)= 3.0_rp/80.0_rp
        weigp(   8)= 3.0_rp/80.0_rp
        weigp(   9)= 3.0_rp/80.0_rp 
        weigp(  10)= 9.0_rp/40.0_rp

! * * * * * * * * * * * * * * * * * * * * * * * 
! *
! *  Regla calculada para ngaus==15, p4 exacto  
! *   
! * * * * * * * * * * * * * * * * * * * * * * *
!    else if(ngaus==15) then
!        posgp(1, 1)= 0.0_rp 
!        posgp(2, 1)= 0.0_rp
!        posgp(1, 2)= 1.0_rp
!        posgp(2, 2)= 0.0_rp
!        posgp(1, 3)= 0.0_rp
!        posgp(2, 3)= 1.0_rp
!        posgp(1, 4)= 1.0_rp/4.0_rp
!        posgp(2, 4)= 0.0_rp
!        posgp(1, 5)= 1.0_rp/2.0_rp
!        posgp(2, 5)= 0.0_rp
!        posgp(1, 6)= 3.0_rp/4.0_rp 
!        posgp(2, 6)= 0.0_rp
!        posgp(1, 7)= 3.0_rp/4.0_rp
!        posgp(2, 7)= 1.0_rp/4.0_rp
!        posgp(1, 8)= 1.0_rp/2.0_rp
!        posgp(2, 8)= 1.0_rp/2.0_rp
!        posgp(1, 9)= 1.0_rp/4.0_rp
!        posgp(2, 9)= 3.0_rp/4.0_rp
!        posgp(1,10)= 0.0_rp
!        posgp(2,10)= 3.0_rp/4.0_rp
!        posgp(1,11)= 0.0_rp
!        posgp(2,11)= 1.0_rp/2.0_rp
!        posgp(1,12)= 0.0_rp
!        posgp(2,12)= 1.0_rp/4.0_rp
!        posgp(1,13)= 1.0_rp/4.0_rp
!        posgp(2,13)= 1.0_rp/4.0_rp
!        posgp(1,14)= 1.0_rp/2.0_rp
!        posgp(2,14)= 1.0_rp/4.0_rp
!        posgp(1,15)= 1.0_rp/4.0_rp
!        posgp(2,15)= 1.0_rp/2.0_rp
!        weigp(1) =  0.0_rp
!        weigp(2) =  0.0_rp
!        weigp(3) =  0.0_rp
!        weigp(4) =  2.0_rp/45.0_rp
!        weigp(5) = -1.0_rp/90.0_rp
!        weigp(6) =  2.0_rp/45.0_rp
!        weigp(7) =  2.0_rp/45.0_rp
!        weigp(8) = -1.0_rp/90.0_rp
!        weigp(9) =  2.0_rp/45.0_rp
!        weigp(10) = 2.0_rp/45.0_rp
!        weigp(11) =-1.0_rp/90.0_rp
!        weigp(12) = 2.0_rp/45.0_rp
!        weigp(13) = 4.0_rp/45.0_rp
!        weigp(14) = 4.0_rp/45.0_rp
!        weigp(15) = 4.0_rp/45.0_rp

! * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *
! *  Regla opcional calculada para ngaus==15, p4 exacto
! *  elemento p4 modificado                   
! *   
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     else if(ngaus==15) then
!        posgp(1, 1)= 0.0_rp 
!        posgp(2, 1)= 0.0_rp
!        posgp(1, 2)= 1.0_rp
!        posgp(2, 2)= 0.0_rp
!        posgp(1, 3)= 0.0_rp
!        posgp(2, 3)= 1.0_rp
!        posgp(1, 4)= 1.0_rp/4.0_rp
!        posgp(2, 4)= 0.0_rp
!        posgp(1, 5)= 1.0_rp/2.0_rp
!        posgp(2, 5)= 0.0_rp
!        posgp(1, 6)= 3.0_rp/4.0_rp 
!        posgp(2, 6)= 0.0_rp
!        posgp(1, 7)= 3.0_rp/4.0_rp
!        posgp(2, 7)= 1.0_rp/4.0_rp
!        posgp(1, 8)= 1.0_rp/2.0_rp
!        posgp(2, 8)= 1.0_rp/2.0_rp
!        posgp(1, 9)= 1.0_rp/4.0_rp
!        posgp(2, 9)= 3.0_rp/4.0_rp
!        posgp(1,10)= 0.0_rp
!        posgp(2,10)= 3.0_rp/4.0_rp
!        posgp(1,11)= 0.0_rp
!        posgp(2,11)= 1.0_rp/2.0_rp
!        posgp(1,12)= 0.0_rp
!        posgp(2,12)= 1.0_rp/4.0_rp
!        posgp(1,13)= 1.0_rp/6.0_rp
!        posgp(2,13)= 1.0_rp/6.0_rp
!        posgp(1,14)= 2.0_rp/3.0_rp
!        posgp(2,14)= 1.0_rp/6.0_rp
!        posgp(1,15)= 1.0_rp/6.0_rp
!        posgp(2,15)= 2.0_rp/3.0_rp
!        weigp(   1)= 2.0_rp/135.0_rp
!        weigp(   2)= 2.0_rp/135.0_rp
!        weigp(   3)= 2.0_rp/135.0_rp
!        weigp(   4)=-4.0_rp/135.0_rp
!        weigp(   5)= 11.0_rp/180.0_rp
!        weigp(   6)=-4.0_rp/135.0_rp
!        weigp(   7)=-4.0_rp/135.0_rp
!        weigp(   8)= 11.0_rp/180.0_rp
!        weigp(   9)=-4.0_rp/135.0_rp
!        weigp(  10)=-4.0_rp/135.0_rp
!        weigp(  11)= 11.0_rp/180.0_rp
!        weigp(  12)=-4.0_rp/135.0_rp
!        weigp(  13)= 3.0_rp/20.0_rp 
!        weigp(  14)= 3.0_rp/20.0_rp 
!        weigp(  15)= 3.0_rp/20.0_rp 
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *
! *  Regla opcional calculada para ngaus==15, p5 exacto
! *  elemento p4 modificado  
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
     if(ngaus==4) then
        posgp(1,1)= 0.0_rp 
        posgp(2,1)= 0.0_rp
        posgp(3,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(3,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(3,3)= 0.0_rp
        posgp(1,4)= 0.0_rp
        posgp(2,4)= 0.0_rp
        posgp(3,4)= 1.0_rp
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 1.0_rp/24.0_rp
     else if(ngaus==5) then
        posgp(1,1)= 0.0_rp 
        posgp(2,1)= 0.0_rp
        posgp(3,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(3,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(3,3)= 0.0_rp
        posgp(1,4)= 0.0_rp
        posgp(2,4)= 0.0_rp
        posgp(3,4)= 1.0_rp
        posgp(1,5)= 1.0_rp/4.0_rp
        posgp(2,5)= 1.0_rp/4.0_rp
        posgp(3,5)= 1.0_rp/4.0_rp
        weigp(  1)= 1.0_rp/120.0_rp
        weigp(  2)= 1.0_rp/120.0_rp
        weigp(  3)= 1.0_rp/120.0_rp
        weigp(  4)= 1.0_rp/120.0_rp
        weigp(  5)= 2.0_rp/15.0_rp
     else if(ngaus==10) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.5_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 0.5_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.0_rp
        posgp(3,10)= 0.5_rp
        weigp(   1)=-1.0_rp/120.0_rp
        weigp(   2)=-1.0_rp/120.0_rp 
        weigp(   3)=-1.0_rp/120.0_rp 
        weigp(   4)=-1.0_rp/120.0_rp 
        weigp(   5)= 1.0_rp/30.0_rp 
        weigp(   6)= 1.0_rp/30.0_rp 
        weigp(   7)= 1.0_rp/30.0_rp 
        weigp(   8)= 1.0_rp/30.0_rp 
        weigp(   9)= 1.0_rp/30.0_rp 
        weigp(  10)= 1.0_rp/30.0_rp 
     else if(ngaus==11) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.5_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 0.5_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.0_rp
        posgp(3,10)= 0.5_rp
        posgp(1,11)= 1.0_rp/4.0_rp
        posgp(2,11)= 1.0_rp/4.0_rp
        posgp(3,11)= 1.0_rp/4.0_rp
        weigp(   1)= 1.0_rp/360.0_rp
        weigp(   2)= 1.0_rp/360.0_rp
        weigp(   3)= 1.0_rp/360.0_rp
        weigp(   4)= 1.0_rp/360.0_rp
        weigp(   5)= 1.0_rp/90.0_rp
        weigp(   6)= 1.0_rp/90.0_rp
        weigp(   7)= 1.0_rp/90.0_rp
        weigp(   8)= 1.0_rp/90.0_rp
        weigp(   9)= 1.0_rp/90.0_rp
        weigp(  10)= 1.0_rp/90.0_rp
        weigp(  11)= 4.0_rp/45.0_rp
     else if(ngaus==15) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.0_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.5_rp
        posgp(2, 9)= 0.0_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.5_rp
        posgp(3,10)= 0.5_rp
        posgp(1,11)= 1.0_rp/3.0_rp
        posgp(2,11)= 1.0_rp/3.0_rp
        posgp(3,11)= 0.0_rp
        posgp(1,12)= 1.0_rp/3.0_rp
        posgp(2,12)= 0.0_rp
        posgp(3,12)= 1.0_rp/3.0_rp
        posgp(1,13)= 1.0_rp/3.0_rp
        posgp(2,13)= 1.0_rp/3.0_rp
        posgp(3,13)= 1.0_rp/3.0_rp
        posgp(1,14)= 0.0_rp
        posgp(2,14)= 1.0_rp/3.0_rp
        posgp(3,14)= 1.0_rp/3.0_rp
        posgp(1,15)= 1.0_rp/4.0_rp
        posgp(2,15)= 1.0_rp/4.0_rp
        posgp(3,15)= 1.0_rp/4.0_rp
        weigp(   1)= 17.0_rp/5040.0_rp
        weigp(   2)= 17.0_rp/5040.0_rp
        weigp(   3)= 17.0_rp/5040.0_rp
        weigp(   4)= 17.0_rp/5040.0_rp
        weigp(   5)=  2.0_rp/315.0_rp
        weigp(   6)=  2.0_rp/315.0_rp
        weigp(   7)=  2.0_rp/315.0_rp
        weigp(   8)=  2.0_rp/315.0_rp
        weigp(   9)=  2.0_rp/315.0_rp
        weigp(  10)=  2.0_rp/315.0_rp
        weigp(  11)=  9.0_rp/840.0_rp
        weigp(  12)=  9.0_rp/840.0_rp
        weigp(  13)=  9.0_rp/840.0_rp
        weigp(  14)=  9.0_rp/840.0_rp
        weigp(  15)= 16.0_rp/315.0_rp
     else if(ngaus==20) then                           ! Integrates P2
        posgp(1, 1)= 0.0_rp                                ! and quartic mono-
        posgp(2, 1)= 0.0_rp                                ! mials to avoid zero 
        posgp(3, 1)= 0.0_rp                                ! weights at the 
        posgp(1, 2)= 1.0_rp                                ! edges
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 1.0_rp/3.0_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 2.0_rp/3.0_rp
        posgp(2, 6)= 0.0_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 2.0_rp/3.0_rp
        posgp(2, 7)= 1.0_rp/3.0_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 1.0_rp/3.0_rp
        posgp(2, 8)= 2.0_rp/3.0_rp
        posgp(3, 8)= 0.0_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 2.0_rp/3.0_rp
        posgp(3, 9)= 0.0_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 1.0_rp/3.0_rp
        posgp(3,10)= 0.0_rp
        posgp(1,11)= 0.0_rp
        posgp(2,11)= 0.0_rp
        posgp(3,11)= 1.0_rp/3.0_rp
        posgp(1,12)= 2.0_rp/3.0_rp
        posgp(2,12)= 0.0_rp
        posgp(3,12)= 1.0_rp/3.0_rp
        posgp(1,13)= 0.0_rp
        posgp(2,13)= 2.0_rp/3.0_rp
        posgp(3,13)= 1.0_rp/3.0_rp
        posgp(1,14)= 0.0_rp
        posgp(2,14)= 0.0_rp
        posgp(3,14)= 2.0_rp/3.0_rp
        posgp(1,15)= 1.0_rp/3.0_rp
        posgp(2,15)= 0.0_rp
        posgp(3,15)= 2.0_rp/3.0_rp
        posgp(1,16)= 0.0_rp
        posgp(2,16)= 1.0_rp/3.0_rp
        posgp(3,16)= 2.0_rp/3.0_rp
        posgp(1,17)= 1.0_rp/3.0_rp
        posgp(2,17)= 1.0_rp/3.0_rp
        posgp(3,17)= 0.0_rp
        posgp(1,18)= 1.0_rp/3.0_rp
        posgp(2,18)= 0.0_rp
        posgp(3,18)= 1.0_rp/3.0_rp
        posgp(1,19)= 1.0_rp/3.0_rp
        posgp(2,19)= 1.0_rp/3.0_rp
        posgp(3,19)= 1.0_rp/3.0_rp
        posgp(1,20)= 0.0_rp
        posgp(2,20)= 1.0_rp/3.0_rp
        posgp(3,20)= 1.0_rp/3.0_rp
        w1 = 383.0_rp/2400.0_rp
        w2 = 1.0_rp/240.0_rp - w1
        w3 = 3.0_rp/80.0_rp  - w2
        weigp(   1)= w1
        weigp(   2)= w1
        weigp(   3)= w1
        weigp(   4)= w1
        weigp(   5)= w2
        weigp(   6)= w2
        weigp(   7)= w2
        weigp(   8)= w2
        weigp(   9)= w2
        weigp(  10)= w2
        weigp(  11)= w2
        weigp(  12)= w2
        weigp(  13)= w2
        weigp(  14)= w2
        weigp(  15)= w2
        weigp(  16)= w2
        weigp(  17)= w3
        weigp(  18)= w3
        weigp(  19)= w3
        weigp(  20)= w3
     else
        istop=1
     end if
  end if
!
! Errors
!
  if(istop==1) call runend('RUTCLO: NOT AVAILABLE INTEGRATION RULE')

end subroutine rutclo

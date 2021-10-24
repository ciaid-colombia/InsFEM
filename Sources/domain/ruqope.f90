subroutine ruqope(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
!
!     This routine sets up the integration constants of open
!     integration rules for brick elements:
! 
!          NDIME = 1             NDIME = 2             NDIME = 3
! 
!      NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!      -----  ----------     -----  ----------     -----  ----------
!        1      q1           1 x 1     q1          1x1x1     q1	
!        2      q3           2 x 2     q3          2x2x2     q3   
!        3      q5           3 x 3     q5          3x3x3     q5
!        4      q7           4 x 4     q7          4x4x4     q7
!        5      q9           5 x 5     q9          5x5x5     q9
!        6      q11          6 x 6     q11         6x6x6     q11
!        7      q13          7 x 7     q13         7x7x7     q13
!        8      q15          8 x 8     q15         8x8x8     q15
!       16      q31         16 x 16    q31        16x16x16   q31
! 
!-----------------------------------------------------------------------
  
  use     typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: posgl(16),weigl(16)
  integer(ip)              :: nlocs,igaus,ilocs,jlocs,klocs
  
  if(ndime==1) then
     nlocs=ngaus
  else if(ndime==2) then
     nlocs=nint(sqrt(real(ngaus,rp)))
  else
     nlocs=nint(real(ngaus,rp)**(1.0_rp/3.0_rp))
  end if

  if(nlocs==1) then
     posgl(1)=0.0_rp
     weigl(1)=2.0_rp
  else if(nlocs==2) then
     posgl(1)=-0.577350269189626_rp
     posgl(2)= 0.577350269189626_rp
     weigl(1)= 1.0_rp
     weigl(2)= 1.0_rp
  else if(nlocs==3) then
     posgl(1)=-0.774596669241483_rp
     posgl(2)= 0.0_rp
     posgl(3)= 0.774596669241483_rp
     weigl(1)= 0.555555555555556_rp
     weigl(2)= 0.888888888888889_rp
     weigl(3)= 0.555555555555556_rp
  else if(nlocs==4)  then
     posgl(1)=-0.861136311594053_rp
     posgl(2)=-0.339981043584856_rp
     posgl(3)= 0.339981043584856_rp
     posgl(4)= 0.861136311594053_rp
     weigl(1)= 0.347854845137454_rp
     weigl(2)= 0.652145154862546_rp
     weigl(3)= 0.652145154862546_rp
     weigl(4)= 0.347854845137454_rp
  else if(nlocs==5)  then
            posgl(1) = -0.906179845938664_rp
            posgl(2) = -0.538469310105683_rp
            posgl(3) =  0.0_rp
            posgl(4) =  0.538469310105683_rp
            posgl(5) =  0.906179845938664_rp
            weigl(1) =  0.236926885056189_rp
            weigl(2) =  0.478628670499366_rp
            weigl(3) =  0.568888888888889_rp
            weigl(4) =  0.478628670499366_rp
            weigl(5) =  0.236926885056189_rp
  else if(nlocs==6)  then
            posgl(1) = -0.932469514203152_rp
            posgl(2) = -0.661209386466265_rp
            posgl(3) = -0.238619186083197_rp
            posgl(4) =  0.238619186083197_rp
            posgl(5) =  0.661209386466265_rp
            posgl(6) =  0.932469514203152_rp
            weigl(1) =  0.171324492379170_rp
            weigl(2) =  0.360761573048139_rp
            weigl(3) =  0.467913934572691_rp
            weigl(4) =  0.467913934572691_rp
            weigl(5) =  0.360761573048139_rp
            weigl(6) =  0.171324492379170_rp
  else if(nlocs==7)  then
            posgl(1) = -0.949107912342759_rp
            posgl(2) = -0.741531185599394_rp
            posgl(3) = -0.405845151377397_rp
            posgl(4) =  0.0_rp
            posgl(5) =  0.405845151377397_rp
            posgl(6) =  0.741531185599394_rp
            posgl(7) =  0.949107912342759_rp
            weigl(1) =  0.129484966168870_rp
            weigl(2) =  0.279705391489277_rp
            weigl(3) =  0.381830050505119_rp
            weigl(4) =  0.417959183673469_rp
            weigl(5) =  0.381830050505119_rp
            weigl(6) =  0.279705391489277_rp
            weigl(7) =  0.129484966168870_rp
  else if(nlocs==8)  then
            posgl(1) = -0.960289856497536_rp
            posgl(2) = -0.796666477413627_rp
            posgl(3) = -0.525532409916329_rp
            posgl(4) = -0.183434642495650_rp
            posgl(5) =  0.183434642495650_rp
            posgl(6) =  0.525532409916329_rp
            posgl(7) =  0.796666477413627_rp
            posgl(8) =  0.960289856497536_rp

            weigl(1) =  0.101228536290376_rp
            weigl(2) =  0.222381034453374_rp
            weigl(3) =  0.313706645877887_rp
            weigl(4) =  0.362683783378362_rp
            weigl(5) =  0.362683783378362_rp
            weigl(6) =  0.313706645877887_rp
            weigl(7) =  0.222381034453374_rp
            weigl(8) =  0.101228536290376_rp
  else if(nlocs==16)  then
            posgl( 1) =-0.98940093499165_rp
            posgl( 2) =-0.94457502307323_rp
            posgl( 3) =-0.86563120238783_rp
            posgl( 4) =-0.75540440835500_rp
            posgl( 5) =-0.61787624440264_rp
            posgl( 6) =-0.45801677765723_rp
            posgl( 7) =-0.28160355077926_rp
            posgl( 8) =-0.09501250983764_rp
            posgl( 9) = 0.09501250983764_rp
            posgl(10) = 0.28160355077926_rp
            posgl(11) = 0.45801677765723_rp
            posgl(12) = 0.61787624440264_rp
            posgl(13) = 0.75540440835500_rp
            posgl(14) = 0.86563120238783_rp
            posgl(15) = 0.94457502307323_rp
            posgl(16) = 0.98940093499165_rp

            weigl( 1) =  0.02715245941175_rp
            weigl( 2) =  0.06225352393865_rp
            weigl( 3) =  0.09515851168249_rp
            weigl( 4) =  0.12462897125553_rp
            weigl( 5) =  0.14959598881658_rp
            weigl( 6) =  0.16915651939500_rp
            weigl( 7) =  0.18260341504492_rp
            weigl( 8) =  0.18945061045507_rp
            weigl( 9) =  0.18945061045507_rp
            weigl(10) =  0.18260341504492_rp
            weigl(11) =  0.16915651939500_rp
            weigl(12) =  0.14959598881658_rp
            weigl(13) =  0.12462897125553_rp
            weigl(14) =  0.09515851168249_rp
            weigl(15) =  0.06225352393865_rp
            weigl(16) =  0.02715245941175_rp
  end if

  if(ndime==1) then
     igaus=0
     do ilocs=1,nlocs
        igaus=igaus+1
        weigp(  igaus)=weigl(ilocs)
        posgp(1,igaus)=posgl(ilocs)
     end do
  else if(ndime==2) then
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           igaus=igaus+1
           weigp(  igaus)=weigl(ilocs)*weigl(jlocs)
           posgp(1,igaus)=posgl(ilocs)
           posgp(2,igaus)=posgl(jlocs)
        end do
     end do
  else if(ndime==3) then
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           do klocs=1,nlocs
              igaus=igaus+1
              weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
              posgp(1,igaus)=posgl(ilocs)
              posgp(2,igaus)=posgl(jlocs)
              posgp(3,igaus)=posgl(klocs)
           end do
        end do
     end do
  end if

end subroutine ruqope


subroutine ruqclo(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
!
!    This routine sets up the integration constants of closed 
!    integration rules for brick elements:
!
!         NDIME = 1             NDIME = 2             NDIME = 3
!
!     NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!     -----  ----------     -----  ----------     -----  ----------  
!       2      q1           2 x 2     q1          2x2x2     q1   
!       3      q2           3 x 3     q2          3x3x3     q2
!       4      q3           4 x 4     q3          4x4x4     q3
!       5      q4           5 x 5     q4          5x5x5     q4
!                             8       p2            20      p2 
!
!-----------------------------------------------------------------------
  
  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: posgl(5), weigl(5)
  integer(ip)              :: inoga(27),igaus,nlocs,ilocs,jlocs,klocs
!
! 2D 8-point integration
!
  if(ndime==2.and.ngaus==8) then
     do igaus=5,8
        weigp(igaus-4)=-1.0_rp/3.0_rp
        weigp(igaus  )= 4.0_rp/3.0_rp
     end do
! Añadido 
     posgp(1,1)= 0.0_rp 
     posgp(1,2)= 0.0_rp 
     posgp(1,3)= 0.0_rp 
     posgp(1,4)= 0.0_rp
     posgp(2,1)= 0.0_rp 
     posgp(2,2)= 0.0_rp 
     posgp(2,3)= 0.0_rp 
     posgp(2,4)= 0.0_rp 
! Hasta aqui
     posgp(1,5)= 0.0_rp 
     posgp(1,6)= 1.0_rp
     posgp(1,7)= 0.0_rp
     posgp(1,8)=-1.0_rp
     posgp(2,5)=-1.0_rp
     posgp(2,6)= 0.0_rp
     posgp(2,7)= 1.0_rp
     posgp(2,8)= 0.0_rp
     return
  end if
!
! 3D 20-point integration
!
  if(ndime==3.and.ngaus==20) then
     nlocs=nint(dfloat(ngaus)**(1.0_rp/3.0_rp))
     do igaus=1,8
        weigp(igaus)=-1.0_rp
     end do
     do igaus=9,20
        weigp(igaus)= 4.0_rp/3.0_rp
     end do
     posgl(1)=-1.0_rp
     posgl(2)= 0.0_rp
     posgl(3)= 1.0_rp
     call chaord(inoga,27)
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           do klocs=1,nlocs
              igaus=igaus+1
              if(inoga(igaus).le.20) then
                 posgp(1,inoga(igaus))=posgl(ilocs)
                 posgp(2,inoga(igaus))=posgl(jlocs)
                 posgp(3,inoga(igaus))=posgl(klocs)
              end if
           end do
        end do
     end do
     return
  end if
!
! Rules obtained from one-dimensional integration
!
  if(ndime==1) then
     nlocs=ngaus
  else if(ndime==2) then
     nlocs=nint(sqrt(dfloat(ngaus)))
  else
     nlocs=nint(dfloat(ngaus)**(1.0_rp/3.0_rp))
  end if

  call chaord(inoga,nlocs**ndime)

  if(nlocs==2) then
     posgl(1)=-1.0_rp
     posgl(2)= 1.0_rp
     weigl(1)= 1.0_rp
     weigl(2)= 1.0_rp
  else if(nlocs==3) then
     posgl(1)=-1.0_rp
     posgl(2)= 0.0_rp
     posgl(3)= 1.0_rp
     weigl(1)= 1.0_rp/3.0_rp
     weigl(2)= 4.0_rp/3.0_rp
     weigl(3)= 1.0_rp/3.0_rp
  else if(nlocs==4) then
     posgl(1)=-1.0_rp
     posgl(2)=-1.0_rp/3.0_rp
     posgl(3)= 1.0_rp/3.0_rp
     posgl(4)= 1.0_rp
     weigl(1)= 1.0_rp/4.0_rp
     weigl(2)= 3.0_rp/4.0_rp
     weigl(3)= 3.0_rp/4.0_rp
     weigl(4)= 1.0_rp/4.0_rp
!* * * * * * * * * * * * * * * * * * 
!* 
!* Valores calculados para ngaus==5
!*
!* * * * * * * * * * * * * * * * * * 
  else if(nlocs==5) then
     posgl(1)=-1.0_rp
     posgl(2)=-1.0_rp/2.0_rp
     posgl(3)= 0.0_rp
     posgl(4)= 1.0_rp/2.0_rp
     posgl(5)= 1.0_rp
     weigl(1)=  7.0_rp/45.0_rp
     weigl(2)= 32.0_rp/45.0_rp
     weigl(3)=  4.0_rp/15.0_rp
     weigl(4)= 32.0_rp/45.0_rp
     weigl(5)=  7.0_rp/45.0_rp
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
           weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
           posgp(1,inoga(igaus))=posgl(ilocs)
           posgp(2,inoga(igaus))=posgl(jlocs)
        end do
     end do
     if(ngaus==5) then                                    !  4       3
        do igaus=1,4                                      !
           weigp(igaus)=1.0_rp/3.0_rp                     !      5
        end do                                            !
        weigp(  5)=8.0_rp/3.0_rp                          !  1       2
        posgp(1,5)=0.0_rp
        posgp(2,5)=0.0_rp
     end if
  else if(ndime==3) then
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           do klocs=1,nlocs
              igaus=igaus+1
              weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
              posgp(1,inoga(igaus))=posgl(ilocs)
              posgp(2,inoga(igaus))=posgl(jlocs)
              posgp(3,inoga(igaus))=posgl(klocs)
           end do
        end do
     end do
     if(ngaus==9) then                                 !    8------7
        do igaus=1,8                                   !   /|     /|
           weigp(igaus)=weigp(igaus)/3.0_rp            !  5------6 |
        end do                                         !  | | 9  | |
        weigp(  9)=16.0_rp/3.0_rp                      !  | 4----|-3
        posgp(1,9)=0.0_rp                              !  |/     |/
        posgp(2,9)=0.0_rp                              !  1------2
        posgp(3,9)=0.0_rp
     end if
  end if

end subroutine ruqclo

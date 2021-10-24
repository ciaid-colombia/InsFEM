subroutine rupclo(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
! 
!     This routine sets up the integration constants of closed rules for
!     PRISMS
! 
!             NDIME = 3    
! 
!          NGAUS  EXACT POL.
!          -----  ----------  
!            3       p1       
! 
!-----------------------------------------------------------------------
  
  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  integer(ip)              :: istop
  
!
! Area integral
!
  if(ndime==3) then
     if(ngaus==6) then
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
        posgp(1,5)= 1.0_rp
        posgp(2,5)= 0.0_rp
        posgp(3,5)= 1.0_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 1.0_rp
        posgp(3,6)= 1.0_rp
        weigp(  1)= 1.0_rp/12.0_rp
        weigp(  2)= 1.0_rp/12.0_rp
        weigp(  3)= 1.0_rp/12.0_rp
        weigp(  4)= 1.0_rp/12.0_rp
        weigp(  5)= 1.0_rp/12.0_rp
        weigp(  6)= 1.0_rp/12.0_rp
     else
        istop=1
     end if
  else if(ndime==2 .and. ngaus==0) then
  else
     istop=1
  end if
!
! Errors
!
  if(istop==1) call runend('RUTCLO: NOT AVAILABLE INTEGRATION RULE')

end subroutine rupclo

subroutine rupope(ndime,ngaus,posgp,weigp)

!-----------------------------------------------------------------------
! 
!     This routine sets up the integration constants of open rules for
!     PRYSMS
! 
!             NDIME = 3    
! 
!          NGAUS  EXACT POL.
!          -----  ----------  
!            1       p1       
!            6       p2       
!       
!-----------------------------------------------------------------------
  
  use      typre
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  integer(ip)              :: istop
  
  istop=0
!
! Volume integral (triangles)
!
  if(ndime==3) then
     if(ngaus==1) then
        posgp(1,1)= 1.0_rp/3.0_rp
        posgp(2,1)= 1.0_rp/3.0_rp
        posgp(3,1)= 1.0_rp/2.0_rp
        weigp(  1)= 1.0_rp/2.0_rp
     else if(ngaus==6) then
        posgp(1,1)= 2.0_rp/3.0_rp
        posgp(2,1)= 1.0_rp/6.0_rp
        posgp(3,1)= 0.21132486540518711774542560974902_rp
        posgp(1,2)= 1.0_rp/6.0_rp
        posgp(2,2)= 2.0_rp/3.0_rp
        posgp(3,2)= 0.21132486540518711774542560974902_rp
        posgp(1,3)= 1.0_rp/6.0_rp
        posgp(2,3)= 1.0_rp/6.0_rp
        posgp(3,3)= 0.21132486540518711774542560974902_rp
        posgp(1,4)= 2.0_rp/3.0_rp
        posgp(2,4)= 1.0_rp/6.0_rp
        posgp(3,4)= 0.78867513459481288225457439025098_rp
        posgp(1,5)= 1.0_rp/6.0_rp
        posgp(2,5)= 2.0_rp/3.0_rp
        posgp(3,5)= 0.78867513459481288225457439025098_rp
        posgp(1,6)= 1.0_rp/6.0_rp
        posgp(2,6)= 1.0_rp/6.0_rp
        posgp(3,6)= 0.78867513459481288225457439025098_rp
        weigp(  1)= 1.0_rp/12.0_rp
        weigp(  2)= 1.0_rp/12.0_rp
        weigp(  3)= 1.0_rp/12.0_rp
        weigp(  4)= 1.0_rp/12.0_rp
        weigp(  5)= 1.0_rp/12.0_rp
        weigp(  6)= 1.0_rp/12.0_rp
     else
        istop=1
     end if
  else
!     istop=1   ! WARNING: Lo he puesto pq viene aqui por boundary shape functions
  end if
!      
! Errors
!
  if(istop==1) call runend('RUTOPE: NOT AVAILABLE QUADRATURE')

end subroutine rupope

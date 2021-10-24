subroutine elemdetjm(nnode,ndime,elcod,detjm)
!-----------------------------------------------------------------------
  !****f* Domain/elmdel
! NAME
!    elmdel
! DESCRIPTION
!    This routine calculates the Cartesian derivatives cartd and detjm.
! USES
!    invmtx
! USED BY
!    nsm_elmope
!    tem_elmope
!    cdr_elmope
!    tem_exaerr
!    extnor
! SOURCE
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: nnode,ndime
  real(rp),    intent(in)  :: elcod(ndime,nnode)
  real(rp),    intent(out) :: detjm
  real(rp)                 :: xjacm(ndime,ndime)
  real(rp)                 :: xjaci(ndime,ndime)
  
  if(ndime==2.and.nnode==3) then
     
     xjacm(1,1) = elcod(1,2)-elcod(1,1)
     xjacm(2,1) = elcod(2,2)-elcod(2,1)
     xjacm(1,2) = elcod(1,3)-elcod(1,1)
     xjacm(2,2) = elcod(2,3)-elcod(2,1)
     
     detjm=xjacm(1,1)*xjacm(2,2)-xjacm(2,1)*xjacm(1,2)


  elseif (ndime==3.and.nnode==4) then
     
     xjacm(1,1) = elcod(1,2)-elcod(1,1)
     xjacm(2,1) = elcod(2,2)-elcod(2,1)
     xjacm(3,1) = elcod(3,2)-elcod(3,1)
     xjacm(1,2) = elcod(1,3)-elcod(1,1)
     xjacm(2,2) = elcod(2,3)-elcod(2,1)
     xjacm(3,2) = elcod(3,3)-elcod(3,1)
     xjacm(1,3) = elcod(1,4)-elcod(1,1)
     xjacm(2,3) = elcod(2,4)-elcod(2,1)
     xjacm(3,3) = elcod(3,4)-elcod(3,1)
     
     xjaci(1,1)  = xjacm(2,2)*xjacm(3,3) - xjacm(3,2)*xjacm(2,3)
     xjaci(2,1)  =-xjacm(2,1)*xjacm(3,3) + xjacm(3,1)*xjacm(2,3)
     xjaci(3,1)  = xjacm(2,1)*xjacm(3,2) - xjacm(3,1)*xjacm(2,2)
     xjaci(2,2) =  xjacm(1,1)*xjacm(3,3) - xjacm(3,1)*xjacm(1,3)
     xjaci(3,2) = -xjacm(1,1)*xjacm(3,2) + xjacm(1,2)*xjacm(3,1)
     xjaci(3,3) =  xjacm(1,1)*xjacm(2,2) - xjacm(2,1)*xjacm(1,2)
     xjaci(1,2) = -xjacm(1,2)*xjacm(3,3) + xjacm(3,2)*xjacm(1,3)
     xjaci(1,3) =  xjacm(1,2)*xjacm(2,3) - xjacm(2,2)*xjacm(1,3)
     xjaci(2,3) = -xjacm(1,1)*xjacm(2,3) + xjacm(2,1)*xjacm(1,3)
     
     detjm = xjacm(1,1)*xjaci(1,1) + xjacm(1,2)*xjaci(2,1)+ xjacm(1,3)*xjaci(3,1)
     
  end if
  
end subroutine elemdetjm
!-----------------------------------------------------------------------
! NOTES
!
! P1 Element in 2D
! ----------------
!
! detjm=(-x1+x2)*(-y1+y3)-(-y1+y2)*(-x1+x3)
!                _                       _    
!               | -y3+y2  -y1+y3   y1-y2  |   
! cartd=1/detjm |                         |   
!               |_ x3-x2   x1-x3  -x1+x2 _|
!
!***
!-----------------------------------------------------------------------

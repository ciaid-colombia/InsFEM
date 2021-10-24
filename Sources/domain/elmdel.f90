subroutine elmdel(nnode,ndime,elcod,cartd,detjm,linea,xjacm,xjaci)
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
  integer(ip), intent(out) :: linea
  real(rp),    intent(out) :: detjm,cartd(ndime,nnode)
  real(rp),    intent(out) :: xjacm(ndime,ndime)
  real(rp),    intent(out) :: xjaci(ndime,ndime)
  real(rp)                 :: denom
  
  if(ndime==2.and.nnode==3) then
     linea=1
     
     xjacm(1,1) = elcod(1,2)-elcod(1,1)
     xjacm(2,1) = elcod(2,2)-elcod(2,1)
     xjacm(1,2) = elcod(1,3)-elcod(1,1)
     xjacm(2,2) = elcod(2,3)-elcod(2,1)
     
     detjm=xjacm(1,1)*xjacm(2,2)-xjacm(2,1)*xjacm(1,2)
     denom=1.0_rp/detjm
     xjaci(1,1) = xjacm(2,2)*denom
     xjaci(2,2) = xjacm(1,1)*denom
     xjaci(2,1) =-xjacm(2,1)*denom
     xjaci(1,2) =-xjacm(1,2)*denom  
     
     cartd(1,1) = -xjaci(1,1) - xjaci(2,1) 
     cartd(1,2) =  xjaci(1,1)
     cartd(1,3) =  xjaci(2,1)
     cartd(2,1) = -xjaci(1,2) - xjaci(2,2) 
     cartd(2,2) =  xjaci(1,2)
     cartd(2,3) =  xjaci(2,2)

     
     !detjm=1.0_rp/(&
     !     &  (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3))&
     !     & -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3)))
     !cartd(1,1)=-elcod(2,3)+elcod(2,2)
     !cartd(1,2)=-elcod(2,1)+elcod(2,3)
     !cartd(1,3)= elcod(2,1)-elcod(2,2)
     !cartd(2,1)= elcod(1,3)-elcod(1,2)
     !cartd(2,2)= elcod(1,1)-elcod(1,3)
     !cartd(2,3)=-elcod(1,1)+elcod(1,2)
     !cartd     = cartd*detjm
     !detjm     = 1.0_rp/detjm
  elseif (ndime==3.and.nnode==4) then
     linea=1
     
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
     denom=1.0_rp/detjm
     xjaci = xjaci*denom
     
     cartd(1,1) = -xjaci(1,1) - xjaci(2,1) - xjaci(3,1)
     cartd(1,2) =  xjaci(1,1)
     cartd(1,3) =  xjaci(2,1)
     cartd(1,4) =  xjaci(3,1)
     cartd(2,1) = -xjaci(1,2) - xjaci(2,2) - xjaci(3,2)
     cartd(2,2) =  xjaci(1,2)
     cartd(2,3) =  xjaci(2,2)
     cartd(2,4) =  xjaci(3,2)
     cartd(3,1) = -xjaci(1,3) - xjaci(2,3) - xjaci(3,3)
     cartd(3,2) =  xjaci(1,3)
     cartd(3,3) =  xjaci(2,3)
     cartd(3,4) =  xjaci(3,3)
     
     !detjm=1.0_rp/(&
     !     & -(elcod(1,1)-elcod(1,2))*(elcod(2,2)-elcod(2,3))*(elcod(3,3)-elcod(3,4)) &
     !     & -(elcod(2,1)-elcod(2,2))*(elcod(3,2)-elcod(3,3))*(elcod(1,3)-elcod(1,4)) &
     !     & -(elcod(3,1)-elcod(3,2))*(elcod(1,2)-elcod(1,3))*(elcod(2,3)-elcod(2,4)) &
     !     & +(elcod(3,1)-elcod(3,2))*(elcod(2,2)-elcod(2,3))*(elcod(1,3)-elcod(1,4)) &
     !     & +(elcod(1,1)-elcod(1,2))*(elcod(3,2)-elcod(3,3))*(elcod(2,3)-elcod(2,4)) &
     !     & +(elcod(2,1)-elcod(2,2))*(elcod(1,2)-elcod(1,3))*(elcod(3,3)-elcod(3,4)) )
     !cartd(1,1)=-elcod(2,3)*elcod(3,4)-elcod(2,2)*elcod(3,3)-elcod(2,4)*elcod(3,2)+elcod(3,3)*elcod(2,4)+elcod(3,2)*elcod(2,3)+elcod(3,4)*elcod(2,2)
     !cartd(1,2)=+elcod(2,4)*elcod(3,1)+elcod(2,3)*elcod(3,4)+elcod(2,1)*elcod(3,3)-elcod(3,4)*elcod(2,1)-elcod(3,3)*elcod(2,4)-elcod(3,1)*elcod(2,3)
     !cartd(1,3)=-elcod(2,1)*elcod(3,2)-elcod(2,4)*elcod(3,1)-elcod(2,2)*elcod(3,4)+elcod(3,1)*elcod(2,2)+elcod(3,4)*elcod(2,1)+elcod(3,2)*elcod(2,4)
     !cartd(1,4)=+elcod(2,2)*elcod(3,3)+elcod(2,1)*elcod(3,2)+elcod(2,3)*elcod(3,1)-elcod(3,2)*elcod(2,3)-elcod(3,1)*elcod(2,2)-elcod(3,3)*elcod(2,1)
     !cartd(2,1)=-elcod(3,3)*elcod(1,4)-elcod(3,2)*elcod(1,3)-elcod(3,4)*elcod(1,2)+elcod(1,3)*elcod(3,4)+elcod(1,2)*elcod(3,3)+elcod(1,4)*elcod(3,2)
     !cartd(2,2)=+elcod(3,4)*elcod(1,1)+elcod(3,3)*elcod(1,4)+elcod(3,1)*elcod(1,3)-elcod(1,4)*elcod(3,1)-elcod(1,3)*elcod(3,4)-elcod(1,1)*elcod(3,3)
     !cartd(2,3)=-elcod(3,1)*elcod(1,2)-elcod(3,4)*elcod(1,1)-elcod(3,2)*elcod(1,4)+elcod(1,1)*elcod(3,2)+elcod(1,4)*elcod(3,1)+elcod(1,2)*elcod(3,4)
     !cartd(2,4)=+elcod(3,2)*elcod(1,3)+elcod(3,1)*elcod(1,2)+elcod(3,3)*elcod(1,1)-elcod(1,2)*elcod(3,3)-elcod(1,1)*elcod(3,2)-elcod(1,3)*elcod(3,1)
     !cartd(3,1)=-elcod(1,3)*elcod(2,4)-elcod(1,2)*elcod(2,3)-elcod(1,4)*elcod(2,2)+elcod(2,3)*elcod(1,4)+elcod(2,2)*elcod(1,3)+elcod(2,4)*elcod(1,2)
     !cartd(3,2)=+elcod(1,4)*elcod(2,1)+elcod(1,3)*elcod(2,4)+elcod(1,1)*elcod(2,3)-elcod(2,4)*elcod(1,1)-elcod(2,3)*elcod(1,4)-elcod(2,1)*elcod(1,3)
     !cartd(3,3)=-elcod(1,1)*elcod(2,2)-elcod(1,4)*elcod(2,1)-elcod(1,2)*elcod(2,4)+elcod(2,1)*elcod(1,2)+elcod(2,4)*elcod(1,1)+elcod(2,2)*elcod(1,4)
     !cartd(3,4)=+elcod(1,2)*elcod(2,3)+elcod(1,1)*elcod(2,2)+elcod(1,3)*elcod(2,1)-elcod(2,2)*elcod(1,3)-elcod(2,1)*elcod(1,2)-elcod(2,3)*elcod(1,1)
     !cartd     = cartd*detjm
     !detjm     = 1.0_rp/detjm
  else
     linea=0
  end if
  
end subroutine elmdel
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

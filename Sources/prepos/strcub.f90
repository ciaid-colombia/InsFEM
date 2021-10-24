subroutine strcub(pnode,npoin,kpoin,lpoin,lnods,strea)

!-----------------------------------------------------------------------
!****f* Mathru/strcub
! NAME 
!    strcub
! DESCRIPTION
!    This routine computes the stream-function at the interior points of
!    the 16-noded cubic quadrilateral and the 10-noded cubic triangle
! USES
! USED BY
!    strloc
!***
!-----------------------------------------------------------------------
!  use      def_kintyp
  use   typre
  implicit none
  integer(ip), intent(in)    :: pnode,npoin
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(inout) :: lpoin(npoin),kpoin
  real(rp),    intent(inout) :: strea(npoin)
  integer(ip)                :: inode,ipoin,ipo10
  real(rp)                   :: contr
  real(rp)                   :: strlo(12)
!
! 16-noded element
!
  if(pnode==16) then
     do inode=1,12
        ipoin=lnods(inode)
        strlo(inode)=strea(ipoin)/real(lpoin(ipoin))
     end do
     strea(lnods(13))=strlo(12)+strlo( 5)-strlo(1)
     strea(lnods(14))=strlo( 6)+strlo( 7)-strlo(2)      
     strea(lnods(15))=strlo( 8)+strlo( 9)-strlo(3)
     strea(lnods(16))=strlo(10)+strlo(11)-strlo(4)      
     lpoin(lnods(13))=lpoin(lnods(13))+1
     lpoin(lnods(14))=lpoin(lnods(14))+1
     lpoin(lnods(15))=lpoin(lnods(15))+1
     lpoin(lnods(16))=lpoin(lnods(16))+1
     kpoin = kpoin + 4
!
! 10-noded element
!
  else if(pnode==10) then
     ipo10=lnods(10)
     strea(ipo10)=0.0
     do inode=1,3
        ipoin=lnods(inode)
        contr=strea(ipoin)/float(lpoin(ipoin))
        strea(ipo10)=strea(ipo10)-contr
     end do
     do inode=4,9
        ipoin=lnods(inode)
        contr=strea(ipoin)/float(lpoin(ipoin))
        strea(ipo10)=strea(ipo10)+contr
     end do
     strea(ipo10)=strea(ipo10)/3.0
     lpoin(ipo10)=lpoin(ipo10)+1
     kpoin = kpoin + 1
  end if
  
end subroutine strcub

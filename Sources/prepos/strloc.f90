subroutine strloc(pnode,ndime,npoin,&
     lnods,lpoin,kpoin,elcod,elvel,strea,&
     nocha,noinv,westr)

!-----------------------------------------------------------------------
!****f* Mathru/strloc
! NAME 
!    strloc
! DESCRIPTION
!    This routine computes the stream-function for each element
! USES
! USED BY
!    stream 
!***
!-----------------------------------------------------------------------
!  use      def_kintyp
  use   typre  
  implicit none
  integer(ip), intent(in)    :: pnode,ndime,npoin
  real(rp),    intent(in)    :: elcod(ndime,pnode),elvel(ndime*pnode)
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(in)    :: nocha(pnode),noinv(pnode)
  real(rp),    intent(in)    :: westr(pnode-1)
  integer(ip), intent(inout) :: lpoin(npoin)
  real(rp),    intent(inout) :: strea(npoin)
  integer(ip)                :: kpoin,inoco,itouc,ipoco,mnodb,mnoda
  integer(ip)                :: icoun,ivxin,ipoi1,ipoi2
  integer(ip)                :: ivxfi,jpoin,ipini,ipfin,inode,lapoi
  integer(ip)                :: ipoi3,ipoi4,ipo13,ipo14,ipo15,ipo16
  real(rp)                   :: xcomp,ycomp,unori,unorf,strfa,preno
  real(rp)                   :: strep,w1,w2,w3,w4
  logical(lg)                :: lofin
!
! Identify first node where the streamfunction is known
!
  inoco=0
  itouc=0
  do while(itouc==0)
     inoco=inoco+1
     ipoco=lnods(inoco)
     if(lpoin(ipoco)>=1) itouc=1
  end do
!
! Compute the streamfuntion at the nodes on the edges of the element
!
  inoco=nocha(inoco)
  do icoun=1,pnode
     mnodb=inoco+icoun
     lofin=.false.
     do while((.not.lofin).and.(mnodb<=pnode))
        if (nocha(mnodb)==0) then
           mnodb=mnodb+1
        else
           lofin=.true.
        end if
     end do
     if(mnodb>pnode) mnodb=mod(mnodb,pnode)
     if(mnodb==0) mnodb=pnode
     mnoda=mnodb-1
     if(mnoda==0) then
        mnoda=pnode
        do while(nocha(mnoda)==0)
           mnoda=mnoda-1
        end do
     end if
     mnodb=noinv(mnodb)
     mnoda=noinv(mnoda)
     xcomp=elcod(1,mnodb)-elcod(1,mnoda)
     ycomp=elcod(2,mnodb)-elcod(2,mnoda)
     ivxin=(mnoda-1)*2+1
     ivxfi=(mnodb-1)*2+1
     unori=-elvel(ivxin)*ycomp+elvel(ivxin+1)*xcomp
     unorf=-elvel(ivxfi)*ycomp+elvel(ivxfi+1)*xcomp
     ipini=lnods(mnoda)
     ipfin=lnods(mnodb)
     strfa=strea(ipini)/real(lpoin(ipini))
     strea(ipfin)=strea(ipfin)+strfa-0.5*(unori+unorf)
     if(lpoin(ipfin).eq.0) kpoin=kpoin+1
     lpoin(ipfin)=lpoin(ipfin)+1
  end do
! 
! Compute the streamfunction at the interior nodes. For the cubic elements
! P3 and Q3 (10 and 16 nodes), a special routine is called. If the same
! algorithm as for the rest of elements is to be used, only the corner
! nodes are employed to compute the streamfunction (see routine CHANUM)
!
  if((pnode==10).or.(pnode==16)) then
     call strcub(pnode,npoin,kpoin,lpoin,lnods,strea)
  else 
     if(nocha(pnode)==0) then
        preno=0.0
        do inode=1,min(pnode-1,8)
           jpoin=lnods(inode)
           strep=strea(jpoin)/real(lpoin(jpoin))
           preno=preno+westr(inode)*strep
        end do
        lapoi=lnods(pnode)
        strea(lapoi)=preno  
        lpoin(lapoi)=lpoin(lapoi)+1
        kpoin=kpoin+1
     end if
  end if
  
end subroutine strloc

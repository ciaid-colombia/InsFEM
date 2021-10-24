module Mod_Setext
  use typre
  use Mod_Mesh
  implicit none
  
contains
subroutine setext(a,zerod,exwor)

!-----------------------------------------------------------------------
!
!    This routine computes the local basis at a boundary node and
!    constructs the vector LPOTY
!
!-----------------------------------------------------------------------
  implicit none
  class(FemMesh) :: a
  real(rp)    :: zerod
  real(rp)    :: exwor(a%ndime,a%npoin)
  
  
  integer(ip) :: ibopo,kdime
  integer(ip) :: kboty,kelty,nlboe,idime,ipoin
  real(rp)    :: xnorm,xnor1,xnor2,xnor3,xnoma
  
!
! Loop over the number of nodal points.
!
  a%exnor=0.0_rp
  ibopo=0
  do ipoin=1,a%npoin
     if(a%lpoty(ipoin)==1) then
        xnorm=0.0_rp
        do idime=1,a%ndime
           xnorm=xnorm+exwor(idime,ipoin)*exwor(idime,ipoin)
        end do
        a%lpoty(ipoin)=ibopo+1
        ibopo=ibopo+1
        ! Normalize exterior vectors.
        do idime=1,a%ndime
           a%exnor(idime,1,ibopo)=exwor(idime,ipoin)/sqrt(xnorm)
        end do
     end if
  end do
!
! Look for additional boundary nodes not detected by exwor and
! explicitely declared by lnodb: OJO check
!      
!  if(ibopo/=nbopo) then 
!     do iboun=1,nboun
!        nlboe=size(lnodb(iboun)%ive)+1
!        kboty=ltypb(iboun)
!        kelty=ltype(lboel(iboun)%ive(nlboe))
!        do inodb=1,nnodb(kboty)
!           ipoin=lnodb(iboun)%ive(inodb)
!           if(a%lpoty(ipoin)==0) then
!              ibopo=ibopo+1
!              a%lpoty(ipoin)=ibopo
!              call extcog(&
!                   nnodb(kboty),a%ndime,a%npoin,nnode(kelty),ltopo,ntens,&
!                   lnodb(iboun)%ive,coord,exwor(1,ipoin))
!           end if
!        end do
!     end do
!  end if
 
  nodes: do ipoin=1,a%npoin

     ibopo=a%lpoty(ipoin)
     if(ibopo/=0) then
!
! Construct the tangent vector in 2D.
!
        if(a%ndime==2) then
           a%exnor(1,2,ibopo)=-a%exnor(2,1,ibopo)
           a%exnor(2,2,ibopo)= a%exnor(1,1,ibopo)
!
! Look for e_k such that n x e_k is maximum.
!
        else if(a%ndime==3) then
           xnor1=a%exnor(1,1,ibopo)*a%exnor(1,1,ibopo)
           xnor2=a%exnor(2,1,ibopo)*a%exnor(2,1,ibopo)
           xnor3=a%exnor(3,1,ibopo)*a%exnor(3,1,ibopo)
           exwor(1,ipoin)=xnor2+xnor3
           exwor(2,ipoin)=xnor1+xnor3
           exwor(3,ipoin)=xnor1+xnor2
           xnoma=0.0_rp
           do idime=1,3
              if(exwor(idime,ipoin).gt.xnoma+1.0e-8_rp) then
                 xnoma=exwor(idime,ipoin)
                 kdime=idime
              end if
           end do
           xnoma=sqrt(xnoma)
!
! Set t_1 = e_k x n, first tangent vector.
!
           if(kdime==1) then
              a%exnor(1,2,ibopo)= 0.0_rp
              a%exnor(2,2,ibopo)=-a%exnor(3,1,ibopo)/xnoma
              a%exnor(3,2,ibopo)= a%exnor(2,1,ibopo)/xnoma
           else if(kdime==2) then
              a%exnor(1,2,ibopo)= a%exnor(3,1,ibopo)/xnoma
              a%exnor(2,2,ibopo)= 0.0_rp
              a%exnor(3,2,ibopo)=-a%exnor(1,1,ibopo)/xnoma
           else if(kdime==3) then
              a%exnor(1,2,ibopo)=-a%exnor(2,1,ibopo)/xnoma
              a%exnor(2,2,ibopo)= a%exnor(1,1,ibopo)/xnoma
              a%exnor(3,2,ibopo)= 0.0_rp
           end if
!            
! Set t_2 = n x t_1, second tangent vector.            
!
           a%exnor(1,3,ibopo)=a%exnor(2,1,ibopo)*a%exnor(3,2,ibopo)&
                -a%exnor(3,1,ibopo)*a%exnor(2,2,ibopo)
           a%exnor(2,3,ibopo)=a%exnor(3,1,ibopo)*a%exnor(1,2,ibopo)&
                -a%exnor(1,1,ibopo)*a%exnor(3,2,ibopo)
           a%exnor(3,3,ibopo)=a%exnor(1,1,ibopo)*a%exnor(2,2,ibopo)&
                -a%exnor(2,1,ibopo)*a%exnor(1,2,ibopo)
        end if
     end if
  end do nodes
  
end subroutine setext
end module
 

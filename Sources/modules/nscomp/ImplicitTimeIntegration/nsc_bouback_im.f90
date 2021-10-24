module Mod_nsc_bouback_im
contains
subroutine nsc_bouback_im(e,gpmom,tresh,wmatr)
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_bouback
! NAME 
!    nsc_bouback
! DESCRIPTION
!    This routine computes a penalization for the NS momentum equation at
!    a given integration point of a boundary received by argument
!    due to the use of a reverse flow filter. The algorithm is:
!    - Compute the normal momentum m at the integration point x.
!    - Evaluate m.n 
!      if(m.n>tresh) then
!        epsilon=0
!      else
!        epsilon= 0.5 * tresh
!      end if
!    - Include epsilon(m.n) contribution on the boundary
! USES
!    vecnor
!    momentum
!    treshold
! USED BY
!    nsc_bouope
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp) :: gpmom(e%ndime),tresh
   real(rp),    intent(inout) :: wmatr(e%ndime+2,e%mnode,e%ndime+2,e%mnode)

   integer(ip)                :: idime
   integer(ip)                :: inodb,jnode,inode,jnodb

   real(rp)                   :: epsi   
   real(rp)                   :: nomom


   !normal momentum 
   nomom = dot_product(e%baloc(:,e%ndime),gpmom)

   if(nomom<=tresh) then
      return
   else
      epsi = 0.5_rp*tresh
   end if

   ! Compute prescribed traction
   do inodb = 1,e%pnodb
      inode = e%lboel(inodb)
      do jnodb = 1,e%pnodb
         jnode = e%lboel(jnodb)
         do idime = 1,e%ndime
            wmatr(1+idime,inode,1+idime,jnode) = wmatr(1+idime,inode,1+idime,jnode) + &
                  epsi*e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)*e%baloc(idime,e%ndime)
         enddo
      enddo
   enddo

end subroutine nsc_bouback_im
end module

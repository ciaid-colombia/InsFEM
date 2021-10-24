subroutine nsc_bouback_ex(e,gpmom,tresh,tract)
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
   real(rp) :: tract(e%ndime+2)

   integer(ip)                :: ldime,idime,ldimb

   real(rp)                   :: epsi   
   real(rp)                   :: nomom


   !normal momentum 
   nomom = dot_product(e%baloc(:,e%ndime),gpmom)

   if(nomom<=tresh) then
      return
   else
      epsi = -0.5_rp*tresh
   end if

   ! Compute prescribed traction
   tract(2:e%ndime+1) = epsi*gpmom(1:e%ndime)*e%baloc(1:e%ndime,e%ndime)

end subroutine nsc_bouback_ex

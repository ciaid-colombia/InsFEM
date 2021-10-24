module Mod_nsc_bouwal
contains
subroutine nsc_bouwal(e,gpvel,visac,denac,delta,tract)

!-----------------------------------------------------------------------
!****f* Nscomp/nsc_bouwal
! NAME 
!    nsm_bouwal
! DESCRIPTION
!    This routine computes the surface traction for the NS equations at
!    a given integration point of a boundary IBOUN received by argument
!    due to the use of a turbulent wall law. The algorithm is:
!    - Compute the tangential velocity u at the integration point x.
!      In fact, u is not necessarily tangential at x. For example:
!      -> u    -> u      
!      o---x---o---x---o
!                      |\ u
!                      | 
!    - Given the distance to the wall y, compute (U*) (friction velocity)
!    - Compute y+=(y)(U*)/(nu)
!      if(y+>5) then
!        t=-(rho)*(U*)^2*(u/|u|)
!      else
!        u+=y+ => U*^2=u*nu/y so that
!        t=-mu*u/y
!      end if
! USES
!    vecnor
!    frivel
! USED BY
!    nsc_bouope
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp) :: gpvel(e%ndime),visac,denac,delta
   real(rp) :: tract(e%ndime+2)

   real(rp),parameter         :: zensi = epsilon(1.0_rp)

   integer(ip)                :: ldime,idime,ldimb

   real(rp)                   :: vikin,velfr,yplus    ! nu, U*, y+, y
   real(rp)                   :: tveno,tvelo(e%ndime),ovelo(e%ndime)             ! |u|, u

   if(delta<zensi) return

   !normal velocity, projection of a over b: a_b = (a . b) b
   ovelo = dot_product(e%baloc(:,e%ndime),gpvel)*e%baloc(:,e%ndime)
   
   tvelo = gpvel - ovelo
   call vecnor(tvelo,e%ndime,tveno,2)                       ! |u|
   vikin=visac/denac                                      ! nu
   if(tveno<=zensi) then
      velfr=0.0_rp
      tract=0.0_rp
      return
   else
      call frivel(delta,tveno,vikin,velfr)                  ! U*
   end if

   ! Compute prescribed traction
 
   yplus=delta*velfr/vikin
   if(yplus<5.0_rp) then
      tract(2:e%ndime+1) = -denac*vikin*tvelo/delta                  ! t=-mu*u/y
   else
      tract(2:e%ndime+1) = -denac*velfr*velfr*tvelo/tveno            ! t=-rho*U*^2*(u/|u|)
   end if

end subroutine nsc_bouwal
end module

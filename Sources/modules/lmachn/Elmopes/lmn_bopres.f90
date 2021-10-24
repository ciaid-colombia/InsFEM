subroutine lmn_bopres(e,wmatr)
!-----------------------------------------------------------------------
!****f* lmachn/lmn_bopres
! NAME 
!    lmn_bopres
! DESCRIPTION
!    This routine computes the contribution to WMATR for the Navier-
!    Stokes equations due to integration on the boundary of v*p*n.
!    It is the same idea as open boundaries but only the pressure
!    part and it is done when there are walls.
!
! USES
! USED BY
!    lmn_bouope 
!***
!-----------------------------------------------------------------------
  use typre
  use Mod_Element
  implicit none
  
  class(FiniteElement) :: e
  real(rp)            :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)

  integer(ip)       :: ldime,inodb,jnodb,inode,jnode
  real(rp)          :: prod1,xmuit


    do ldime=1,e%ndime
        prod1=e%baloc(ldime,e%ndime)
        do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            xmuit=e%shapb(inodb,e%igaub)*prod1
            do jnodb=1,e%pnodb
                jnode = e%lboel(jnodb)
                wmatr(ldime,inode,e%ndime+1,jnode) = &
                wmatr(ldime,inode,e%ndime+1,jnode)   &
                    +xmuit*e%shapb(jnodb,e%igaub)
            end do
            
        end do
    end do

 
end subroutine lmn_bopres

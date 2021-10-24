subroutine lmn_bouopb(a,e,wmatr,acvis,fvins)
!-----------------------------------------------------------------------
!    This routine computes the contribution to WMATR for the Navier-
!    Stokes equations due to open boundaries. In this case, the term
!    S.n (S = -p I + 2 mu Sym(grad(u)) being the Cauchy stress tensor) 
!    is considered unknown.
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_LowMach
   implicit none
   class(LowMachProblem)      :: a
   class(FiniteElement)        :: e
   real(rp),    intent(inout) :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode)
   real(rp),    intent(in)    :: acvis,fvins

   integer(ip)                :: inodb,ldime,jnode,kdime,inode,jnodb
   real(rp)                   :: xmuit
   real(rp)                   :: prod1
   !
   ! Contribution from the viscous term: 2 mu Sym(grad(u).n
   !
   do inodb=1,e%pnodb
      inode = e%lboel(inodb)
      xmuit=e%shapb(inodb,e%igaub)*acvis
      do ldime=1,e%ndime
         do jnode=1,e%pnode
            wmatr(ldime,inode,ldime,jnode) = wmatr(ldime,inode,ldime,jnode)  &
             -xmuit*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))  
         end do
         if (fvins>zelmn) then
            do jnode=1,e%pnode
               do kdime=1,e%ndime
                  wmatr(ldime,inode,kdime,jnode)=wmatr(ldime,inode,kdime,jnode)&
                        -xmuit*e%cartb(ldime,jnode)*e%baloc(kdime,e%ndime)
               end do
            end do
         end if
      end do
   end do
   
   !Contribution from the pressure term.     
   do ldime=1,e%ndime
      prod1=e%baloc(ldime,e%ndime)
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         xmuit=e%shapb(inodb,e%igaub)*prod1
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            wmatr(ldime,inode,e%ndime+2,jnode)= wmatr(ldime,inode,e%ndime+2,jnode) &
               +xmuit*e%shapb(jnodb,e%igaub)
         end do
      end do
   end do
   
  
end subroutine lmn_bouopb

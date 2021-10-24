subroutine lmn_boucon(a,e,wmatr,acden)
   use typre
   use Mod_Element
   use Mod_LowMach
   implicit none
   class(LowMachProblem)      :: a
   class(FiniteElement)        :: e
   real(rp),    intent(inout) :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode)
   real(rp),    intent(in)    :: acden
   integer(ip)                :: inodb,jnode,inode,jnodb
   
   do inodb=1,e%pnodb
      inode = e%lboel(inodb)
      do jnodb=1,e%pnodb
         jnode = e%lboel(jnodb)
         wmatr(e%ndime+2,inode,1:e%ndime,jnode) = acden* e%shapb(inodb,e%igaub) &
         * e%shapb(jnodb,e%igaub)* e%baloc(1:e%ndime,e%ndime) &
         + wmatr(e%ndime+2,inode,1:e%ndime,jnode) 
      end do
   end do
  
end subroutine lmn_boucon

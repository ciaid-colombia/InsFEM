subroutine sup_bouopb_p(e,wmatr,lawvi,LawviParam,fvins)
!-----------------------------------------------------------------------
!    This routine computes the contribution to WMATR for the Navier-
!    Stokes equations due to open boundaries. In this case, the term
!    S.n (S = -p I + Sigma) being the Cauchy stress tensor) 
!    is considered unknown.
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement)        :: e
   real(rp),    intent(inout) :: wmatr(e%ndime+(e%ndime -1)*(e%ndime-1)+2+1,e%mnode,e%ndime+(e%ndime -1)*(e%ndime-1)+2+1,e%mnode)
   real(rp), parameter        :: zensi = epsilon(1.0_rp)
   integer(ip), intent(in)    :: lawvi
   real(rp),  intent(in)      :: LawviParam(*),fvins 
   integer(ip)                :: inodb,ldime,jnode,kdime,inode,jnodb,auxtens
   real(rp)                   :: xmuit,auxvis
   real(rp)                   :: prod1
   
   auxtens=(e%ndime-1)*(e%ndime-1)+2
   
   if(lawvi<0)then
      auxvis=LawviParam(1)*LawviParam(2) !mu*beta      
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         xmuit=e%shapb(inodb,e%igaub)*auxvis
         do ldime=1,e%ndime
            do jnode=1,e%pnode
               wmatr(auxtens+ldime,inode,auxtens+ldime,jnode) = wmatr(auxtens+ldime,auxtens+inode,ldime,jnode)  &
               -xmuit*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))  
            end do
            if (fvins>zensi) then
               do jnode=1,e%pnode
                  do kdime=1,e%ndime
                     wmatr(auxtens+ldime,inode,auxtens+kdime,jnode)=wmatr(auxtens+ldime,inode,auxtens+kdime,jnode)&
                           -xmuit*e%cartb(ldime,jnode)*e%baloc(kdime,e%ndime)
                  end do
               end do
            end if
         end do
      end do   
   end if  

   ! Contribution from the viscous term: 2 mu Sym(grad(u).n => (in three field  sigma.n)     
   if(e%ndime==2)then
      do inodb=1,e%pnodb
         inode=e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode=e%lboel(jnodb)
            wmatr(auxtens+1,inode,1,jnode)=wmatr(auxtens+1,inode,1,jnode) + e%shapb(inodb,e%igaub)*e%baloc(1,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+1,inode,2,jnode)=wmatr(auxtens+1,inode,2,jnode) + 0.0_rp 
            wmatr(auxtens+1,inode,3,jnode)=wmatr(auxtens+1,inode,3,jnode) + e%shapb(inodb,e%igaub)*e%baloc(2,e%ndime)*e%shapb(jnodb,e%igaub)
            
            wmatr(auxtens+2,inode,1,jnode)=wmatr(auxtens+2,inode,1,jnode) + 0.0_rp 
            wmatr(auxtens+2,inode,2,jnode)=wmatr(auxtens+2,inode,2,jnode) + e%shapb(inodb,e%igaub)*e%baloc(2,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+2,inode,3,jnode)=wmatr(auxtens+2,inode,3,jnode) + e%shapb(inodb,e%igaub)*e%baloc(1,e%ndime)*e%shapb(jnodb,e%igaub)
         end do
      end do
    elseif(e%ndime==3)then
      do inodb=1,e%pnodb
         inode=e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode=e%lboel(jnodb)            
            wmatr(auxtens+1,inode,1,jnode)=wmatr(auxtens+1,inode,1,jnode) + e%shapb(inodb,e%igaub)*e%baloc(1,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+1,inode,2,jnode)=wmatr(auxtens+1,inode,2,jnode) + 0.0_rp 
            wmatr(auxtens+1,inode,3,jnode)=wmatr(auxtens+1,inode,3,jnode) + 0.0_rp
            wmatr(auxtens+1,inode,4,jnode)=wmatr(auxtens+1,inode,4,jnode) + 0.0_rp
            wmatr(auxtens+1,inode,5,jnode)=wmatr(auxtens+1,inode,5,jnode) + e%shapb(inodb,e%igaub)*e%baloc(3,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+1,inode,6,jnode)=wmatr(auxtens+1,inode,6,jnode) + e%shapb(inodb,e%igaub)*e%baloc(2,e%ndime)*e%shapb(jnodb,e%igaub)
          
            
            wmatr(auxtens+2,inode,1,jnode)=wmatr(auxtens+2,inode,1,jnode) + 0.0_rp 
            wmatr(auxtens+2,inode,2,jnode)=wmatr(auxtens+2,inode,2,jnode) + e%shapb(inodb,e%igaub)*e%baloc(2,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+2,inode,3,jnode)=wmatr(auxtens+2,inode,3,jnode) + 0.0_rp
            wmatr(auxtens+2,inode,4,jnode)=wmatr(auxtens+2,inode,4,jnode) + e%shapb(inodb,e%igaub)*e%baloc(3,e%ndime)*e%shapb(jnodb,e%igaub) 
            wmatr(auxtens+2,inode,5,jnode)=wmatr(auxtens+2,inode,5,jnode) + 0.0_rp
            wmatr(auxtens+2,inode,6,jnode)=wmatr(auxtens+2,inode,6,jnode) + e%shapb(inodb,e%igaub)*e%baloc(1,e%ndime)*e%shapb(jnodb,e%igaub)
            
            wmatr(auxtens+3,inode,1,jnode)=wmatr(auxtens+3,inode,1,jnode) + 0.0_rp 
            wmatr(auxtens+3,inode,2,jnode)=wmatr(auxtens+3,inode,2,jnode) + 0.0_rp
            wmatr(auxtens+3,inode,3,jnode)=wmatr(auxtens+3,inode,3,jnode) + e%shapb(inodb,e%igaub)*e%baloc(3,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+3,inode,4,jnode)=wmatr(auxtens+3,inode,4,jnode) + e%shapb(inodb,e%igaub)*e%baloc(2,e%ndime)*e%shapb(jnodb,e%igaub) 
            wmatr(auxtens+3,inode,5,jnode)=wmatr(auxtens+3,inode,5,jnode) + e%shapb(inodb,e%igaub)*e%baloc(1,e%ndime)*e%shapb(jnodb,e%igaub)
            wmatr(auxtens+3,inode,6,jnode)=wmatr(auxtens+3,inode,6,jnode) + 0.0_rp
         end do
      end do    
   end if     
   
  
end subroutine sup_bouopb_p

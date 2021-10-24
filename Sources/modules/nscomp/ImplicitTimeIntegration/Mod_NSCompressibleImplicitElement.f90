module Mod_NSCompressibleImplicitElement
   use typre
   use Mod_Element
   
contains
   
   subroutine nsc_elmbdq(e,dvolu,dtinv,Add,Amd,LTdd,LTdm,elmdq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block rho,q 
      !    (q, rho/dt) + (q, Add·rho) + Tau_d(q·L*dd, Add·rho) + Tau_m(q·L*dm, Amd·rho) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Add(e%mnode)
      real(rp),    intent(in)    :: Amd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm),Amd)
      
      do jnode=1,e%pnode
         elmdq(1:e%pnode,jnode) = elmdq(1:e%pnode,jnode) + &
            (e%shape(1:e%pnode,e%igaus)*(e%shape(jnode,e%igaus)*dtinv+Add(jnode)) + &
             LTdd(1:e%pnode)*Add(jnode) + aux(1:e%pnode,jnode))*dvolu
      end do

   end subroutine nsc_elmbdq

   subroutine nsc_elmbmq(e,dvolu,dtinv,Adm,Amm,LTdd,LTdm,elmmq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block mom,q 
      !  +(q, Adm·mom) + Tau_d(q·L*dd, Adm·mom) + Tau_m(q·L*dm, Amm·mom) 
      !  +Tau_m(q·L*dm, mom/dt) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Adm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Amm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmmq(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmq(inode,:,jnode) = elmmq(inode,:,jnode) + &
             ((e%shape(inode,e%igaus)+LTdd(inode))*Adm(:,jnode) + &
             matmul(LTdm(:,inode),Amm(:,:,jnode)) + &
             LTdm(:,inode)*e%shape(jnode,e%igaus)*dtinv)*dvolu 
         end do
      end do

   end subroutine nsc_elmbmq
  
   subroutine nsc_elmbeq(e,dvolu,Ame,LTdm,elmeq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block ene,q 
      !  + Tau_m(q·L*dm, Ame·ene) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ame(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeq(e%mnode,e%mnode)

      integer(ip)                :: inode,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm),Ame)
      forall (inode=1:e%pnode, jnode=1:e%pnode)
         elmeq(inode,jnode) = elmeq(inode,jnode) + aux(inode,jnode)*dvolu 
      end forall

   end subroutine nsc_elmbeq

   subroutine nsc_elmbdn(e,dvolu,dtinv,Add,Amd,Aed,LTmd,LTmm,LTme,elmdn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block rho,n 
      !  + (n, Amd·rho) + Tau_d(n·L*md, Add·rho) + Tau_d(n·L*md, rho/dt) 
      !  + Tau_m(n·L*mm, Amd·rho) + Tau_e(n·L*me, Aed·rho)  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Add(e%mnode)
      real(rp),    intent(in)    :: Amd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aed(e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmdn(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Amd)
         do jnode=1,e%pnode
            elmdn(idime,1:e%pnode,jnode) = elmdn(idime,1:e%pnode,jnode) + &
            (e%shape(1:e%pnode,e%igaus)*Amd(idime,jnode) + &
             LTmd(idime,1:e%pnode)*(Add(jnode)+e%shape(jnode,e%igaus)*dtinv)  + &
             aux(1:e%pnode,jnode) + LTme(idime,1:e%pnode)*Aed(jnode))*dvolu 
         end do
      end do

   end subroutine nsc_elmbdn

   subroutine nsc_elmbmn(e,dvolu,dtinv,Adm,Amm,Aem,LTmd,LTmm,LTme,elmmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block mom,n 
      !  (n, mom/dt) + (n, Amm·mom) + Tau_m(n·L*mm, mom/dt) 
      !  + Tau_d(n·L*md, Adm·mom) + Tau_m(n·L*mm, Amm·mom) + Tau_e(n·L*me, Aem·mom)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Adm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Amm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmmn(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: idime,inode,jdime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      do idime=1,e%ndime
         do inode=1,e%pnode
            elmmn(idime,1:e%pnode,idime,inode) = elmmn(idime,1:e%pnode,idime,inode) + & 
            e%shape(1:e%pnode,e%igaus)*e%shape(inode,e%igaus)*dtinv*dvolu
            do jdime=1,e%ndime
               do jnode=1,e%pnode
                  elmmn(idime,inode,jdime,jnode) = elmmn(idime,inode,jdime,jnode) + &
                (e%shape(inode,e%igaus)*Amm(idime,jdime,jnode)+ &
                   LTmd(idime,inode)*Adm(jdime,jnode) + &
                   dot_product(LTmm(idime,:,inode),Amm(:,jdime,jnode)) +&
                   LTmm(idime,jdime,inode)*e%shape(jnode,e%igaus)*dtinv + &                
                   LTme(idime,inode)*Aem(jdime,jnode))*dvolu 
               end do
            end do
         end do
      end do

   end subroutine nsc_elmbmn

   subroutine nsc_elmben(e,dvolu,dtinv,Ame,Aee,LTmm,LTme,elmen)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block ene,n 
      !  +(n, Ame·ene) + Tau_m(n·L*mm, Ame·ene) + Tau_e(n·L*me, Aee·ene)
      !  + Tau_e(n·L*me, ene/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ame(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aee(e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmen(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Ame)
         do jnode=1,e%pnode
            elmen(idime,1:e%pnode,jnode) = elmen(idime,1:e%pnode,jnode) + &
            (e%shape(1:e%pnode,e%igaus)*Ame(idime,jnode) + aux(1:e%pnode,jnode) + &
             LTme(idime,1:e%pnode)*(Aee(jnode)+e%shape(jnode,e%igaus)*dtinv))*dvolu 
         end do
      end do

   end subroutine nsc_elmben

   subroutine nsc_elmbdg(e,dvolu,dtinv,Add,Amd,Aed,LTed,LTem,LTee,elmdg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block rho,g 
      !  +(g, Aed·rho) + Tau_d(g·L*ed, rho/dt) + Tau_d(g·L*ed, Add·rho) 
      !  + Tau_m(g·L*em, Amd·rho) + Tau_e(g·L*ee, Aed·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Add(e%mnode)
      real(rp),    intent(in)    :: Amd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aed(e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmdg(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      integer(ip)                :: jnode
      
      aux = matmul(transpose(LTem),Amd)      
      do jnode=1,e%pnode
         elmdg(1:e%pnode,jnode) = elmdg(1:e%pnode,jnode) + &
         ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*Aed(jnode) + &
          LTed(1:e%pnode)*(e%shape(jnode,e%igaus)*dtinv+Add(jnode)) + &
          aux(1:e%pnode,jnode))*dvolu 
      end do

   end subroutine nsc_elmbdg

   subroutine nsc_elmbmg(e,dvolu,dtinv,Adm,Amm,Aem,LTed,LTem,LTee,elmmg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block mom,g 
      !  +(g, Aem·mom) + Tau_d(g·L*ed, Adm·mom) + Tau_m(g·L*em, Amm·mom) + Tau_e(g·L*ee, Aem·mom)
      !  + Tau_m(g·L*em, mom/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Adm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Amm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmmg(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmg(inode,:,jnode) = elmmg(inode,:,jnode) + &
            ((e%shape(inode,e%igaus)+LTee(inode))*Aem(:,jnode) + &
              LTed(inode)*Adm(:,jnode) + &
              matmul(LTem(:,inode),Amm(:,:,jnode)) + &
              LTem(:,inode)*e%shape(jnode,e%igaus)*dtinv)*dvolu 
         end do
      end do

   end subroutine nsc_elmbmg

   subroutine nsc_elmbeg(e,dvolu,dtinv,Ame,Aee,LTem,LTee,elmeg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block ene,g 
      !  (g, ene/dt) + (g, Aee·ene) + Tau_m(g·L*em, Ame·ene) + Tau_e(g·L*ee, Aee·ene)
      !  + Tau_e(g·L*ee, ene/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ame(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Aee(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elmeg(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)
      
      aux = matmul(transpose(LTem),Ame)      
      do jnode=1,e%pnode
         elmeg(1:e%pnode,jnode) = elmeg(1:e%pnode,jnode) + &
         ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*(e%shape(jnode,e%igaus)*dtinv+Aee(jnode)) + &
           aux(1:e%pnode,jnode))*dvolu 
      end do

   end subroutine nsc_elmbeg

   subroutine nsc_elmrhd(e,dvolu,LTdm,elexd,elexm,eltemd,eltemm,elrhd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for density
      !    (q, fd) + Tau_m(q.L*dm, fm)
      !    + (q, rho_n/dt) + Tau_m(q·L*dm, mom_n/dt) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: elexd(1),elexm(e%ndime)
      real(rp),    intent(in)    :: eltemd(1),eltemm(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhd(e%mnode)

      real(rp)                   :: aux(e%mnode)
      real(rp)                   :: auxd,auxm(e%ndime)
     
      auxd=elexd(1)+eltemd(1)
      auxm=elexm+eltemm
      
      aux = matmul(auxm,LTdm)

      elrhd(1:e%pnode) = elrhd(1:e%pnode) + &
       (e%shape(1:e%pnode,e%igaus)*auxd+aux(1:e%pnode))*dvolu 

   end subroutine nsc_elmrhd

   subroutine nsc_elmrhm(e,dvolu,LTmd,LTmm,LTme,elexd,elexm,elexe,eltemd,eltemm,elteme,elrhm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for density
      !    (n, fm) + (n, mom_n/dt)
      !    + Tau_d(n·L*md, rho_n/dt) + Tau_d(n·L*md, fd)
      !    + Tau_m(n·L*mm, mom_n/d)  + Tau_m(n·L*mm, fm)
      !    + Tau_e(n·L*me, ene_n/dt) + Tau_e(n·L*me, fe) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: elexd(1),elexm(e%ndime),elexe(1)
      real(rp),    intent(in)    :: eltemd(1),eltemm(e%ndime),elteme(1)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhm(e%ndime,e%mnode)

      integer(ip)                :: idime
      real(rp)                   :: aux(e%mnode)
      real(rp)                   :: auxd,auxm(e%ndime),auxe
     
      auxd=elexd(1)+eltemd(1)
      auxm=elexm+eltemm
      auxe=elexe(1)+elteme(1)

      do idime=1,e%ndime
         aux(:) =  matmul(transpose(LTmm(idime,:,:)),auxm) 
         elrhm(idime,1:e%pnode) = elrhm(idime,1:e%pnode) + &
          (e%shape(1:e%pnode,e%igaus)*auxm(idime) + &
          LTmd(idime,1:e%pnode)*auxd + LTme(idime,1:e%pnode)*auxe + &
           aux(1:e%pnode))*dvolu 
      end do

   end subroutine nsc_elmrhm

   subroutine nsc_elmrhe(e,dvolu,LTed,LTem,LTee,elexd,elexm,elexe,eltemd,eltemm,elteme,elrhe)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for density
      !    (g, fe) + (g, ene_n/dt)
      !    + Tau_d(g·L*ed, rho_n/dt) + Tau_d(g·L*ed, fd)
      !    + Tau_m(g·L*em, mom_n/d)  + Tau_m(g·L*em, fm)
      !    + Tau_e(g·L*ee, ene_n/dt) + Tau_e(g·L*ee, fe) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: elexd(1),elexm(e%ndime),elexe(1)
      real(rp),    intent(in)    :: eltemd(1),eltemm(e%ndime),elteme(1)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhe(e%mnode)

      integer(ip)                :: inode
      real(rp)                   :: aux(e%mnode)
      real(rp)                   :: auxd,auxm(e%ndime),auxe
     
      auxd=elexd(1)+eltemd(1)
      auxm=elexm+eltemm
      auxe=elexe(1)+elteme(1)

      aux =  matmul(auxm,LTem) 

      elrhe(1:e%pnode) = elrhe(1:e%pnode) + &
       ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*auxe + &
       LTed(1:e%pnode)*auxd + aux(1:e%pnode))*dvolu 

   end subroutine nsc_elmrhe

   subroutine nsc_elmbdq_diff(e,dvolu,Kmd,Ked,LTdm,LTde,elmdq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion lhs terms for ASGS for block rho,q
      !    - Tau_m(q·L*dm, Kmd·rho) - Tau_e(q·L*de, Ked·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Kmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Ked(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm(:,:)),Kmd)

      do jnode=1,e%pnode
            elmdq(1:e%pnode,jnode) = elmdq(1:e%pnode,jnode) - &
            (aux(1:e%pnode,jnode) + LTde(1:e%pnode)*Ked(jnode))*dvolu
      end do

   end subroutine nsc_elmbdq_diff

   subroutine nsc_elmbmq_diff(e,dvolu,Kmm,Kem,LTdm,LTde,elmmq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion lhs terms for ASGS for block mom,q 
      !  - Tau_m(q·L*dm, Kmm·mom) - Tau_e(q·L*de, Kem·mom) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Kmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmq(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmq(inode,:,jnode) = elmmq(inode,:,jnode) - &
             (matmul(LTdm(:,inode),Kmm(:,:,jnode)) + &
             LTde(inode)*Kem(:,jnode))*dvolu 
         end do
      end do

   end subroutine nsc_elmbmq_diff
   
   subroutine nsc_elmbeq_diff(e,dvolu,Kee,LTde,elmeq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion lhs terms for ASGS for block tem,q
      !    - Tau_e(q·L*de, Kee·tem)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Kee(e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeq(e%mnode,e%mnode)

      integer(ip)                :: jnode

      do jnode=1,e%pnode
            elmeq(1:e%pnode,jnode) = elmeq(1:e%pnode,jnode) -LTde(1:e%pnode)*Kee(jnode)*dvolu
      end do

   end subroutine nsc_elmbeq_diff

   subroutine nsc_elmbdn_diff(e,dvolu,KGmd,Kmd,Ked,LTmm,LTme,elmdn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion lhs terms for ASGS for block rho,n 
      !  +(d_k n, Kmd_k·rho) - Tau_m(n·L*mm, Kmd·rho) - Tau_e(n·L*me, Ked·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGmd(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Ked(e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdn(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)
      real(rp)                   :: aux1(e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Kmd)
         do jnode=1,e%pnode
             aux1(:) = matmul(transpose(e%cartd(:,:)),KGmd(:,idime,jnode))
            elmdn(idime,1:e%pnode,jnode) = elmdn(idime,1:e%pnode,jnode) + &
             (aux1(1:e%pnode)-aux(1:e%pnode,jnode)-LTme(idime,1:e%pnode)*Ked(jnode))*dvolu
         end do
      end do

   end subroutine nsc_elmbdn_diff

   subroutine nsc_elmbmn_diff(e,dvolu,KGmm,Kmm,Kem,LTmm,LTme,elmmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block mom,n 
      !  (d_k n, Kmm_k·mom) - Tau_m(n·L*mm, Kmm·mom) - Tau_e(n·L*me, Kem·mom)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGmm(e%ndime,e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmn(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: idime,jdime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)
      real(rp)                   :: aux1(e%mnode)

      do idime=1,e%ndime
         do jdime=1,e%ndime
            aux = matmul(transpose(LTmm(idime,:,:)),Kmm(:,jdime,:))
            do jnode=1,e%pnode
               aux1(:) = matmul(transpose(e%cartd(:,:)),KGmm(:,idime,jdime,jnode))
               elmmn(idime,1:e%pnode,jdime,jnode) = elmmn(idime,1:e%pnode,jdime,jnode) + &
                  (aux1(1:e%pnode)-aux(1:e%pnode,jnode)-LTme(idime,1:e%pnode)*Kem(jdime,jnode))*dvolu 
            end do
         end do
      end do
      
   end subroutine nsc_elmbmn_diff

   subroutine nsc_elmben_diff(e,dvolu,Kee,LTme,elmen)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion lhs terms for ASGS for block ene,n 
      !  - Tau_e(n·L*me, Kee·ene)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Kee(e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmen(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode

      
      do idime=1,e%ndime
         do jnode=1,e%pnode
            elmen(idime,1:e%pnode,jnode) = elmen(idime,1:e%pnode,jnode) - &
             LTme(idime,1:e%pnode)*Kee(jnode)*dvolu 
         end do
      end do

   end subroutine nsc_elmben_diff

   subroutine nsc_elmbdg_diff(e,dvolu,KGed,Kmd,Ked,LTem,LTee,elmdg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block rho,g 
      !  +(d_k g, Ked_k·rho) - Tau_m(g·L*em, Kmd·rho) - Tau_e(g·L*ee,Ked·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGed(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Ked(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdg(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      real(rp)                   :: aux1(e%mnode,e%mnode)
      integer(ip)                :: jnode
      
      aux = matmul(transpose(e%cartd),KGed) 
      aux1 = matmul(transpose(LTem),Kmd)      
      do jnode=1,e%pnode
         elmdg(1:e%pnode,jnode) = elmdg(1:e%pnode,jnode) + &
         (aux(1:e%pnode,jnode) - aux1(1:e%pnode,jnode) - LTee(1:e%pnode)*Ked(jnode))*dvolu 
      end do

   end subroutine nsc_elmbdg_diff

   subroutine nsc_elmbmg_diff(e,dvolu,KGem,Kmm,Kem,LTem,LTee,elmmg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusive lhs terms for ASGS for block mom,g 
      !  +(d_k g, Kem_k·mom) - Tau_m(g·L*em, Kmm·mom) - Tau_e(g·L*ee, Kem·mom)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGem(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmg(e%mnode,e%ndime,e%mnode)

      real(rp)                   :: aux(e%mnode)
      real(rp)                   :: aux1
      integer(ip)                :: inode,jnode,jdime
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            do jdime=1,e%ndime
               aux(:) = matmul(transpose(e%cartd(:,:)),KGem(:,jdime,jnode))
               aux1 = dot_product(LTem(:,inode),Kmm(:,jdime,jnode))
               elmmg(inode,jdime,jnode) = elmmg(inode,jdime,jnode) + &
                (aux(inode) - aux1 - LTee(inode)*Kem(jdime,jnode))*dvolu
            end do
         end do
      end do

   end subroutine nsc_elmbmg_diff

   subroutine nsc_elmbeg_diff(e,dvolu,KGee,Kee,LTee,elmeg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusive lhs terms for ASGS for block ene,g 
      !  + (d_k g, Kee_k·ene) - Tau_e(g·L*ee, Kee·ene)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGee(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Kee(e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeg(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      integer(ip)                :: jnode

      aux = matmul(transpose(e%cartd),KGee) 
      do jnode = 1,e%pnode
         elmeg(1:e%pnode,jnode) = elmeg(1:e%pnode,jnode) + &
           (aux(1:e%pnode,jnode) - LTee(1:e%pnode)*Kee(jnode))*dvolu 
      end do

   end subroutine nsc_elmbeg_diff

   subroutine nsc_elmbdq_reac(e,dvolu,Sdd,Smd,LTdd,LTdm,elmdq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,q
      !   -(q, Sdd·rho) - Tau_d(q·L*dd, Sdd·rho) - Tau_m(q·L*dm, Smd·rho) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Sdd(e%mnode)
      real(rp),    intent(in)    :: Smd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm(:,:)),Smd)

      do jnode=1,e%pnode
            elmdq(1:e%pnode,jnode) = elmdq(1:e%pnode,jnode)&
                                  -((e%shape(1:e%pnode,e%igaus)+LTdd(1:e%pnode))*Sdd(jnode)&
                                    + aux(1:e%pnode,jnode))*dvolu
      end do

   end subroutine nsc_elmbdq_reac

   subroutine nsc_elmbmq_reac(e,dvolu,Smm,LTdm,elmmq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,q
      !   - Tau_m(q·L*dm, Smm·mom) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Smm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmq(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmq(inode,:,jnode) = elmmq(inode,:,jnode) - &
             matmul(LTdm(:,inode),Smm(:,:,jnode))*dvolu 
         end do
      end do

   end subroutine nsc_elmbmq_reac

   subroutine nsc_elmbeq_reac(e,dvolu,See,LTde,elmeq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,q
      !  - Tau_e(q·L*de, See·ene) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: See(e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeq(e%mnode,e%mnode)

      integer(ip)                :: jnode

      do jnode=1,e%pnode
            elmeq(1:e%pnode,jnode) = elmeq(1:e%pnode,jnode)&
                                  -LTde(1:e%pnode)*See(jnode)*dvolu 
      end do

   end subroutine nsc_elmbeq_reac

   subroutine nsc_elmbdn_reac(e,dvolu,Sdd,Smd,Sed,LTmd,LTmm,LTme,elmdn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,n 
      !  -(n, Smd·rho) - Tau_m(n·L*md, Sdd·rho) - Tau_m(n·L*mm, Smd·rho) - Tau_e(n·L*me, Sed·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Sdd(e%mnode)
      real(rp),    intent(in)    :: Smd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Sed(e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdn(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Smd)
         do jnode=1,e%pnode
            elmdn(idime,1:e%pnode,jnode) = elmdn(idime,1:e%pnode,jnode) - &
             (e%shape(1:e%pnode,e%igaus)*Smd(idime,jnode) + &
              LTmd(idime,1:e%pnode)*Sdd(jnode) +&
              aux(1:e%pnode,jnode)+LTme(idime,1:e%pnode)*Sed(jnode))*dvolu
         end do
      end do

   end subroutine nsc_elmbdn_reac

   subroutine nsc_elmbmn_reac(e,dvolu,Smm,Sem,LTmm,LTme,elmmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block mom,n 
      !  -(n, Smm·mom) - Tau_m(n·L*mm, Smm·mom) - Tau_e(n·L*me, Sem·mom)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Smm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Sem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmn(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: idime,jdime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      do idime=1,e%ndime
         do jdime=1,e%ndime
            aux = matmul(transpose(LTmm(idime,:,:)),Smm(:,jdime,:))
            do jnode=1,e%pnode
               elmmn(idime,1:e%pnode,jdime,jnode) = elmmn(idime,1:e%pnode,jdime,jnode) - &
                (aux(1:e%pnode,jnode)+e%shape(1:e%pnode,e%igaus)*Smm(idime,jdime,jnode) + &
                LTme(idime,1:e%pnode)*Sem(jdime,jnode))*dvolu 
            end do
         end do
      end do
      
   end subroutine nsc_elmbmn_reac

   subroutine nsc_elmben_reac(e,dvolu,See,LTme,elmen)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,n 
      !  - Tau_e(n·L*me, See·ene)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: See(e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmen(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         do jnode=1,e%pnode
            elmen(idime,1:e%pnode,jnode) = elmen(idime,1:e%pnode,jnode) - &
              LTme(idime,1:e%pnode)*See(jnode)*dvolu
         end do
      end do

   end subroutine nsc_elmben_reac

   subroutine nsc_elmbdg_reac(e,dvolu,Sdd,Smd,Sed,LTed,LTem,LTee,elmdg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,g 
      !  -(g, Sed·rho) - Tau_d(g·L*ed, Sdd·rho) - Tau_m(g·L*em, Smd·rho) - Tau_e(g·L*ee,Sed·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Sdd(e%mnode)
      real(rp),    intent(in)    :: Smd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Sed(e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdg(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      integer(ip)                :: jnode
      
      aux = matmul(transpose(LTem),Smd)      
      do jnode=1,e%pnode
         elmdg(1:e%pnode,jnode) = elmdg(1:e%pnode,jnode) - &
         ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*Sed(jnode) + &
          aux(1:e%pnode,jnode)+LTed(1:e%pnode)*Sdd(jnode))*dvolu 
      end do

   end subroutine nsc_elmbdg_reac

   subroutine nsc_elmbmg_reac(e,dvolu,Smm,Sem,LTem,LTee,elmmg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusive lhs terms for ASGS for block mom,g 
      !  -(g, Sem·mom) - Tau_m(g·L*em, Smm·mom)- Tau_e(g·L*ee, Sem·mom)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Smm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Sem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmg(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmg(inode,:,jnode) = elmmg(inode,:,jnode) - &
             (matmul(LTem(:,inode),Smm(:,:,jnode)) + &
             (e%shape(inode,e%igaus)+LTee(inode))*Sem(:,jnode))*dvolu
         end do
      end do

   end subroutine nsc_elmbmg_reac

   subroutine nsc_elmbeg_reac(e,dvolu,See,LTee,elmeg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the reaction lhs terms for ASGS for block rho,g 
      !  -(g, See·ene) - Tau_e(g·L*ee,See·ene)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: See(e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeg(e%mnode,e%mnode)

      integer(ip)                :: jnode
      
      do jnode=1,e%pnode
         elmeg(1:e%pnode,jnode) = elmeg(1:e%pnode,jnode) - &
         (e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*See(jnode)*dvolu 
      end do

   end subroutine nsc_elmbeg_reac

end module


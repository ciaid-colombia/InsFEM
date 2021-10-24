module Mod_NSCompressiblePrimitiveElement
   use typre
   use Mod_Element
   
   implicit none

   real(rp), parameter :: zeroc = epsilon(0.0_rp)
contains
   
  subroutine nsc_pr_elmbdq(e,dvolu,Edd,Emd,Eed,LTdd,LTdm,LTde,elmdq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block pre,q 
      !    (q, Atdd·pre/dt + Add·pre) + Tau_d(q·L*dd, Atdd·pre/dt + Add·pre) 
      !    + Tau_m(q·L*dm, Atmd·pre/dt + Amd·pre) 
      !    + Tau_e(q·L*de, Ated·pre/dt + Aed·pre) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edd(e%mnode)
      real(rp),    intent(in)    :: Emd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eed(e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm),Emd)
      
      do jnode=1,e%pnode
         elmdq(1:e%pnode,jnode) = elmdq(1:e%pnode,jnode) + &
            ((e%shape(1:e%pnode,e%igaus)+LTdd(1:e%pnode))*Edd(jnode) + &
              aux(1:e%pnode,jnode) + LTde(1:e%pnode)*Eed(jnode))*dvolu
      end do

   end subroutine nsc_pr_elmbdq
   
   subroutine nsc_pr_elmbmq(e,dvolu,Edm,Emm,Eem,LTdd,LTdm,LTde,elmmq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block vel,q 
      !  + (q, Adm·vel) + Tau_d(q·L*dd, Adm·vel) + Tau_m(q·L*dm, Atmm·vel/dt + Amm·vel) 
      !  + Tau_e(q·L*de, Atem·vel/dt + Aem·vel) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Emm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmq(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmq(inode,:,jnode) = elmmq(inode,:,jnode) + &
             ((e%shape(inode,e%igaus)+LTdd(inode))*Edm(:,jnode) + &
             matmul(LTdm(:,inode),Emm(:,:,jnode)) + &
             LTde(inode)*Eem(:,jnode))*dvolu 
         end do
      end do

   end subroutine nsc_pr_elmbmq
  
   subroutine nsc_pr_elmbeq(e,dvolu,Ede,Eme,Eee,LTdd,LTdm,LTde,elmeq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block tem,q 
      !   (q, Atde·tem/dt + Ade·tem) + Tau_d(q·L*dd, Atde·tem/dt + Ade·tem) 
      !    + Tau_m(q·L*dm, Atme·tem/dt + Ame·tem) 
      !    + Tau_e(q·L*de, Atee·tem/dt + Aee·tem) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ede(e%mnode)
      real(rp),    intent(in)    :: Eme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eee(e%mnode)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeq(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      aux = matmul(transpose(LTdm),Eme)

      do jnode=1,e%pnode
         elmeq(1:e%pnode,jnode) = elmeq(1:e%pnode,jnode) + &
            ((e%shape(1:e%pnode,e%igaus)+LTdd(1:e%pnode))*Ede(jnode) + &
              aux(1:e%pnode,jnode) + LTde(1:e%pnode)*Eee(jnode))*dvolu
      end do

   end subroutine nsc_pr_elmbeq

   subroutine nsc_pr_elmbdn(e,dvolu,Edd,Emd,Eed,LTmd,LTmm,LTme,elmdn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block pre,n 
      !  (n, Atmd·pre/dt + Amd·pre) + Tau_d(n·L*md, Atdd·pre/dt + Add·pre) 
      !  + Tau_m(n·L*mm, Atmd·pre/dt + Amd·pre) 
      !  + Tau_e(n·L*me, Ated·pre/dt + Aed·pre)  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edd(e%mnode)
      real(rp),    intent(in)    :: Emd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eed(e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdn(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Emd)
         do jnode=1,e%pnode
            elmdn(idime,1:e%pnode,jnode) = elmdn(idime,1:e%pnode,jnode) + &
            (e%shape(1:e%pnode,e%igaus)*Emd(idime,jnode) + &
             LTmd(idime,1:e%pnode)*Edd(jnode)  + &
             aux(1:e%pnode,jnode) + LTme(idime,1:e%pnode)*Eed(jnode))*dvolu 
         end do
      end do

   end subroutine nsc_pr_elmbdn

   subroutine nsc_pr_elmbmn(e,dvolu,Edm,Emm,Eem,LTmd,LTmm,LTme,elmmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block vel,n 
      !  (n, Atmm·vel/dt + Amm·vel) + Tau_d(n·L*md, Adm·vel) 
      !  + Tau_m(n·L*mm, Atmm·vel/dt + Amm·vel) 
      !  + Tau_e(n·L*me, Atem·vel/dt + Aem·vel)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Emm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmn(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: idime,inode,jdime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      do idime=1,e%ndime
         do inode=1,e%pnode
            do jdime=1,e%ndime
               do jnode=1,e%pnode
                  elmmn(idime,inode,jdime,jnode) = elmmn(idime,inode,jdime,jnode) + &
                  (e%shape(inode,e%igaus)*Emm(idime,jdime,jnode)+ &
                   dot_product(LTmm(idime,:,inode),Emm(:,jdime,jnode)) +&
                   LTmd(idime,inode)*Edm(jdime,jnode) + &
                   LTme(idime,inode)*Eem(jdime,jnode))*dvolu 
               end do
            end do
         end do
      end do

   end subroutine nsc_pr_elmbmn

   subroutine nsc_pr_elmben(e,dvolu,Ede,Eme,Eee,LTmd,LTmm,LTme,elmen)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block tem,n 
      !  (n, Atme·tem/dt + Ame·tem) + Tau_d(n·L*md, Atde·tem/dt + Ade·tem) 
      !  + Tau_m(n·L*mm, Atme·tem/dt + Ame·tem) 
      !  + Tau_e(n·L*me, Atee·tem/dt + Aee·tem)  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ede(e%mnode)
      real(rp),    intent(in)    :: Eme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eee(e%mnode)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmen(e%ndime,e%mnode,e%mnode)

      integer(ip)                :: idime,jnode
      real(rp)                   :: aux(e%mnode,e%mnode)

      
      do idime=1,e%ndime
         aux = matmul(transpose(LTmm(idime,:,:)),Eme)
         do jnode=1,e%pnode
            elmen(idime,1:e%pnode,jnode) = elmen(idime,1:e%pnode,jnode) + &
            (e%shape(1:e%pnode,e%igaus)*Eme(idime,jnode) + aux(1:e%pnode,jnode) + &
             LTmd(idime,1:e%pnode)*Ede(jnode) + LTme(idime,1:e%pnode)*Eee(jnode))*dvolu 
         end do
      end do

   end subroutine nsc_pr_elmben

   subroutine nsc_pr_elmbdg(e,dvolu,Edd,Emd,Eed,LTed,LTem,LTee,elmdg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block pre,g 
      !  (g, Ated·pre/dt + Aed·pre) + Tau_d(g·L*ed, Atdd·pre/dt + Add·pre) 
      !  + Tau_m(g·L*em, Atmd·pre/dt + Amd·pre)
      !  + Tau_e(g·L*ee, Ated·pre/dt + Aed·pre)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edd(e%mnode)
      real(rp),    intent(in)    :: Emd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eed(e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdg(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      integer(ip)                :: jnode
      
      aux = matmul(transpose(LTem),Emd)      
      do jnode=1,e%pnode
         elmdg(1:e%pnode,jnode) = elmdg(1:e%pnode,jnode) + &
         ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*Eed(jnode) + &
          LTed(1:e%pnode)*Edd(jnode) + aux(1:e%pnode,jnode))*dvolu 
      end do

   end subroutine nsc_pr_elmbdg

   subroutine nsc_pr_elmbmg(e,dvolu,Edm,Emm,Eem,LTed,LTem,LTee,elmmg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block vel,g 
      !  (g, Atem·vel/dt + Aem·vel) + Tau_d(g·L*ed, Adm·vel) 
      !  + Tau_m(g·L*em, Atmm·vel/dt + Amm·vel)
      !  + Tau_e(g·L*ee, Atem·vel/dt + Aem·vel)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Edm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Emm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmmg(e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmmg(inode,:,jnode) = elmmg(inode,:,jnode) + &
            ((e%shape(inode,e%igaus)+LTee(inode))*Eem(:,jnode) + &
              LTed(inode)*Edm(:,jnode)+matmul(LTem(:,inode),Emm(:,:,jnode)))*dvolu 
         end do
      end do

   end subroutine nsc_pr_elmbmg

   subroutine nsc_pr_elmbeg(e,dvolu,Ede,Eme,Eee,LTed,LTem,LTee,elmeg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block tem,g 
      !  (g, Atee·tem/dt + Aee·tem) + Tau_d(g·L*ed, Atde·tem/dt + Ade·tem) 
      !  + Tau_m(g·L*em, Atme·tem/dt + Ame·tem)
      !  + Tau_e(g·L*ee, Atee·tem/dt + Aee·tem)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: Ede(e%mnode)
      real(rp),    intent(in)    :: Eme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: Eee(e%mnode)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmeg(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux(e%mnode,e%mnode)
      
      aux = matmul(transpose(LTem),Eme)      
      do jnode=1,e%pnode
         elmeg(1:e%pnode,jnode) = elmeg(1:e%pnode,jnode) + &
         ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*Eee(jnode) + &
           LTed(1:e%pnode)*Ede(jnode) + aux(1:e%pnode,jnode))*dvolu 
      end do

   end subroutine nsc_pr_elmbeg

   subroutine nsc_pr_elmrhd(e,dvolu,Erhd,Erhm,Erhe,LTdd,LTdm,LTde,elrhd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for mass equation
      !    + (q, Atdd·pre_n/dt + Atde·tem_n/dt) + (q, fd) 
      !    + Tau_d(q·L*dd, Atdd·pre_n/dt + Atde·tem_n/dt) + Tau_d(q·L*dd, fd) 
      !    + Tau_m(q·L*dm, Atmd·pre_n/dt + Atmm·vel_n/dt + Atme·tem_n/dt) + Tau_m(q.L*dm, fm)
      !    + Tau_e(q·L*de, Ated·pre_n/dt + Atem·vel_n/dt + Atee·tem_n/dt) + Tau_e(q.L*de, fe)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: Erhd,Erhe
      real(rp),    intent(in)    :: Erhm(e%ndime)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhd(e%mnode)

      real(rp)                   :: aux(e%mnode)

      aux = matmul(Erhm,LTdm)

      elrhd(1:e%pnode) = elrhd(1:e%pnode) + &
       ((e%shape(1:e%pnode,e%igaus)+LTdd(1:e%pnode))*Erhd + &
        aux(1:e%pnode) + LTde(1:e%pnode)*Erhe)*dvolu 

   end subroutine nsc_pr_elmrhd

   subroutine nsc_pr_elmrhm(e,dvolu,Erhd,Erhm,Erhe,LTmd,LTmm,LTme,elrhm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for momentum equation
      !    + (n, Atmd·pre_n/dt + Atmm·vel_n/dt + Atme·tem_n/dt) + (n, fm) 
      !    + Tau_d(n·L*md, Atdd·pre_n/dt + Atde·tem_n/dt) + Tau_d(n·L*md, fd) 
      !    + Tau_m(n·L*mm, Atmd·pre_n/dt + Atmm·vel_n/dt + Atme·tem_n/dt) + Tau_m(n.L*mm, fm)
      !    + Tau_e(n·L*me, Ated·pre_n/dt + Atem·vel_n/dt + Atee·tem_n/dt) + Tau_e(n.L*me, fe)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: Erhd,Erhe
      real(rp),    intent(in)    :: Erhm(e%ndime)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhm(e%ndime,e%mnode)

      integer(ip)                :: idime
      real(rp)                   :: aux(e%mnode)

      do idime=1,e%ndime
         aux(:) =  matmul(transpose(LTmm(idime,:,:)),Erhm) 
         elrhm(idime,1:e%pnode) = elrhm(idime,1:e%pnode) + &
          (e%shape(1:e%pnode,e%igaus)*Erhm(idime) + &
          LTmd(idime,1:e%pnode)*Erhd + LTme(idime,1:e%pnode)*Erhe + &
           aux(1:e%pnode))*dvolu 
      end do

   end subroutine nsc_pr_elmrhm

   subroutine nsc_pr_elmrhe(e,dvolu,Erhd,Erhm,Erhe,LTed,LTem,LTee,elrhe)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for energy
      !    + (g, Ated·pre_n/dt + Atem·vel_n/dt + Atee·tem_n/dt) + (g, fe) 
      !    + Tau_d(g·L*ed, Atdd·pre_n/dt + Atde·tem_n/dt) + Tau_d(g·L*ed, fd) 
      !    + Tau_m(g·L*em, Atmd·pre_n/dt + Atmm·vel_n/dt + Atme·tem_n/dt) + Tau_m(g.L*em, fm)
      !    + Tau_e(g·L*ee, Ated·pre_n/dt + Atem·vel_n/dt + Atee·tem_n/dt) + Tau_e(g.L*ee, fe)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: Erhd,Erhe
      real(rp),    intent(in)    :: Erhm(e%ndime)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhe(e%mnode)

      real(rp)                   :: aux(e%mnode)

      aux =  matmul(Erhm,LTem) 

      elrhe(1:e%pnode) = elrhe(1:e%pnode) + &
       ((e%shape(1:e%pnode,e%igaus)+LTee(1:e%pnode))*Erhe + &
       LTed(1:e%pnode)*Erhd + aux(1:e%pnode))*dvolu 

   end subroutine nsc_pr_elmrhe

!Radial Damping Diffusion
   subroutine nsc_pr_elmbdq_diff(e,dvolu,KGdd,elmdq)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusive lhs terms for ASGS for block ene,g 
      !  + (d_k g, Kdd_k·rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      type(FiniteElement) :: e
      
      real(rp),    intent(in)    :: KGdd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      real(rp)                   :: aux(e%mnode,e%mnode)
      integer(ip)                :: jnode

      aux = matmul(transpose(e%cartd),KGdd) 
      do jnode = 1,e%pnode
         elmdq(1:e%pnode,jnode) = elmdq(1:e%pnode,jnode) + &
                    aux(1:e%pnode,jnode)*dvolu 
      end do

   end subroutine nsc_pr_elmbdq_diff
end module


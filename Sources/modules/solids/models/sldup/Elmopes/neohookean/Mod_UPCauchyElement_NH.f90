submodule(Mod_CauchyElement) UPCauchyElement_NH

   implicit none

contains

   !-------------LHS------------

   module subroutine up_elmVU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !  -tau_s*(symgrad(v_h)_ij,
      !                          1/2(F_jK*du_i/dX_K + F_iK*du_j/dX_K) 
      !                          - 1/nd*F_lK*du_l/dX_K*ii_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),             intent(in)    :: Fmat(nd,nd),Finv(nd,nd)
      real(rp),             intent(inout) :: elmat(nd,e%pnode,nd,e%pnode)
      integer(ip)           :: inode,jnode,idime
      real(rp)              :: B(nd,nd),B1(nd,nd)
      real(rp)              :: BB(nd,nd),Bt(nd,nd),Bts(nd,nd)
      real(rp)              :: aux,aux2,aux3

      aux  = (1.0_rp/(J*4.0_rp))
      aux2 = (1.0_rp/(J*nd*2.0_rp))
      aux3 = (1.0_rp/(4.0_rp*mu))

      do inode=1,e%pnode
          do jnode=1,e%pnode

              B   = 0.0_rp
              B1  = 0.0_rp
              Bt  = 0.0_rp
              Bts = 0.0_rp

              !1/4(dv_j/dx_i+dv_i/dx_j,F_iK*du_j/dX_K + F_jK*du_i/dX_K)
              call sld_calculate_2Fdu_tensor(e,nd,inode,jnode,Fmat,B1)

              !1/2(dv_j/dx_i+dv_i/dx_j,1/nd*F_lK*du_l/dX_K*ii_ij)
              call sld_calculate_traceFdu_tensor(e,nd,inode,jnode,Fmat,Bt)

              B   = aux*B1 - aux2*Bt

              !1/2(dv_j/dx_i+dv_i/dx_j,(J/2)F^-1_Kl*du_l/dX_K*S_ij)
              call sld_calculate_traceFduS_tensor(e,nd,tn,inode,jnode,Finv,gdev,Bts)
              B = B - aux3*Bts

              elmat(1:nd,inode,1:nd,jnode) = elmat(1:nd,inode,1:nd,jnode)+B(:,:)*dvol

          end do
      end do

   end subroutine up_elmVU_NH

   module subroutine up_elmVP_NH(e,nd,dvol,mu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                (d(v_i)/dx_i),p) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,mu
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,1,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: G(nd),GN(nd),aux

      aux = 1.0_rp/(2.0_rp*mu)
      do inode=1,e%pnode

          call sld_calculateB_inditial(e,nd,inode,G)

          do jnode=1,e%pnode

              GN = 0.0_rp
              GN = aux*G*e%shape(jnode,e%igaus)*dvol

              elmat(1:nd,inode,1,jnode) = elmat(1:nd,inode,1,jnode) + GN(:)

          end do
      end do

   end subroutine up_elmVP_NH

end submodule UPCauchyElement_NH

submodule(Mod_CauchyElement) SUPCauchyElement_NH_SGS
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental stabilization for
      ! the three field neo-hookean formulation. This will be 
      ! un-understadable without prior reading and understanding
      ! of the formulation. Go read! Seriously, its freaking crazy! OMG.
      !
      !-----------------------------------------------------------------------

   implicit none

contains

   !-------------------Tau_u matrix calculation----------

   module subroutine sup_elm_usgs_VU_NH_dynsgs(e,nd,tn,dvol,tautau,densi,dtinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental matrix 
      !            (v_h,i*(1-tau_k^-1*tau_t), rho*U,i^n+1/dt2)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau,densi,dtinv
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,nd,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: BG(nd,nd)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              call sld_calculateVTauTau_U(e,nd,tn,inode,jnode,tautau,densi,dtinv,BG)
              
              elmat(1:nd,inode,1:nd,jnode) = elmat(1:nd,inode,1:nd,jnode) + BG(:,:)*dvol

          end do
      end do

   end subroutine sup_elm_usgs_VU_NH_dynsgs

   module subroutine sup_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol,tautau,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental matrix 
      !            (v_h,i*(1-tau_k^-1*tau_t), -div(S_ij),i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,tn,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: BG(nd,tn)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              call sld_calculateVTauTau_S(e,nd,tn,inode,jnode,tautau,BG)
              
              elmat(1:nd,inode,1:tn,jnode) = elmat(1:nd,inode,1:tn,jnode) + BG(:,:)*dvol

          end do
      end do

   end subroutine sup_elm_usgs_VS_NH_dynsgs

   module subroutine sup_elm_usgs_VP_NH_dynsgs(e,nd,tn,dvol,tautau,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental matrix 
      !            (v_h,i*(1-tau_k^-1*tau_t), -grad(P),i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,1,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: BG(nd)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              call sld_calculateVTauTau_P(e,nd,tn,inode,jnode,tautau,BG)
              
              elmat(1:nd,inode,1,jnode) = elmat(1:nd,inode,1,jnode) + BG(:)*dvol

          end do
      end do

   end subroutine sup_elm_usgs_VP_NH_dynsgs

   module subroutine sup_elm_usgs_EU_NH_dyn(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,densi,dtinv2,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*(2(F_ik*d(xi_ij)/dx_k - 1/nd*F_lk*d(xi_ll)/dx_k)
      !       -(J/mu)F^-1_lk*d(xi_ij*S_ij)/dx_k),rho*U_i*dtinv2)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd),dtinv2
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu,densi
      real(rp),    intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: Bi(tn,nd),B(tn,nd),U(nd)
      real(rp)                   :: BG(tn),aux

      aux = densi*dtinv2*dvol

      do inode=1,e%pnode

          B  = 0.0_rp
          Bi = 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - (J/(2.0_rp*mu))*B

          Bi = aux*Bi

          do jnode=1,e%pnode

              !du/dt
              call sld_calculateU_inditial(e,nd,jnode,U)
              
              BG  = matmul(Bi,U)

              elmat(1:tn,inode,1,jnode) = elmat(1:tn,inode,1,jnode) + BG(:)

          end do
      end do

   end subroutine sup_elm_usgs_EU_NH_dyn

   module subroutine sup_elm_usgs_QU_NH_dyn(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,densi,dtinv2,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !         -tau_u*((2mu/lam*nd)*F_iK*dq/dX_K 
      !     -d/dX_K((Jp/lam)-1)q*F^(-1)_Ki,rho*U_i*dtinv2)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd),densi,dtinv2
      real(rp),   intent(inout)  :: elmat(1,e%pnode,nd,e%pnode)
      integer(ip)                :: idime,jdime,inode,jnode
      real(rp)                   :: GG,G(nd),U(nd),C(nd),aux

      aux=densi*dtinv2*dvol

      do inode=1,e%pnode

          G = 0.0_rp

          !(2mu/lam*nd)*F_ik*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = ((2.0*mu)/(lam*nd))*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gpress(1),Finv,C)
          G = G - C

          G = aux*G

          do jnode=1,e%pnode

              !du/dt
              call sld_calculateU_inditial(e,nd,jnode,U)

              GG = dot_product(G,U)

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + GG

          end do
      end do

   end subroutine sup_elm_usgs_QU_NH_dyn

   module subroutine sup_elm_usgs_ES_NH(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !       -(J/2mu)F^-1_Ki*d(xi_lm*S_ml)/dX_K),-div(S_iq)_i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu
      real(rp),   intent(inout)  :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: BB(tn,tn),Bt(tn,tn),Bts(tn,tn),aux,aux2
      
      aux  = 1.0_rp/nd
      aux2 = J/(2.0_rp*mu)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              BB  = 0.0_rp
              Bt  = 0.0_rp
              Bts = 0.0_rp

              call sld_calculate_FdN_tensor(e,nd,tn,inode,jnode,Fmat,BB)
              call sld_calculate_FdN_tensor_trace(e,nd,tn,inode,jnode,Fmat,Bt)

              !Negative sign because of -div(S_ij)
              BB = -BB + aux*Bt

              call sld_calculate_trace_dNS_F_tensor(e,nd,tn,inode,jnode,Finv,grsig,gdev,Bts)
              !Positive sign because of -div(S_ij)
              BB = BB + aux2*Bts

              BB = BB*dvol

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode)+BB(:,:)

          end do

      end do

   end subroutine sup_elm_usgs_ES_NH

   module subroutine sup_elm_usgs_EP_NH(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*(2(F_ik*d(xi_ij)/dx_k - 1/nd*F_lk*d(xi_ll)/dx_k)
      !       -(J/mu)F^-1_lk*d(xi_ij*S_ij)/dx_k),-grad(p)_i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu
      real(rp),    intent(inout) :: elmat(tn,e%pnode,1,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: Bi(tn,nd),B(tn,nd),G(nd)
      real(rp)                   :: BG(tn)

      do inode=1,e%pnode

          B  = 0.0_rp
          Bi = 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - (J/(2.0_rp*mu))*B

          Bi = Bi*dvol

          do jnode=1,e%pnode

              !-grad(p)
              call sld_calculateB_inditial(e,nd,jnode,G)

              BG  = -matmul(Bi,G)

              elmat(1:tn,inode,1,jnode) = elmat(1:tn,inode,1,jnode) + BG(:)

          end do
      end do

   end subroutine sup_elm_usgs_EP_NH

   module subroutine sup_elm_usgs_QS_NH(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !     -tau_u*((2mu/lam*nd)*F_iK*d(q_h)/dX_K
      !     - d/dX_K(((Jp/lam)-1)*F^-1_Ki*q_h),-div(S_ij)_i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,tn,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: C(nd),G(nd),B(nd,tn)
      real(rp)                   :: B_aux(tn,nd)
      real(rp)                   :: GB(tn),GBt(tn)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              GB = 0.0_rp

              !(2.0*mu/lam*nd)*F_iK*d(q_h)/dX_K
              call sld_calculate_dqF_vector(e,nd,tn,inode,jnode,Fmat,GBt)
              GB = ((2.0*mu)/(lam*nd))*GBt

              !-d/dX_K((Jp/lam)-1)*F^-1_Ki*q_h)
              call sld_calculate_dqFinv_vector(e,nd,tn,inode,jnode,J,lam,gradp,gpress(1),Finv,GBt)
              GB = GB - GBt

              !Negative sign because of -div(S_ij)
              elmat(1,inode,1:tn,jnode) = elmat(1,inode,1:tn,jnode) - GB*dvol

          end do
      end do

   end subroutine sup_elm_usgs_QS_NH

   module subroutine sup_elm_usgs_QP_NH(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !         -tau_u*((2mu/lam*nd)*F_iK*dq/dX_K 
      !     -d/dX_K((Jp/lam)-1)q*F^(-1)_Ki,-grad(p_h))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                :: idime,jdime,inode,jnode
      real(rp)                   :: GG,G(nd),G_j(nd),C(nd)

      do inode=1,e%pnode

          G = 0.0_rp

          !(2mu/lam*nd)*F_ik*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = ((2.0*mu)/(lam*nd))*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gpress(1),Finv,C)
          G = G - C

          G = G*dvol

          do jnode=1,e%pnode

              !-grad(p_h)
              call sld_calculateB_inditial(e,nd,jnode,G_j)

              GG = -dot_product(G,G_j)

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + GG

          end do
      end do

   end subroutine sup_elm_usgs_QP_NH

   !-------------------Tau_s matrix calculation----------

   module subroutine sup_elm_ssgs_VU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !  -tau_s*(symgrad(v_h)_ij,
      !            2(-F_jK*du_i/dX_K+1/nd*F_lK*du_l/dX_K*ii_ij)
      !            +(J/mu)F^-1_lK*du_l/dX_K*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),             intent(inout) :: elmat(nd,e%pnode,nd,e%pnode)
      integer(ip)           :: inode,jnode,idime
      real(rp)              :: B(nd,nd),B1(nd,nd)
      real(rp)              :: BB(nd,nd),Bt(nd,nd),Bts(nd,nd)
      real(rp)              :: aux,aux2

      aux  = (1.0_rp/(2.0_rp*nd))
      aux2 = (J/(4.0_rp*mu))

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

              B   = -0.25_rp*B1 + aux*Bt

              !1/2(dv_j/dx_i+dv_i/dx_j,(J/2mu)F^-1_Kl*du_l/dX_K*S_ij)
              call sld_calculate_traceFduS_tensor(e,nd,tn,inode,jnode,Finv,gdev,Bts)
              B = B + aux2*Bts

              elmat(1:nd,inode,1:nd,jnode) = elmat(1:nd,inode,1:nd,jnode)+B(:,:)*dvol

          end do
      end do

   end subroutine sup_elm_ssgs_VU_NH

   module subroutine sup_elm_ssgs_VS_NH(e,nd,tn,dvol,J,mu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !               -tau_s(symgrad(v_h)_ij,(J/2mu)*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu
      real(rp),             intent(inout) :: elmat(nd,e%pnode,tn,e%pnode)
      integer(ip)                         :: inode,jnode
      real(rp)                            :: B(nd,tn)
      real(rp)                            :: BN(nd,tn),aux

      aux = J/(2.0_rp*mu)

      do inode=1,e%pnode

          !grad(v_h)_ij
          call sld_calculateBt(e,nd,tn,inode,B)

          do jnode=1,e%pnode

              BN  = 0.0_rp

              !(grad(v_h)_ij,(J/2mu)*S_ij)
              BN  = aux*B*e%shape(jnode,e%igaus)*dvol

              elmat(1:nd,inode,1:tn,jnode) = elmat(1:nd,inode,1:tn,jnode) + BN(1:nd,1:tn)

          end do
      end do

   end subroutine sup_elm_ssgs_VS_NH

   module subroutine sup_elm_ssgs_EU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_s*((J/2mu)*xi_ij,2(-F_iK*du_j/dX_K+1/nd*F_lK*du_l/dX_K*ii_ij)
      !                       +(J/mu)F^-1_lK*du_l/dX_K*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),    intent(in)    :: Fmat(nd,nd),Finv(nd,nd)
      real(rp),    intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: ii(tn,tn)
      real(rp)                   :: B(tn,nd),B1(tn,nd),Bts(tn,nd)
      real(rp)                   :: Bt(tn,nd)
      real(rp)                   :: N,aux1,aux2

      aux1 = 1.0_rp/nd
      aux2 = J/(2.0_rp*mu)

      do inode=1,e%pnode

          ! (J/2mu)*xi_ij
          N   = aux2*e%shape(inode,e%igaus)*dvol

          do jnode=1,e%pnode

              B   = 0.0_rp
              B1  = 0.0_rp
              Bt  = 0.0_rp
              Bts = 0.0_rp

              !  -0.5*(F_jK*du_i/dx_K+F_iK*du_j/dx_K)
              call sld_calculate_2Fdu(e,nd,tn,jnode,Fmat,B1)

              !  (1/nd*F_lK*du_l/dX_K*ii_ij)
              call sld_calculate_traceFdu(e,nd,tn,jnode,Fmat,Bt)

              B = (-0.5_rp*B1+aux1*Bt)*N

              !(J/2mu)F^-1_lK*du_l/dX_K*S_ij
              call sld_calculate_traceFduS(e,nd,tn,jnode,Finv,gdev,Bts)
              B = B + aux2*N*Bts

              elmat(1:tn,inode,1:nd,jnode) = elmat(1:tn,inode,1:nd,jnode) + B(:,:)

          end do
      end do

   end subroutine sup_elm_ssgs_EU_NH

   module subroutine sup_elm_ssgs_ES_NH(e,nd,tn,dvol,J,mu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                 -tau_s*(J/mu*xi_ij,J/mu*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip) :: inode,jnode
      real(rp)    :: ii(tn,tn)
      real(rp)    :: Ni,Nj,aux
      real(rp)    :: NN(tn,tn)

      call getDiagTensor2Txy(nd,tn,ii)
      aux = (J/(2.0_rp*mu))
      ii = aux*aux*ii*dvol

      do inode=1,e%pnode

          ! J/2mu*xi_ij
          Ni   = e%shape(inode,e%igaus)

          do jnode=1,e%pnode

              Nj  = e%shape(jnode,e%igaus)

              ! (J/2mu*xi_ij,J/2mu*S_ij)
              NN  = Ni*Nj*ii

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode) + NN(1:tn,1:tn)

          end do
      end do

   end subroutine sup_elm_ssgs_ES_NH

   !-------------------Tau_p matrix calculation----------

   module subroutine sup_elm_psgs_VU_NH(e,nd,dvol,J,mu,lam,gpress,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !             -tau_p*(div(v_h),((Jp/lam)-1)*Finv_Kl-(2mu/lam*nd)*Fmat_lK)*du_l/dx_K)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,gpress(1),J,mu,lam
      real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,nd,e%pnode)
      integer(ip)                      :: inode,jnode,idime,jdime
      real(rp)                         :: GG(nd,nd),Gi(nd),B(nd)
      real(rp)                         :: mat(nd,nd),Baux(nd)

      do inode=1,e%pnode

          !div(v_h)
          call sld_calculateB_inditial(e,nd,inode,Gi)

          do jnode=1,e%pnode

              mat = 0.0_rp
              B   = 0.0_rp

              mat = (J/lam)*gpress(1)*Finv
              mat = mat - Finv
              call sld_calculate_tracePressFinv(e,nd,jnode,mat,Baux)

              B = Baux

              mat = 0.0_rp
              mat = -((2.0*mu)/(lam*nd))*Fmat
              call sld_calculate_tracePressFmat(e,nd,jnode,mat,Baux)
              B = B + Baux

              !(div(v_h),mat_lk*du_l/dx_k)
              GG = 0.0_rp
              do idime=1,nd
                  do jdime=1,nd
                      GG(idime,jdime) = Gi(idime)*B(jdime)
                  end do
              end do

              elmat(1:nd,inode,1:nd,jnode) = elmat(1:nd,inode,1:nd,jnode)+GG*dvol

          end do
      end do

   end subroutine sup_elm_psgs_VU_NH

   module subroutine sup_elm_psgs_VP_NH(e,nd,dvol,J,lam,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                -tau_p*(div(v_h),(J/lam)*p_h)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,lam,J
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,1,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: G(nd),GN(nd)

      G  = 0.0_rp
      GN = 0.0_rp

      do inode=1,e%pnode

          !div(v_h)
          call sld_calculateB_inditial(e,nd,inode,G)

          do jnode=1,e%pnode

              !(div(v_h),(J/lam)*p_h)
              GN = G*e%shape(jnode,e%igaus)*dvol

              elmat(1:nd,inode,1,jnode) = elmat(1:nd,inode,1,jnode) + (J/lam)*GN(:)

          end do
      end do

   end subroutine sup_elm_psgs_VP_NH

   module subroutine sup_elm_psgs_QU_NH(e,nd,dvol,J,mu,lam,gpress,Fmat,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      !This routine computes the elemental material matrix 
      !-tau_p*((J/lam)*q_h,((Jp/lam)-1)*Finv_Kl-(2mu/lam*nd)*Fmat_lK)*du_l/dX_K)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,J,mu,lam,gpress(1)
      real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),    intent(inout)       :: elmat(1,e%pnode,nd,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: G(nd),NG(nd),N,mat(nd,nd),Gaux(nd)


      do inode=1,e%pnode

          !(J/lam)*q_h
          N = 0.0_rp
          N = (J/lam)*e%shape(inode,e%igaus)*dvol

          do jnode=1,e%pnode

              G   = 0.0_rp
              NG  = 0.0_rp
              mat = 0.0_rp

              mat= (J/lam)*gpress(1)*Finv
              mat= mat - Finv
              call sld_calculate_tracePressFinv(e,nd,jnode,mat,Gaux)

              G = Gaux

              mat= 0.0_rp
              mat= -((2.0*mu)/(lam*nd))*Fmat
              call sld_calculate_tracePressFmat(e,nd,jnode,mat,Gaux)
              G = G + Gaux

              !(J/lam)*q_h*mat_lK*du_l/dX_K
              NG = N*G

              elmat(1,inode,1:nd,jnode) = elmat(1,inode,1:nd,jnode) + NG(1:nd)

          end do
      end do

   end subroutine sup_elm_psgs_QU_NH

   module subroutine sup_elm_psgs_QP_NH(e,nd,dvol,J,lam,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !             -tau_p*((J/lam)*q_h,(J/lam)*p_h)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,J,lam
      real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: NN,Ni,Nj

      do inode=1,e%pnode

          !(J/lam)*q_h
          Ni = e%shape(inode,e%igaus)*(J/lam)*dvol

          do jnode=1,e%pnode

              !(J/lam)*p_h
              Nj = e%shape(jnode,e%igaus)*(J/lam)
              NN = Ni*Nj

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + NN

          end do
      end do

   end subroutine sup_elm_psgs_QP_NH

   !---------------------------Higher order derivatives----------------------

   module subroutine sup_elm_usgs_EU_NH_dyn_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,densi,dtinv2,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*(2(F_ik*d(xi_ij)/dx_k - 1/nd*F_lk*d(xi_ll)/dx_k)
      !       -(J/mu)F^-1_lk*d(xi_ij*S_ij)/dx_k),rho*U_i*dtinv2)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd),dtinv2
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J,densi
      real(rp),    intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: Bi(tn,nd),B(tn,nd),U(nd)
      real(rp)                   :: BG(tn),aux

      aux = densi*dtinv2*dvol

      do inode=1,e%pnode

          B  = 0.0_rp
          Bi = 0.0_rp

          call sld_calculate_NdF(e,nd,tn,inode,divFmat,B)
          Bi = B

          !    xi_lm*S_lm*d(F^-1_Ki)/dX_K
          call sld_calculate_trace_NS_dF(e,nd,tn,inode,J,gdev,Jder,divFinv,B)
          Bi = Bi - (1.0_rp/(2.0_rp*mu))*B

          Bi = aux*Bi

          do jnode=1,e%pnode

              !du/dt
              call sld_calculateU_inditial(e,nd,jnode,U)
              
              BG  = matmul(Bi,U)

              elmat(1:tn,inode,1,jnode) = elmat(1:tn,inode,1,jnode) + BG(:)

          end do
      end do

   end subroutine sup_elm_usgs_EU_NH_dyn_nonli

   module subroutine sup_elm_usgs_QU_NH_dyn_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,densi,dtinv2,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !         -tau_u*((2mu/lam*nd)*F_iK*dq/dX_K 
      !     -d/dX_K((Jp/lam)-1)q*F^(-1)_Ki,rho*U_i*dtinv2)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd),densi,dtinv2
      real(rp),   intent(inout)  :: elmat(1,e%pnode,nd,e%pnode)
      integer(ip)                :: idime,jdime,inode,jnode
      real(rp)                   :: GG,G(nd),U(nd),C(nd),aux

      aux=densi*dtinv2*dvol

      do inode=1,e%pnode

          G = 0.0_rp

          !(2mu/lam*nd)*F_ik*dq/dX_K
          call sld_calculate_qdF(e,nd,tn,inode,divFmat,C)
          G = ((2.0*mu)/(lam*nd))*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_qdFinv(e,nd,tn,inode,J,lam,gpress(1),Jder,divFinv,C)
          G = G - C

          G = aux*G

          do jnode=1,e%pnode

              !du/dt
              call sld_calculateU_inditial(e,nd,jnode,U)

              GG = dot_product(G,U)

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + GG

          end do
      end do

   end subroutine sup_elm_usgs_QU_NH_dyn_nonli

   module subroutine sup_elm_usgs_ES_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*((xi_ij*d(F_jK)/dX_K-1/nd*xi_ll*d(F_iK)/dX_K)
      !       -(1/2mu)*d(J*F^-1_Ki)/dX_K*xi_lm*S_ml),-div(S_iq)_i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J
      real(rp),   intent(inout)  :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: BB(tn,tn),Bt(tn,tn),Bts(tn,tn),aux,aux2
      
      aux  = 1.0_rp/nd
      aux2 = J/(2.0_rp*mu)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              BB  = 0.0_rp
              Bt  = 0.0_rp
              Bts = 0.0_rp

              call sld_calculate_NdF_tensor(e,nd,tn,inode,jnode,divFmat,BB)
              call sld_calculate_NdF_tensor_trace(e,nd,tn,inode,jnode,divFmat,Bt)

              !Negative sign because of -div(S_ij)
              BB = -BB + aux*Bt

              call sld_calculate_trace_NS_dF_tensor(e,nd,tn,inode,jnode,gdev,J,Jder,divFinv,Bts)
              !Positive sign because of -div(S_ij)
              BB = BB + aux2*Bts

              BB = BB*dvol

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode)+BB(:,:)

          end do

      end do

   end subroutine sup_elm_usgs_ES_NH_nonli

   module subroutine sup_elm_usgs_EP_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*(2(xi_ij*d(F_iK)/dx_K - 1/nd*xi_ll*d(F_lK)/dx_K)
      !       -(J/mu)xi_ij*S_ij*d(F^-1_lK)/dx_K),-grad(p)_i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J
      real(rp),    intent(inout) :: elmat(tn,e%pnode,1,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: Bi(tn,nd),B(tn,nd),G(nd)
      real(rp)                   :: BG(tn)

      do inode=1,e%pnode

          B  = 0.0_rp
          Bi = 0.0_rp

          call sld_calculate_NdF(e,nd,tn,inode,divFmat,B)
          Bi = B

          !    xi_lm*S_lm*d(F^-1_Ki)/dX_K
          call sld_calculate_trace_NS_dF(e,nd,tn,inode,J,gdev,Jder,divFinv,B)
          Bi = Bi - (1.0_rp/(2.0_rp*mu))*B

          Bi = Bi*dvol

          do jnode=1,e%pnode

              !-grad(p)
              call sld_calculateB_inditial(e,nd,jnode,G)

              BG  = -matmul(Bi,G)

              elmat(1:tn,inode,1,jnode) = elmat(1:tn,inode,1,jnode) + BG(:)

          end do
      end do

   end subroutine sup_elm_usgs_EP_NH_nonli

   module subroutine sup_elm_usgs_QS_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !     -tau_u*((2mu/lam*nd)*F_iK*d(q_h)/dX_K
      !     - d/dX_K(((Jp/lam)-1)*F^-1_Ki*q_h),-div(S_ij)_i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,tn,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: C(nd),G(nd),B(nd,tn)
      real(rp)                   :: B_aux(tn,nd)
      real(rp)                   :: GB(tn),GBt(tn)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              GB = 0.0_rp

              !(2.0*mu/lam*nd)*F_iK*d(q_h)/dX_K
              call sld_calculate_qdF_vector(e,nd,tn,inode,jnode,divFmat,GBt)
              GB = ((2.0*mu)/(lam*nd))*GBt

              !-d/dX_K((Jp/lam)-1)*F^-1_Ki*q_h)
              call sld_calculate_qdFinv_vector(e,nd,tn,inode,jnode,J,lam,gpress(1),Jder,divFinv,GBt)
              GB = GB - GBt

              !Negative sign because of -div(S_ij)
              elmat(1,inode,1:tn,jnode) = elmat(1,inode,1:tn,jnode) - GB*dvol

          end do
      end do

   end subroutine sup_elm_usgs_QS_NH_nonli

   module subroutine sup_elm_usgs_QP_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !         -tau_u*((2mu/lam*nd)*q*dF_iK/dX_K 
      !     -(Jp/lam)-1)q*d/dX_K(F^(-1)_Ki),-grad(p_h))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: GG,G(nd),G_j(nd),C(nd)

      do inode=1,e%pnode

          G = 0.0_rp

          !(2mu/lam*nd)*F_ik*dq/dX_K
          call sld_calculate_qdF(e,nd,tn,inode,divFmat,C)
          G = ((2.0*mu)/(lam*nd))*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_qdFinv(e,nd,tn,inode,J,lam,gpress(1),Jder,divFinv,C)
          G = G - C

          G = G*dvol

          do jnode=1,e%pnode

              !-grad(p_h)
              call sld_calculateB_inditial(e,nd,jnode,G_j)

              GG = -dot_product(G,G_j)

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + GG

          end do
      end do

   end subroutine sup_elm_usgs_QP_NH_nonli

end submodule SUPCauchyElement_NH_SGS

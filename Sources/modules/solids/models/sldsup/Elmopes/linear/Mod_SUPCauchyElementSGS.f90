submodule(Mod_CauchyElement) SUPCauchyElementSGS
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental stabilization for
      ! the three field linear formulation. 
      ! Taken from 'A mixed three-field FE formulation for stress accurate 
      !             analysis including the incompressible limit' 
      !             Authors: Chiumenti, Cervera, Codina
      !             doi    : 10.1016/j.cma.2014.08.004
      !
      !-----------------------------------------------------------------------

   implicit none

contains

    !-------------------Tau_u matrix calculation----------

   module subroutine sup_elm_usgs_ES(e,ndime,tn,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B*B_t,   (div y,div s)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: dvol
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode,i,j
      real(rp)    :: B_i(tn,ndime),B_j(ndime,tn)
      real(rp)    :: B_aux(tn,ndime)
      real(rp)    :: BB(tn,tn)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              BB  = 0.0_rp

              call sld_calculate_FdN_tensor_lin(e,ndime,tn,inode,jnode,BB)

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode)+BB(:,:)*dvol
          end do
      end do

   end subroutine sup_elm_usgs_ES

   module subroutine sup_elm_usgs_EP(e,ndime,tn,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B*G, (div y, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime,tn
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(tn,e%pnode,1,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: B(tn,ndime),G(ndime)
      real(rp)                         :: BG(tn)

      do inode=1,e%pnode

          call sld_calculateB(e,ndime,tn,inode,B)

          do jnode=1,e%pnode

              call sld_calculateB_inditial(e,ndime,jnode,G)

              BG = matmul(B,G)

              elmat(1:tn,inode,1,jnode) = elmat(1:tn,inode,1,jnode) + BG(:)*dvol

          end do
      end do

   end subroutine sup_elm_usgs_EP

   module subroutine sup_elm_usgs_QS(e,ndime,tn,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    G*B
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime,tn
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(1,e%pnode,tn,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: G(ndime),B(ndime,tn)
      real(rp)                         :: GB(tn)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              GB = 0.0_rp

              !d/dX_i(q_h)
              call sld_calculate_dqFinv_vector_lin(e,ndime,tn,inode,jnode,GB)

              !Negative sign because of -div(S_ij)
              elmat(1,inode,1:tn,jnode) = elmat(1,inode,1:tn,jnode) - GB*dvol

          end do
      end do


   end subroutine sup_elm_usgs_QS

   module subroutine sup_elm_usgs_QP(e,ndime,tn,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    G_t*G
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime,tn
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: GG,G_i(ndime),G_j(ndime)

      do inode=1,e%pnode

          call sld_calculateB_inditial(e,ndime,inode,G_i)
          do jnode=1,e%pnode

              call sld_calculateB_inditial(e,ndime,jnode,G_j)

              GG = dot_product(G_i,G_j)*dvol

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + GG

          end do
      end do

   end subroutine sup_elm_usgs_QP

    !-------------------Tau_s matrix calculation----------

   module subroutine sup_elm_ssgs_VU(e,ndime,tn,dvol,c_dev,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B_t*P*Ñ_s 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: dvol
      real(rp),             intent(in)    :: c_dev(tn,tn)
      real(rp),             intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp)    :: CB(tn,ndime)
      real(rp)    :: B_j(tn,ndime),B_i(ndime,tn)
      real(rp)    :: B_aux(tn,ndime)
      real(rp)    :: BCB(ndime,ndime)

      do inode=1,e%pnode

          call sld_calculateB(e,ndime,tn,inode,B_aux)
          B_i   = transpose(B_aux)

          do jnode=1,e%pnode

              call sld_calculateB(e,ndime,tn,jnode,B_j)

              CB  = matmul(c_dev,B_j)
              BCB = matmul(B_i,CB)*dvol

              !Notice the minus sign
              elmat(1:ndime,inode,1:ndime,jnode) = elmat(1:ndime,inode,1:ndime,jnode) - BCB(:,:)

          end do
      end do

   end subroutine sup_elm_ssgs_VU

   module subroutine sup_elm_ssgs_VS(e,ndime,tn,dvol,P,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B_t*P*Ñ_s 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: P(tn,tn),dvol
      real(rp),             intent(inout) :: elmat(ndime,e%pnode,tn,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp) :: B(ndime,tn)
      real(rp) :: auxB(tn,ndime)
      real(rp) :: BPN(ndime,tn)
      real(rp) :: ii(tn,tn),N(tn,tn)
      real(rp) :: PN(tn,tn)

      call getDiagTensor(ndime,tn,ii)

      do inode=1,e%pnode

          call sld_calculateB(e,ndime,tn,inode,auxB)
          B   = transpose(auxB)

          do jnode=1,e%pnode

              N = e%shape(jnode,e%igaus)*ii

              PN  = matmul(P,N)
              BPN = matmul(B,PN)*dvol

              elmat(1:ndime,inode,1:tn,jnode) = elmat(1:ndime,inode,1:tn,jnode) + BPN(1:ndime,1:tn)

          end do
      end do

   end subroutine sup_elm_ssgs_VS

   module subroutine sup_elm_ssgs_EU(e,ndime,tn,dvol,P,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    N_t*P*B,    (y,P:grad_s u) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: P(tn,tn),dvol
      real(rp),             intent(inout) :: elmat(tn,e%pnode,ndime,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp)    :: ii(tn,tn)
      real(rp)    :: N(tn,tn)
      real(rp)    :: B(tn,ndime)
      real(rp)    :: PB(tn,ndime)
      real(rp)    :: NPB(tn,ndime)

      call getDiagTensor(ndime,tn,ii)

      do inode=1,e%pnode

          N    = e%shape(inode,e%igaus)*(ii)

          do jnode=1,e%pnode

              call sld_calculateB(e,ndime,tn,jnode,B)

              PB  = matmul(P,B)
              NPB = matmul(N,PB)*dvol

              elmat(1:tn,inode,1:ndime,jnode) = elmat(1:tn,inode,1:ndime,jnode) + NPB(:,:)

          end do
      end do

   end subroutine sup_elm_ssgs_EU

   module subroutine sup_elm_ssgs_ES(e,ndime,tn,dvol,D,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    -Ñ_s_t*D*Ñ_s ,     (y,D:s)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: D(tn,tn),dvol
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp)    :: ii(tn,tn)
      real(rp)    :: N_i(tn,tn) ,N_j(tn,tn)
      real(rp)    :: B(tn,ndime)
      real(rp)    :: NN(tn,tn)

      !call getDiagTensor(ndime,tn,ii)

      !do inode=1,e%pnode

      !    N_i   = e%shape(inode,e%igaus)*(ii)

      !    do jnode=1,e%pnode

      !        N_j  = e%shape(jnode,e%igaus)*ii

      !        DN   = matmul(D,N_j)
      !        NDN  = matmul(N_i,DN)*dvol

      !        !Notice the minus sign
      !        elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode) - NDN(1:tn,1:tn)

      !    end do
      !end do
      
      call getDiagTensor(ndime,tn,ii)

      do inode=1,e%pnode

          ! xi_ij
          N_i = e%shape(inode,e%igaus)*ii

          do jnode=1,e%pnode

              N_j  = e%shape(jnode,e%igaus)*ii

              ! (xi_ij,S_ij)
              NN = matmul(D,N_j)
              NN = matmul(N_i,NN)*dvol

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode) - NN(1:tn,1:tn)

          end do
      end do

   end subroutine sup_elm_ssgs_ES

   !-------------------Tau_p matrix calculation----------

   module subroutine sup_elm_psgs_VU(e,ndime,dvol,K,G,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    G_t*G
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime
      real(rp),    intent(in)          :: dvol,K,G
      real(rp),    intent(inout)       :: elmat(ndime,e%pnode,ndime,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode,i,j
      real(rp)                         :: GG(ndime,ndime),G_i(ndime),G_j(ndime),C_vol

      C_vol=min(K,(2.0_rp/3.0_rp)*G)

      do inode=1,e%pnode

          call sld_calculateB_inditial(e,ndime,inode,G_i)
          do jnode=1,e%pnode

              call sld_calculateB_inditial(e,ndime,jnode,G_j)

              GG = 0.0_rp
              do i=1,ndime
                  do j=1,ndime
                      GG(i,j) = G_i(i)*G_j(j)
                  end do
              end do

              elmat(1:ndime,inode,1:ndime,jnode) = elmat(1:ndime,inode,1:ndime,jnode) - C_vol*GG*dvol

          end do
      end do

   end subroutine sup_elm_psgs_VU

   module subroutine sup_elm_psgs_VP(e,ndime,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    GN_p
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(ndime,e%pnode,1,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: G(ndime),GN(ndime)

      G  = 0.0_rp
      GN = 0.0_rp

      do inode=1,e%pnode

          call sld_calculateB_inditial(e,ndime,inode,G)

          do jnode=1,e%pnode

              GN = G*e%shape(jnode,e%igaus)*dvol

              elmat(1:ndime,inode,1,jnode) = elmat(1:ndime,inode,1,jnode) + GN(1:ndime)

          end do
      end do

   end subroutine sup_elm_psgs_VP

   module subroutine sup_elm_psgs_QU(e,ndime,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    N_iG_j
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(1,e%pnode,ndime,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: G(ndime),NG(ndime),N

      G  = 0.0_rp
      NG = 0.0_rp
      N  = 0.0_rp

      do inode=1,e%pnode
          N = e%shape(inode,e%igaus)
          do jnode=1,e%pnode

              call sld_calculateB_inditial(e,ndime,jnode,G)
              NG = N*G*dvol

              elmat(1,inode,1:ndime,jnode) = elmat(1,inode,1:ndime,jnode) + NG(1:ndime)

          end do
      end do

   end subroutine sup_elm_psgs_QU

   module subroutine sup_elm_psgs_QP(e,ndime,dvol,D_vol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    -D_vol*N_j_t*N_i
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: ndime
      real(rp),    intent(in)          :: dvol,D_vol
      real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                      :: idime,jdime,inode,jnode
      real(rp)                         :: DNN,N_i,N_j

      do inode=1,e%pnode
          N_i = e%shape(inode,e%igaus)
          do jnode=1,e%pnode

              N_j = e%shape(jnode,e%igaus)
              !Notice the minus sign
              DNN = -D_vol*N_i*N_j*dvol

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + DNN

          end do
      end do

   end subroutine sup_elm_psgs_QP

end submodule SUPCauchyElementSGS

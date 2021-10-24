submodule(Mod_CauchyElement) SUPCauchyElement

   implicit none

contains

   module subroutine sup_elmVS(e,ndime,tn,dvol,P,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B_t*P*Ñ_s 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: dvol,P(tn,tn)
      real(rp),             intent(inout) :: elmat(ndime,e%pnode,tn,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp) :: B(ndime,tn)
      real(rp) :: auxB(tn,ndime)
      real(rp) :: BPN(ndime,tn)
      real(rp) :: ii(tn,tn),N(tn,tn)
      real(rp) :: PN(tn,tn)

      call getDiagTensor2Txy(ndime,tn,ii)

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

   end subroutine sup_elmVS

   module subroutine sup_elmVP(e,ndime,dvol,elmat)
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

   end subroutine sup_elmVP

   module subroutine sup_elmEU(e,ndime,tn,dvol,P,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    N_t*P*B,   (y,P:grad_s(u))
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

      call getDiagTensor2Txy(ndime,tn,ii)

      do inode=1,e%pnode

          N    = e%shape(inode,e%igaus)*(ii)

          do jnode=1,e%pnode

              call sld_calculateB(e,ndime,tn,jnode,B)

              PB  = matmul(P,B)
              NPB = matmul(N,PB)*dvol

              elmat(1:tn,inode,1:ndime,jnode) = elmat(1:tn,inode,1:ndime,jnode) + NPB(:,:)

          end do
      end do

   end subroutine sup_elmEU

   module subroutine sup_elmES(e,ndime,tn,dvol,D,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    -N_s_t*D*N_s,    (y,D:s)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: ndime,tn
      real(rp),             intent(in)    :: dvol,D(tn,tn)
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip) :: idime,jdime,inode,jnode
      real(rp)    :: ii(tn,tn)
      real(rp)    :: N_i(tn,tn),N_j(tn,tn)
      real(rp)    :: B(tn,ndime)
      real(rp)    :: DN(tn,tn)
      real(rp)    :: NDN(tn,tn)

      call getDiagTensor2Txy(ndime,tn,ii)

      do inode=1,e%pnode

          N_i    = e%shape(inode,e%igaus)*(ii)

          do jnode=1,e%pnode

              N_j  = e%shape(jnode,e%igaus)*ii

              DN   = matmul(D,N_j)
              NDN  = matmul(N_i,DN)*dvol

              !Notice the minus sign
              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode) - NDN(1:tn,1:tn)

          end do
      end do

   end subroutine sup_elmES

   module subroutine sup_elmQU(e,ndime,dvol,elmat)
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
      real(rp)                         :: G(ndime),NG(ndime)

      G  = 0.0_rp
      NG = 0.0_rp

      do inode=1,e%pnode
          do jnode=1,e%pnode

              call sld_calculateB_inditial(e,ndime,jnode,G)
              NG = G*e%shape(inode,e%igaus)*dvol

              elmat(1,inode,1:ndime,jnode) = elmat(1,inode,1:ndime,jnode) + NG(1:ndime)

          end do
      end do

   end subroutine sup_elmQU

   module subroutine sup_elmQP(e,ndime,dvol,D_vol,elmat)
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
      real(rp)                         :: DNN

      do inode=1,e%pnode
          do jnode=1,e%pnode

              !Notice the minus sign
              DNN = -D_vol*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + DNN

          end do
      end do

   end subroutine sup_elmQP

end submodule SUPCauchyElement

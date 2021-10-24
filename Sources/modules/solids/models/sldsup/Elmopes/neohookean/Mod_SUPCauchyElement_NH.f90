submodule(Mod_CauchyElement) SUPCauchyElement_NH

   implicit none

contains

    !-------------RHS------------
   module subroutine sup_rhs_U_NH(e,nd,tn,dvol,gdev,gpress,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms 
      !    -0.5*(d(v_i)/dx_j+d(v_j)/dx_i),S_ij) - (d(v_i)/dx_i),p)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement)       :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: dvol,gdev(tn),gpress(1)
      real(rp),    intent(inout) :: elrhs(nd,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: Bi(nd,tn)
      real(rp)                   :: rhs(nd),G(nd)

      do inode=1,e%pnode

          rhs = 0.0_rp

          !-grad(v_i)(S_ij)
          call sld_calculateBt(e,nd,tn,inode,Bi)
          rhs = rhs - matmul(Bi,gdev)

          !-div(v_i)(P)
          call sld_calculateB_inditial(e,nd,inode,G)
          rhs = rhs - G*gpress(1)

          !Assembly
          elrhs(:,inode) = elrhs(:,inode) + rhs*dvol
      end do

   end subroutine sup_rhs_U_NH

   module subroutine sup_rhs_S_NH(e,nd,tn,dvol,mu,J,b,gdev,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !       (xi_ij,0.5*(b-(b_ii/nd)*ii)) - (xi_ij,(J/2mu)S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,gdev(tn),b(nd,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode,jnode
      real(rp)    :: ii(tn,tn),i(nd,nd),b_ii
      real(rp)    :: N(tn),C(nd,nd),AU2(nd,nd),Cvoigt(tn)
      real(rp)    :: AU1(tn),rhs(tn),rhs2(tn)

      call getDiagTensor2Txy(nd,tn,ii)
      call get2ndIITensor(nd,i)

      call trace(nd,b,b_ii)

      !Stress tensor term: (b-(b_ii/nd)*ii)
      C = (b - (b_ii/nd)*i)
      !C = b

      call getVoigtStress(tn,nd,C,Cvoigt)

      !-A(U): -(J/2mu)S_ij)+0.5(b-(b_ii/nd)*ii)
      AU1 = 0.5_rp*(Cvoigt-(J/mu)*gdev)*dvol
      N   = matmul(ii,AU1)

      do inode=1,e%pnode

          rhs  = 0.0_rp

          !Testing against shape function
          rhs = e%shape(inode,e%igaus)*N

          !Assembly
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_S_NH

   module subroutine sup_rhs_P_NH(e,nd,tn,dvol,mu,lam,J,b,gpress,elrhs,divdisp)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !   (q_h,-(J/lam)*p + log(J) + (mu/lam)*((b_ii/nd)-1.0_rp))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,lam,J,dvol,gpress(1),b(nd,nd),divdisp
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                         :: inode,jnode
      real(rp)                            :: N,b_ii,AU,rhs

      call trace(nd,b,b_ii)

      !-A(U)
      AU = -(J/lam)*gpress(1)+(mu/lam)*((b_ii/nd)-1.0_rp)
      AU = AU + log(J) !+ divdisp

      do inode=1,e%pnode

          rhs = 0.0_rp
          rhs = e%shape(inode,e%igaus)*AU*dvol

          !Assembly
          elrhs(1,inode) = elrhs(1,inode) + rhs

      end do

   end subroutine sup_rhs_P_NH

   !-------------LHS------------

   module subroutine sup_elmVS_NH(e,nd,tn,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                (d(v_i)/dx_j),S_ij) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol
      real(rp),             intent(inout) :: elmat(nd,e%pnode,tn,e%pnode)
      integer(ip)                         :: inode,jnode
      real(rp)                            :: Baux(nd,tn),B(nd,tn)

      do inode=1,e%pnode

          call sld_calculateBt(e,nd,tn,inode,B)

          do jnode=1,e%pnode

              Baux = 0.0_rp

              Baux = B*e%shape(jnode,e%igaus)*dvol

              elmat(1:nd,inode,1:tn,jnode) = elmat(1:nd,inode,1:tn,jnode) + Baux(1:nd,1:tn)

          end do
      end do

   end subroutine sup_elmVS_NH

   module subroutine sup_elmVP_NH(e,nd,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                (d(v_i)/dx_i),p) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,1,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: G(nd),GN(nd)

      do inode=1,e%pnode

          call sld_calculateB_inditial(e,nd,inode,G)

          do jnode=1,e%pnode

              GN = 0.0_rp
              GN = G*e%shape(jnode,e%igaus)*dvol

              elmat(1:nd,inode,1,jnode) = elmat(1:nd,inode,1,jnode) + GN(:)

          end do
      end do

   end subroutine sup_elmVP_NH

   module subroutine sup_elmEU_NH(e,nd,tn,dvol,J,mu,Finv,Fmat,gdev,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !  (xi_ij,(-F_iK*du_j/dx_K+1/nd*F_lK*du_l/dx_K*ii_ij)
      !       +(J/mu)F^-1_lK*du_l/dx_K*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,J,dvol,Finv(nd,nd),Fmat(nd,nd),gdev(tn)
      real(rp),             intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
      integer(ip) :: inode,jnode
      real(rp)    :: B(tn,nd),B1(tn,nd),Bt(tn,nd)
      real(rp)    :: Bts(tn,nd),N,aux1,aux2

      aux1 = 1.0_rp/nd
      aux2 = J/(2.0_rp*mu)

      do inode=1,e%pnode

          ! xi_ij
          N  = e%shape(inode,e%igaus)*dvol

          do jnode=1,e%pnode

              B   = 0.0_rp
              B1  = 0.0_rp
              Bt  = 0.0_rp
              Bts = 0.0_rp

              !  -(F_jK*du_i/dx_K+F_iK*du_j/dx_K)
              call sld_calculate_2Fdu(e,nd,tn,jnode,Fmat,B1)

              !  (1/nd*F_lK*du_l/dx_K*ii_ij)
              call sld_calculate_traceFdu(e,nd,tn,jnode,Fmat,Bt)

              B = (-0.5_rp*B1 + aux1*Bt)*N

              !(J/2mu)F^-1_Kl*du_l/dX_K*S_ij
              call sld_calculate_traceFduS(e,nd,tn,jnode,Finv,gdev,Bts)
              B = B + aux2*N*Bts

              elmat(1:tn,inode,1:nd,jnode) = elmat(1:tn,inode,1:nd,jnode) + B(:,:)

          end do
      end do

   end subroutine sup_elmEU_NH

   module subroutine sup_elmES_NH(e,nd,tn,dvol,J,mu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !         (xi_ij,(J/2mu)*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
      integer(ip) :: inode,jnode
      real(rp)    :: ii(tn,tn)
      real(rp)    :: Ni,Nj
      real(rp)    :: NN(tn,tn),aux

      aux = J/(2.0_rp*mu)

      call getDiagTensor2Txy(nd,tn,ii)
      ii = aux*ii*dvol

      do inode=1,e%pnode

          Ni = e%shape(inode,e%igaus)

          do jnode=1,e%pnode

              Nj = e%shape(jnode,e%igaus)

              NN = Ni*Nj*ii

              elmat(1:tn,inode,1:tn,jnode) = elmat(1:tn,inode,1:tn,jnode) + NN(1:tn,1:tn)

          end do
      end do

   end subroutine sup_elmES_NH

   module subroutine sup_elmQU_NH(e,nd,dvol,J,mu,lam,Finv,Fmat,gpress,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !   (q_h,((Jp/lam)-1)*F^-1_Kl-(2mu/lam*nd)*F_lK)*du_l/dx_K)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,mu,lam,J,gpress(1)
      real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),    intent(inout)       :: elmat(1,e%pnode,nd,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: Ni(nd),mat(nd,nd)
      real(rp)                         :: B(nd),Baux(nd)

      do inode=1,e%pnode

          Ni = 0.0_rp
          Ni = e%shape(inode,e%igaus)*dvol

          do jnode=1,e%pnode

              B = 0.0_rp

              mat = (J/lam)*gpress(1)*Finv
              mat = mat - Finv
              call sld_calculate_tracePressFinv(e,nd,jnode,mat,Baux)

              B = Baux

              mat = 0.0_rp
              mat=-((2.0*mu)/(lam*nd))*Fmat
              call sld_calculate_tracePressFmat(e,nd,jnode,mat,Baux)
              B = B + Baux

              elmat(1,inode,1:nd,jnode) = elmat(1,inode,1:nd,jnode) + Ni*B(1:nd)

          end do
      end do

   end subroutine sup_elmQU_NH

   module subroutine sup_elmQP_NH(e,nd,dvol,J,lam,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                     (q_h,(J/lam)*p)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,J,lam
      real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
      integer(ip)                      :: inode,jnode
      real(rp)                         :: Ni,Nj

      do inode=1,e%pnode

          Ni = e%shape(inode,e%igaus)

          do jnode=1,e%pnode

              Nj = e%shape(jnode,e%igaus)

              elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + (J/lam)*Ni*Nj*dvol

          end do
      end do

   end subroutine sup_elmQP_NH

end submodule SUPCauchyElement_NH

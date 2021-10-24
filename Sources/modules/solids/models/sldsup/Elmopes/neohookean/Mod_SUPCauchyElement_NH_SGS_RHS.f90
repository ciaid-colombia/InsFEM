submodule(Mod_CauchyElement) SUPCauchyElement_NH_SGS_RHS
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
    !-------------ADJOINT EXTERNAL FORCES------------

   module subroutine sup_rhs_usgs_S_NH_extforces(e,nd,tn,dvol,mu,J,tau,Fmat,Finv,grsig,gdev,elrhs,elext)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !           -(J/2mu)F^-1_Ki*d/dX_K(xi_lm*S_ml),
      !           - F_ext,i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,tau,grsig(tn,nd),gdev(tn)
      real(rp),             intent(in)    :: elext(nd),Fmat(nd,nd),Finv(nd,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd)
      real(rp)    :: B(tn,nd),aux(nd),BS(tn),aux2
      real(rp)    :: AU(nd),rhs(tn),fext(tn)

      !external forces
      aux = tau*elext(:)*dvol 
      aux2= J/(2.0_rp*mu)

      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          fext= 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - aux2*B

          fext = matmul(Bi,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) - fext(:)

      end do

   end subroutine sup_rhs_usgs_S_NH_extforces

   module subroutine sup_rhs_usgs_P_NH_extforces(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,elrhs,elext)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,
      !              - F_ext,i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)           :: nd,tn
      real(rp),   intent(in)           :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)           :: lam,J,mu,tau
      real(rp),   intent(in)           :: gpress(1),elext(nd),gradp(nd)
      real(rp),   intent(inout)        :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux,aux1(nd)
      real(rp)                  :: fext

      !external forces
      aux1 = tau*elext(:)*dvol 

      !Constant part of adjoint
      aux   = (2.0_rp*mu)/(lam*nd)

      do inode=1,e%pnode

          G   = 0.0_rp
          fext= 0.0_rp

          !-----For the remaining terms----

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = aux*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gpress(1),Finv,C)
          G = G-C

          fext = dot_product(G,aux1)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) - fext

      end do

   end subroutine sup_rhs_usgs_P_NH_extforces

    !-------------ADJOINT RHS------------

   module subroutine sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    tau_s*(symgrad(v_h)_ij,
      !                    -1/2(b - (b_ll/nd)*i)_ij+(J/2mu)*S_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau
      real(rp),             intent(in)    :: AU(tn)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bt(nd,tn),aux(tn)
      real(rp)    :: rhs(nd)

      aux  = AU*tau*dvol

      do inode=1,e%pnode

          rhs  = 0.0_rp

          !grad(v_h)_ij
          call sld_calculateBt(e,nd,tn,inode,Bt)

          !Testing against adjoint
          rhs = matmul(Bt,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_ssgs_U_NH

   module subroutine sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !             tau_p*(div(v_h),
      !    (J/lam)*P - ln(J) - (mu/lam)*(b_ll/nd-1))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau,AU
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      real(rp)                            :: N(nd),rhs(nd),aux
      integer(ip)                         :: inode

      !Pre multiplying constant part of adjoint
      aux = AU*tau*dvol

      do inode=1,e%pnode

          rhs = 0.0_rp

          !div(v_h)
          call sld_calculateB_inditial(e,nd,inode,N)

          !Testing against adjoint
          rhs = N*aux

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_psgs_U_NH

   module subroutine sup_rhs_usgs_S_NH(e,nd,tn,dvol,mu,J,tau,grsig,gdev,Fmat,Finv,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !           -(J/2mu)F^-1_Ki*d(xi_lm*S_ml)/dX_K),
      !           -div(S_ij),i - grad(P),i - F_ext,i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),             intent(in)    :: grsig(tn,nd),AU(nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd),G(nd),aux1,aux(nd)
      real(rp)    :: B(tn,nd),DS(nd),BS(tn)
      real(rp)    :: rhs(tn)

      !A(U)
      aux1= J/(2.0_rp*mu)
      aux = tau*AU*dvol

      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          rhs = 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - aux1*B

          !Testing against adjoint
          rhs = matmul(Bi,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_usgs_S_NH

   module subroutine sup_rhs_ssgs_S_NH(e,nd,tn,dvol,mu,J,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    tau_s*(J/2mu*xi_ij,
      !                     J/2mu*S_ij-0.5*(b-(b_ll/nd)*ii)_ij)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: AU(tn),dvol,tau,J,mu
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: N(tn),rhs(tn),aux(tn),ii(tn,tn)

      call getDiagTensor2Txy(nd,tn,ii)

      !A(U): (J/2mu*xi_ij,(J/2mu)S_ij-1/2(b-(b_ll/nd)*ii)_ij)
      aux = 0.5_rp*(J/mu)*AU*tau*dvol

      !Pre multiplying constant part of adjoint
      N  = matmul(ii,aux)

      do inode=1,e%pnode

          rhs  = 0.0_rp

          ! xi_ij
          !Testing against adjoint
          rhs = e%shape(inode,e%igaus)*N

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_ssgs_S_NH

   module subroutine sup_rhs_usgs_P_NH(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,
      !              -div(S_ij)_i - grad(p)_i - F_ext,i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau,AU(nd)
      real(rp),   intent(in)    :: gpress(1),gradp(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux1,aux(nd)
      real(rp)                  :: rhs

      !Constant part of adjoint
      aux1   = (2.0_rp*mu)/(lam*nd)

      !A(U)
      aux= tau*AU*dvol

      do inode=1,e%pnode

          G   = 0.0_rp
          rhs = 0.0_rp

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = aux1*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gpress(1),Finv,C)
          G = (G-C)

          !Testing against adjoint
          rhs = dot_product(G,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) + rhs

      end do

   end subroutine sup_rhs_usgs_P_NH

   module subroutine sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,J,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    -tau_p((J/lam)q_h,
      !       (J/lam)*gpress(1)-ln(J)-(mu/lam)*(b_ll/nd-1))
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: lam,J,dvol
      real(rp),             intent(in)    :: tau,AU
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                         :: inode
      real(rp)                            :: aux,rhs

      !A(U)
      !Pre multiplying constant part of adjoint
      aux = tau*(J/lam)*AU*dvol

      do inode=1,e%pnode

          rhs = 0.0_rp

          !Testing against adjoint
          rhs = e%shape(inode,e%igaus)*aux

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) + rhs

      end do

   end subroutine sup_rhs_psgs_P_NH

   !-------------ADJOINT RHS Dynamic SGS------------

   module subroutine sup_rhs_usgs_U_NH_AU_dynsgs(e,nd,tn,dvol,tautau,elext,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                     (Vh*(1-tau^-1*tau_t),F_i - AU)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tautau
      real(rp),             intent(in)    :: elext(nd),AU(nd)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      integer(ip) :: inode
      real(rp)    :: B(nd),aux(nd)
      real(rp)    :: rhs(nd)

      aux = (elext-AU)*(1.0_rp-tautau)*dvol

      do inode=1,e%pnode

          rhs = 0.0_rp

          call sld_calculateVTauTauRHS_U(e,nd,inode,aux,rhs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) - rhs(:)

      end do

   end subroutine 

   module subroutine sup_rhs_usgs_U_NH_SGSACCEL_dynsgs(e,nd,tn,dvol,dtinv2,densi,tautau,vesgs,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                     (Vh*tau^-1*tau_t,rho*Usgs_i^n/dtinv2)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tautau
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      integer(ip) :: inode
      real(rp)    :: usgs(nd),rhs(nd)

      usgs = densi*dtinv2*dvol*vesgs*tautau

      do inode=1,e%pnode

          rhs = 0.0_rp

          call sld_calculateVTauTauRHS_U(e,nd,inode,usgs,rhs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine 

   module subroutine sup_rhs_usgs_S_NH_dynsgs(e,nd,tn,dvol,mu,J,tau,grsig,gdev,Fmat,Finv,dtinv2,densi,vesgs,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !           -(J/2mu)F^-1_Ki*d(xi_lm*S_ml)/dX_K),
      !                      Usgs_i^n)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd),grsig(tn,nd)
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd),G(nd),usgs(nd)
      real(rp)    :: B(tn,nd),DS(nd),BS(tn)
      real(rp)    :: AU(nd),rhs(tn)

      usgs = densi*dtinv2*dvol*vesgs

      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          rhs = 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - (J/(2.0_rp*mu))*B

          !Testing against adjoint
          rhs = matmul(Bi,usgs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) - tau*rhs(:)

      end do

   end subroutine 
  
   module subroutine sup_rhs_usgs_P_NH_dynsgs(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,dtinv2,densi,vesgs,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,
      !                                            Usgs_i^n)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: gpress(1),gradp(nd)
      real(rp),   intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux,usgs(nd)
      real(rp)                  :: AU(nd),rhs

      !Constant part of adjoint
      aux   = (2.0_rp*mu)/(lam*nd)

      usgs = densi*dtinv2*dvol*vesgs

      do inode=1,e%pnode

          G   = 0.0_rp
          rhs = 0.0_rp

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = aux*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gpress(1),Finv,C)
          G = (G-C)

          !Testing against adjoint
          rhs = dot_product(G,usgs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) - tau*rhs

      end do

   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_dynsgs_nonli(e,nd,tn,dvol,mu,J,tau,gdev,Jder,divFmat,divFinv,dtinv2,densi,vesgs,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !           -(J/2mu)F^-1_Ki*d(xi_lm*S_ml)/dX_K),
      !                      Usgs_i^n)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: divFmat(nd),divFinv(nd)
      real(rp),             intent(in)    :: mu,dvol,J,tau,Jder(nd)
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd),G(nd),usgs(nd)
      real(rp)    :: B(tn,nd),DS(nd),BS(tn)
      real(rp)    :: AU(nd),rhs(tn),aux1

      usgs = densi*dtinv2*dvol*vesgs
      aux1= 1.0_rp/(2.0_rp*mu)

      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          rhs = 0.0_rp

          call sld_calculate_NdF(e,nd,tn,inode,divFmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_NS_dF(e,nd,tn,inode,J,gdev,Jder,divFinv,B)
          Bi = Bi - aux1*B

          !Testing against adjoint
          rhs = matmul(Bi,usgs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) - tau*rhs(:)

      end do

   end subroutine 
  
   module subroutine sup_rhs_usgs_P_NH_dynsgs_nonli(e,nd,tn,dvol,mu,lam,J,tau,gpress,Jder,divFmat,divFinv,dtinv2,densi,vesgs,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,
      !                                            Usgs_i^n)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: Jder(nd),gpress(1)
      real(rp),   intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux,usgs(nd)
      real(rp)                  :: AU(nd),rhs

      !Constant part of adjoint
      aux   = (2.0_rp*mu)/(lam*nd)

      usgs = densi*dtinv2*dvol*vesgs

      do inode=1,e%pnode

          G   = 0.0_rp
          rhs = 0.0_rp

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_qdF(e,nd,tn,inode,divFmat,C)
          G = aux*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_qdFinv(e,nd,tn,inode,J,lam,gpress(1),Jder,divFinv,C)
          G = G-C

          !Testing against adjoint
          rhs = dot_product(G,usgs)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) - tau*rhs

      end do

   end subroutine  

   !----------------------ADJOINT TESTING OF PROJECTION ----------------------

   module subroutine sup_rhs_ssgs_U_NH_repro(e,nd,tn,dvol,tau,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    tau_s*(symgrad(v_h)_ij,repro_S))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau,gpres(tn)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bt(nd,tn),aux
      real(rp)    :: rhs(nd)

      aux = dvol*tau
      do inode=1,e%pnode

          rhs  = 0.0_rp

          !grad(v_h)_ij
          call sld_calculateBt(e,nd,tn,inode,Bt)

          !Testing against adjoint
          rhs = matmul(Bt,gpres)

          !Assembly: RHS - tau*L*(Vh)*gpres)
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)*aux

      end do

   end subroutine sup_rhs_ssgs_U_NH_repro

   module subroutine sup_rhs_psgs_U_NH_repro(e,nd,tn,dvol,tau,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                  tau_p*(div(v_h),gpres_P)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol
      real(rp),             intent(in)    :: tau,gpres
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
      real(rp)                            :: N(nd),aux,rhs(nd)
      integer(ip)                         :: inode

      aux = tau*dvol

      do inode=1,e%pnode

          rhs = 0.0_rp

          !div(v_h)
          call sld_calculateB_inditial(e,nd,inode,N)

          !Testing against adjoint
          rhs = N*gpres*aux

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_psgs_U_NH_repro

   module subroutine sup_rhs_usgs_S_NH_repro(e,nd,tn,dvol,mu,J,tau,grsig,Fmat,Finv,gdev,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u*((F_jK*d(xi_ij)/dX_K-1/nd*F_iK*d(xi_ll)/dX_K)
      !           -(J/2mu)F^-1_Ki*d(xi_lm*S_ml)/dX_K),gpres_U)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau,gpres(nd),gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd),grsig(tn,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd)
      real(rp)    :: B(tn,nd)
      real(rp)    :: aux,rhs(tn)

      aux = dvol*tau
      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          rhs = 0.0_rp

          call sld_calculate_FdN(e,nd,tn,inode,Fmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_dNS_F(e,nd,tn,inode,Finv,grsig,gdev,B)
          Bi = Bi - (J/(2.0_rp*mu))*B

          !Testing against adjoint
          rhs = matmul(Bi,gpres)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + aux*rhs(:)

      end do

   end subroutine sup_rhs_usgs_S_NH_repro

   module subroutine sup_rhs_ssgs_S_NH_repro(e,nd,tn,dvol,mu,J,tau,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                tau_s*(J/2mu*xi_ij,gpres_S)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau,gpres(tn)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: ii(tn,tn),i(tn)
      real(rp)    :: N(tn),C(tn)
      real(rp)    :: rhs(tn)

      call getDiagTensor2Txy(nd,tn,ii)

      !Pre multiplying constant part of adjoint
      N  = 0.5_rp*(J/mu)*matmul(ii,gpres)*dvol

      do inode=1,e%pnode

          rhs  = 0.0_rp

          ! xi_ij
          !Testing against adjoint
          rhs = e%shape(inode,e%igaus)*N

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_ssgs_S_NH_repro

   module subroutine sup_rhs_usgs_P_NH_repro(e,nd,tn,dvol,mu,lam,J,tau,gradp,gppress,Fmat,Finv,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,gpres_U)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: gradp(nd),gppress(nd),gpres(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux,rhs

      !Constant part of adjoint
      aux   = (2.0_rp*mu)/(lam*nd)

      do inode=1,e%pnode

          G   = 0.0_rp
          rhs = 0.0_rp

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_dqF(e,nd,inode,Fmat,C)
          G = aux*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_dqFinv(e,nd,inode,J,lam,gradp,gppress(1),Finv,C)
          G = (G-C)*dvol

          !Testing against adjoint
          rhs = dot_product(G,gpres)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) + tau*rhs

      end do

   end subroutine sup_rhs_usgs_P_NH_repro

   module subroutine sup_rhs_psgs_P_NH_repro(e,nd,tn,dvol,mu,lam,J,b,tau,gpres,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                     -tau_p((J/lam)q_h,gpres_P)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,lam,J,dvol,b(nd,nd)
      real(rp),             intent(in)    :: tau,gpres
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                         :: inode
      real(rp)                            :: AU,b_ii,rhs

      AU = tau*(J/lam)*dvol*gpres

      do inode=1,e%pnode

          rhs = 0.0_rp

          !Testing against adjoint
          rhs = e%shape(inode,e%igaus)*AU

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) + rhs

      end do

   end subroutine sup_rhs_psgs_P_NH_repro

    !-----------------------------Higher order derivative RHS------------

   module subroutine sup_rhs_usgs_S_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! -tau_u*(2(xi_ij*d(F_iK)/dx_K - 1/nd*xi_ll*d(F_lK)/dx_K)
      !       -(J/mu)xi_ij*S_ij*d(F^-1_lK)/dx_K),-grad(p)_i)
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: Jder(nd),gdev(tn),mu,tau,J
      real(rp),   intent(in)    :: AU(nd)
      real(rp),   intent(inout) :: elrhs(tn,e%pnode)
      integer(ip) :: inode
      real(rp)    :: Bi(tn,nd),G(nd),aux1,aux(nd)
      real(rp)    :: B(tn,nd),DS(nd),BS(tn)
      real(rp)    :: rhs(tn)

      !A(U)
      aux = tau*AU*dvol
      aux1= 1.0_rp/(2.0_rp*mu)

      do inode=1,e%pnode

          B   = 0.0_rp
          Bi  = 0.0_rp
          rhs = 0.0_rp

          call sld_calculate_NdF(e,nd,tn,inode,divFmat,B)
          Bi = B

          !    F^-1_Ki*d(xi_lm*S_lm)/dX_K
          call sld_calculate_trace_NS_dF(e,nd,tn,inode,J,gdev,Jder,divFinv,B)
          Bi = Bi - aux1*B

          !Testing against adjoint
          rhs = matmul(Bi,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(:,inode) = elrhs(:,inode) + rhs(:)

      end do

   end subroutine sup_rhs_usgs_S_NH_nonli

   module subroutine sup_rhs_usgs_P_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,tau,AU,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      ! tau_u((2mu/lam*nd)*dqF-((Jp/lam)-1)*dqF^-1)_i,
      !              -div(S_ij)_i - grad(p)_i - F_ext,i)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: Jder(nd),gpress(1),mu,J,lam,tau
      real(rp),   intent(in)    :: AU(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
      integer(ip)               :: inode
      real(rp)                  :: C(nd),G(nd)
      real(rp)                  :: aux1,aux(nd)
      real(rp)                  :: rhs

      !Constant part of adjoint
      aux1   = (2.0_rp*mu)/(lam*nd)

      !A(U)
      aux= tau*AU*dvol

      do inode=1,e%pnode

          G   = 0.0_rp
          rhs = 0.0_rp

          !(2mu/lam*nd)*F_iK*dq/dX_K
          call sld_calculate_qdF(e,nd,tn,inode,divFmat,C)
          G = aux1*C

          !-d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
          call sld_calculate_qdFinv(e,nd,tn,inode,J,lam,gpress(1),Jder,divFinv,C)
          G = G-C

          !Testing against adjoint
          rhs = dot_product(G,aux)

          !Assembly: RHS + tau*L*(Vh)(-f + A(Uh))
          elrhs(1,inode) = elrhs(1,inode) + rhs

      end do

   end subroutine sup_rhs_usgs_P_NH_nonli

end submodule SUPCauchyElement_NH_SGS_RHS

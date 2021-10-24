submodule(Mod_CauchyElement) SUPCauchyElement_NH_Residual
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

    !-------------ADJOINT RHS------------

   module subroutine sup_residualU_dyn(e,nd,accel,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                         -(du/dt2)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: accel(nd)
      real(rp),             intent(inout) :: res(nd)

      res = res + accel

   end subroutine sup_residualU_dyn


   module subroutine sup_residualU(e,nd,gradp,divstr,fext,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !               -(-div(S_ij),i - grad(P),i) + F_ext,i
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: gradp(nd),divstr(nd),fext(nd)
      real(rp),             intent(inout) :: res(nd)

      res = res + fext + gradp + divstr

   end subroutine sup_residualU

   module subroutine sup_residualS(e,nd,tn,mu,J,b,gdev,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                -[(J/2mu)*gdev - 1/2(b - (b_ll/nd)*i)_ij]
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: b(nd,nd),gdev(tn),mu,J
      real(rp),             intent(inout) :: res(tn)
      real(rp)    :: i(nd,nd),b_ii
      real(rp)    :: C(nd,nd),Cvoigt(tn)
      real(rp)    :: AU(tn)

      call get2ndIITensor(nd,i)
      call trace(nd,b,b_ii)

      !Stress tensor term: (b - (b_ll/nd)*i)_ij
      C = (b - (b_ii/nd)*i)

      call getVoigtStress(tn,nd,C,Cvoigt)

      !A(U): (J/2mu)*gdev-0.5(b-(b_ll/nd)*i)_ij
      AU  = 0.5_rp*((J/mu)*gdev-Cvoigt)

      res = -AU

   end subroutine sup_residualS

   module subroutine sup_residualP(e,nd,mu,lam,J,b,gpress,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !            -[ (J/lam)*press - ln(J) - (mu/lam)*(b_ll/nd-1) ]
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd),gpress(1)
      real(rp),             intent(inout) :: res
      real(rp)                            :: b_ii,AU

      call trace(nd,b,b_ii)

      !A(U):
      AU =(J/lam)*gpress(1)-(mu/lam)*((b_ii/nd)-1.0_rp)
      AU = AU - log(J)

      res = -AU

   end subroutine sup_residualP

   !--------------------Orthogonal Residual---------------

   module subroutine sup_residualS_OSS(e,nd,tn,b,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                 -1/2(b - (b_ll/nd)*i)_ij
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: b(nd,nd)
      real(rp),             intent(inout) :: res(tn)
      real(rp)    :: i(nd,nd),b_ii
      real(rp)    :: C(nd,nd),Cvoigt(tn)
      real(rp)    :: AU(tn)

      call get2ndIITensor(nd,i)
      call trace(nd,b,b_ii)

      !Stress tensor term: (b - (b_ll/nd)*i)_ij
      C = (b - (b_ii/nd)*i)

      call getVoigtStress(tn,nd,C,Cvoigt)

      !A(U): -0.5(b-(b_ll/nd)*i)_ij
      AU  = 0.5_rp*Cvoigt

      res = -AU

   end subroutine sup_residualS_OSS


   module subroutine sup_residualP_OSS(e,nd,mu,lam,J,b,res)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !                  - ln(J) - (mu/lam)*(b_ll/nd-1)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd)
      real(rp),             intent(inout) :: res
      real(rp)                            :: b_ii,AU

      call trace(nd,b,b_ii)

      !A(U):
      AU =-(mu/lam)*((b_ii/nd)-1.0_rp)
      AU = AU -log(J)

      res = -AU

   end subroutine sup_residualP_OSS

end submodule SUPCauchyElement_NH_Residual

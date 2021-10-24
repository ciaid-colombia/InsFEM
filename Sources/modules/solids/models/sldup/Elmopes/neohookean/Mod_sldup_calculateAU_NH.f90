module Mod_sldup_calculateAU
  use typre
  implicit none   
   
  contains

  subroutine getAU_U(nd,gradp,divstr,AU)
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: gradp(nd),divstr(nd)
      real(rp),             intent(inout) :: AU(nd)

      !A(U)
      AU = - gradp - divstr

  end subroutine

  subroutine getAU_U_dyn(nd,densi,accel,gradp,divstr,AU)
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: gradp(nd),divstr(nd),accel(nd),densi
      real(rp),             intent(inout) :: AU(nd)

      !A(U)
      AU = densi*accel - gradp - divstr

  end subroutine

  subroutine getAU_P(nd,mu,lam,J,b,gpress,AU)
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd)
      real(rp),             intent(in)    :: gpress(1)
      real(rp),             intent(inout) :: AU
      real(rp)                            :: b_ii

      call trace(nd,b,b_ii)
      !A(U):
      AU = (J/lam)*gpress(1)-(mu/lam)*((b_ii/nd)-1.0_rp)
      AU = AU - log(J)

  end subroutine

  subroutine getAU_P_OSGS(nd,mu,lam,J,b,gpress,AU)
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd)
      real(rp),             intent(in)    :: gpress(1)
      real(rp),             intent(inout) :: AU
      real(rp)                            :: b_ii

      call trace(nd,b,b_ii)
      !A(U):
      AU = -(mu/lam)*((b_ii/nd)-1.0_rp)
      AU = AU - log(J)

  end subroutine

  end module

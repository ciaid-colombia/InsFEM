module Mod_sldsup_calculateAU
  use typre
  implicit none   
   
  contains

  subroutine getAU_S(nd,tn,mu,J,b,gdev,AU)
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,J,gdev(tn),b(nd,nd)
      real(rp),             intent(inout) :: AU(tn)
      real(rp)    :: i(tn),b_ii
      real(rp)    :: C(tn),bvoigt(tn)

      call get2ndIITensor_voigt(nd,tn,i)
      call trace(nd,b,b_ii)
      call getVoigtStress(tn,nd,b,bvoigt)

      !Stress tensor term: (b - (b_ll/nd)*i)_ij
      C = (bvoigt - (b_ii/nd)*i)

      !A(U): (J/2mu)*S_ij -0.5(b-(b_ll/nd)*i)_ij
      AU  = 0.5_rp*((J/mu)*gdev-C)

  end subroutine

  subroutine getAU_S_OSGS(nd,tn,mu,b,gdev,AU)
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,gdev(tn),b(nd,nd)
      real(rp),             intent(inout) :: AU(tn)
      real(rp)    :: i(tn),b_ii
      real(rp)    :: C(tn),bvoigt(tn)

      call get2ndIITensor_voigt(nd,tn,i)
      call trace(nd,b,b_ii)
      call getVoigtStress(tn,nd,b,bvoigt)

      !Stress tensor term: (b - (b_ll/nd)*i)_ij
      C = (bvoigt - (b_ii/nd)*i)

      !A(U): (J/2mu)*S_ij -0.5(b-(b_ll/nd)*i)_ij
      AU  = -0.5_rp*C

  end subroutine

  end module

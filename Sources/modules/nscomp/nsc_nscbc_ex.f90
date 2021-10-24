module Mod_nsc_nscbc_ex
contains
subroutine nsc_nscbc_ex(e,accph,accvh,gpden,gpmom,gpene,grden,grmom,grene,coef,tract)
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_nscbc
! NAME 
!    nsc_nscbc
! DESCRIPTION
!    This routine computes a non-reflecting boundary condition 
!    for subsonic flows, based on the NSCBC approach.
!    Higher order derivatives (stress and heat flux) are neglected.
!    It recalculates the essential boundary conditions.
!    Transverse, and pressure terms are neglected (commented).
! USED BY
!    nsc_Endbouope
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp) :: accph,accvh
   real(rp) :: gpden,gpmom(e%ndime),gpene
   real(rp) :: grden(e%ndime),grmom(e%ndime,e%ndime),grene(e%ndime)
   real(rp) :: coef(3)
   real(rp),    intent(inout) :: tract(e%ndime+2)

   integer(ip)                :: idime
   integer(ip)                :: inodb,jnode,inode,jnodb

   real(rp) :: cgamma,aux_k,gptem 
   real(rp) :: lomom(e%ndime),auxlogrm(e%ndime,e%ndime)
   real(rp) :: logrd(e%ndime),logrm(e%ndime,e%ndime),logre(e%ndime)
   real(rp) :: logrv(e%ndime,e%ndime),logrp(e%ndime)
   real(rp) :: gppre,lospd,dens,veloc(e%ndime),cartvel(e%ndime),tempe 
   real(rp) :: L1,L2,L5,d(4) 

   cgamma = accph/accvh
   aux_k = dot_product(gpmom,gpmom)/(gpden*gpden)

   !baloc is (t1,t2,n)
   !therefore (ndime) is the normal component.
   !(1) is a transversal for 2d and (2) transversal for 3d
   !local density derivative
   logrd = matmul(grden,e%baloc)
  
   !local momentum and derivatives 
   lomom = matmul(gpmom,e%baloc)
   ! d_l m_l = T_lc (d_c m_c) T_cl
   auxlogrm = matmul(transpose(grmom),e%baloc)
   logrm = matmul(transpose(auxlogrm),e%baloc)

   !local energy derivative
   logre = matmul(grene,e%baloc)

   !Calulate local pressure and gradient
   gppre = (cgamma-1.0_rp)*(gpene-gpden*aux_k/2.0_rp)

   !d_j(p)=(gamma-1)*(dj(etot) - m/rho dj(m) + m·m/rho·rho dj(rho))
   logrp = 0.0_rp
   logrp = logrp + (cgamma-1.0_rp) * (logre + aux_k*logrd)
   logrp = logrp - (cgamma-1.0_rp) * matmul(lomom,logrm)/gpden

   !Calulate local temperature
   gptem = gppre/(gpden*(accph-accvh))
   
   !Calculate speed of sound
   lospd = sqrt(cgamma*gppre/gpden)

   !Calculate local velocity derivatives 
   !dj(ui) = dj(mi)/rho - ui/rho (dj rho)
   do idime=1,e%ndime
      logrv(idime,1:e%ndime) = logrm(idime,1:e%ndime)/gpden -&
                               lomom(idime)*logrd(1:e%ndime)/(gpden*gpden) 
   end do

   !Calculate amplitudes of characteristic waves

   if(e%ndime == 3) then 
        veloc(2) = veloc(2) + lomom(3)*logrv(2,3)/gpden
   end if

   !L1 Linear Relaxation Method
   L1 = coef(1)*abs(1.0_rp-aux_k/(lospd*lospd))*(gppre-coef(3))/(sqrt(2.0_rp)*gpden*coef(2))
!   L1 = (-lospd+lomom(e%ndime)/gpden)*(-gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))
   !L2 and L5
   L2 = (lomom(e%ndime)/gpden)*(lospd*lospd*logrd(e%ndime)-logrp(e%ndime))
   L5 = (lospd+lomom(e%ndime)/gpden)*(gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))

   !Calculate vector of derivatives  
   !d1, d2, ...,d5
   d = 0.0_rp
   d(2) = (L1+L5)/2.0_rp 
   d(1) = (L2+d(2))/(lospd*lospd) 
   d(3) = (L5-L1)/(2.0_rp*gpden*lospd)
   d(4) = (lomom(e%ndime)/gpden)*logrv(1,e%ndime)

   !Density equation
   dens = 0.0_rp
   dens = dens + d(1) 

   !Temperature equation
   tempe = 0.0_rp
   tempe = tempe - gptem*L2/(gpden*lospd*lospd)
   tempe = tempe + gptem*(cgamma-1.0_rp)*d(2)/(gpden*lospd*lospd) 

   veloc = 0.0_rp 
   veloc(e%ndime) = veloc(e%ndime) + d(3) 
   veloc(1) = veloc(1) + d(4) 

   cartvel = matmul(e%baloc,veloc)

   ! Compute prescribed traction
   tract(1)           = dens
   tract(2:e%ndime+1) = cartvel(1:e%ndime)
   tract(e%ndime+2)   = tempe

end subroutine nsc_nscbc_ex
end module

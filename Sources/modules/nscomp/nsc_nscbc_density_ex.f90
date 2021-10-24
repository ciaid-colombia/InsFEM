subroutine nsc_nscbc_density_ex(e,dtinv,accph,accvh,gpden,gpmom,gpene,gpden0,gpmom0,gpene0,grden,grmom,grene,incre)
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
   real(rp) :: dtinv,accph,accvh
   real(rp) :: gpden,gpmom(e%ndime),gpene,gpden0,gpmom0(e%ndime),gpene0
   real(rp) :: grden(e%ndime),grmom(e%ndime,e%ndime),grene(e%ndime)
   real(rp),    intent(inout) :: incre(e%ndime+2)

   integer(ip)                :: idime
   integer(ip)                :: inodb,jnode,inode,jnodb

   real(rp) :: cgamma,aux_k,gptem,aux_k0 
   real(rp) :: lomom(e%ndime),auxlogrm(e%ndime,e%ndime),lomom0(e%ndime)
   real(rp) :: logrd(e%ndime),logrm(e%ndime,e%ndime),logre(e%ndime)
   real(rp) :: logrv(e%ndime,e%ndime),logrp(e%ndime)
   real(rp) :: gppre,lospd,dens,veloc(e%ndime),cartvel(e%ndime),gppre0,gptem0 
   real(rp) :: L1,L2,L5,d(2) 

   cgamma = accph/accvh
   aux_k = dot_product(gpmom,gpmom)/(gpden*gpden)
   aux_k0 = dot_product(gpmom0,gpmom0)/(gpden0*gpden0)

   !baloc is (t1,t2,n)
   !therefore (ndime) is the normal component.
   !(1) is a transversal for 2d and (2) transversal for 3d
   !local density derivative
   logrd = matmul(grden,e%baloc)
  
   !local momentum and derivatives 
   lomom = matmul(gpmom,e%baloc)
   lomom0 = matmul(gpmom0,e%baloc)
   ! d_l m_l = T_lc (d_c m_c) T_cl
   auxlogrm = matmul(transpose(grmom),e%baloc)
   logrm = matmul(transpose(auxlogrm),e%baloc)

   !local energy derivative
   logre = matmul(grene,e%baloc)

   !Calulate local pressure and gradient
   gppre = (cgamma-1.0_rp)*(gpene-gpden*aux_k/2.0_rp)
   gppre0 = (cgamma-1.0_rp)*(gpene0-gpden0*aux_k0/2.0_rp)

   !d_j(p)=(gamma-1)*(dj(etot) - m/rho dj(m) + m·m/rho·rho dj(rho))
   logrp = 0.0_rp
   logrp = logrp + (cgamma-1.0_rp) * (logre + aux_k*logrd)
   logrp = logrp - (cgamma-1.0_rp) * matmul(lomom,logrm)/gpden

   !Calulate local temperature
   gptem = gppre/(gpden*(accph-accvh))
   gptem0 = gppre0/(gpden0*(accph-accvh))
   
   !Calculate speed of sound
   lospd = sqrt(cgamma*gppre/gpden)

   !Calculate local velocity derivatives 
   !dj(ui) = dj(mi)/rho - ui/rho (dj rho)
   do idime=1,e%ndime
      logrv(idime,1:e%ndime) = logrm(idime,1:e%ndime)/gpden -&
                               lomom(idime)*logrd(1:e%ndime)/(gpden*gpden) 
   end do


   !Calculate amplitude of outgoing characteristic wave (L1) in terms of interior points
   !Subsonic outgoing wave lambda_1 = (+) soundspeed (-) velocity
   !           <---|      (+n) boundary             
   !           <---|----  (+c)             
   !               |             
   !           ----|--->  (-u1)             
   L1 = (lospd-lomom(e%ndime)/gpden)*(-gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))
   !L5 in terms of L1 and velocity derivative
   L5 = L1 - 2*lospd*(lomom(e%ndime)-gpden*lomom0(e%ndime)/gpden0)*dtinv

   !Calculate vector of derivatives  
   d = 0.0_rp
   d(2) = (L1+L5)/2.0_rp 
   L2 = (cgamma-1.0_rp)*d(2)+ gpden*lospd*lospd*(1.0_rp-gptem0/gptem)*dtinv
   d(1) = (L2+d(2))/(lospd*lospd) 


   !Density equation
   dens = 0.0_rp
   dens = dens + d(1) 

   ! Compute prescribed traction
   incre(1)           = dens

end subroutine nsc_nscbc_density_ex

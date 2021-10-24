module Mod_NSCompressiblePrimitiveNBC
   use typre
   use Mod_Element
   
   implicit none

contains
   
subroutine nsc_pr_WeakElementalAssembly(e,coeff,norder,Abmat,Sbmat,Rhmat,wmatr,wrhsi)
!-----------------------------------------------------------------------
! NAME 
!    nsc_pr_WeakElementalAssembly
! DESCRIPTION
!    This routine computes a penalization for the NS momentum equation at
!    a given integration point of a boundary received by argument
!    due to the weak imposition of conditions. 
! USED BY
!    nsc_pr_bouope
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp) :: coeff,Rhmat(e%ndime+2),norder(e%mnodb)
   real(rp),    intent(in) :: Abmat(e%ndime+2,e%ndime+2)
   real(rp),    intent(in) :: Sbmat(e%ndime+2,e%ndime+2)
   real(rp),    intent(inout) :: wmatr(e%ndime+2,e%mnode,e%ndime+2,e%mnode)
   real(rp),    intent(inout) :: wrhsi(e%ndime+2,e%mnode)

   integer(ip)                :: idofn,jdofn
   integer(ip)                :: inodb,jnode,inode,jnodb

   do inodb = 1,e%pnodb
      inode = e%lboel(inodb)
      do idofn = 1,e%ndime+2
         do jdofn = 1,e%ndime+2
            do jnodb = 1,e%pnodb
               jnode = e%lboel(jnodb)
                  wmatr(idofn,inode,jdofn,jnode) = wmatr(idofn,inode,jdofn,jnode) + &
                        coeff*e%shapb(inodb,e%igaub)*Abmat(idofn,jdofn)*norder(jnodb) +&
                        coeff*e%shapb(inodb,e%igaub)*Sbmat(idofn,jdofn)*e%shapb(jnodb,e%igaub)
            enddo
         enddo
         wrhsi(idofn,inode) = wrhsi(idofn,inode) +&
               coeff*Rhmat(idofn)*e%shapb(inodb,e%igaub)
      enddo
   enddo
end subroutine nsc_pr_WeakElementalAssembly

subroutine nsc_pr_nscbc_ex(e,accph,accvh,gppre,gpvel,gptem,grpre,grvel,grtem,relpre,reltem,coef,tract)
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
   real(rp) :: gppre,gpvel(e%ndime),gptem
   real(rp) :: grpre(e%ndime),grvel(e%ndime,e%ndime),grtem(e%ndime)
   real(rp) :: relpre,reltem,coef(3)
   real(rp),    intent(inout) :: tract(e%ndime+2)

   real(rp) :: cgamma,lomach
   real(rp) :: lovel(e%ndime),auxlogrv(e%ndime,e%ndime)
   real(rp) :: logrp(e%ndime),logrv(e%ndime,e%ndime),logrt(e%ndime)
   real(rp) :: gpden,logrd(e%ndime)
   real(rp) :: lospd,tempe,veloc(e%ndime),cartvel(e%ndime),press 
   real(rp) :: L1,L2,L5,d(4) 

   cgamma = accph/accvh

   !baloc is (t1,t2,n)
   !therefore (ndime) is the normal component.
   !(1) is a transversal for 2d and (2) transversal for 3d

   !local pressure derivative
   logrp = matmul(grpre,e%baloc)
  
   !local velocity and derivatives 
   lovel = matmul(gpvel,e%baloc)
   ! d_l v_l = T_lc (d_c v_c) T_cl
   auxlogrv = matmul(transpose(grvel),e%baloc)
   logrv = matmul(transpose(auxlogrv),e%baloc)

   !local temperature derivative
   logrt = matmul(grtem,e%baloc)

   !Calulate local density and gradient
   gpden = (gppre+relpre)/((accph-accvh)*(gptem+reltem))

   !d_j(rho)=1/(R)*(dj(pre)/gptem - pre/gptem*gptem (dj(tem)))
   logrd = (logrp - (gppre+relpre)*logrt/(gptem+reltem))/((accph-accvh)*(gptem+reltem))

   !Calculate speed of sound
   lospd = sqrt(cgamma*(gppre+relpre)/gpden)
   !Local Mach number
   lomach = sqrt(dot_product(gpvel,gpvel))/lospd

   !Calculate amplitudes of characteristic waves

   if(e%ndime == 3) then 
      veloc(2) = lovel(3)*logrv(2,3)
   end if


   !L1 Linear Relaxation Method
   L1 = coef(1)*lospd*abs(1.0_rp-lomach*lomach)*(gppre-coef(3))/coef(2)
!   L1 = (-lospd+lovel(e%ndime))*(-gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))
 
   !L2 and L5
   L2 = lovel(e%ndime)*(lospd*lospd*logrd(e%ndime)-logrp(e%ndime))
   L5 = (lospd+lovel(e%ndime))*(gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))

   !Calculate vector of derivatives  
   !d1, d2, ...,d5
   d = 0.0_rp
   d(2) = (L1+L5)/2.0_rp 
   d(1) = (L2+d(2))/(lospd*lospd) 
   d(3) = (L5-L1)/(2.0_rp*gpden*lospd)
   d(4) = lovel(e%ndime)*logrv(1,e%ndime)

   !Pressure equation
   press = 0.0_rp
   press = press + d(2) 

   !Temperature equation
   tempe = 0.0_rp
   tempe = tempe - (gptem+reltem)*L2/(gpden*lospd*lospd)
   tempe = tempe + (gptem+reltem)*(cgamma-1.0_rp)*d(2)/(gpden*lospd*lospd) 

   veloc = 0.0_rp 
   veloc(e%ndime) = veloc(e%ndime) + d(3) 
   veloc(1) = veloc(1) + d(4) 
  
   cartvel = matmul(e%baloc,veloc)

   ! Compute prescribed traction
   tract(1)           = press
   tract(2:e%ndime+1) = cartvel(1:e%ndime)
   tract(e%ndime+2)   = tempe

end subroutine nsc_pr_nscbc_ex

subroutine nsc_pr_nscbc_press_ex(e,dtinv,accph,accvh,gppre,gpvel,gptem,gpvel0,grpre,grvel,relpre,reltem,incre)
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_nscbc
! NAME 
!    nsc_nscbc
! DESCRIPTION
!    This routine computes a non-reflecting inlet boundary condition 
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
   real(rp) :: gppre,gpvel(e%ndime),gptem,gpvel0(e%ndime)
   real(rp) :: grpre(e%ndime),grvel(e%ndime,e%ndime)
   real(rp) :: relpre,reltem
   real(rp),    intent(inout) :: incre(e%ndime+2)

   real(rp) :: cgamma
   real(rp) :: lovel(e%ndime),auxlogrv(e%ndime,e%ndime),lovel0(e%ndime)
   real(rp) :: logrp(e%ndime),logrv(e%ndime,e%ndime)
   real(rp) :: gpden
   real(rp) :: lospd,press 
   real(rp) :: L1,L5

   cgamma = accph/accvh

   !baloc is (t1,t2,n)
   !therefore (ndime) is the normal component.
   !(1) is a transversal for 2d and (2) transversal for 3d
   !local pressure derivative
   logrp = matmul(grpre,e%baloc)
  
   !local velocity and derivatives 
   lovel = matmul(gpvel,e%baloc)
   lovel0 = matmul(gpvel0,e%baloc)
   ! d_l v_l = T_lc (d_c v_c) T_cl
   auxlogrv = matmul(transpose(grvel),e%baloc)
   logrv = matmul(transpose(auxlogrv),e%baloc)

   !Calulate local density and gradient
   gpden = (gppre+relpre)/((accph-accvh)*(gptem+reltem))

   !Calculate speed of sound
   lospd = sqrt(cgamma*(gppre+relpre)/gpden)

   !Calculate amplitude of outgoing characteristic wave (L1) in terms of interior points
   !Subsonic outgoing wave lambda_1 = (+) soundspeed (-) velocity
   !           <---|      (+n) boundary             
   !           <---|----  (+c)             
   !               |             
   !           ----|--->  (-u1)             
   L1 = (lospd-lovel(e%ndime))*(-gpden*lospd*logrv(e%ndime,e%ndime)+logrp(e%ndime))
   !L5 in terms of L1 and velocity derivative
   L5 = L1 - 2*lospd*gpden*(lovel(e%ndime)-lovel0(e%ndime))*dtinv

   !Pressure equation
   press = 0.0_rp
   press = press + (L1+L5)/2.0_rp 

   incre(1) = press

end subroutine nsc_pr_nscbc_press_ex

end module


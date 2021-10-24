module Mod_LowMachElement_Residual
   use typre
   use Mod_Element
   
contains
   
   !Pressure subscale residual
   subroutine lmn_elmrfe_pre(e,Integrator,dtinv,actex,veadv,degau,grate,divve,resim)
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)     :: e
      type(TimeIntegratorDt1) :: Integrator

      real(rp) :: veadv(e%ndime),degau(*),grate(e%ndime),divve,resim,dtinv,actex,gprhs(1)
      
      call Integrator%GetRHS(1,degau(2),gprhs)
      !Temporal derivative LHS
!      resim = resim - (degau(1) - degau(2))*dtinv
      resim = resim - (degau(1)*Integrator%Coefficients(1)-gprhs(1))*dtinv
      
      !Gradient temperature
      resim = resim + degau(1)*actex*dot_product(veadv,grate)

      !Divergence velocity
      resim = resim - degau(1)*divve

   end subroutine lmn_elmrfe_pre

   !Velocity subscale residual
   subroutine lmn_elmrfe_mom(e,dtinv,acden,vegau,grapr,elext,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: dtinv,acden,vegau(e%ndime),grapr(e%ndime),resim(e%ndime)
      real(rp) :: elext(e%ndime)
      
      integer(ip) :: idime,jdime
      
      !Temporal derivative LHS
      resim(1:e%ndime) = resim(1:e%ndime) - acden*dtinv*vegau(1:e%ndime)
      
      !RHS (Temporal derivative + External forces)
      resim(1:e%ndime) = resim(1:e%ndime) + elext(1:e%ndime)
      
      !Pressure Gradient
      resim(1:e%ndime) = resim(1:e%ndime) - grapr(1:e%ndime)

   end subroutine lmn_elmrfe_mom

   subroutine lmn_elmrfe_mom_trm(e,grapr,elext,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: grapr(e%ndime),resim(e%ndime)
      real(rp) :: elext(e%ndime)
      
      !RHS (Temporal derivative + External forces)
      resim(1:e%ndime) = resim(1:e%ndime) + elext(1:e%ndime)
      
      !Pressure Gradient
      resim(1:e%ndime) = resim(1:e%ndime) - grapr(1:e%ndime)

   end subroutine lmn_elmrfe_mom_trm
   
   subroutine lmn_elmrfe_mom_adv(e,acden,veadv,grave,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: acden,veadv(e%ndime),grave(e%ndime,e%ndime),resim(e%ndime)
      real(rp) :: elext(e%ndime)
      
      integer(ip) :: idime,jdime

      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(1:e%ndime) = resim(1:e%ndime) - acden*veadv(jdime)*grave(1:e%ndime,jdime)
      end do

   end subroutine lmn_elmrfe_mom_adv

   subroutine lmn_elmrfe_mom_nonlinear(e,acvis,elvel,gprep)
      implicit none
      class(FiniteElement) :: e
      real(rp)             :: acvis,elvel(e%ndime,e%pnode),gprep(e%ndime)
      integer(ip)          :: idime,inode,pos(e%ndime,e%ndime)
      real(rp) :: aux1(e%pnode)

      select case (e%ndime)
         case (2)
            pos(1,:)=(/1,3/)
            pos(2,:)=(/3,2/)
         case (3)
            pos(1,:)=(/1,4,5/)
            pos(2,:)=(/4,2,6/)
            pos(3,:)=(/5,6,3/)
         case default
      end select     
     

      do inode = 1,e%pnode         ! Contribution from the laplacian term
         aux1(inode) = sum(e%hessi(1:e%ndime,inode)) 
      end do   
      do idime = 1,e%ndime
         do inode = 1,e%pnode
            gprep(idime) = gprep(idime) + acvis*(aux1(inode)*elvel(idime,inode) + &
               dot_product(e%hessi(pos(idime,1:e%ndime),inode),elvel(1:e%ndime,inode)))
         end do   
      end do

   end subroutine

   !Temperature subscale residual
   subroutine lmn_elmrfe_ene(e,Integrator,LHSdtinv,dtinv,acden,actex,tegau,tegau_lin,pther,elext,resim)
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)     :: e
      type(TimeIntegratorDt1) :: Integrator
      real(rp) :: LHSdtinv,dtinv,acden,actex,tegau,tegau_lin,resim
      real(rp) :: elext,pther(*),gprhs(1)
      
      call Integrator%GetRHS(1,pther(3),gprhs)

      !Temporal derivative LHS
      resim = resim - acden*LHSdtinv*tegau_lin
      
      !RHS (Temporal derivative + External forces)
      resim = resim + elext
      
      !Thermodynamic pressure
      resim = resim + actex*tegau*dtinv*(pther(1)*Integrator%Coefficients(1)-gprhs(1)) !(dpth/dt)
!      resim = resim + actex*tegau*dtinv*(pther(1)-pther(3)) !(dpth/dt)

      !Contribution from the reaction term, which is in the finite element space
!      resim = resim + acrcp*tegau

   end subroutine lmn_elmrfe_ene

   subroutine lmn_elmrfe_ene_trm(e,Integrator,dtinv,actex,tegau,pther,elext,resim)
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)     :: e
      type(TimeIntegratorDt1) :: Integrator
      real(rp) :: dtinv,actex,tegau,resim
      real(rp) :: elext,pther(*),gprhs(1)
      
      call Integrator%GetRHS(1,pther(3),gprhs)

      !RHS (Temporal derivative + External forces)
      resim = resim + elext
      
      !Thermodynamic pressure
      resim = resim + actex*tegau*dtinv*(pther(1)*Integrator%Coefficients(1)-gprhs(1)) !(dpth/dt)
!      resim = resim + actex*tegau*dtinv*(pther(1)-pther(3)) !(dpth/dt)

   end subroutine lmn_elmrfe_ene_trm
   
   subroutine lmn_elmrfe_ene_adv(e,acden,veadv,grate,resim)
      implicit none
      class(FiniteElement)     :: e
      real(rp) :: acden,veadv(e%ndime),grate(e%ndime),resim
      integer(ip) :: idime,jdime
      
      !Contribution from the convective term
      resim = resim - acden*dot_product(veadv(1:e%ndime),grate(1:e%ndime))

   end subroutine lmn_elmrfe_ene_adv

   subroutine lmn_elmrfe_ene_nonlinear(e,actco,eltem,gprep)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: actco,eltem(e%pnode),gprep
      integer(ip)          :: inode
      real(rp) :: aux1(e%pnode)
     
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo   
         gprep = gprep + actco*dot_product(aux1,eltem(1:e%pnode))
   
   end subroutine

end module

module Mod_TemperatureElement
   use typre
   use Mod_Element
   implicit none

contains

   subroutine tem_ComputeTestf(e,acden,timom,AgradV,acrcp,testf)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: acden,AGradV(e%pnode), testf(e%pnode),timom,acrcp
      
      testf = timom*(acden*AGradV - acrcp*e%shape(1:e%pnode,e%igaus))
      
   end subroutine   
      
   subroutine tem_ComputeExternalForces(e,acsou,acsph,elext)
      use typre

      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: acden, acsou,acsph,elext
      
      elext = elext + acsou/acsph

   end subroutine
   
   subroutine tem_TimeIntegrationToElext(e,Integrator,acden,dtinv,gptem,elext)
      use typre
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)      :: Integrator
      real(rp)                   :: gptem(*),elext
      real(rp)                   :: acden,dtinv
      
      real(rp)                   :: gprhs(1)
      
      !Time integration
      call Integrator%GetRHS(1,gptem(2),gprhs)
      elext = elext + acden*dtinv*gprhs(1)
          
   end subroutine

   subroutine tem_elmrhu(e,dvolu,testf,elext,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: testf(e%pnode),elext
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%pnode)

      integer(ip)                :: inode
      real(rp)                   :: aux,tmp

      do inode=1,e%pnode
         aux = (e%shape(inode,e%igaus) + testf(inode))*dvolu
         elrhs(inode) = elext*aux + elrhs(inode)
      end do

   end subroutine tem_elmrhu
   
   subroutine tem_elmbuv(e,dvolu,denac,acrcp,dtinv,vtemp,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: vtemp(:)
      real(rp),    intent(in)    :: testf(:)
      real(rp),    intent(in)    :: dvolu,denac,dtinv,acrcp
      real(rp),    intent(inout) :: elmat(:,:)

      integer(ip)                :: jnode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         !(rho/dt+s)*Mass + rho*a·grad
         aux1 = (denac*dtinv+acrcp)*(e%shape(jnode,e%igaus)) + denac*vtemp(jnode)
         elmat(1:e%pnode,jnode) = (e%shape(1:e%pnode,e%igaus) + testf(1:e%pnode))*aux1*dvolu + elmat(1:e%pnode,jnode)
      end do

   end subroutine tem_elmbuv
   
   subroutine tem_elmrhu_rep(e,dvolu,testf,elext,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: testf(e%pnode),elext
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%pnode)

      integer(ip)                :: inode
      real(rp)                   :: aux,tmp

      do inode=1,e%pnode
         aux = (e%shape(inode,e%igaus) )*dvolu
         elrhs(inode) = elext*aux + elrhs(inode)
      end do

   end subroutine tem_elmrhu_rep
   
   subroutine tem_elmbuv_rep(e,dvolu,denac,acrcp,dtinv,vtemp,testf,elmat)
      !This subroutine does not include the time derivative and reaction in the
      !residual
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, a*grad u) + (v, u/dt) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: vtemp(:)
      real(rp),    intent(in)    :: testf(:)
      real(rp),    intent(in)    :: dvolu,denac,dtinv,acrcp
      real(rp),    intent(inout) :: elmat(:,:)

      integer(ip)                :: jnode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         !(rho/dt+s)*Mass + rho*a·grad
         aux1 = (denac*dtinv+acrcp)*(e%shape(jnode,e%igaus)) + denac*vtemp(jnode)
         elmat(1:e%pnode,jnode) = (e%shape(1:e%pnode,e%igaus))*aux1*dvolu + testf(1:e%pnode)*denac*vtemp(jnode)*dvolu + elmat(1:e%pnode,jnode)
      end do

   end subroutine tem_elmbuv_rep

   subroutine tem_ComputeTestfNonLinear(e,acvis,timom,kfl_stabm,testf)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: acvis, timom, testf(e%pnode)
      integer(ip)         :: kfl_stabm
      
      integer(ip) :: inode
      real(rp)    :: aux1
      
      !Stabilization terms : -tau visco *lapla(v)
      aux1 = real(kfl_stabm)*timom*acvis
      do inode = 1,e%pnode
         testf(inode) = testf(inode) + aux1*sum(e%hessi(1:e%ndime,inode))
      enddo
   end subroutine
   
   subroutine tem_elmbuv_lap(e,dvolu,visac,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    - tau1*(L*v, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: visac
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)

      integer(ip)                :: jnode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         aux1 = -visac*sum(e%hessi(1:e%ndime,jnode))
         elmat(1:e%pnode,jnode) = (testf(1:e%pnode))*aux1*dvolu + elmat(1:e%pnode,jnode)
      end do

   end subroutine tem_elmbuv_lap

   
   !-----------------------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   
   !This suborutine computes the residual at a Gauss point
   subroutine tem_elmrfe_oto(e,dtinv,acden,acsph,acrcp,acsou,veadv, &
         tegau,grate,elext,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: dtinv,acden,acsph,acrcp,acsou,veadv(e%ndime),tegau,grate(e%ndime),resim
      real(rp) :: elext
      
      integer(ip) :: idime,jdime
      
      resim = 0.0_rp
      
      !Temporal derivative LHS is in the finite element space
      !resim = resim + acden*dtinv*tegau
      
      !RHS (Temporal derivative + External forces) is in the finite element space
      !resim = resim - elext
      
      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim = resim + acden*veadv(jdime)*grate(jdime) 
      end do
      
      !Contribution from the reaction term is in the finite element space
      !resim = resim + acrcp*tegau
      
   end subroutine tem_elmrfe_oto
   
   !This suborutine computes the FULL residual at a Gauss point (For dissipation)
   subroutine tem_elmrfe_osa(e,dtinv,acden,acsph,acrcp,acsou,veadv, &
         tegau,grate,elext,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: dtinv,acden,acsph,acrcp,acsou,veadv(e%ndime),tegau,grate(e%ndime),resim
      real(rp) :: elext
      
      integer(ip) :: idime,jdime
      
      resim = 0.0_rp
      
      !Temporal derivative LHS, which is in the finite element space
      resim = resim + acden*dtinv*tegau
      
      !RHS (Temporal derivative + External forces), which is in the finite element space
      resim = resim - elext
      
      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim = resim + acden*veadv(jdime)*grate(jdime) 
      end do
      
      !Contribution from the reaction term, which is in the finite element space
      resim = resim + acrcp*tegau
      
   end subroutine tem_elmrfe_osa
   
   subroutine tem_elmrfe_oto_nonlinear(e,acvis,eltem,gprep)
      class(FiniteElement) :: e
      real(rp) :: acvis, eltem(e%pnode), gprep
      
      integer(ip) :: idime,inode
      real(rp)    :: aux1(e%pnode)
      
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo   
      gprep = gprep -  acvis*dot_product(aux1,eltem(1:e%pnode))
   end subroutine
   
   
   subroutine tem_elmrep(e,dvol,gprep,elrep)
      !Contribution to the residual rhs
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: dvol
      real(rp) :: gprep
      real(rp) :: elrep(*)
      
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elrep(inode) = elrep(inode) + e%shape(inode,e%igaus)*gprep*dvol
      enddo
      
   end subroutine  
   
   subroutine tem_elmrhs_oss(e,dvol,testf,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau*(L*v, proj(Res))
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none

    class(FiniteElement)        :: e
    real(rp),    intent(in)    :: gprep
    real(rp),    intent(in)    :: testf(e%pnode)
    real(rp),    intent(in)    :: dvol
    real(rp),    intent(inout) :: elrhs(e%pnode)

    integer(ip)                :: inode,idime
    real(rp)                   :: aux2,tmp2

    
    do inode=1,e%pnode
       elrhs(inode) = gprep*testf(inode)*dvol + elrhs(inode)
    end do

  end subroutine tem_elmrhs_oss
  
  
  !For dynamic subgrid scales
  subroutine tem_elmrhs_dss(e,dvol,acden,dtinv,testf,tesgs,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau*(L*v, proj(Res))
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none

    class(FiniteElement)        :: e
    real(rp),    intent(in)    :: dvol,acden,dtinv,tesgs
    real(rp),    intent(in)    :: testf(e%pnode)
    real(rp),    intent(inout) :: elrhs(e%pnode)

    integer(ip)                :: inode,idime
    real(rp)                   :: aux

    aux = acden*dtinv*dvol*tesgs
    elrhs(1:e%pnode) = elrhs(1:e%pnode) + aux*testf(1:e%pnode)
    
  end subroutine tem_elmrhs_dss
   
end module

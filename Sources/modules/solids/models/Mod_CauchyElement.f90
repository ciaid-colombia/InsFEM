module Mod_CauchyElement
   use typre
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_CauchyElement_tools
   use Mod_SUPCauchyElement_NH_tools
   use Mod_SUPCauchyElement_lin_tools
   use Mod_UPCauchyElement_NH_tools
   implicit none
 
   interface 

   !-------------------Irreducible formulation functions--------
   module subroutine sld_TimeInt2extF(e,nsteps,Integrator,dtinv2,gpdisp,elext,densi)
      implicit none
      class(FiniteElement), intent(in)    :: e
      type(TimeIntegratorDt2), intent(in):: Integrator
      integer(ip) , intent(in)    :: nsteps
      real(rp)    , intent(in)    :: gpdisp(e%ndime,nsteps-1)
      real(rp)    , intent(in)    :: densi,dtinv2
      real(rp)    , intent(inout) :: elext(e%ndime)
    end subroutine  

   module subroutine sld_elmrhu(e,dvolu,elext,elrhs)
      implicit none
      class(FiniteElement), intent(in) :: e
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
    end subroutine  

   module subroutine sld_ComputeExternalForces(e,densi,grnor,gravi,traction,elext,dvolu,force_factor)
      implicit none
      class(FiniteElement), intent(in)    :: e
      real(rp), intent(in) :: densi,grnor,gravi(e%ndime)
      real(rp), intent(in) :: traction(e%ndime),dvolu,force_factor
      real(rp), intent(inout) :: elext(e%ndime)
    end subroutine  

   module subroutine sld_buildMassMatrix(e,dvol,ndime,densi,elmat,dtinv)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip) , intent(in)    :: ndime
      real(rp)    , intent(in)    :: densi,dvol,dtinv
      real(rp)    , intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
    end subroutine  

   module subroutine sld_elmk(e,ndime,tn,c_elas,dvol,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)    :: tn,ndime
      real(rp),    intent(in)    :: c_elas(tn,tn),dvol
      real(rp),    intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
    end subroutine  

   module subroutine sld_elmgeo(e,ndime,tn,stress,dvol,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)    :: tn,ndime
      real(rp)   , intent(in)    :: stress(tn),dvol
      real(rp)   , intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
    end subroutine  

   module subroutine getInternalForce(e,ndime,tn,stress,dvol,force)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in) :: tn,ndime
      real(rp)   , intent(in) :: stress(tn),dvol
      real(rp)   , intent(inout) :: force(ndime,e%pnode)
    end subroutine  
 
    !This subroutines are found in the sldsup linear submodule
    module subroutine sup_elmVS(e,ndime,tn,dvol,P,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn 
       real(rp),             intent(in)    :: dvol,P(tn,tn)
       real(rp),             intent(inout) :: elmat(ndime,e%pnode,tn,e%pnode)
    end subroutine  
    
    module subroutine sup_elmVP(e,ndime,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(ndime,e%pnode,1,e%pnode)
    end subroutine  

    module subroutine sup_elmEU(e,ndime,tn,dvol,P,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol,P(tn,tn)
       real(rp),             intent(inout) :: elmat(tn,e%pnode,ndime,e%pnode)
     end subroutine  

    module subroutine sup_elmES(e,ndime,tn,dvol,D,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol,D(tn,tn)
       real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elmQU(e,ndime,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol
       real(rp),    intent(inout)       :: elmat(1,e%pnode,ndime,e%pnode)
    end subroutine  

    module subroutine sup_elmQP(e,ndime,dvol,D_vol,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol,D_vol
       real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
    end subroutine  

     !-------------------------TU---------------------------
   module subroutine sup_elm_usgs_EU_NH_dyn(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,densi,dtinv2,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd),dtinv2
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu,densi
      real(rp),    intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
    end subroutine  

    module subroutine sup_elm_usgs_ES(e,ndime,tn,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elm_usgs_EP(e,ndime,tn,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(tn,e%pnode,1,e%pnode)
    end subroutine  

   module subroutine sup_elm_usgs_QU_NH_dyn(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,densi,dtinv2,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd),densi,dtinv2
      real(rp),   intent(inout)  :: elmat(1,e%pnode,nd,e%pnode)
    end subroutine  

    module subroutine sup_elm_usgs_QS(e,ndime,tn,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(1,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elm_usgs_QP(e,ndime,tn,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(1,e%pnode,1,e%pnode)
    end subroutine  

     !-------------------------TS---------------------------
    module subroutine sup_elm_ssgs_VU(e,ndime,tn,dvol,c_dev,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(in)    :: c_dev(tn,tn)
       real(rp),             intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
    end subroutine  

    module subroutine sup_elm_ssgs_VS(e,ndime,tn,dvol,P,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: P(tn,tn),dvol
       real(rp),             intent(inout) :: elmat(ndime,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elm_ssgs_EU(e,ndime,tn,dvol,P,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: P(tn,tn),dvol
       real(rp),             intent(inout) :: elmat(tn,e%pnode,ndime,e%pnode)
    end subroutine  

    module subroutine sup_elm_ssgs_ES(e,ndime,tn,dvol,D,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: ndime,tn
       real(rp),             intent(in)    :: D(tn,tn),dvol
       real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elm_psgs_VU(e,ndime,dvol,K,G,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol,G,K
       real(rp),    intent(inout)       :: elmat(ndime,e%pnode,ndime,e%pnode)
    end subroutine  

    module subroutine sup_elm_psgs_VP(e,ndime,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol
       real(rp),    intent(inout)       :: elmat(ndime,e%pnode,1,e%pnode)
    end subroutine  

    module subroutine sup_elm_psgs_QU(e,ndime,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol
       real(rp),    intent(inout)       :: elmat(1,e%pnode,ndime,e%pnode)
    end subroutine  

    module subroutine sup_elm_psgs_QP(e,ndime,dvol,D_vol,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: ndime
       real(rp),    intent(in)          :: dvol,D_vol
       real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
    end subroutine  

     !-------------------------Dynamic SGS---------------------------

   module subroutine sup_elm_usgs_VU_NH_dynsgs(e,nd,tn,dvol,tautau,densi,dtinv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau,densi,dtinv
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,nd,e%pnode)
    end subroutine  

   module subroutine sup_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol,tautau,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,tn,e%pnode)
    end subroutine  

   module subroutine sup_elm_usgs_VP_NH_dynsgs(e,nd,tn,dvol,tautau,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,1,e%pnode)
    end subroutine  

    !-------------------------RHS---------------------------
    !This subroutines are found in the sldsup neohookean submodule
   module subroutine sup_rhs_U_NH_dyn(e,nd,tn,dvol,gpaccel,elrhs)
      implicit none 
      class(FiniteElement)       :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: dvol,gpaccel(nd)
      real(rp),    intent(inout) :: elrhs(nd,e%pnode)
   end subroutine

   module subroutine sup_rhs_U_NH(e,nd,tn,dvol,gdev,gpress,elrhs)
      implicit none 
      class(FiniteElement)       :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: dvol,gdev(tn),gpress(1) 
      real(rp),    intent(inout) :: elrhs(nd,e%pnode)
    end subroutine

   module subroutine sup_rhs_S_NH(e,nd,tn,dvol,mu,J,b,gdev,elrhs) 
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,J,dvol,gdev(tn),b(nd,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
    end subroutine  

   module subroutine sup_rhs_P_NH(e,nd,tn,dvol,mu,lam,J,b,gpress,elrhs,divdisp)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,lam,J,dvol,gpress(1),b(nd,nd),divdisp
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
    end subroutine  

    module subroutine sup_elmVS_NH(e,nd,tn,dvol,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: nd,tn
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(nd,e%pnode,tn,e%pnode)
    end subroutine  
    
    module subroutine sup_elmVP_NH(e,nd,dvol,elmat) 
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: nd
       real(rp),             intent(in)    :: dvol
       real(rp),             intent(inout) :: elmat(nd,e%pnode,1,e%pnode)
    end subroutine  

    module subroutine sup_elmEU_NH(e,nd,tn,dvol,J,mu,Finv,Fmat,gdev,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: nd,tn
       real(rp),             intent(in)    :: mu,J,dvol,Finv(nd,nd),Fmat(nd,nd),gdev(tn)
       real(rp),             intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
     end subroutine  

    module subroutine sup_elmES_NH(e,nd,tn,dvol,J,mu,elmat)
       implicit none
       class(FiniteElement), intent(in)    :: e
       integer(ip),          intent(in)    :: nd,tn
       real(rp),             intent(in)    :: dvol,J,mu
       real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
    end subroutine  

    module subroutine sup_elmQU_NH(e,nd,dvol,J,mu,lam,Finv,Fmat,gpress,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: nd
       real(rp),    intent(in)          :: dvol,mu,lam,J,gpress(1)
       real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
       real(rp),    intent(inout)       :: elmat(1,e%pnode,nd,e%pnode)
    end subroutine  

    module subroutine sup_elmQP_NH(e,nd,dvol,J,lam,elmat)
       implicit none
       class(FiniteElement), intent(in) :: e
       integer(ip), intent(in)          :: nd
       real(rp),    intent(in)          :: dvol,J,lam
       real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
    end subroutine  

   module subroutine sup_elm_usgs_ES_NH(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu
      real(rp),   intent(inout)  :: elmat(tn,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_EP_NH(e,nd,tn,dvol,J,mu,grsig,gdev,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: grsig(tn,nd),gdev(tn),J,mu
      real(rp),   intent(inout) :: elmat(tn,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_QS_NH(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd)
      real(rp),    intent(inout) :: elmat(1,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_QP_NH(e,nd,tn,dvol,J,mu,lam,gradp,gpress,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),gradp(nd)
      real(rp),   intent(inout) :: elmat(1,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_elm_ssgs_VU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),             intent(inout) :: elmat(nd,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_elm_ssgs_VS_NH(e,nd,tn,dvol,J,mu,elmat)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu
      real(rp),             intent(inout) :: elmat(nd,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_ssgs_EU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),             intent(in)    :: Fmat(nd,nd),Finv(nd,nd)
      real(rp),             intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_elm_ssgs_ES_NH(e,nd,tn,dvol,J,mu,elmat) 
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu
      real(rp),             intent(inout) :: elmat(tn,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_psgs_VU_NH(e,nd,dvol,J,mu,lam,gpress,Fmat,Finv,elmat)
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,gpress(1),J,mu,lam
      real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_elm_psgs_VP_NH(e,nd,dvol,J,lam,elmat)
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,lam,J
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_elm_psgs_QU_NH(e,nd,dvol,J,mu,lam,gpress,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,J,mu,lam,gpress(1)
      real(rp),    intent(in)          :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),    intent(inout)       :: elmat(1,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_elm_psgs_QP_NH(e,nd,dvol,J,lam,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,J,lam
      real(rp),    intent(inout)       :: elmat(1,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau,AU,elrhs)
      implicit none 
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau
      real(rp),             intent(in)    :: AU(tn)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau,AU
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH(e,nd,tn,dvol,mu,J,tau,grsig,gdev,Fmat,Finv,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd)
      real(rp),             intent(in)    :: grsig(tn,nd),AU(nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_ssgs_S_NH(e,nd,tn,dvol,mu,J,tau,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: AU(tn),dvol,tau,J,mu
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau,AU(nd)
      real(rp),   intent(in)    :: gpress(1),gradp(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,J,tau,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: lam,J,dvol
      real(rp),             intent(in)    :: tau,AU
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_extforces(e,nd,tn,dvol,mu,J,tau,Fmat,Finv,grsig,gdev,elrhs,elext)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,tau,grsig(tn,nd),gdev(tn)
      real(rp),             intent(in)    :: elext(nd),Fmat(nd,nd),Finv(nd,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH_extforces(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,elrhs,elext)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)              :: nd,tn
      real(rp),   intent(in)              :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)              :: lam,J,mu,tau
      real(rp),   intent(in)              :: gpress(1),elext(nd),gradp(nd)
      real(rp),   intent(inout)           :: elrhs(1,e%pnode)
 
   end subroutine  

   !--------------------------RESIDUAL CALCULATION----------------
   module subroutine sup_residualU_dyn(e,nd,accel,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: accel(nd)
      real(rp),             intent(inout) :: res(nd)
   end subroutine  

   module subroutine sup_residualU(e,nd,gradp,divstr,fext,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: gradp(nd),divstr(nd),fext(nd)
      real(rp),             intent(inout) :: res(nd)
   end subroutine  

   module subroutine sup_residualS(e,nd,tn,mu,J,b,gdev,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: b(nd,nd),gdev(tn),mu,J
      real(rp),             intent(inout) :: res(tn)
   end subroutine  

   module subroutine sup_residualP(e,nd,mu,lam,J,b,gpress,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd),gpress(1)
      real(rp),             intent(inout) :: res
   end subroutine  


   module subroutine sup_residualS_OSS(e,nd,tn,b,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: b(nd,nd)
      real(rp),             intent(inout) :: res(tn)
   end subroutine  

   module subroutine sup_residualP_OSS(e,nd,mu,lam,J,b,res)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd
      real(rp),             intent(in)    :: mu,lam,J,b(nd,nd)
      real(rp),             intent(inout) :: res
   end subroutine  

   !---------------------ADJOINT TESTING OF PROJECTION -----------

   module subroutine sup_rhs_ssgs_U_NH_repro(e,nd,tn,dvol,tau,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tau,gpres(tn)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_psgs_U_NH_repro(e,nd,tn,dvol,tau,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol
      real(rp),             intent(in)    :: tau,gpres
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_repro(e,nd,tn,dvol,mu,J,tau,grsig,Fmat,Finv,gdev,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau,gpres(nd),gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd),grsig(tn,nd)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_ssgs_S_NH_repro(e,nd,tn,dvol,mu,J,tau,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau,gpres(tn)
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH_repro(e,nd,tn,dvol,mu,lam,J,tau,gradp,gppress,Fmat,Finv,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: gradp(nd),gppress(nd),gpres(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_psgs_P_NH_repro(e,nd,tn,dvol,mu,lam,J,b,tau,gpres,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,lam,J,dvol,b(nd,nd)
      real(rp),             intent(in)    :: tau,gpres
      real(rp),             intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_U_NH_AU_dynsgs(e,nd,tn,dvol,tautau,elext,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tautau
      real(rp),             intent(in)    :: elext(nd),AU(nd)
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_U_NH_SGSACCEL_dynsgs(e,nd,tn,dvol,dtinv2,densi,tautau,vesgs,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,tautau
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_dynsgs(e,nd,tn,dvol,mu,J,tau,grsig,gdev,Fmat,Finv,dtinv2,densi,vesgs,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: mu,dvol,J,tau
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: Finv(nd,nd),Fmat(nd,nd),grsig(tn,nd)
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH_dynsgs(e,nd,tn,dvol,mu,lam,J,tau,gradp,gpress,Fmat,Finv,dtinv2,densi,vesgs,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,Fmat(nd,nd),Finv(nd,nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: gpress(1),gradp(nd)
      real(rp),   intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   !------------------------------HIGHER ORDER DERIVATIVES-------------

   module subroutine sup_elm_usgs_ES_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J
      real(rp),   intent(inout)  :: elmat(tn,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_EP_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J
      real(rp),    intent(inout) :: elmat(tn,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_QS_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,elmat)
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,tn,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_QP_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd)
      real(rp),   intent(inout)  :: elmat(1,e%pnode,1,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,tau,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: J,Jder(nd),gdev(tn),mu,tau
      real(rp),   intent(in)    :: AU(nd)
      real(rp),   intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,tau,AU,elrhs)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: Jder(nd),gpress(1),mu,J,tau,lam
      real(rp),   intent(in)    :: AU(nd)
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_EU_NH_dyn_nonli(e,nd,tn,dvol,J,mu,gdev,Jder,divFmat,divFinv,densi,dtinv2,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd),dtinv2
      real(rp),   intent(in)     :: Jder(nd),gdev(tn),mu,J,densi
      real(rp),    intent(inout) :: elmat(tn,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_elm_usgs_QU_NH_dyn_nonli(e,nd,tn,dvol,J,mu,lam,gpress,Jder,divFmat,divFinv,densi,dtinv2,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)     :: lam,J,mu,gpress(1),Jder(nd),densi,dtinv2
      real(rp),   intent(inout)  :: elmat(1,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_S_NH_dynsgs_nonli(e,nd,tn,dvol,mu,J,tau,gdev,Jder,divFmat,divFinv,dtinv2,densi,vesgs,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: divFmat(nd),divFinv(nd)
      real(rp),             intent(in)    :: mu,dvol,J,tau,Jder(nd)
      real(rp),             intent(in)    :: gdev(tn)
      real(rp),             intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),             intent(inout) :: elrhs(tn,e%pnode)
   end subroutine  

   module subroutine sup_rhs_usgs_P_NH_dynsgs_nonli(e,nd,tn,dvol,mu,lam,J,tau,gpress,Jder,divFmat,divFinv,dtinv2,densi,vesgs,elrhs)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),intent(in)    :: nd,tn
      real(rp),   intent(in)    :: dvol,divFmat(nd),divFinv(nd)
      real(rp),   intent(in)    :: lam,J,mu,tau
      real(rp),   intent(in)    :: Jder(nd),gpress(1)
      real(rp),   intent(in)    :: dtinv2,vesgs(nd),densi
      real(rp),   intent(inout) :: elrhs(1,e%pnode)
   end subroutine  

   !--------------------------UP FORMULATION------------------------------
   module subroutine up_elmVU_NH(e,nd,tn,dvol,J,mu,gdev,Fmat,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in)    :: e
      integer(ip),          intent(in)    :: nd,tn
      real(rp),             intent(in)    :: dvol,J,mu,gdev(tn)
      real(rp),             intent(in)    :: Fmat(nd,nd),Finv(nd,nd)
      real(rp),             intent(inout) :: elmat(nd,e%pnode,nd,e%pnode)
   end subroutine  

   module subroutine up_elmVP_NH(e,nd,dvol,mu,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)          :: nd
      real(rp),    intent(in)          :: dvol,mu
      real(rp),    intent(inout)       :: elmat(nd,e%pnode,1,e%pnode)
   end subroutine  

   !--------------------------DYNAMIC SGS------------------------------
   module subroutine up_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol,tautau,divstr,Finv,elmat)
      implicit none
      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau,Finv(nd,nd),divstr(nd)
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,nd,e%pnode)
    end subroutine  

   end interface

end module Mod_CauchyElement

#include "standard/Elmopes/Mod_CauchyElement_Irreducible.f90"   
#include "sldsup/Elmopes/linear/Mod_SUPCauchyElement.f90"   
#include "sldsup/Elmopes/linear/Mod_SUPCauchyElementSGS.f90"   
#include "sldsup/Elmopes/neohookean/Mod_SUPCauchyElement_NH.f90"   
#include "sldsup/Elmopes/neohookean/Mod_SUPCauchyElement_NH_SGS.f90"   
#include "sldsup/Elmopes/neohookean/Mod_SUPCauchyElement_NH_SGS_RHS.f90"   
#include "sldsup/Elmopes/neohookean/Mod_SUPCauchyElement_NH_Residual.f90"   
#include "sldup/Elmopes/neohookean/Mod_UPCauchyElement_NH.f90"
#include "sldup/Elmopes/neohookean/Mod_UPCauchyElement_NH_SGS.f90"

module Mod_LowMachElement
   use typre
   use Mod_Element
   use Mod_LowMachElement_OSS
   use Mod_LowMachElement_Newton
   use Mod_LowMachElement_Residual
   
contains

   subroutine ComputeTauCon(e,acden,staco,chale,timom,ticon)
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco,chale(2),acden
      real(rp)                   :: timom,ticon
      
      ticon = chale(2)*chale(2) / (staco*acden*timom*e%npol*e%npol)   
   end subroutine

   subroutine lmn_TimeIntegrationToElext_mom(e,Integrator,acden,dtinv,gpvel,elext)
      use typre
      use Mod_Element
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp)                   :: gpvel(e%ndime,*),elext(e%ndime)
      real(rp)                   :: acden,dtinv
      
      real(rp)                   :: gprhs(e%ndime)
      
      !Time integration
      call Integrator%GetRHS(e%ndime,gpvel(:,2),gprhs)
      elext = elext + acden*dtinv*gprhs(1:e%ndime)
       
   end subroutine
   
   subroutine lmn_TimeIntegrationToElext_ene(e,Integrator,acden,dtinv,gptem,elext)
      use typre
      use Mod_Element
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp)                   :: gptem(*),elext
      real(rp)                   :: acden,dtinv
      
      real(rp)                   :: gprhs(1)
      
      !Time integration
      call Integrator%GetRHS(1,gptem(2),gprhs)
      elext = elext + acden*dtinv*gprhs(1)
      
   end subroutine

   subroutine lmn_ComputeExternalForces_mom(e,acden,grnor,gravi,elext)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: acden, elext(e%ndime),grnor,gravi(e%ndime)
      
      elext = elext + acden*grnor*gravi(1:e%ndime)

   end subroutine
 
   subroutine lmn_ComputeExternalForces_ene(e,kfl_sourc,acsph,Source,gpsou,elext)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: elext,acsph,Source,gpsou
      integer(ip)                :: kfl_sourc 
      elext = elext + Source/acsph
      if (kfl_sourc==2) elext = elext + gpsou/acsph
   end subroutine

!U-V terms  
   subroutine lmn_elmbuv(e,dvolu,acden,dtinv,vtemp,testf_mom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu,acden,dtinv
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: jnode,idime,inode
      real(rp)                   :: aux1

      !rho*Mass/dt + rho*a·grad
      do jnode=1,e%pnode
         aux1 = acden*(e%shape(jnode,e%igaus)*dtinv + vtemp(jnode))
         do idime=1,e%ndime
            do inode=1,e%pnode
               elmat(idime,inode,idime,jnode) = (e%shape(inode,e%igaus) + &
               testf_mom(idime,inode))*aux1*dvolu + elmat(idime,inode,idime,jnode)
            end do
         end do
      end do

   end subroutine lmn_elmbuv
   
   subroutine lmn_elmbuv_nln(e,dvolu,visac,testf_mom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    - tau1*(L*v, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: visac
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: jnode,inode,jdime,pos(e%ndime,e%ndime)
      real(rp)                   :: aux1

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
     
      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            do jdime = 1,e%ndime
               elmat(jdime,inode,jdime,jnode) = elmat(jdime,inode,jdime,jnode) - &
               visac*testf_mom(jdime,inode)*dvolu* (sum(e%hessi(1:e%ndime,jnode)) + &
               0.333_rp*sum(e%hessi(pos(jdime,1:e%ndime),jnode)))
            enddo
         enddo
      enddo
   end subroutine lmn_elmbuv_nln

   subroutine lmn_elmvis_div(e,dvolt0,visac,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for viscosity for block U,V 
    !     mu*div(e·grad u)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement)   :: e
    real(rp)                   :: dvolt0
    real(rp),    intent(in)    :: visac
    real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

    integer(ip)                :: inode,jnode,idime,jdime
    real(rp)                   :: tmp

    tmp = dvolt0*visac
    forall (idime=1:e%ndime, inode=1:e%pnode, jdime=1:e%ndime, jnode=1:e%pnode)
       elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) + &
       tmp*e%cartd(idime,jnode)*e%cartd(jdime,inode)
    end forall

   end subroutine lmn_elmvis_div
 
   subroutine lmn_elmbuv_div(e,acvis,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !     -2/3*mu* (grad v,I div u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e   
      real(rp)                   :: dvol,acvis
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime

      forall (idime=1:e%ndime,inode=1:e%pnode,jdime=1:e%ndime, jnode=1:e%pnode)
         elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) - &
         0.666_rp*acvis*dvol*sum(e%cartd(:,inode))*e%cartd(jdime,jnode)
      end forall

   end subroutine lmn_elmbuv_div

   subroutine lmn_elmbuv_con(e,acden,ticon,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    rho*ticon*(div v, div u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e   
      real(rp)                   :: dvol, acden, ticon
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime

      forall (idime=1:e%ndime,inode=1:e%pnode,jdime=1:e%ndime, jnode=1:e%pnode)
         elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) + &
               acden*ticon*dvol*e%cartd(idime,inode)*e%cartd(jdime,jnode)
      end forall

   end subroutine lmn_elmbuv_con

! U-Q terms 
   subroutine lmn_elmbuq(e,timom,dvolu,acden,dtinv,vtemp,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     (div u, q) + tau1*(grad q, Lu) + tau1*(grad q, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: dvolu,timom,acden,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      real(rp)                   :: tmp,aux1

      tmp = acden*dtinv

      do jnode=1,e%pnode
         aux1 = (vtemp(jnode)*acden+tmp*e%shape(jnode,e%igaus))*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = elmat(1,inode,:,jnode) + &
            (-e%shape(jnode,e%igaus)*e%cartd(:,inode) + e%cartd(:,inode)*aux1)*dvolu*acden
         end do
      end do

   end subroutine lmn_elmbuq

   subroutine lmn_elmbuq_nln(e,dvolu,visac,acden,timom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q 
      !    - tau1*(grad q, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: timom
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: visac, acden
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: jnode,inode,jdime,pos(e%ndime,e%ndime)
      real(rp)                   :: aux1

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
     
      aux1 = acden*visac*timom*dvolu
      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            do jdime = 1,e%ndime
               elmat(1,inode,jdime,jnode) = elmat(1,inode,jdime,jnode) - & 
               aux1*e%cartd(jdime,inode)*(sum(e%hessi(1:e%ndime,jnode)) + &
               0.333_rp*sum(e%hessi(pos(jdime,1:e%ndime),jnode)))
            enddo
         enddo
      enddo 
   end subroutine lmn_elmbuq_nln

   !T-V term
   subroutine lmn_elmbtv(e,acden,ticon,dvol,aGradN,actex,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block T,V 
      !    rho*ticon*alpha*(a.gradT,divv)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)        :: e  
      real(rp),    intent(in)    :: aGradN(e%pnode)
      real(rp)                   :: dvol, acden, ticon, actex
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime

      do jnode=1,e%pnode
         do inode=1,e%pnode
            elmat(1:e%ndime,inode,1,jnode) =  elmat(1:e%ndime,inode,1,jnode) - &
               acden*actex*ticon*aGradN(jnode)*e%cartd(1:e%ndime,inode)*dvol
         end do
      end do
   end subroutine lmn_elmbtv

!P-V term
   subroutine lmn_elmbpv(e,dvolu,testf_mom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,V 
      !     - (div v, p) + tau1*(L*v, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)

      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         do inode=1,e%pnode
            elmat(1:e%ndime,inode,1,jnode) = -e%cartd(1:e%ndime,inode)*e%shape(jnode,e%igaus)*dvolu + &
               e%cartd(1:e%ndime,jnode)*testf_mom(1:e%ndime,inode)*dvolu + &
               elmat(1:e%ndime,inode,1,jnode)
         end do
      end do

   end subroutine lmn_elmbpv

!P-Q term
   subroutine lmn_elmbpq(e,dvolu,acden,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,Q
      !    (grad q, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,acden
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)

      integer(ip)                :: inode,jnode

      forall(inode=1:e%pnode,jnode=1:e%pnode)
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + dot_product(e%cartd(:,inode),e%cartd(:,jnode))*dvolu*acden
      end forall

   end subroutine lmn_elmbpq

   subroutine lmn_elmbpq_cnf(e,epspe,acden,acvis,dvol,elmat)

      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,Q
      !   eps*(q,p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: acden,acvis,dvol,epspe
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)

      integer(ip)                :: inode,jnode

      forall(inode=1:e%pnode,jnode=1:e%pnode)
            elmat(1,inode,1,jnode) = epspe*acden/acvis*dvol*e%shape(jnode,e%igaus)*e%shape(inode,e%igaus) + elmat(1,inode,1,jnode)
      end forall

   end subroutine lmn_elmbpq_cnf

!T-W terms
   subroutine lmn_elmbtw(e,dvol,acden,dtinv,vtemp,testf_ene,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block T,W
      !   (w, a·grad T) - tau1*(L*w, LT) + (W, T/dt) - tau1*(L*W, T/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf_ene(e%pnode)
      real(rp),    intent(in)    :: dvol,acden,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)

      integer(ip)                :: jnode,inode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         !rho*Mass/dt + rho*a·grad
         aux1 = acden*(e%shape(jnode,e%igaus)*dtinv + vtemp(jnode))
!         aux1 = (acden*dtinv+acrcp)*(e%shape(jnode,e%igaus)) + acden*vtemp(jnode)
         elmat(1,1:e%pnode,1,jnode) = elmat(1,1:e%pnode,1,jnode) + aux1*dvol*(e%shape(1:e%pnode,e%igaus) + testf_ene(1:e%pnode))
      end do

   end subroutine lmn_elmbtw

   subroutine lmn_elmbtw_nln(e,dvolu,actco,testf_ene,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    - tau1*(L*v, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: testf_ene(e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: actco
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)

      integer(ip)                :: jnode,inode
 
      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) - &
            actco*testf_ene(inode)*dvolu* (sum(e%hessi(1:e%ndime,jnode))) 
         enddo
      enddo
   end subroutine lmn_elmbtw_nln

!U-RHS term  
   subroutine lmn_elmrhu(e,Integrator,dvolu,dtinv,ticon,gpden,gpvel,testf_mom,elext,elrhs)
      use Mod_TimeIntegrator
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode),elext(e%ndime)
      real(rp),    intent(in)    :: dvolu,ticon,dtinv,gpden(*),gpvel(e%ndime,*)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)

      integer(ip)                :: inode,idime
      real(rp)                   :: aux(e%ndime),tmp,gprhs(1)

      call Integrator%GetRHS(1,gpden(2),gprhs(1))
      do inode=1,e%pnode
         do idime=1,e%ndime
            aux(idime) = (e%shape(inode,e%igaus) + testf_mom(idime,inode))*dvolu
            elrhs(idime,inode) = elext(idime)*aux(idime) + elrhs(idime,inode) -dtinv* &
            dvolu*ticon*(gpden(1)*Integrator%Coefficients(1) - gprhs(1)) *e%cartd(idime,inode)
         end do
      end do

   end subroutine lmn_elmrhu

!P-RHS terms
   subroutine lmn_elmrhp(e,Integrator,timom,dvolu,dtinv,gpden,elext,elrhs)
      use Mod_TimeIntegrator
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    tau1*(f, grad q) + tau1*(grad q, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(in)    :: timom,dvolu,gpden(*),dtinv
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

      integer(ip)                :: inode
      real(rp)                   :: tmp1,aux(e%ndime),gprhs(1)

      call Integrator%GetRHS(1,gpden(2),gprhs(1))
      tmp1 = dvolu*timom*gpden(1)
      aux = elext*tmp1
      do inode=1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + dot_product(e%cartd(:,inode),aux) -dtinv* &
         e%shape(inode,e%igaus) * (gpden(1)*Integrator%Coefficients(1)-gprhs(1)) *dvolu
      end do

   end subroutine lmn_elmrhp

   subroutine lmn_elmrhp_cnf(e,epspe,acden,acvis,dvolu,gppre,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    tau1*(f, grad q) + tau1*(grad q, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gppre,acden,acvis,epspe,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

      integer(ip)                :: inode
      
      do inode=1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + epspe*acden/acvis*gppre*e%shape(inode,e%igaus)*dvolu
      end do

   end subroutine lmn_elmrhp_cnf
   
   
!T-RHS term 
   subroutine lmn_elmrht(e,Integrator,dvolu,dtinv,testf_ene,elext,actex,acden,gptem,gtemp,gppth,elrhs)
      use Mod_TimeIntegrator
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (w, Q) + (Pe, Q) + (w+Pe, alpha*T*pth/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp),    intent(in)    :: testf_ene(e%pnode),elext
      real(rp),    intent(in)    :: dvolu,dtinv,actex,acden
      real(rp),    intent(in)    :: gtemp,gppth(*),gptem(*)
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

      integer(ip)                :: inode
      real(rp)                   :: aux,tmp,gprhs(1)

      call Integrator%GetRHS(1,gppth(3),gprhs(1))
      tmp = dtinv*actex*gtemp* (gppth(1)*Integrator%Coefficients(1) - gprhs(1))
      do inode=1,e%pnode
         aux = dvolu*(e%shape(inode,e%igaus) + testf_ene(inode))
         elrhs(1,inode) = elrhs(1,inode) + (elext + tmp) * aux
      end do

   end subroutine lmn_elmrht
   
   !-----------------------------------------------------------------------------------
   !Compute adjoint (L*)
    
   subroutine lmn_ComputeTestf_mom(e,acden,timom,AGradV,testf_mom)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: acden, AGradV(e%pnode), testf_mom(e%ndime,e%pnode), timom
      integer(ip) :: inode
      
      do inode=1,e%pnode
         testf_mom(1:e%ndime,inode) = acden*AGradV(inode)*timom
      end do

   end subroutine

   subroutine lmn_ComputeTestf_ene(e,acden,tiene,AGradT,testf_ene)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: acden,AGradT(e%pnode),testf_ene(e%pnode),tiene
      integer(ip) :: inode
      
      do inode=1,e%pnode
         testf_ene(inode) = (acden*AGradT(inode))*tiene !-acrcp*e%shape(1:e%pnode,e%igaus))*tiene
      end do
  
   end subroutine

   subroutine lmn_ComputeTestfNonLinear_mom(e,acvis,timom,kfl_stabm,testf_mom)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: acvis, timom, testf_mom(e%ndime,e%pnode)
      integer(ip)         :: kfl_stabm
      
      integer(ip) :: inode, idime, pos(e%ndime,e%ndime)
      real(rp)    :: aux1 

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

      !Stabilization terms : +tau visco *lapla(v)
      aux1 = real(kfl_stabm)*timom*acvis
      do inode = 1,e%pnode
         do idime = 1,e%ndime
            testf_mom(idime,inode) = testf_mom(idime,inode) + aux1 * &
            (sum(e%hessi(1:e%ndime,inode)) + 0.333_rp*sum(e%hessi(pos(idime,1:e%ndime),inode)))
         enddo
      enddo

   end subroutine

   subroutine lmn_ComputeTestfNonLinear_ene(e,actco,tiene,kfl_stabm,testf_ene)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp)            :: actco, tiene, testf_ene(e%pnode)
      integer(ip)         :: kfl_stabm
      
      integer(ip) :: inode
      real(rp)    :: aux1 

      !Stabilization terms : +tau actco *lapla(w)
      aux1 = real(kfl_stabm)*tiene*actco
      do inode = 1,e%pnode
         testf_ene(inode) = testf_ene(inode) + aux1 * sum(e%hessi(1:e%ndime,inode))
      enddo

   end subroutine

   !For dynamic subscales
   subroutine lmn_ComputeTestfDSS_mom(e,acden,timom,dtinv,testf_mom)
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp)            :: acden,dtinv,testf_mom(e%ndime,e%pnode),timom
      real(rp)            :: gprhs(1)
      integer(ip)         :: inode
     
      do inode=1,e%pnode 
         testf_mom(1:e%ndime,inode) = testf_mom(1:e%ndime,inode) - timom*dtinv* &
         e%shape(inode,e%igaus)*acden
      end do

   end subroutine

   subroutine lmn_ComputeTestfDSS_ene(e,acden,tiene,dtinv,testf_ene)
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp)            :: acden,dtinv,testf_ene(e%pnode),tiene,gprhs(1)
      
      testf_ene(1:e%pnode) = testf_ene(1:e%pnode) - dtinv*tiene* &
      e%shape(1:e%pnode,e%igaus)*acden

   end subroutine

   !-----------------------------------------------------------------------------------
   !Dynamic subscales RHS
   subroutine lmn_elmrhu_dss(e,dvol,dtinv,acden,testf,vesgs,elrhu)
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: dvol,dtinv,acden
      real(rp),    intent(in)    :: testf(e%ndime,e%pnode)
      real(rp),    intent(in)    :: vesgs(e%ndime)
      real(rp),    intent(inout) :: elrhu(e%ndime,e%pnode)
      
      integer(ip) :: inode,idime
      real(rp) :: aux(e%ndime)
      
      aux(1:e%ndime) = acden*dtinv*dvol*vesgs(1:e%ndime)
      !Compute contributions to RHS : Block U
      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhu(idime,inode) = elrhu(idime,inode) + aux(idime)*testf(idime,inode) 
         end do
      end do
   end subroutine 
  
   subroutine lmn_elmrhp_dss(e,dvol,dtinv,acden1,acden2,timom,vesgs,elrhp)
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: dvol,dtinv,timom,acden1,acden2
      real(rp),    intent(in)    :: vesgs(e%ndime)
      real(rp),    intent(inout) :: elrhp(1,e%pnode)
      
      integer(ip) :: inode
      real(rp) :: aux
 
      !Compute contributions to RHS : Block P
      aux = acden1*acden2*dtinv*timom*dvol
      do inode=1,e%pnode
         elrhp(1,inode) = elrhp(1,inode) + aux*dot_product(e%cartd(1:e%ndime,inode),vesgs(1:e%ndime)) 
      end do
   end subroutine  
  
   subroutine lmn_elmrht_dss(e,dvol,dtinv,acden,testf_ene,tesgs,elrht)
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: dvol,dtinv,timom,acden
      real(rp),    intent(in)    :: tesgs, testf_ene(e%pnode)
      real(rp),    intent(inout) :: elrht(1,e%pnode)
      
      integer(ip) :: inode
      real(rp) :: aux
 
      !Compute contributions to RHS : Block T
      aux = acden*dtinv*tesgs*dvol
      do inode=1,e%pnode
         elrht(1,inode) = elrht(1,inode) + aux*testf_ene(inode) 
      end do
   end subroutine  
  
end module

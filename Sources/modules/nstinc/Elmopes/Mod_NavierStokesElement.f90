module Mod_NavierStokesElement
   use typre
   use Mod_Element
   use Mod_NavierStokesElement_Residuals
   use Mod_NavierStokesElement_OSS
   
contains
   
   subroutine ComputeTauDiv(e,staco,chale,timom,tidiv)
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: staco,chale(2),timom
      real(rp),    intent(inout) :: tidiv
      real(rp)                   :: appop

      appop = staco*chale(2)*chale(2)/(e%npol*e%npol)                 ! Alg. PPO inverse (c*h^2)
      tidiv = appop/timom
   end subroutine
   
   subroutine nsi_TimeIntegrationToElTemp(e,Integrator,acden,dtinv,gpvel,eltemp)
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)       :: e
      type(TimeIntegratorDt1)    :: Integrator
      real(rp),    intent(inout) :: eltemp(e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime,*),acden,dtinv
      real(rp)                   :: gprhs(e%ndime)
      
      !Time integration
      call Integrator%GetRHS(e%ndime,gpvel(:,2),gprhs)
      eltemp = eltemp + acden*dtinv*gprhs(1:e%ndime)
   end subroutine
   
   subroutine nsi_ComputeExternalForces(e,acden,grnor,gravi,elext)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: acden,grnor,gravi(e%ndime)
      real(rp),    intent(inout) :: elext(e%ndime)
      
      elext = elext + acden*grnor*gravi(1:e%ndime)
   end subroutine
   
   subroutine nsm_elmrhu(e,dvolu,testf,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: testf(e%pnode),elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux,tmp

      do inode=1,e%pnode
         aux = (e%shape(inode,e%igaus) + testf(inode))*dvolu
         elrhs(:,inode) = (elext(:)+eltemp(:))*aux + elrhs(:,inode)
      end do
   end subroutine nsm_elmrhu   
  
   subroutine nsm_elmbuv(e,dvolu,denac,dtinv,vtemp,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu,denac,dtinv
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
      integer(ip)                :: jnode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         !rho*Mass/dt + rho*a·grad
         aux1 = denac*(e%shape(jnode,e%igaus)*dtinv + vtemp(jnode))
         elmat(1:e%pnode,jnode) = (e%shape(1:e%pnode,e%igaus) + testf(1:e%pnode))*aux1*dvolu + elmat(1:e%pnode,jnode)
      end do
   end subroutine nsm_elmbuv
   
   subroutine nsm_elmbuv_lap(e,dvolu,visac,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    - tau1*(L*v, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: visac
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
      integer(ip)                :: jnode,inode
      real(rp)                   :: aux1

      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            elmat(inode,jnode) = elmat(inode,jnode) - testf(inode)*visac*sum(e%hessi(1:e%ndime,jnode))*dvolu
         enddo
      enddo
   end subroutine nsm_elmbuv_lap
   
   subroutine nsm_elmbuq_lap(e,dvolu,visac,timom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q 
      !    - tau1*(grad q, visco*lapla(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: timom
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: visac
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: jnode,inode,jdime
      real(rp)                   :: aux1

      aux1 = visac*timom*dvolu
      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            do jdime = 1,e%ndime
               elmat(1,inode,jdime,jnode) = elmat(1,inode,jdime,jnode) - aux1*(e%cartd(jdime,inode)*sum(e%hessi(1:e%ndime,jnode)))
            enddo
         enddo
      enddo
   end subroutine nsm_elmbuq_lap
   
   subroutine nsm_elmrhp(e,timom,dvolu,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    tau1*(f, grad q) + tau1*(grad q, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: tmp1,aux(e%ndime)

      tmp1 = dvolu*timom
      aux = (elext+eltemp)*tmp1
      do inode=1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + dot_product(e%cartd(:,inode),aux)
      end do
   end subroutine nsm_elmrhp
  
   subroutine nsm_elmbpv(e,dvolu,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,V 
      !     - (div v, p) + tau1*(L*v, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         aux1 = e%shape(jnode,e%igaus)*dvolu
         do inode=1,e%pnode
            aux2 = testf(inode)*dvolu
            elmat(:,inode,1,jnode) = -e%cartd(:,inode)*aux1 + &
                    e%cartd(:,jnode)*aux2 + elmat(:,inode,1,jnode)
         end do
      end do
   end subroutine nsm_elmbpv

   subroutine nsm_elmbuq(e,timom,dvolu,denac,dtinv,vtemp,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     (div u, q) + tau1*(grad q, Lu) + tau1*(grad q, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: dvolu,timom,denac,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: tmp,aux1

      tmp = denac*dtinv

      do jnode=1,e%pnode
         aux1 = (vtemp(jnode)*denac+tmp*e%shape(jnode,e%igaus))*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = (e%shape(inode,e%igaus)*e%cartd(:,jnode) +&
               e%cartd(:,inode)*aux1)*dvolu + elmat(1,inode,:,jnode)
         end do
      end do
   end subroutine nsm_elmbuq

   subroutine nsm_elmbuq_slip(e,timom,dvolu,denac,dtinv,vtemp,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     (div u, q) + tau1*(grad q, Lu) + tau1*(grad q, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: dvolu,timom,denac,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: tmp,aux1

      tmp = denac*dtinv
      do jnode=1,e%pnode
         aux1 = (vtemp(jnode)*denac+tmp*e%shape(jnode,e%igaus))*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = (-e%shape(jnode,e%igaus)*e%cartd(:,inode) +&
               e%cartd(:,inode)*aux1)*dvolu + elmat(1,inode,:,jnode)
         end do
      end do
   end subroutine nsm_elmbuq_slip

   subroutine nsm_elmbpv_slip(e,dvolu,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,V 
      !     - (div v, p) + tau1*(L*v, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         do inode=1,e%pnode
            aux1 = e%shape(inode,e%igaus)*dvolu
            aux2 = testf(inode)*dvolu
            elmat(:,inode,1,jnode) = e%cartd(:,jnode)*aux1 +&
               e%cartd(:,jnode)*aux2 + elmat(:,inode,1,jnode)
         end do
      end do
   end subroutine nsm_elmbpv_slip
   
   subroutine nsm_elmbpq(e,dvolu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,Q
      !    (grad q, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode

      forall(inode=1:e%pnode,jnode=1:e%pnode)
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) +&
               dot_product(e%cartd(:,inode),e%cartd(:,jnode))*dvolu
      end forall
   end subroutine nsm_elmbpq
   
   subroutine nsm_elmbpq_pen(e,acvis,dvolu,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,Q
      !    (grad q, grad p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dvolu,acvis
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      real(rp)                   :: penalty
      integer(ip)                :: inode,jnode,k
      
      k=7
      penalty = (10.0_rp**(-k)*acvis)*dvolu
      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmat(1,inode,1,inode) = elmat(1,inode,1,inode) + &
               e%shape(inode,e%igaus)*penalty*e%shape(jnode,e%igaus)
         end do
      end do
   end subroutine nsm_elmbpq_pen

   subroutine nsm_elmdiv(e,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !     (div v, div u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e   
      real(rp)                   :: dvol
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime

      forall (idime=1:e%ndime,inode=1:e%pnode,jdime=1:e%ndime,jnode=1:e%pnode)
         elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) + &
            dvol*e%cartd(idime,inode)*e%cartd(jdime,jnode)
      end forall
   end subroutine nsm_elmdiv
   
   subroutine nsm_elmdiv_explicit(e,dvol,elvel,elmat,elrhu)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !     (div v, div u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e   
      real(rp)                   :: dvol
      real(rp),    intent(inout) :: elrhu(e%ndime,e%mnode)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      real(rp)                   :: elvel(e%ndime,e%pnode)
      real(rp)                   :: divu
      integer(ip)                :: inode,jnode,idime,jdime

      do jnode = 1,e%pnode
         do jdime = 1,e%ndime
            do inode = 1,e%pnode
               do idime = 1,e%ndime
                  if (idime == jdime) then
                     elmat(idime,inode,jdime,jnode) =  &
                     elmat(idime,inode,jdime,jnode) + &
                        dvol*e%cartd(idime,inode)*e%cartd(jdime,jnode)
                  else
                     elrhu(idime,inode) = elrhu(idime,inode) - &
                        dvol*e%cartd(idime,inode)*e%cartd(jdime,jnode)*elvel(jdime,jnode)
                  endif
               enddo
            enddo
         enddo
      enddo
   end subroutine nsm_elmdiv_explicit
   
   subroutine nsm_elmvis_div(e,dvolt0,visac,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for viscosity for block U,V 
    !     mu*div(e·grad u)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement)       :: e
    real(rp)                   :: dvolt0
    real(rp),    intent(in)    :: visac
    real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
    integer(ip)                :: inode,jnode,idime,jdime
    real(rp)                   :: tmp

    tmp = dvolt0*visac
    forall (idime=1:e%ndime, inode=1:e%pnode, jdime=1:e%ndime, jnode=1:e%pnode)
       elmat(idime,inode,jdime,jnode) =    &
            elmat(idime,inode,jdime,jnode) + tmp*e%cartd(idime,jnode)*e%cartd(jdime,inode)
    end forall
   end subroutine nsm_elmvis_div
  
   subroutine nsm_ComputeTestf(e,acden,timom,AGradV,testf)
      implicit none
      class(FiniteElement)     :: e
      real(rp),  intent(in)    :: acden,AGradV(e%pnode),timom
      real(rp),  intent(inout) :: testf(e%pnode)
      
      testf = acden*AGradV*timom
   end subroutine
   
   subroutine nsm_ComputeTestfNonLinear(e,acvis,timom,kfl_stabm,testf)
      implicit none
      class(FiniteElement)     :: e
      real(rp),  intent(in)    :: acvis,timom
      real(rp),  intent(inout) :: testf(e%pnode)
      integer(ip), intent(in)  :: kfl_stabm
      integer(ip)              :: inode
      real(rp)                 :: aux1
      
      !Stabilization terms : +tau visco *lapla(v)
      aux1 = real(kfl_stabm)*timom*acvis
      do inode = 1,e%pnode
         testf(inode) = testf(inode) + aux1*sum(e%hessi(1:e%ndime,inode))
      enddo
   end subroutine 
  
end module

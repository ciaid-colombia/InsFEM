module Mod_LowMachElement_OSS

   use typre
   use Mod_Element
   
contains

! U-V terms 
   subroutine lmn_elmbuv_trm(e,dvolu,acden,dtinv,vtemp,testf_mom,elmat)
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu,acden,dtinv
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: jnode,idime,inode
      real(rp)                   :: aux1,aux2

      !rho*Mass/dt + rho*a·grad
      do jnode=1,e%pnode
         aux1 = acden*e%shape(jnode,e%igaus)*dtinv*dvolu
         aux2 = acden*vtemp(jnode)*dvolu
         do idime=1,e%ndime
            do inode=1,e%pnode
               elmat(idime,inode,idime,jnode) = (e%shape(inode,e%igaus) + &
               testf_mom(idime,inode))*aux2 + aux1*e%shape(inode,e%igaus) + elmat(idime,inode,idime,jnode)
            end do
         end do
      end do

   end subroutine lmn_elmbuv_trm

! U-Q terms 
   subroutine lmn_elmbuq_trm(e,timom,dvolu,acden,dtinv,vtemp,elmat)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: dvolu,timom,acden,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1

      do jnode=1,e%pnode
         aux1 = vtemp(jnode)*acden*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = elmat(1,inode,:,jnode) + &
            (-e%shape(jnode,e%igaus)*e%cartd(:,inode) + e%cartd(:,inode)*aux1)*dvolu*acden
         end do
      end do

   end subroutine lmn_elmbuq_trm

!T-W terms
   subroutine lmn_elmbtw_trm(e,dvol,acden,dtinv,vtemp,testf_ene,elmat)
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf_ene(e%pnode)
      real(rp),    intent(in)    :: dvol,acden,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)

      integer(ip)                :: jnode,inode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         aux1 = acden*e%shape(jnode,e%igaus)*dtinv*dvol
!         aux1 = (acden*dtinv+acrcp)*(e%shape(jnode,e%igaus))*dvol
         aux2 = acden*vtemp(jnode)*dvol
         elmat(1,1:e%pnode,1,jnode) = elmat(1,1:e%pnode,1,jnode) + aux1*e%shape(1:e%pnode,e%igaus) + aux2*(e%shape(1:e%pnode,e%igaus) + testf_ene(1:e%pnode))
      end do

   end subroutine lmn_elmbtw_trm

!U-RHS term  
   subroutine lmn_elmrhu_trm(e,Integrator,dvolu,dtinv,ticon,gpden,gpvel,testf_mom,elext,elrhs)
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
      real(rp)                   :: aux(e%ndime),tmp,gprhs_vel(e%ndime),gprhs_den(1)

      call Integrator%GetRHS(e%ndime,gpvel(:,2),gprhs_vel)
      call Integrator%GetRHS(1,gpden(2),gprhs_den(1))

      do inode=1,e%pnode
         do idime=1,e%ndime
            aux(idime) = (e%shape(inode,e%igaus) + testf_mom(idime,inode))*dvolu
            elrhs(idime,inode) = elext(idime)*aux(idime) + &
            gpden(1)*dtinv*gprhs_vel(idime) * e%shape(inode,e%igaus)*dvolu - &
            dtinv*dvolu*ticon* (gpden(1)*Integrator%Coefficients(1) - gprhs_den(1))* &
            e%cartd(idime,inode) + elrhs(idime,inode) 
         end do
      end do

   end subroutine lmn_elmrhu_trm

!T-RHS term 
   subroutine lmn_elmrht_trm(e,Integrator,dvolu,dtinv,testf_ene,elext,actex,acden,gptem,gtemp,gppth,elrhs)
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
      real(rp)                   :: aux,tmp,gprhs_tem(1),gprhs_pth(1)

      call Integrator%GetRHS(1,gptem(2),gprhs_tem)
      call Integrator%GetRHS(1,gppth(3),gprhs_pth(1))
      tmp = dtinv*actex*gtemp* (gppth(1)*Integrator%Coefficients(1) - gprhs_pth(1))
      do inode=1,e%pnode
         aux = dvolu*(e%shape(inode,e%igaus) + testf_ene(inode))
         elrhs(1,inode) = elrhs(1,inode) + (elext + tmp) * aux + acden*dtinv*gprhs_tem(1)*e%shape(inode,e%igaus)*dvolu
      end do

   end subroutine lmn_elmrht_trm
   
   subroutine lmn_elmrhu_oss(e,ticon,dvol,testf_mom,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau1*(L* v, proj(Res(u))) + tau2*(grad v, proj(Res(p)))
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none
    class(FiniteElement)       :: e
    real(rp),    intent(in)    :: gprep(e%ndime+1)
    real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
    real(rp),    intent(in)    :: ticon,dvol
    real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
    integer(ip)                :: inode,idime

    do inode=1,e%pnode
       do idime=1,e%ndime
          elrhs(idime,inode) = - gprep(idime)*testf_mom(idime,inode)*dvol &
               + gprep(e%ndime+1)*e%cartd(idime,inode)*ticon*dvol + elrhs(idime,inode)
       end do
    end do

  end subroutine lmn_elmrhu_oss
   
   subroutine lmn_elmrhp_oss(e,timom,dvol,acden,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau1*(grad q, proj(Res(u)))
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none
    class(FiniteElement)       :: e
    real(rp),    intent(in)    :: gprep(e%ndime)
    real(rp),    intent(in)    :: timom,dvol,acden
    real(rp),    intent(inout) :: elrhs(1,e%pnode)
    integer(ip)                :: inode,idime

    do inode=1,e%pnode
       do idime=1,e%ndime
          elrhs(1,inode) = - e%cartd(idime,inode)*gprep(idime)*dvol*timom*acden + elrhs(1,inode)
       end do
    end do

   end subroutine lmn_elmrhp_oss
  
   subroutine lmn_elmrht_oss(e,dvol,testf_ene,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau3*(test w, proj(Res))
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none
    class(FiniteElement)       :: e
    real(rp),    intent(in)    :: gprep
    real(rp),    intent(in)    :: dvol,testf_ene(e%pnode)
    real(rp),    intent(inout) :: elrhs(1,e%pnode)
    integer(ip)                :: inode

    do inode=1,e%pnode
       elrhs(1,inode) = - testf_ene(inode)*gprep*dvol + elrhs(1,inode)
    end do

   end subroutine lmn_elmrht_oss
  
end module


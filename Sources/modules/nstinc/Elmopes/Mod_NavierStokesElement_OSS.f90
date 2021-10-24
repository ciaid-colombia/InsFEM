module Mod_NavierStokesElement_OSS
   use typre
   use Mod_Element
   
contains
   
   subroutine nsm_elmrhu_dss(e,dvol,dtinv,acden,testf,vesgs,elrhu)
      use typre
      implicit none
      class(FiniteElement)       :: e
      real(rp)                   :: dvol,dtinv,acden
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: vesgs(*)
      real(rp),    intent(inout) :: elrhu(e%ndime,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux(e%ndime)
      
      aux = acden*vesgs(1:e%ndime)*dtinv
      
      !Compute contributions to RHS : Block U
      do inode=1,e%pnode
         elrhu(1:e%ndime,inode) = aux*testf(inode)*dvol + elrhu(1:e%ndime,inode)
      end do
   end subroutine
   
   subroutine nsm_elmrhp_dss(e,dvol,dtinv,acden,timom,vesgs,elrhp)
      use typre
      implicit none
      class(FiniteElement)       :: e
      real(rp)                   :: acden,dvol,dtinv,timom
      real(rp),    intent(in)    :: vesgs(*)
      real(rp),    intent(inout) :: elrhp(1,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux
            
      !Compute contributions to RHS : Block P
      aux = acden*dtinv*timom*dvol
      do inode=1,e%pnode
         elrhp(1,inode) = aux*dot_product(e%cartd(1:e%ndime,inode),vesgs(1:e%ndime)) + elrhp(1,inode)
      end do
   end subroutine

   !We skip the part of the residual which belongs to the FiniteElement space
   subroutine nsm_elmrhu_trm(e,dvolu,testf,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + (v, u_n/dt) 
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
         aux = e%shape(inode,e%igaus)*dvolu 
         elrhs(:,inode) = (elext(:)+eltemp(:))*aux + elrhs(:,inode)
      end do
   end subroutine nsm_elmrhu_trm
   
   subroutine nsm_elmbuv_trm(e,dvolu,denac,dtinv,vtemp,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) 
      ! The finite element part of the residual is skiped (OSS skip)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu,denac,dtinv
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
      integer(ip)                :: jnode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         !rho*Mass/dt + rho*a·grad
         aux1 = denac*(e%shape(jnode,e%igaus)*dtinv + vtemp(jnode))
         elmat(1:e%pnode,jnode) = (e%shape(1:e%pnode,e%igaus)*aux1 + testf(1:e%pnode)*denac*vtemp(jnode))*dvolu + elmat(1:e%pnode,jnode)
      end do
   end subroutine nsm_elmbuv_trm 
    
   subroutine nsm_elmrhp_trm(e,timom,dvolu,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for OSS, 
      ! if the finite elemen part of the solution is skiped
      ! Nothing needs to be done  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

   end subroutine nsm_elmrhp_trm
   
   subroutine nsm_elmbuq_trm(e,timom,dvolu,denac,dtinv,vtemp,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,Q
      !  when the finite element part of the residual is skiped
      !     (div u, q) + tau1*(grad q, Lu) 
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
         aux1 = vtemp(jnode)*denac*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = (e%shape(inode,e%igaus)*e%cartd(:,jnode) + &
               e%cartd(:,inode)*aux1)*dvolu + elmat(1,inode,:,jnode)
         end do
      end do
   end subroutine nsm_elmbuq_trm
   
   subroutine nsm_elmrhu_oss(e,tidiv,dvol,testf,gprep,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for OSS
      !    tau*(L*v, proj(Res))
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: gprep(e%ndime+1)
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: tidiv,dvol
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux2,tmp2

      tmp2 = tidiv*dvol
      do inode=1,e%pnode
         aux2 = testf(inode)*dvol
         do idime=1,e%ndime
            elrhs(idime,inode) = elrhs(idime,inode) + gprep(idime)*aux2 &
               + gprep(e%ndime+1)*e%cartd(idime,inode)*tmp2
         end do
      end do
   end subroutine nsm_elmrhu_oss
   
   subroutine nsm_elmrhp_oss(e,timom,dvol,gprep,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for OSS
      !    tau1*(grad q, proj(Res))
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: gprep(e%ndime)
      real(rp),    intent(in)    :: timom,dvol
      real(rp),    intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: tmp1

      tmp1 = dvol*timom
      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = e%cartd(idime,inode)*gprep(idime)*tmp1 + elrhs(1,inode)
         end do
      end do
   end subroutine nsm_elmrhp_oss
  
   subroutine nsm_elmrhu_split(e,dvolu,testf,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Split OSS
      !    (v, f) + (v, u_n/dt) 
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
         aux = e%shape(inode,e%igaus)*dvolu 
         elrhs(:,inode) = (elext(:)+eltemp(:))*aux + elrhs(:,inode)
      end do
   end subroutine nsm_elmrhu_split
   
   subroutine nsm_elmbuv_split(e,dvolu,denac,dtinv,vtemp,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for split OSS for block U,V 
      !   (v, a·grad u) + s*(v,u) + tau1*(L*v, Lu) + (v, u/dt) 
      ! The finite element part of the residual is skiped (OSS skip)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: vtemp(e%pnode)
      real(rp),    intent(in)    :: testf(e%pnode)
      real(rp),    intent(in)    :: dvolu,denac,dtinv
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
      integer(ip)                :: jnode
      real(rp)                   :: aux1,aux2

      do jnode=1,e%pnode
         !rho*Mass/dt + rho*a·grad
         aux1 = denac*(e%shape(jnode,e%igaus)*dtinv + vtemp(jnode))
         elmat(1:e%pnode,jnode) = (e%shape(1:e%pnode,e%igaus)*aux1 + testf(1:e%pnode)*denac*vtemp(jnode))*dvolu + elmat(1:e%pnode,jnode)
      end do
   end subroutine nsm_elmbuv_split
   
   subroutine nsm_elmrhp_split(e,timom,dvolu,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for OSS, 
      ! if the finite elemen part of the solution is skiped
      ! Nothing needs to be done  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

   end subroutine nsm_elmrhp_split
   
   subroutine nsm_elmrhp_split2(e,timom,dvolu,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for split OSS, 
      ! if the external force term is NOT skiped in stabilization term
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%pnode)
      real(rp)                   :: aux(e%ndime),tmp1
      integer(ip)                :: inode
      
      tmp1 = dvolu*timom
      aux = elext(1:e%ndime)*tmp1
      do inode=1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + dot_product(e%cartd(:,inode),aux)
      end do
   end subroutine nsm_elmrhp_split2
   
   subroutine nsm_elmrhu_split2(e,dvolu,testf,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + (v, u_n/dt) 
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
         aux = e%shape(inode,e%igaus)*dvolu 
         elrhs(:,inode) = (elext(:)+eltemp(:))*aux + elrhs(:,inode)
      end do
   end subroutine nsm_elmrhu_split2
   
   subroutine nsm_elmbuq_split(e,timom,dvolu,denac,dtinv,vtemp,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for split OSS for block U,Q
      !        !     (div u, q) 
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
         !aux1 = vtemp(jnode)*denac*timom
         do inode=1,e%pnode
            elmat(1,inode,:,jnode) = (e%shape(inode,e%igaus)*e%cartd(:,jnode))*dvolu + elmat(1,inode,:,jnode)
         end do
      end do
   end subroutine nsm_elmbuq_split
   
   subroutine nsm_elmbpv_split(e,dvolu,testf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for split OSS for block P,V 
      !     - (div v, p) 
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
            elmat(:,inode,1,jnode) = -e%cartd(:,inode)*aux1  + elmat(:,inode,1,jnode)
         end do
      end do
   end subroutine nsm_elmbpv_split
  
   subroutine nsm_elmrhu_skipt(e,dvolu,testf,elext,eltemp,elrhs)
      !-----------------------------------------------------------------------
      !
      ! Only the external forces in the residual, not the temporal terms
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: testf(e%pnode),elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux1,tmp,aux2

      do inode=1,e%pnode
         aux1 = (e%shape(inode,e%igaus))*dvolu
         aux2 = (testf(inode))*dvolu
         elrhs(:,inode) = (elext(:)+eltemp(:))*aux1 + elext(:)*aux2 + elrhs(:,inode)
      end do
   end subroutine nsm_elmrhu_skipt
  
   subroutine nsm_elmrhu_oss_coriolis(e,dvol,denac,timom,gprep,AxesAngularVeloc,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for OSS
      !    tau*(Omega X v, proj(Res))
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: gprep(e%ndime+1)
      real(rp),    intent(in)    :: denac
      real(rp),    intent(in)    :: timom,dvol
      real(rp),    intent(in)    :: AxesAngularVeloc(e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux1,aux2

      aux1 = denac*timom*dvol
      do inode=1,e%pnode
         aux2 = aux1*e%shape(inode,e%igaus)
         elrhs(1,inode) = elrhs(1,inode) + aux2*(AxesAngularVeloc(3)*gprep(2)-AxesAngularVeloc(2)*gprep(3))
         elrhs(2,inode) = elrhs(2,inode) + aux2*(AxesAngularVeloc(1)*gprep(3)-AxesAngularVeloc(3)*gprep(1))
         elrhs(3,inode) = elrhs(3,inode) + aux2*(AxesAngularVeloc(2)*gprep(1)-AxesAngularVeloc(1)*gprep(2))
      end do
   end subroutine nsm_elmrhu_oss_coriolis
  
end module


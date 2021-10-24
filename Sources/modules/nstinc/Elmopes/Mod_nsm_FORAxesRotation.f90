 module Mod_nsm_FORAxesRotation
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersFORAxesRotation

   type, extends(PointerSetter) :: SPFORAxesRotation
contains
      procedure :: SpecificSet => SpecificSetFORAxesRotation
   end type
   type(SPFORAxesRotation) :: SetPointersFORAxesRotation

contains

   !----------------------------------------------------------------------------
   !Setting Pointers

   subroutine SpecificSetFORAxesRotation(d)
      implicit none
      class(SPFORAxesRotation) :: d
      integer(ip) :: ndime
      integer(ip) :: kfl_nonlinear

      if (a%kfl_FORAxesRotation == 1) then

         call a%Mesh%GetNdime(ndime)
         if(ndime==2)then
            call runend('nsi_Elmope: FOR Axis rotation not possible for 2d case')   
         end if
  
         !Matrices
         call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsCoriolis)

         !Non-linear elements
         call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
         if (kfl_nonlinear == 1) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsNonLinearCoriolis)
         endif
      endif
            
   end subroutine

   !---------------------------------------------------------------------------

   subroutine InGaussElmatsCoriolis
      implicit none
      
      !Compute contributions to elemental matrix : Block U,V
      call nsm_elmbuv_coriolis(e,timom,dvol,acden,a%FORAxesAngularVeloc,elmuv)
      
      !Compute contributions to elemental matrix : Block U,Q
      call nsm_elmbuq_coriolis(e,timom,dvol,acden,a%FORAxesAngularVeloc,elmuq)     

   end subroutine

   subroutine InGaussElmatsNonLinearCoriolis
      implicit none
      
      call nsm_elmbuv_lap_coriolis(e,timom,dvol,acden,acvis,a%FORAxesAngularVeloc,elmuv)
      
   end subroutine

   !---------------------------------------------------------------------------
   subroutine nsm_elmbuv_coriolis(e,timom,dvolu,denac,AxesAngularVeloc,elmuv)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   -(Omega X v, u) + tau1*(Omega X v, Omega X u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(in)    :: dvolu,denac,timom
      real(rp), intent(in)    :: AxesAngularVeloc(e%ndime)
      real(rp), intent(inout) :: elmuv(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip) :: jnode

      ! ([v]'·rho[Skew mat]',-[I]+ tau*rho[Skew mat]·[u])
      do jnode=1,e%pnode
         elmuv(1,1:e%pnode,1,jnode) = elmuv(1,1:e%pnode,1,jnode) + &
                 timom*acden*acden*(AxesAngularVeloc(3)**2+AxesAngularVeloc(2)**2)*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(1,1:e%pnode,2,jnode) = elmuv(1,1:e%pnode,2,jnode) + &
                  (-acden*AxesAngularVeloc(3)-timom*acden*acden*AxesAngularVeloc(1)*AxesAngularVeloc(2))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(1,1:e%pnode,3,jnode) = elmuv(1,1:e%pnode,3,jnode) + &
                  (acden*AxesAngularVeloc(2)-timom*acden*acden*AxesAngularVeloc(1)*AxesAngularVeloc(3))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(2,1:e%pnode,1,jnode) = elmuv(2,1:e%pnode,1,jnode) + &
                  (acden*AxesAngularVeloc(3)-timom*acden*acden*AxesAngularVeloc(1)*AxesAngularVeloc(2))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(2,1:e%pnode,2,jnode) = elmuv(2,1:e%pnode,2,jnode) + &
                  timom*acden*acden*(AxesAngularVeloc(1)**2+AxesAngularVeloc(3)**2)*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(2,1:e%pnode,3,jnode) = elmuv(2,1:e%pnode,3,jnode) + &
                  (-acden*AxesAngularVeloc(1)-timom*acden*acden*AxesAngularVeloc(2)*AxesAngularVeloc(3))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(3,1:e%pnode,1,jnode) = elmuv(3,1:e%pnode,1,jnode) + &
                  (-acden*AxesAngularVeloc(2)-timom*acden*acden*AxesAngularVeloc(1)*AxesAngularVeloc(3))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(3,1:e%pnode,2,jnode) = elmuv(3,1:e%pnode,2,jnode) + &
                  (acden*AxesAngularVeloc(1)-timom*acden*acden*AxesAngularVeloc(2)*AxesAngularVeloc(3))*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
         elmuv(3,1:e%pnode,3,jnode) = elmuv(3,1:e%pnode,3,jnode) + &
                  timom*acden*acden*(AxesAngularVeloc(2)**2+AxesAngularVeloc(1)**2)*e%shape(1:e%pnode,e%igaus)*e%shape(jnode,e%igaus)
      end do
   end subroutine nsm_elmbuv_coriolis
   
   subroutine nsm_elmbuq_coriolis(e,timom,dvolu,denac,AxesAngularVeloc,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     + tau1*(grad q, Omega X u) 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(in)    :: dvolu,timom,denac
      real(rp), intent(in)    :: AxesAngularVeloc(e%ndime)
      real(rp), intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip) :: jnode
      real(rp)    :: aux1

      ! ([grad q],tau*rho[Skew mat]·[u])
      aux1 = denac*timom*dvolu
      do jnode=1,e%pnode
         elmat(1,1:e%pnode,1,jnode) = elmat(1,1:e%pnode,1,jnode) + &
                (e%cartd(2,1:e%pnode)*AxesAngularVeloc(3)-e%cartd(3,1:e%pnode)*AxesAngularVeloc(2))*e%shape(jnode,e%igaus)*aux1
         elmat(1,1:e%pnode,2,jnode) = elmat(1,1:e%pnode,2,jnode) + &
                (e%cartd(3,1:e%pnode)*AxesAngularVeloc(1)-e%cartd(1,1:e%pnode)*AxesAngularVeloc(2))*e%shape(jnode,e%igaus)*aux1
         elmat(1,1:e%pnode,3,jnode) = elmat(1,1:e%pnode,3,jnode) + &
                (e%cartd(1,1:e%pnode)*AxesAngularVeloc(2)-e%cartd(2,1:e%pnode)*AxesAngularVeloc(1))*e%shape(jnode,e%igaus)*aux1
      end do

   end subroutine nsm_elmbuq_coriolis

   subroutine nsm_elmbuv_lap_coriolis(e,timom,dvolu,denac,visac,AxesAngularVeloc,elmuv)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !    + tau1*(visco*lapla(v), Omega X u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(in)    :: dvolu,timom,denac
      real(rp), intent(in)    :: visac,AxesAngularVeloc(e%ndime)
      real(rp), intent(inout) :: elmuv(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip) :: jnode,inode
      real(rp)    :: aux1

      aux1 = denac*visac*timom*dvolu
      forall(inode=1:e%pnode,jnode=1:e%pnode)
             elmuv(1,inode,2,jnode) = elmuv(1,inode,2,jnode) &
                          - aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(3)*e%shape(jnode,e%igaus)
             elmuv(1,inode,3,jnode) = elmuv(1,inode,3,jnode) &
                          + aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(2)*e%shape(jnode,e%igaus)
             elmuv(2,inode,1,jnode) = elmuv(2,inode,1,jnode) &
                          + aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(3)*e%shape(jnode,e%igaus)
             elmuv(2,inode,3,jnode) = elmuv(2,inode,3,jnode) &
                          - aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(1)*e%shape(jnode,e%igaus)
             elmuv(3,inode,1,jnode) = elmuv(3,inode,1,jnode) &
                          - aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(2)*e%shape(jnode,e%igaus)
             elmuv(3,inode,2,jnode) = elmuv(3,inode,2,jnode) &
                             + aux1*sum(e%hessi(1:e%ndime,inode))*AxesAngularVeloc(1)*e%shape(jnode,e%igaus)
      end forall
   end subroutine nsm_elmbuv_lap_coriolis

end module

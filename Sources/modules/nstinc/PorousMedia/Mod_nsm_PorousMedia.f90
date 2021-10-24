module Mod_nsm_PorousMedia
   use typre
   use Mod_Mesh
   use Mod_Element
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersPorousMedia

   type, extends(PointerSetter) :: SPPorousMedia
contains
      procedure :: SpecificSet => SpecificSetPorousMedia
   end type
   type(SPPorousMedia) :: SetPointersPorousMedia

 contains
 
   subroutine SpecificSetPorousMedia(d)
      implicit none
      class(SPPorousMedia) :: d
 
      if (a%kfl_PorousMedia == 1) then
         call ConcatenateProcedures(ProcHook_InGaussElmats,PorousBox)
      endif
   end subroutine 
   
!----------------------------------------------------------------------------
   subroutine PorousBox
      implicit none
      real(rp) :: gpcod(3)
      real(rp) ::  midPoint
      
      midPoint = a%PME%minX + a%PME%Dis       !Divided to 2 box
   
      call e%interpg(e%ndime,e%elcod,gpcod)
      if (a%PME%minX <= gpcod(1) .and. gpcod(1) < midPoint  .and. a%PME%minY <= gpcod(2) .and. gpcod(2) <= a%PME%maxY ) then      !Box 1
         call elm_porous(e,dvol,a%PME%perm1,wrmat1)
      elseif(midPoint <= gpcod(1) .and. gpcod(1)<= a%PME%maxX .and. a%PME%minY <= gpcod(2) .and. gpcod(2) <= a%PME%maxY ) then    !Box 2
         call elm_porous(e,dvol,a%PME%perm2,wrmat1)
      endif
   end subroutine
   
   subroutine elm_porous(e,dvol,perm,myelmat)
         implicit none
         type(FiniteElement)     :: e
         real(rp), intent(in)    :: dvol
         real(rp), intent(inout) :: myelmat(e%mnode,e%mnode)
         integer(ip) :: inode,jnode
         real(rp)    :: tmp, perm,invperm,acvis

         invperm = 1.0_rp/perm
         do inode=1,e%pnode
            do jnode=1,e%pnode
               myelmat(inode,jnode) = e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*invperm*dvol+ myelmat(inode,jnode)
            end do
         end do
   end subroutine  
   
 end module
 

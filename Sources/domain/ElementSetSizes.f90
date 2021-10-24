subroutine ElementSetSizes(a,e)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   class(FiniteElement) :: e

   e%ndime = a%ndime
   e%mnode = a%mnode
   e%mnodb = a%mnodb
   e%mface = a%mface
   e%ntens = a%ntens
   e%mlapl = a%mlapl
   
   !Initializations for e%nelty = 1
   e%ielty0 = 1
   e%ielty  = 1
   e%iblty  = 1
   e%iblty0 = 1
   e%pnode = a%nnode(1)
   e%pnodb = a%nnodb(1)
   
   e%kfl_perio = .false.
   if (a%nslave > 0 .and. a%kfl_perio == 1) e%kfl_perio = .true.
   
      if(e%ndime==2) then
         select case(e%pnode)
            case(3,4)
               e%npol=1
            case(6,9)
               e%npol=2
            case(10,16)
               e%npol=3
            case(15,25)
               e%npol=4
         end select
      endif
      if(e%ndime==3) then
         select case(e%pnode)
            case(4,8)
               e%npol=1
            case(10,27)
               e%npol=2
            case(20,64)
               e%npol=3
            case(35,125)
               e%npol=4
         end select
      endif
      e%npol4 = e%npol*e%npol*e%npol*e%npol
     
end subroutine


subroutine ElementLoad(a,ielem,e)
   use typre
   use Mod_Element
   use Mod_ElementWithDataStructures
   use Mod_Mesh
   implicit none
   class(FemMesh), target :: a
   integer(ip)            :: ielem,inode,ipoin
   class(FiniteElement)    :: e
   
   integer(ip) :: ispos

   
!   e%pface = a%nface(e%ielty)
   !Nothing is computed
   e%kfl_IsAlreadyComputedLinearDerivatives = .false.
   
   ispos = a%nnode(1)*(ielem-1)  
   
   !I have done something, reset everything
   if (e%ielty0 == -1) then
      if (a%nelty == 1) then
         e%pnode = a%nnode(1)
         call a%ElementSetPointers(e)
      endif
   endif
      
   !------------------------------------------------------------------------------
   if (a%nelty > 1) then  !Otherwise it has been done in the allocation
      e%pnode = a%pnods(ielem+1)-a%pnods(ielem)
      e%ielty = a%ltype(e%pnode)
      ispos = a%pnods(ielem)-1
      if (e%ielty /= e%ielty0) then
         call a%ElementSetPointers(e)
         e%ielty0 = e%ielty
      endif   
   endif
  !------------------------------------------------------------------------------
  
   e%lnods(1:e%pnode) = a%lnods(ispos+1:ispos+e%pnode)
   call e%gather(e%ndime,e%elcod,a%coord)

   !For periodic boundary conditions, modify so that now we work with the Masters
   if (e%kfl_perio .eqv. .true.) e%lnods(1:e%pnode) = a%MasterSlave(a%lnods(ispos+1:ispos+e%pnode))

   !For mesh movement
   if (a%kfl_alemov==1) then
      e%elcod(:,1:e%pnode) = e%elcod(:,1:e%pnode) + a%displ(:,e%lnods(1:e%pnode),1)
   end if
   
   !If using ElementDataStructures
   if (a%kfl_UseElementDataStructures) then
      select type (e)
         type is (FiniteElementWDS)
         if (e%kfl_ClosedRule == 0) then
            e%EDS => a%ElementDataStructure(ielem)
         elseif (e%kfl_ClosedRule == 1) then
            e%EDS => a%ElementDataStructureClosedRule(ielem)
         endif
      end select
   endif   
   

end subroutine

subroutine FaceLoad(a,iface,ielem,e)
   use typre
   use Mod_Mesh
   use Mod_Element
   use Mod_ElementWithDataStructures   
   implicit none
   class(FemMesh), target       :: a
   class(FiniteElement), target :: e
   integer(ip)                  :: ielem,iface,inode,inodb,icount

   
   call a%ElementLoad(ielem,e)
   
   if (a%nelty > 1) then 
      e%iblty = a%nnodb(e%ielty)
      if (e%iblty /= e%iblty0) then
         call a%BoundarySetPointers(e)
         e%iblty0 = e%iblty
      endif   
   endif

   e%lnodb(1:e%pnodb) = e%lnods(a%cfael(1:e%pnodb,iface,e%ielty))
   e%bocod(:,1:e%pnodb) = e%elcod(:,a%cfael(1:e%pnodb,iface,e%ielty))
   e%lboel => a%cfael(1:e%pnodb,iface,e%ielty)

   if (a%kfl_UseElementDataStructures) then
       select type (e)
          type is (FiniteElementWDS)
          if (e%kfl_ClosedRule == 0) then
             e%EDS => a%ElementDataStructure(ielem)
          elseif (e%kfl_ClosedRule == 1) then
             e%EDS => a%ElementDataStructureClosedRule(ielem)
          endif
       end select
   endif   
   
   
end subroutine


subroutine BoundaryLoad(a,iboun,e)
   use typre
   use Mod_Element
   use Mod_ElementWithDataStructures   
   use Mod_Mesh
   implicit none
   class(FemMesh), target :: a
   integer(ip)            :: iboun
   class(FiniteElement)    :: e
   
   integer(ip) :: ispos2,ielem
   
   ! Element properties and dimensions
   ispos2 = (a%nnodb(1)+1)*(iboun-1)                               
   
   !------------------------------------------------------------------------
   if (a%nelty > 1) then 
      ispos2 = a%pboel(iboun)-1
      e%iblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      if (e%iblty /= e%iblty0) then
         call a%BoundarySetPointers(e)
         e%iblty0 = e%iblty
      endif   
   endif
   !------------------------------------------------------------------------------
   e%lboel => a%lboel(ispos2+1:ispos2+e%pnodb+1)
   ielem = e%lboel(e%pnodb+1) 
   call a%ElementLoad(ielem,e)
   e%lnodb(1:e%pnodb) = e%lnods(e%lboel(1:e%pnodb))
   e%bocod(:,1:e%pnodb) = e%elcod(:,e%lboel(1:e%pnodb))
   
   !If using ElementDataStructures
   if (a%kfl_UseElementDataStructures) then
      select type (e)
         type is (FiniteElementWDS)
         if (e%kfl_ClosedRule == 0) then
            e%BDS => a%BoundaryDataStructure(iboun)
         elseif (e%kfl_ClosedRule == 1) then
            e%BDS => a%BoundaryDataStructureClosedRule(iboun)
         endif
      end select
   endif   
   
end subroutine


subroutine ElementSetPointers(a,e)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh),target :: a
   class(FiniteElement)   :: e
   
   if ((e%ndime ==2 .and. e%pnode == 3) .or. (e%ndime==3.and.e%pnode==4)) then
      e%linea = 1
   else
      e%linea = 0
   endif
   
   if (e%kfl_ClosedRule == 0) then
      e%pgaus = a%ngaus(e%ielty)
      e%shape => a%shape(:,:,e%ielty)
      e%deriv => a%deriv(:,:,:,e%ielty)
      e%weigp => a%weigp(:,e%ielty)
      if (e%linea == 0) e%heslo => a%heslo(:,:,:,e%ielty)
   else
      e%pgaus = a%nnode(e%ielty)
      e%shape => a%shape_clo(:,:,e%ielty)
      e%deriv => a%deriv_clo(:,:,:,e%ielty)
      e%weigp => a%weigp_clo(:,e%ielty)
      if (e%linea == 0) e%heslo => a%heslo_clo(:,:,:,e%ielty)
   endif
   
   e%shacg => a%shacg(:,e%ielty)
   e%dercg => a%dercg(:,:,e%ielty)
   e%weicg = a%weicg(e%ielty)
   
   e%pnode = a%nnode(e%ielty)
   e%hnatu = a%hnatu(e%ielty)
   e%shaga => a%shaga(:,:,e%ielty)
   e%llapl = a%llapl(e%ielty)
  
   !Boundaries of all elements 
   e%pface = a%nface(e%ielty)

   e%shapb => a%shapb(:,:,e%ielty)
   e%derib => a%derib(:,:,:,e%ielty)
   e%weigb => a%weigb(:,e%ielty)
   e%shagb => a%shagb(:,:,e%ielty)
   
   e%pnodb = a%nnodb(e%ielty)
   e%pgaub = a%ngaub(e%ielty)
 
end subroutine


subroutine BoundarySetPointers(a,e)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh),target :: a
   class(FiniteElement)   :: e
   
   e%shapb => a%shapb(:,:,e%iblty)
   e%derib => a%derib(:,:,:,e%iblty)
   e%weigb => a%weigb(:,e%iblty)
   e%shagb => a%shagb(:,:,e%iblty)
   
   e%pnodb = a%nnodb(e%iblty)
   e%pgaub = a%ngaub(e%iblty)

end subroutine


subroutine ElementAlloc(a,e,Memor,rule,outstr)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_ElementWithDataStructures
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   class(FiniteElement), pointer :: e
   type(MemoryMan)     :: Memor
   character(*)        :: outstr
   character(6)        :: rule
   
   if (a%kfl_UseElementDataStructures .eqv. .false.) then
      allocate(FiniteElement::e)
      
   else
      allocate(FiniteElementWDS::e)
   endif
      
      
   !A closed integration rule needs to be enforced for 
   !the smoothing process (both LHS (vmass) and RHS)
   if (rule .eq. 'ForceC') then
      call e%ForceClosedRule
   elseif (rule .eq. 'Defaul') then
      call e%DefaultRule
   else
      call runend('ElementAlloc: wrong rule')
   endif   
   
   call a%ElementSetSizes(e)
   call e%alloc(Memor,outstr)
   call a%ElementSetPointers(e)
   call a%BoundarySetPointers(e)
   call a%GetNpoinLocal(e%npoinLocal)
   
end subroutine

subroutine ElementDealloc(a,e,Memor,rule,outstr)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   class(FiniteElement), pointer :: e 
   type(MemoryMan)     :: Memor
   character(*)       :: outstr
   character(6)        :: rule

   call e%dealloc(Memor,outstr)
   deallocate(e)

end subroutine

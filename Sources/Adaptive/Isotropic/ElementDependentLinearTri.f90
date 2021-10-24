module Mod_AdaptiveTriangle
   use Mod_AdaptiveInterpolationList
   implicit none
   


contains

   subroutine LinearTriangle_GetElementVariation(elcod,ivariation)
      use typre
      implicit none
      real(rp) :: elcod(:,:)
      integer(ip) :: ivariation
      
      ivariation = 1
      
   end subroutine

   subroutine LinearTriangle_ElementGetDimensions(pnodb,nchild,nface,nchildface)
      use typre
      implicit none
      integer(ip), intent(out) :: pnodb,nchild,nface,nchildface
   
      pnodb = 2
      nchild = 4
      nface = 3
      nchildface = 2
   end subroutine
   
   subroutine LinearTriangle_Add_node_to_child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
      !This subroutine adds a new node to a child element when this is built
      use typre
      implicit none
      integer(ip), intent(in) :: nchild,inode,jnode,ipoin
      integer(ip), intent(inout) :: pnods(*),lnods(*),PointerToChildren(*)
      
      integer(ip) :: knode
      integer(ip) :: ielem,ispos

      
      if (inode == jnode) then
         knode = inode
      elseif ((inode == 1 .and. jnode == 2) .or. (inode == 2 .and. jnode == 1)) then
         knode = 4
      elseif ((inode == 1 .and. jnode == 3) .or. (inode == 3 .and. jnode == 1)) then
         knode = 5
      elseif ((inode == 3 .and. jnode == 2) .or. (inode == 2 .and. jnode == 3)) then
         knode = 6
      endif    
      
      select case (knode)
      case(1:3)
         ielem = PointerToChildren(knode)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      case(4)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
      case(5)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      case(6)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
      end select
      
   end subroutine
   
   subroutine LinearTriangle_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
      use typre
      implicit none
      integer(ip), intent(in) :: pnods(*),lnods(*),nchild,PointerToChildren(nchild)
      integer(ip), intent(out) :: nedge,ledger(3,*)
      
      integer(ip) :: ichild1,ispos1,ichild2,ispos2,ichild3,ispos3

      nedge = 3
      
      ichild1 = PointerToChildren(1)
      ispos1 = pnods(ichild1)-1
      
      ichild2 = PointerToChildren(2)
      ispos2 = pnods(ichild2)-1
      
      ichild3 = PointerToChildren(3)
      ispos3 = pnods(ichild3)-1
      
      !first edge
      ledger(1,1) = lnods(ispos1+1)
      ledger(2,1) = lnods(ispos2+1)
      ledger(3,1) = lnods(ispos1+2)
      
      !second edge
      ledger(1,2) = lnods(ispos2+1)
      ledger(2,2) = lnods(ispos3+1)
      ledger(3,2) = lnods(ispos2+2)
      
      !third edge
      ledger(1,3) = lnods(ispos3+1)
      ledger(2,3) = lnods(ispos1+1)
      ledger(3,3) = lnods(ispos3+2)

   end subroutine  

   subroutine LinearTriangle_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
      use typre
      implicit none

      integer(ip) :: iface,nface,pnodb
      integer(ip) :: facnod(pnodb,nface),lnods(*)
      integer(ip) :: sort

      if (iface == 0) then
         facnod(1,1) = lnods(1)
         facnod(2,1) = lnods(2)
         facnod(1,2) = lnods(3)
         facnod(2,2) = lnods(1)
         facnod(1,3) = lnods(2)
         facnod(2,3) = lnods(3)
         if (sort == 1) then
               call sortb(2,facnod(1,2))
               call sortb(2,facnod(1,3))
         endif
      elseif (iface == 1) then
         facnod(1,1) = lnods(1)
         facnod(2,1) = lnods(2)
      elseif (iface == 2) then
         facnod(1,1) = lnods(3)
         facnod(2,1) = lnods(1)
      elseif (iface == 3) then
         facnod(1,1) = lnods(2)
         facnod(2,1) = lnods(3)
      endif   
      if (sort == 1) call sortb(2,facnod(1,1))     

   
   end subroutine

   subroutine LinearTriangle_up_faces(nface,ison,jface)
      !this subroutine gives the new neighbour face in the parent when an element is gone
      use typre
      implicit none
      integer(ip):: ielty,nface,jface,ison
      select case (ison)
         case (1)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 2) then
               jface = 2
            endif
         case (2)
            if (jface == 1) then 
               jface = 3 
            elseif (jface == 2) then
               jface = 1
            endif
         case (3)
            if (jface == 2) then
               jface = 3    
            elseif (jface == 1) then
               jface = 2 
            endif   
         !case (4) !this one will never happen, interior faces of the refined element 
         end select
      
   end subroutine

   subroutine LinearTriangle_down_faces(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1
         facchild(2,1) = 1 !first son, first face
         facchild(1,2) = 2
         facchild(2,2) = 2 !second son, second face
      case (2)
         facchild(1,1) = 1
         facchild(2,1) = 2 !first son, second face
         facchild(1,2) = 3
         facchild(2,2) = 1 !third son, first face
      case (3)
         facchild(1,1) = 2
         facchild(2,1) = 1 !second son, first face
         facchild(1,2) = 3
         facchild(2,2) = 2 !third son, second face
      end select

   end subroutine

   subroutine LinearTriangle_nodes_to_face(pnodb,lnodb,iface)
      use typre
      implicit none
      integer(ip) :: ndime,pnode,pnodb,lnodb(pnodb),iface

      if ((lnodb(1) == 1 .and. lnodb(2) == 2) .or. (lnodb(1) == 2 .and. lnodb(2) == 1)) then
         iface = 1
      elseif ((lnodb(1) == 1 .and. lnodb(2) == 3) .or. (lnodb(1) == 3 .and. lnodb(2) == 1)) then
         iface = 2
      elseif ((lnodb(1) == 2 .and. lnodb(2) == 3) .or. (lnodb(1) == 3 .and. lnodb(2) == 2)) then
         iface = 3
      endif

   end subroutine

   subroutine LinearTriangle_ConnectInteriorFaces(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      lfacs(1,jfpos(1)+3) = son(4)
      lfacs(2,jfpos(1)+3) = 2
      
      lfacs(1,jfpos(4)+2) = son(1)
      lfacs(2,jfpos(4)+2) = 3
      
      lfacs(1,jfpos(2)+3) = son(4)
      lfacs(2,jfpos(2)+3) = 1
      
      lfacs(1,jfpos(4)+1) = son(2)
      lfacs(2,jfpos(4)+1) = 3
      
      lfacs(1,jfpos(3)+3) = son(4)
      lfacs(2,jfpos(3)+3) = 3   
      
      lfacs(1,jfpos(4)+3) = son(3)
      lfacs(2,jfpos(4)+3) = 3
      
   end subroutine


   subroutine LinearTriangle_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      !there are 3 children with faces shared with the parent element
      inchild = 3
      !These elements are child 1, 2 and 3
      jjchild(1) = 1
      jjchild(2) = 2
      jjchild(3) = 3
      !Each child has 2 external faces
      jjnface(1:3)  = 2
      !The external faces are the 2 first ones for all the children
      jjface(1,1:3) = 1
      jjface(2,1:3) = 2
   end subroutine
   
   subroutine LinearTriangle_GetInterpolationList(ivariation,InterpList)
      use typre
      implicit none
      integer(ip) :: ivariation
      type(InterpolationList) :: InterpList
      
      InterpList%n = 0
   end subroutine
   
   subroutine LinearTriangle_GetEdgeNodesForInode(inode,nnode,nodelist)
      use typre
      implicit none
      integer(ip) :: inode, nnode,nodelist(2)
      
      nnode = 2
      
      if (inode == 1) then
         nodelist(1) = 2
         nodelist(2) = 3
      elseif (inode == 2) then
         nodelist(1) = 1
         nodelist(2) = 3
      elseif (inode == 3) then
         nodelist(1) = 1
         nodelist(2) = 2
      endif
   end subroutine
   
   subroutine LinearTriangle_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
      use typre
      implicit none
      integer(ip) :: nInteriorNodes
      
      nInteriorNodes = 0
      
   end subroutine
   
   subroutine LinearTriangle_GetInteriorPointFromParent(inode,jchild,jnode)
      use typre
      implicit none
      integer(ip) :: inode,jchild,jnode
      
      
      call runend('Wrong number of interior nodes')
   end subroutine   
   
   subroutine LinearTriangle_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,iadd,newpo,PointerToChildren(*),pnods(*),lnods(*)
      
      call runend('Should not add interior node in triangle')
   end subroutine    
   
   subroutine LinearTriangle_Up_Edges(ison,childedge,parentedge)
      use typre
      implicit none
      integer(ip) :: ielty,ison, childedge,parentedge
      
      call runend('Linear Triangles do not have edges')
   end subroutine
   
   subroutine LinearTriangle_GetEdges(iedge,lnodedge)
      use typre
      implicit none
      integer(ip) :: ielty,iedge,lnodedge(*)
      
      call runend('Linear Triangles do not have edges')
   end subroutine
   
   subroutine LinearTriangle_GetNedge(nedge)
      use typre
      implicit none
      integer(ip) :: nedge
      
      nedge = 0
   end subroutine
   
   subroutine  LinearTriangle_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
      use typre
      implicit none
      integer(ip) :: nInTheFaceNodes
      
      nInTheFaceNodes = 0
   end subroutine
   
   subroutine LinearTriangle_GetFaceForInTheFaceNode(inode,iface)
      use typre
      implicit none
      integer(ip) :: inode,iface
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine   
   
   subroutine LinearTriangle_GetInTheFacePointFromParent(iface,kchild,knode)
      use typre
      implicit none
      integer(ip) :: iface,kchild,knode
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   subroutine LinearTriangle_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,inode,kpoin,PointerToChildren(*),pnods(*),lnods(*)
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   
    subroutine LinearTriangle_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,3),NewNodalArray(ndofn,3)
      integer(ip) :: ielty,ivariation,ichild
   
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
      elseif (ichild == 2) then
         NewNodalArray(:,1) = OldNodalArray(:,2)
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
      elseif (ichild == 3) then
         NewNodalArray(:,1) = OldNodalArray(:,3)
         NewNodalArray(:,2) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
      elseif (ichild == 4) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
      endif
   
   
   end subroutine
   
   subroutine LinearTriangle_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,3),NewNodalArray(ndofn,3)
      integer(ip) :: ivariation,ichild
      
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
      elseif (ichild == 2) then
         NewNodalArray(:,2) = OldNodalArray(:,1)
      elseif (ichild == 3) then
         NewNodalArray(:,3) = OldNodalArray(:,1)
      
      endif
   
   end subroutine
   
end module
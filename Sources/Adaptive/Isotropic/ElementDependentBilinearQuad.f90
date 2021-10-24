module Mod_AdaptiveQuads
   use Mod_AdaptiveInterpolationList
   implicit none
   


contains

   subroutine BilinearQuad_GetElementVariation(elcod,ivariation)
      use typre
      implicit none
      real(rp) :: elcod(:,:)
      integer(ip) :: ivariation
      
      ivariation = 1
      
   end subroutine

   subroutine BilinearQuad_ElementGetDimensions(pnodb,nchild,nface,nchildface)
      use typre
      implicit none
      integer(ip), intent(out) :: pnodb,nchild,nface,nchildface
   
      pnodb = 2
      nchild = 4
      nface = 4
      nchildface = 2
   end subroutine
   
   subroutine BilinearQuad_Add_node_to_child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
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
         knode = 5
      elseif ((inode == 2 .and. jnode == 3) .or. (inode == 3 .and. jnode == 2)) then
         knode = 6
      elseif ((inode == 3 .and. jnode == 4) .or. (inode == 4 .and. jnode == 3)) then
         knode = 7
      elseif ((inode == 4 .and. jnode == 1) .or. (inode == 1 .and. jnode == 4)) then
         knode = 8   
      endif    
      
      select case (knode)
      case(1:4)
         ielem = PointerToChildren(knode)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      
      case(5)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
      
      case(6)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
      case(7)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin  
      
      case(8)
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
      end select
      
   end subroutine
   
   subroutine BilinearQuad_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
      use typre
      implicit none
      integer(ip), intent(in) :: pnods(*),lnods(*),nchild,PointerToChildren(nchild)
      integer(ip), intent(out) :: nedge,ledger(3,*)
      
      integer(ip) :: ichild(4),ispos(4)

      nedge = 4
      
      ichild(1) = PointerToChildren(1)
      ispos(1) = pnods(ichild(1))-1
      
      ichild(2) = PointerToChildren(2)
      ispos(2) = pnods(ichild(2))-1
      
      ichild(3) = PointerToChildren(3)
      ispos(3) = pnods(ichild(3))-1
      
      ichild(4) = PointerToChildren(4)
      ispos(4) = pnods(ichild(4))-1
      
      !first edge
      ledger(1,1) = lnods(ispos(1)+1)
      ledger(2,1) = lnods(ispos(2)+1)
      ledger(3,1) = lnods(ispos(1)+2)
      
      !second edge
      ledger(1,2) = lnods(ispos(2)+1)
      ledger(2,2) = lnods(ispos(3)+1)
      ledger(3,2) = lnods(ispos(2)+2)
      
      !third edge
      ledger(1,3) = lnods(ispos(3)+1)
      ledger(2,3) = lnods(ispos(4)+1)
      ledger(3,3) = lnods(ispos(3)+2)
      
      !fourth edge
      ledger(1,4) = lnods(ispos(4)+1)
      ledger(2,4) = lnods(ispos(1)+1)
      ledger(3,4) = lnods(ispos(4)+2)

   end subroutine  

   subroutine BilinearQuad_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
      use typre
      implicit none

      integer(ip) :: iface,nface,pnodb
      integer(ip) :: facnod(pnodb,nface),lnods(*)
      integer(ip) :: sort

      if (iface == 0) then
         facnod(1,1) = lnods(1)
         facnod(2,1) = lnods(2)
         facnod(1,2) = lnods(2)
         facnod(2,2) = lnods(3)
         facnod(1,3) = lnods(3)
         facnod(2,3) = lnods(4)
         facnod(1,4) = lnods(4)
         facnod(2,4) = lnods(1)
         if (sort == 1) then
               call sortb(2,facnod(1,2))
               call sortb(2,facnod(1,3))
               call sortb(2,facnod(1,4))
         endif
      elseif (iface == 1) then
         facnod(1,1) = lnods(1)
         facnod(2,1) = lnods(2)
      elseif (iface == 2) then
         facnod(1,1) = lnods(2)
         facnod(2,1) = lnods(3)
      elseif (iface == 3) then
         facnod(1,1) = lnods(3)
         facnod(2,1) = lnods(4)
      elseif (iface == 4) then
         facnod(1,1) = lnods(4)
         facnod(2,1) = lnods(1)   
      endif   
      if (sort == 1) call sortb(2,facnod(1,1))     
   end subroutine

   subroutine BilinearQuad_up_faces(nface,ison,jface)
      !this subroutine gives the new neighbour face in the parent when an element is gone
      use typre
      implicit none
      integer(ip):: ielty,nface,jface,ison
      select case (ison)
         case (1)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 4) then 
               jface = 4
            endif
         case (2)
            if (jface == 1) then
               jface = 2 
            elseif (jface == 4) then
               jface = 1
            endif
         case (3)
            if (jface == 1) then
               jface = 3    
            elseif (jface == 4) then
               jface = 2 
            endif
         case (4)
            if (jface == 1) then
               jface = 4    
            elseif (jface == 4) then
               jface = 3 
            endif
         end select
      
   end subroutine

   subroutine BilinearQuad_down_faces(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1
         facchild(2,1) = 1 !first son, first face
         facchild(1,2) = 2
         facchild(2,2) = 4 !second son, fourth face
      case (2)
         facchild(1,1) = 2
         facchild(2,1) = 1 !first son, second face
         facchild(1,2) = 3
         facchild(2,2) = 4 !third son, first face
      case (3)
         facchild(1,1) = 3
         facchild(2,1) = 1 !second son, first face
         facchild(1,2) = 4
         facchild(2,2) = 4 !third son, second face
      case (4)
         facchild(1,1) = 4
         facchild(2,1) = 1 !second son, first face
         facchild(1,2) = 1
         facchild(2,2) = 4 !third son, second face   
         
      end select

   end subroutine

   subroutine BilinearQuad_nodes_to_face(pnodb,lnodb,iface)
      use typre
      implicit none
      integer(ip) :: ndime,pnode,pnodb,lnodb(pnodb),iface

      if ((lnodb(1) == 1 .and. lnodb(2) == 2) .or. (lnodb(1) == 2 .and. lnodb(2) == 1)) then
         iface = 1
      elseif ((lnodb(1) == 2 .and. lnodb(2) == 3) .or. (lnodb(1) == 3 .and. lnodb(2) == 2)) then
         iface = 2
      elseif ((lnodb(1) == 3 .and. lnodb(2) == 4) .or. (lnodb(1) == 4 .and. lnodb(2) == 3)) then
         iface = 3
      elseif ((lnodb(1) == 4 .and. lnodb(2) == 1) .or. (lnodb(1) == 1 .and. lnodb(2) == 4)) then
         iface = 4   
      endif

   end subroutine

   subroutine BilinearQuad_ConnectInteriorFaces(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      lfacs(1,jfpos(1)+2) = son(2)
      lfacs(2,jfpos(1)+2) = 3
      
      lfacs(1,jfpos(2)+3) = son(1)
      lfacs(2,jfpos(2)+3) = 2
      
      !2
      lfacs(1,jfpos(2)+2) = son(3)
      lfacs(2,jfpos(2)+2) = 3
      
      lfacs(1,jfpos(3)+3) = son(2)
      lfacs(2,jfpos(3)+3) = 2
      
      !3
      lfacs(1,jfpos(3)+2) = son(4)
      lfacs(2,jfpos(3)+2) = 3
      
      lfacs(1,jfpos(4)+3) = son(3)
      lfacs(2,jfpos(4)+3) = 2
      
      !4
      lfacs(1,jfpos(4)+2) = son(1)
      lfacs(2,jfpos(4)+2) = 3
      
      lfacs(1,jfpos(1)+3) = son(4)
      lfacs(2,jfpos(1)+3) = 2
      
   end subroutine


   subroutine BilinearQuad_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      !there are 3 children with faces shared with the parent element
      inchild = 4
      !These elements are child 1, 2 and 3
      jjchild(1) = 1
      jjchild(2) = 2
      jjchild(3) = 3
      jjchild(4) = 4
      !Each child has 2 external faces
      jjnface(1:4)  = 2
      !The external faces are face 1 and 4 for all the children
      jjface(1,1:4) = 1
      jjface(2,1:4) = 4
   end subroutine
   
   subroutine BilinearQuad_GetInterpolationList(ivariation,InterpList)
      use typre
      implicit none
      integer(ip) :: ivariation
      type(InterpolationList) :: InterpList
      
      InterpList%n = 1
      InterpList%childElementList(1) = 1
      InterpList%childObjectiveNodeList(1) = 3
      InterpList%nInterpolatingNodes(1) = 4
      InterpList%InterpolatingNodesList(1:4,1) = (/1_ip, 2_ip, 3_ip, 4_ip /)
      InterpList%InterpolatingCoefficientsList(1:4,1) = 0.25_rp
   end subroutine
   
   subroutine BilinearQuad_GetEdgeNodesForInode(inode,nnode,nodelist)
      use typre
      implicit none
      integer(ip) :: inode, nnode,nodelist(2)
      
      nnode = 2
      
      if (inode == 1) then
         nodelist(1) = 2
         nodelist(2) = 4
      elseif (inode == 2) then
         nodelist(1) = 1
         nodelist(2) = 3
      elseif (inode == 3) then
         nodelist(1) = 2
         nodelist(2) = 4
      elseif (inode == 4) then
         nodelist(1) = 3
         nodelist(2) = 1   
      endif
   end subroutine
   
   subroutine BilinearQuad_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
      use typre
      implicit none
      integer(ip) :: nInteriorNodes
      
      nInteriorNodes = 1
      
   end subroutine
   
   subroutine BilinearQuad_GetInteriorPointFromParent(inode,jchild,jnode)
      use typre
      implicit none
      integer(ip) :: inode,jchild,jnode
      
      if (inode == 1) then
         jchild = 1
         jnode = 3
      else
         call runend('Wrong number of interior nodes')
      endif
   end subroutine
   
   subroutine BilinearQuad_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,iadd,newpo,PointerToChildren(*),pnods(*),lnods(*)
      integer(ip) :: ispos(4),ichild
      
      !Only one interior node needs to be added
      if (iadd == 1) then
         do ichild = 1,4
            ispos(ichild) = pnods(PointerToChildren(ichild))-1
            lnods(ispos(ichild)+3) = newpo
         enddo
      endif
   end subroutine    
   
   subroutine BilinearQuad_Up_Edges(ison,childedge,parentedge)
      use typre
      implicit none
      integer(ip) :: ielty,ison, childedge,parentedge
      
      call runend('BilinearQuads do not have edges')
   end subroutine
   
   subroutine BilinearQuad_GetEdges(iedge,lnodedge)
      use typre
      implicit none
      integer(ip) :: ielty,iedge,lnodedge(*)
      
      call runend('BilinearQuads do not have edges')
   end subroutine
   
   subroutine BilinearQuad_GetNedge(nedge)
      use typre
      implicit none
      integer(ip) :: nedge
      
      nedge = 0
   end subroutine
   
   subroutine  BilinearQuad_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
      use typre
      implicit none
      integer(ip) :: nInTheFaceNodes
      
      nInTheFaceNodes = 0
   end subroutine
   
   subroutine BilinearQuad_GetFaceForInTheFaceNode(inode,iface)
      use typre
      implicit none
      integer(ip) :: inode,iface
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine  
   
   subroutine BilinearQuad_GetInTheFacePointFromParent(iface,kchild,knode)
      use typre
      implicit none
      integer(ip) :: iface,kchild,knode
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   subroutine BilinearQuad_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,inode,kpoin,PointerToChildren(*),pnods(*),lnods(*)
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   subroutine BilinearQuad_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,4),NewNodalArray(ndofn,4)
      integer(ip) :: ielty,ivariation,ichild
   
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
      elseif (ichild == 2) then
         NewNodalArray(:,1) =  OldNodalArray(:,2)
         NewNodalArray(:,2) =  (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) =  (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,4) =  (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
      elseif (ichild == 3) then
         NewNodalArray(:,1) =  OldNodalArray(:,3)
         NewNodalArray(:,2) =  (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
         NewNodalArray(:,3) =  (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,4) =  (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
      elseif (ichild == 4) then
         NewNodalArray(:,1) =  OldNodalArray(:,4)
         NewNodalArray(:,2) =  (OldNodalArray(:,4)+OldNodalArray(:,1))*0.5
         NewNodalArray(:,3) =  (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,4) =  (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
      endif
   
   end subroutine
   
   subroutine BilinearQuad_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,4),NewNodalArray(ndofn,4)
      integer(ip) :: ivariation,ichild
      
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
      elseif (ichild == 2) then
         NewNodalArray(:,2) = OldNodalArray(:,1)
      elseif (ichild == 3) then
         NewNodalArray(:,3) = OldNodalArray(:,1)
      elseif (ichild == 4) then
         NewNodalArray(:,4) = OldNodalArray(:,1)
      endif
   
   end subroutine

   
      

   
end module
module Mod_AdaptiveElementDependentSubroutines
   use typre
   use Mod_AdaptiveTriangle
   use Mod_AdaptiveTetra 
   use Mod_AdaptiveQuads
   use Mod_AdaptiveHexa
   use Mod_AdaptiveInterpolationList
   implicit none
   
   

contains

   subroutine AdaptiveGetElementType(ndime,pnode,ielty)
      use typre
      implicit none
      integer(ip), intent(in) :: ndime,pnode
      integer(ip), intent(out) :: ielty
      
      if (ndime == 2 .and. pnode == 3) then
         ielty = 0      !Triangles
      elseif (ndime == 3 .and. pnode == 4) then
         ielty = 1      !Tetras
      elseif (ndime == 2 .and. pnode == 4) then
         ielty = 2
      elseif (ndime == 3 .and. pnode == 8) then
         ielty = 3
      endif
   end subroutine
   
   !Element types are defined as:
   ! 0: Linear triangles
   ! 1: Linear tetras
   
   ! 2: Bilinear Quads
   ! 3: Trilinear Hecas
   
   !An array of coordinates is required in order to classify tetras into the three varieties

   
   subroutine AdaptiveGetElementVariation(ielty,elcod,ivariation)
      use typre
      implicit none
      integer(ip), intent(in) ::ielty
      real(rp) :: elcod(:,:)
      integer(ip) :: ivariation
      
      select case (ielty)
         case(0)
            call LinearTriangle_GetElementVariation(elcod,ivariation)
         case(1)   
            call LinearTetra_GetElementVariation(elcod,ivariation)
         case(2)   
            call BilinearQuad_GetElementVariation(elcod,ivariation)   
         case(3)   
            call TrilinearHexa_GetElementVariation(elcod,ivariation)      
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
   

   subroutine ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
      use typre
      implicit none
      integer(ip), intent(in) ::ielty
      integer(ip), intent(out) :: pnodb,nchild,nface,nchildface
      select case (ielty)
         case(0)
            call LinearTriangle_ElementGetDimensions(pnodb,nchild,nface,nchildface)
         case(1)   
            call LinearTetra_ElementGetDimensions(pnodb,nchild,nface,nchildface)
         case(2)   
            call BilinearQuad_ElementGetDimensions(pnodb,nchild,nface,nchildface)
         case(3)   
            call TrilinearHexa_ElementGetDimensions(pnodb,nchild,nface,nchildface)      
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
   subroutine AdaptiveGetNedge(ielty,nedge)
      use typre
      implicit none
      integer(ip) :: ielty,nedge
      
      select case(ielty)
         case(0)
            call LinearTriangle_GetNedge(nedge)
         case(1)
            call LinearTetra_GetNedge(nedge)
         case(2)
            call BilinearQuad_GetNedge(nedge)   
         case(3)
            call TrilinearHexa_GetNedge(nedge)   
         case default
            call runend('GetNedge not ready for this element type')
      end select
   end subroutine
   
   subroutine Adaptive_Up_Edges(ielty,ison,childedge,parentedge)
      use typre
      implicit none
      integer(ip) :: ielty,ison, childedge,parentedge
      
      select case(ielty)
         case(0)
            call LinearTriangle_Up_Edges(ison,childedge,parentedge)
         case(1)
            call LinearTetra_Up_Edges(ison,childedge,parentedge)
         case(2)
            call BilinearQuad_Up_Edges(ison,childedge,parentedge)   
         case(3)
            call TrilinearHexa_Up_Edges(ison,childedge,parentedge)   
         case default
            call runend('Up_Edges not ready for this element type')
      end select
   end subroutine
   
   subroutine AdaptiveGetEdges(ielty,iedge,lnodedge)
      use typre
      implicit none
      integer(ip) :: ielty,iedge,lnodedge(*)
      
      select case(ielty)
         case(0)
            call LinearTriangle_GetEdges(iedge,lnodedge)
         case(1)
            call LinearTetra_GetEdges(iedge,lnodedge)
         case(2)
            call BilinearQuad_GetEdges(iedge,lnodedge)   
         case(3)
            call TrilinearHexa_GetEdges(iedge,lnodedge)      
         case default
            call runend('GetEdges not ready for this element type')
      end select
   end subroutine
   
   subroutine Adaptive_GetEdgeNodesForInode(ielty,inode,nnode,nodelist)
      use typre
      implicit none
      integer(ip) :: ielty, inode, nnode,nodelist(*)
      
      select case(ielty)
         case(0)
            call LinearTriangle_GetEdgeNodesForInode(inode,nnode,nodelist)
         case(1)
            call LinearTetra_GetEdgeNodesForInode(inode,nnode,nodelist)
         case(2)
            call BilinearQuad_GetEdgeNodesForInode(inode,nnode,nodelist)
         case(3)
            call TrilinearHexa_GetEdgeNodesForInode(inode,nnode,nodelist)      
         case default
            call runend('GetEdges not ready for this element type')
      end select
      
   end subroutine
   
   subroutine Adaptive_GetNumberOfInteriorNodesToAdd(ielty,nInteriorNodes)
      use typre
      implicit none
      integer(ip) :: ielty,nInteriorNodes
      
      select case(ielty)
         case(0)
            call LinearTriangle_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
         case(1)
            call  LinearTetra_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
         case(2)
            call  BilinearQuad_GetNumberOfInteriorNodesToAdd(nInteriorNodes)   
         case(3)
            call  TrilinearHexa_GetNumberOfInteriorNodesToAdd(nInteriorNodes)     
         case default
            call runend('GetEdges not ready for this element type')
      end select
      
   end subroutine
   
   subroutine Adaptive_GetInteriorPointFromParent(ielty,inode,jchild,jnode)
      use typre
      implicit none
      integer(ip) :: ielty,inode,jchild,jnode
      
      select case(ielty)
         case(0)
            call LinearTriangle_GetInteriorPointFromParent(inode,jchild,jnode)
         case(1)
            call LinearTetra_GetInteriorPointFromParent(inode,jchild,jnode)
         case(2)
            call BilinearQuad_GetInteriorPointFromParent(inode,jchild,jnode)
         case(3)
            call TrilinearHexa_GetInteriorPointFromParent(inode,jchild,jnode)   
         case default
            call runend('GetEdges not ready for this element type')
      end select
   end subroutine   
   
   subroutine Adaptive_AddNodeToChildInterior(ielty,ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ielty,ivariation,iadd,newpo,PointerToChildren(*),pnods(*),lnods(*)
      
      select case(ielty)
         case(0)
            call LinearTriangle_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
         case(1)
            call LinearTetra_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
         case(2)
            call BilinearQuad_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods) 
         case(3)
            call TrilinearHexa_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)      
         case default
            call runend('GetEdges not ready for this element type')
      end select
   end subroutine    
   
   subroutine add_node_to_child(ielty,ivariation,nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
      !This subroutine adds a new node to a child element when this is built
      use typre
      implicit none
      integer(ip), intent(in) :: ielty,ivariation,nchild,inode,jnode,ipoin
      integer(ip), intent(inout) :: pnods(*),lnods(*),PointerToChildren(*)
      
      integer(ip) :: knode
      integer(ip) :: ielem,ispos
       
      select case (ielty)
         case(0)
           call LinearTriangle_Add_Node_To_Child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
         case(1)   
           call LinearTetra_Add_Node_To_Child(ivariation,nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
         case(2)   
           call BilinearQuad_Add_Node_To_Child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren) 
         case(3)   
           call TrilinearHexa_Add_Node_To_Child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)    
         case default
            call runend('ElementNotReady')
      end select

      
   end subroutine

   subroutine get_edges_from_child(ielty,pnods,lnods,nchild,PointerToChildren,nedge,ledger)
      use typre
      implicit none
      integer(ip), intent(in) :: ielty,pnods(*),lnods(*),nchild,PointerToChildren(nchild)
      integer(ip), intent(out) :: nedge,ledger(3,*)
      
      integer(ip) :: ichild1,ispos1,ichild2,ispos2,ichild3,ispos3
      
      select case (ielty)
         case(0)
            call LinearTriangle_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
         case(1)   
            call LinearTetra_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
         case(2)   
            call BilinearQuad_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
         case(3)   
            call TrilinearHexa_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)   
         case default
            call runend('ElementNotReady')
      end select
      
   end subroutine   

   subroutine face_from_el(ielty,nface,pnodb,lnods,iface,sort,facnod)
      use typre
      implicit none
      integer(ip) :: iface,ielty,nface,pnodb
      integer(ip) :: facnod(pnodb,nface),lnods(*)
      integer(ip) :: sort

      select case (ielty)
         case(0)
            call LinearTriangle_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
         case(1)   
            call LinearTetra_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
         case(2)   
            call BilinearQuad_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
         case(3)   
            call TrilinearHexa_face_from_el(nface,pnodb,lnods,iface,sort,facnod)   
         case default
            call runend('ElementNotReady')
      end select   
      
   end subroutine

   subroutine face_compare(facnod1,facnod2,pnodb,isequal)
      use typre
      implicit none
      integer(ip) :: pnodb, facnod1(pnodb), facnod2(pnodb)
      integer(ip) :: inodb,isequal
      isequal = 1
      do inodb = 1,pnodb
         if (facnod1(inodb) /= facnod2(inodb)) then
            isequal = 0
            exit
         endif
      enddo
   end subroutine

   subroutine up_faces(ielty,ivariation,nface,ison,jface)
      !this subroutine gives the new neighbour face in the parent when an element is gone
      use typre
      implicit none
      integer(ip):: ielty,ivariation,nface,jface,ison
      
      select case (ielty)
         case(0)
            call LinearTriangle_up_faces(nface,ison,jface)
         case(1)   
            call LinearTetra_up_faces(ivariation,nface,ison,jface)
         case(2)   
            call BilinearQuad_up_faces(nface,ison,jface)   
         case(3)   
            call TrilinearHexa_up_faces(nface,ison,jface)      
         case default
            call runend('ElementNotReady')
      end select   
      
   end subroutine

   subroutine AdaptiveConnectInteriorFaces(ielty,ivariation,jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: ielty,ivariation,jfpos(*), son(*), lfacs(2,*)
      
      select case (ielty)
         case(0)
            call LinearTriangle_ConnectInteriorFaces(jfpos,son,lfacs)
         case(1)   
            call LinearTetra_ConnectInteriorFaces(ivariation,jfpos,son,lfacs)
         case(2)   
            call BilinearQuad_ConnectInteriorFaces(jfpos,son,lfacs)   
         case(3)   
            call TrilinearHexa_ConnectInteriorFaces(jfpos,son,lfacs)      
         case default
            call runend('ElementNotReady')
      end select   
      
      
   end subroutine

   subroutine AdaptiveGetChildrenExternalFacesList(ielty,ivariation,inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: ielty,ivariation,inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      select case (ielty)
         case(0)
            call LinearTriangle_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)
         case(1)   
            call LinearTetra_AdaptiveGetChildrenExternalFacesList(ivariation,inchild,jjchild,jjnface,jjface)
         case(2)   
            call BilinearQuad_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)   
         case(3)   
            call TrilinearHexa_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)      
         case default
            call runend('ElementNotReady')
      end select   
      
   end subroutine

   subroutine down_faces(ielty,ivariation,jface,facchild)
      use typre
      implicit none
      integer(ip):: ielty,ivariation,jface,facchild(2,*)

         select case (ielty)
            case(0)
               call LinearTriangle_down_faces(jface,facchild)
            case(1)   
               call LinearTetra_down_faces(ivariation,jface,facchild)
            case(2)   
               call BilinearQuad_down_faces(jface,facchild)   
            case(3)
               call TrilinearHexa_down_faces(jface,facchild)
            case default
               call runend('ElementNotReady')
         end select   


   end subroutine

   subroutine nodes_to_face(ielty,pnodb,lnodb,iface)
      use typre
      implicit none
      integer(ip) :: ielty,pnodb,lnodb(pnodb),iface
      
      select case (ielty)
         case(0)
            call LinearTriangle_nodes_to_face(pnodb,lnodb,iface)
         case(1)   
            call LinearTetra_nodes_to_face(pnodb,lnodb,iface)
         case(2)   
            call BilinearQuad_nodes_to_face(pnodb,lnodb,iface)   
         case(3)   
            call TrilinearHexa_nodes_to_face(pnodb,lnodb,iface)      
         case default
            call runend('ElementNotReady')
      end select   
   end subroutine
   
   subroutine Adaptive_GetInterpolationList(ielty,ivariation,InterpList)
      use typre
      implicit none
      integer(ip), intent(in) ::ielty
      integer(ip) :: ivariation
      type(InterpolationList) :: InterpList
      
      select case (ielty)
         case(0)
            call LinearTriangle_GetInterpolationList(ivariation,InterpList)
         case(1)   
            call LinearTetra_GetInterpolationList(ivariation,InterpList)
         case(2)   
            call BilinearQuad_GetInterpolationList(ivariation,InterpList)   
         case(3)   
            call TrilinearHexa_GetInterpolationList(ivariation,InterpList)      
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
   subroutine Adaptive_GetNumberofInTheFaceNodesToAdd(ielty,nInTheFaceNodes)
      use typre
      implicit none
      integer(ip) :: ielty,nInTheFaceNodes
      
       select case (ielty)
         case(0)
            call LinearTriangle_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
         case(1)   
            call LinearTetra_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
         case(2)   
            call BilinearQuad_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
         case(3)   
            call TrilinearHexa_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
         case default
            call runend('ElementNotReady')
      end select
      
   end subroutine  
   
   subroutine Adaptive_GetFaceForInTheFaceNode(ielty,inode,iface)
      use typre
      implicit none
      integer(ip) :: ielty,inode,iface
       
      select case (ielty)
         case(0)
            call LinearTriangle_GetFaceForInTheFaceNode(inode,iface)
         case(1)   
            call LinearTetra_GetFaceForInTheFaceNode(inode,iface)
         case(2)   
            call BilinearQuad_GetFaceForInTheFaceNode(inode,iface)
         case(3)   
            call TrilinearHexa_GetFaceForInTheFaceNode(inode,iface)
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
   subroutine Adaptive_GetInTheFacePointFromParent(ielty,iface,kchild,knode)
      use typre
      implicit none
      integer(ip) :: ielty,iface,kchild,knode
      
      select case (ielty)
         case(0)
            call LinearTriangle_GetInTheFacePointFromParent(iface,kchild,knode)
         case(1)   
            call LinearTetra_GetInTheFacePointFromParent(iface,kchild,knode)
         case(2)   
            call BilinearQuad_GetInTheFacePointFromParent(iface,kchild,knode)
         case(3)   
            call TrilinearHexa_GetInTheFacePointFromParent(iface,kchild,knode)
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
   subroutine Adaptive_AddNodeToChildInTheFace(ielty,ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ielty,ivariation,inode,kpoin,PointerToChildren(*), pnods(*), lnods(*)
   
      select case (ielty)
         case(0)
            call LinearTriangle_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
         case(1)   
            call LinearTetra_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
         case(2)   
            call BilinearQuad_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
         case(3)   
            call TrilinearHexa_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
         case default
            call runend('ElementNotReady')
      end select
   end subroutine
   
    subroutine ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ielty,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(:,:),NewNodalArray(:,:)
      integer(ip) :: ielty,ivariation,ichild
     
      select case (ielty)
         case(0)
            call LinearTriangle_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(1)   
            call LinearTetra_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(2)   
            call BilinearQuad_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(3)   
            call TrilinearHexa_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case default
            call runend('ElementNotReady')
      end select
   
   
   
   end subroutine
   
   subroutine ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ielty,ivariation,ichild     )
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(:,:),NewNodalArray(:,:)
      integer(ip) :: ielty,ivariation,ichild
     
      select case (ielty)
         case(0)
            call LinearTriangle_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(1)   
            call LinearTetra_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(2)   
            call BilinearQuad_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case(3)   
            call TrilinearHexa_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
         case default
            call runend('ElementNotReady')
      end select
   
   
   
   end subroutine

   
end module
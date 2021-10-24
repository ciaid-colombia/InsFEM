module Mod_AdaptiveTetra
   use Mod_AdaptiveInterpolationList
   implicit none

contains

   subroutine LinearTetra_GetElementVariation(elcod,ivariation)
      use typre
      implicit none
      real(rp) :: elcod(:,:)
      integer(ip) :: ivariation
      
      real(rp) :: coord(3,10)
      real(rp) :: diag1(3),diag2(3),diag3(3),norm(3)
      
      
      
      !We need to find which the shortest diagonal is, and then set the variations
      coord(:,5) = (elcod(:,1)+elcod(:,2))*0.5_rp
      coord(:,6) = (elcod(:,1)+elcod(:,3))*0.5_rp
      coord(:,7) = (elcod(:,1)+elcod(:,4))*0.5_rp
      coord(:,8) = (elcod(:,2)+elcod(:,3))*0.5_rp
      coord(:,9) = (elcod(:,2)+elcod(:,4))*0.5_rp
      coord(:,10) = (elcod(:,3)+elcod(:,4))*0.5_rp
      
      diag1(:) = coord(:,7)-coord(:,8)
      diag2(:) = coord(:,6)-coord(:,9)
      diag3(:) = coord(:,5)-coord(:,10)
      
      norm(1) = dot_product(diag1,diag1)
      norm(2) = dot_product(diag2,diag2)
      norm(3) = dot_product(diag3,diag3)
      
      ivariation = minloc(norm,1)
      
   end subroutine   
      

   subroutine LinearTetra_ElementGetDimensions(pnodb,nchild,nface,nchildface)
      use typre
      implicit none
      integer(ip), intent(out) :: pnodb,nchild,nface,nchildface
   
      pnodb = 3
      nchild = 8
      nface = 4
      nchildface = 4
   end subroutine
   
   subroutine LinearTetra_nodes_to_face(pnodb,lnodb,iface)
      use typre
      implicit none
      integer(ip) :: ndime,pnode,pnodb,lnodb(pnodb),iface
      
      integer(ip) :: auxlnodb(pnodb)
      
      auxlnodb = lnodb
      call sortb(3,auxlnodb)
      
      if (auxlnodb(1) == 1 .and. auxlnodb(2) == 2 .and. auxlnodb(3) == 3) then
         iface = 1
      elseif (auxlnodb(1) == 1 .and. auxlnodb(2) == 3 .and. auxlnodb(3) == 4) then
         iface = 2
      elseif (auxlnodb(1) == 2 .and. auxlnodb(2) == 3 .and. auxlnodb(3) == 4) then
         iface = 3
      elseif (auxlnodb(1) == 1 .and. auxlnodb(2) == 2 .and. auxlnodb(3) == 4) then
         iface = 4   
      endif   
   end subroutine

   
   subroutine LinearTetra_GetNedge(nedge)
      use typre
      implicit none
      integer(ip) :: nedge
      
      nedge = 6
   end subroutine
   
   subroutine LinearTetra_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
      use typre
      implicit none
      integer(ip), intent(in) :: pnods(*),lnods(*),nchild,PointerToChildren(nchild)
      integer(ip), intent(out) :: nedge,ledger(3,*)
      integer(ip) :: ispos(4)

      nedge = 6
      
      ispos(:) = pnods(PointerToChildren(:))-1
      
      !first edge
      ledger(1,1) = lnods(ispos(1)+1)
      ledger(2,1) = lnods(ispos(2)+1)
      ledger(3,1) = lnods(ispos(1)+2)
      
      !second edge
      ledger(1,2) = lnods(ispos(1)+1)
      ledger(2,2) = lnods(ispos(3)+1)
      ledger(3,2) = lnods(ispos(1)+3)
      
      !third edge
      ledger(1,3) = lnods(ispos(1)+1)
      ledger(2,3) = lnods(ispos(4)+1)
      ledger(3,3) = lnods(ispos(1)+4)
      
      !fourth edge
      ledger(1,4) = lnods(ispos(2)+1)
      ledger(2,4) = lnods(ispos(3)+1)
      ledger(3,4) = lnods(ispos(2)+2)

      !fifth edge
      ledger(1,5) = lnods(ispos(2)+1)
      ledger(2,5) = lnods(ispos(4)+1)
      ledger(3,5) = lnods(ispos(2)+4)
      
      !sixth edge
      ledger(1,6) = lnods(ispos(3)+1)
      ledger(2,6) = lnods(ispos(4)+1)
      ledger(3,6) = lnods(ispos(3)+4)
      
   end subroutine 
   
   subroutine LinearTetra_Up_Edges(ison,childedge,parentedge)
      use typre
      implicit none
      integer(ip) :: ison, childedge,parentedge
      
      if (ison == 1) then
         select case (childedge)
            case (1)
               parentedge = 1
            case(2)
               parentedge = 2
            case(3)
               parentedge = 3
            case default
               parentedge = -1
         end select
      elseif (ison == 2) then
         select case(childedge)
            case(1)
               parentedge = 4
            case(2)
               parentedge = 1
            case(3)
               parentedge = 5
            case default
               parentedge = -1
         end select
      elseif (ison == 3) then
         select case(childedge)
            case(1)
               parentedge = 2
            case(2)
               parentedge = 4
            case(3)
               parentedge = 6
            case default
               parentedge = -1
         end select
       elseif (ison == 4) then
         select case(childedge)
            case(1)
               parentedge = 6
            case(2)
               parentedge = 5
            case(3)
               parentedge = 3
            case default
               parentedge = -1
         end select   
      !The rest of (interior) children have no parent edges, edges are completely new
      else
         parentedge = -1
      endif   
   end subroutine
   
   subroutine LinearTetra_GetEdges(iedge,lnodedge)
      use typre
      implicit none
      integer(ip) :: iedge,lnodedge(2)
      
      select case(iedge)
         case(1)
            lnodedge(1) = 1
            lnodedge(2) = 2
         case(2)
            lnodedge(1) = 1
            lnodedge(2) = 3
         case(3)
            lnodedge(1) = 1
            lnodedge(2) = 4
         case(4)
            lnodedge(1) = 2
            lnodedge(2) = 3
         case(5)
            lnodedge(1) = 2
            lnodedge(2) = 4
         case(6)
            lnodedge(1) = 3
            lnodedge(2) = 4
      end select
   end subroutine
   
   subroutine LinearTetra_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
      use typre
      implicit none

      integer(ip) :: iface,nface,pnodb
      integer(ip) :: facnod(pnodb,nface),lnods(*)
      integer(ip) :: sort,auxfaces(3),faces(3,4)
      integer(ip) :: jface
      
!       faces(:,1) = (/1, 2, 3 /)
!       faces(:,2) = (/1, 3, 4 /)
!       faces(:,3) = (/2, 4, 3 /)
!       faces(:,4) = (/1, 4, 2 /)
      
      faces = reshape((/   1, 2, 3, &
                           1, 3, 4, &
                           2, 4, 3,  &
                           1, 4, 2  &
                           /), shape(faces))

      if (iface == 0) then
         do jface = 1,4
            auxfaces = faces(:,jface)
            facnod(:,jface) = lnods(auxfaces)
         enddo
         if (sort == 1) then
               call sortb(3,facnod(1,2))
               call sortb(3,facnod(1,3))
               call sortb(3,facnod(1,4))
         endif
      else
         auxfaces = faces(:,iface)
         facnod(:,1) = lnods(auxfaces)
      endif   
      if (sort == 1) call sortb(3,facnod(1,1))     
   end subroutine
   
   !-------------------------------------------------------------------
   !*******************************************************************
   !Interface for variations
   subroutine LinearTetra_Add_node_to_child(ivariation,nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
      !This subroutine adds a new node to a child element when this is built
      use typre
      implicit none
      integer(ip), intent(in) :: ivariation,nchild,inode,jnode,ipoin
      integer(ip), intent(inout) :: pnods(*),lnods(*),PointerToChildren(*)
      
      select case (ivariation)
         case(1)
            call LinearTetra_Add_node_to_child1(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
         case(2)
            call LinearTetra_Add_node_to_child2(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
         case(3)
            call LinearTetra_Add_node_to_child3(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
      end select
      
   end subroutine
   
   subroutine LinearTetra_up_faces(ivariation,nface,ison,jface)
      !this subroutine gives the new neighbour face in the parent when an element is gone
      use typre
      implicit none
      integer(ip):: ivariation,ielty,nface,jface,ison

      select case (ivariation)
         case(1)
            call LinearTetra_up_faces1(nface,ison,jface)
         case(2)
            call LinearTetra_up_faces2(nface,ison,jface)
         case(3)
            call LinearTetra_up_faces3(nface,ison,jface)
      end select
      
   end subroutine
   
   subroutine LinearTetra_down_faces(ivariation,jface,facchild)
   use typre
   implicit none
   integer(ip):: ivariation,jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (ivariation)
         case(1)
            call LinearTetra_down_faces1(jface,facchild)
         case(2)
            call LinearTetra_down_faces2(jface,facchild)
         case(3)
            call LinearTetra_down_faces3(jface,facchild)
      end select

   end subroutine
      
   subroutine LinearTetra_ConnectInteriorFaces(ivariation,jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: ivariation,jfpos(*), son(*), lfacs(2,*)
      
      select case (ivariation)
         case(1)
            call LinearTetra_ConnectInteriorFaces1(jfpos,son,lfacs)
         case(2)
            call LinearTetra_ConnectInteriorFaces2(jfpos,son,lfacs)
         case(3)
            call LinearTetra_ConnectInteriorFaces3(jfpos,son,lfacs)
      end select
   end subroutine

   subroutine LinearTetra_AdaptiveGetChildrenExternalFacesList(ivariation,inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: ivariation,inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      select case (ivariation)
         case(1)
            call LinearTetra_AdaptiveGetChildrenExternalFacesList1(inchild,jjchild,jjnface,jjface)
         case(2)
            call LinearTetra_AdaptiveGetChildrenExternalFacesList2(inchild,jjchild,jjnface,jjface)
         case(3)
            call LinearTetra_AdaptiveGetChildrenExternalFacesList3(inchild,jjchild,jjnface,jjface)
      end select
   end subroutine
   
   
   
   
   !-------------------------------------------------------------------
   !*******************************************************************
   !Variation 1
   
   subroutine LinearTetra_Add_node_to_child1(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
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
      elseif ((inode == 1 .and. jnode == 3) .or. (inode == 3 .and. jnode == 1)) then
         knode = 6
      elseif ((inode == 1 .and. jnode == 4) .or. (inode == 4 .and. jnode == 1)) then
         knode = 7    
      elseif ((inode == 2 .and. jnode == 3) .or. (inode == 3 .and. jnode == 2)) then
         knode = 8
      elseif ((inode == 2 .and. jnode == 4) .or. (inode == 4 .and. jnode == 2)) then
         knode = 9    
      elseif ((inode == 3 .and. jnode == 4) .or. (inode == 4 .and. jnode == 3)) then
         knode = 10        
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
         lnods(ispos+3) = ipoin
      
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
      case(6)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      
      case(7)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
      case(8)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
      case(9)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
      case(10)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      end select
   end subroutine
   
   subroutine LinearTetra_up_faces1(nface,ison,jface)
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
            elseif (jface == 4) then
               jface = 4
            else
               call runend('LinearTetra_up_faces')
            endif
         case (2)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then 
               jface = 3
            else
               call runend('LinearTetra_upfaces')
            endif
               
         case (3)
            if (jface == 1) then
               jface = 1   
            elseif (jface == 2) then
               jface = 3
            elseif (jface == 4) then 
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(4)
            if (jface == 1) then 
               jface = 3   
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(5)
            if (jface == 2) then 
               jface = 1   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(6)
            if (jface == 3) then 
               jface = 4   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(7)
            if (jface == 2) then
               jface = 2   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(8)
            if (jface == 3) then
               jface = 3
            else
               call runend('LinearTetra_upfaces')
            endif
      end select
      
   end subroutine
   
   subroutine LinearTetra_down_faces1(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1
         facchild(2,1) = 1 !first son, first face
         facchild(1,2) = 2
         facchild(2,2) = 1 !second son, first face
         facchild(1,3) = 3
         facchild(2,3) = 1
         facchild(1,4) = 5
         facchild(2,4) = 2 
      case (2)
         facchild(1,1) = 1
         facchild(2,1) = 2 
         facchild(1,2) = 3
         facchild(2,2) = 4 
         facchild(1,3) = 4
         facchild(2,3) = 4 
         facchild(1,4) = 7
         facchild(2,4) = 2 
      case (3)
         facchild(1,1) = 2
         facchild(2,1) = 4 
         facchild(1,2) = 3
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 1 
         facchild(1,4) = 8
         facchild(2,4) = 3 
      case(4)
         facchild(1,1) = 1
         facchild(2,1) = 4 
         facchild(1,2) = 2
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 2 
         facchild(1,4) = 6
         facchild(2,4) = 3 
      end select

   end subroutine
      
   subroutine LinearTetra_ConnectInteriorFaces1(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      !Element 1 to 5
      lfacs(1,jfpos(1)+3) = son(5)
      lfacs(2,jfpos(1)+3) = 1
      
      lfacs(1,jfpos(5)+1) = son(1)
      lfacs(2,jfpos(5)+1) = 3
      
      !Element 2 to 6
      lfacs(1,jfpos(2)+3) = son(6)
      lfacs(2,jfpos(2)+3) = 2
      
      lfacs(1,jfpos(6)+2) = son(2)
      lfacs(2,jfpos(6)+2) = 3
      
      !Element 3 to 7
      lfacs(1,jfpos(3)+3) = son(7)
      lfacs(2,jfpos(3)+3) = 1
      
      lfacs(1,jfpos(7)+1) = son(3)
      lfacs(2,jfpos(7)+1) = 3
      
      !Element 4 to 8
      lfacs(1,jfpos(4)+3) = son(8)
      lfacs(2,jfpos(4)+3) = 2
      
      lfacs(1,jfpos(8)+2) = son(4)
      lfacs(2,jfpos(8)+2) = 3
      
      !Element 5 to 6
      lfacs(1,jfpos(5)+3) = son(6)
      lfacs(2,jfpos(5)+3) = 1
      
      lfacs(1,jfpos(6)+1) = son(5)
      lfacs(2,jfpos(6)+1) = 3
      
      !Element 5 to 7
      lfacs(1,jfpos(5)+4) = son(7)
      lfacs(2,jfpos(5)+4) = 4
      
      lfacs(1,jfpos(7)+4) = son(5)
      lfacs(2,jfpos(7)+4) = 4
      
      !Element 6 to 8
      lfacs(1,jfpos(6)+4) = son(8)
      lfacs(2,jfpos(6)+4) = 4
      
      lfacs(1,jfpos(8)+4) = son(6)
      lfacs(2,jfpos(8)+4) = 4
      
      !Element 7 to 8
      lfacs(1,jfpos(7)+3) = son(8)
      lfacs(2,jfpos(7)+3) = 1
      
      lfacs(1,jfpos(8)+1) = son(7)
      lfacs(2,jfpos(8)+1) = 3
      
      
      
      
   end subroutine

   subroutine LinearTetra_AdaptiveGetChildrenExternalFacesList1(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      integer(ip) :: i
      
      !there are 8 children with faces shared with the parent element
      inchild = 8
      !These elements are child 1 to 8
      do i = 1,8
         jjchild(i) = i
      enddo
      
      !First four children (external) have 3 external faces
      jjnface(1:4) = 3
      !Second (internal) children have 1 external face
      jjnface(5:8) = 1
      
      !Which are the external faces for each children
      jjface(1,1) = 1
      jjface(2,1) = 2
      jjface(3,1) = 4
      
      jjface(1,2) = 1
      jjface(2,2) = 2
      jjface(3,2) = 4
      
      jjface(1,3) = 1
      jjface(2,3) = 2
      jjface(3,3) = 4
      
      jjface(1,4) = 1
      jjface(2,4) = 2
      jjface(3,4) = 4
      
      jjface(1,5) = 2
      
      jjface(1,6) = 3
      
      jjface(1,7) = 2
      
      jjface(1,8) = 3
      
      
   end subroutine
   
   !-------------------------------------------------------------------
   !*******************************************************************
   !Variation 2
   
   subroutine LinearTetra_Add_node_to_child2(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
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
      elseif ((inode == 1 .and. jnode == 3) .or. (inode == 3 .and. jnode == 1)) then
         knode = 6
      elseif ((inode == 1 .and. jnode == 4) .or. (inode == 4 .and. jnode == 1)) then
         knode = 7    
      elseif ((inode == 2 .and. jnode == 3) .or. (inode == 3 .and. jnode == 2)) then
         knode = 8
      elseif ((inode == 2 .and. jnode == 4) .or. (inode == 4 .and. jnode == 2)) then
         knode = 9    
      elseif ((inode == 3 .and. jnode == 4) .or. (inode == 4 .and. jnode == 3)) then
         knode = 10        
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
         lnods(ispos+3) = ipoin
      
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
      case(6)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
      
      case(7)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin

         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
      case(8)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
      case(9)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
      case(10)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      end select
   end subroutine
   
   subroutine LinearTetra_up_faces2(nface,ison,jface)
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
            elseif (jface == 4) then
               jface = 4
            else
               call runend('LinearTetra_up_faces')
            endif
         case (2)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then 
               jface = 3
            else
               call runend('LinearTetra_upfaces')
            endif
               
         case (3)
            if (jface == 1) then
               jface = 1   
            elseif (jface == 2) then
               jface = 3
            elseif (jface == 4) then 
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(4)
            if (jface == 1) then 
               jface = 3   
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(5)
            if (jface == 3) then 
               jface = 4   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(6)
            if (jface == 1) then 
               jface = 1   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(7)
            if (jface == 3) then
               jface = 3   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(8)
            if (jface == 1) then
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
      end select
      
   end subroutine
   
   subroutine LinearTetra_down_faces2(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1
         facchild(2,1) = 1 !first son, first face
         facchild(1,2) = 2
         facchild(2,2) = 1 !second son, first face
         facchild(1,3) = 3
         facchild(2,3) = 1
         facchild(1,4) = 6
         facchild(2,4) = 1 
      case (2)
         facchild(1,1) = 1
         facchild(2,1) = 2 
         facchild(1,2) = 3
         facchild(2,2) = 4 
         facchild(1,3) = 4
         facchild(2,3) = 4 
         facchild(1,4) = 8
         facchild(2,4) = 1 
      case (3)
         facchild(1,1) = 2
         facchild(2,1) = 4 
         facchild(1,2) = 3
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 1 
         facchild(1,4) = 7
         facchild(2,4) = 3 
      case(4)
         facchild(1,1) = 1
         facchild(2,1) = 4 
         facchild(1,2) = 2
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 2 
         facchild(1,4) = 5
         facchild(2,4) = 3 
      end select

   end subroutine
  
   subroutine LinearTetra_ConnectInteriorFaces2(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      !Element 1 to 5
      lfacs(1,jfpos(1)+3) = son(5)
      lfacs(2,jfpos(1)+3) = 1
      
      lfacs(1,jfpos(5)+1) = son(1)
      lfacs(2,jfpos(5)+1) = 3
      
      !Element 2 to 6
      lfacs(1,jfpos(2)+3) = son(6)
      lfacs(2,jfpos(2)+3) = 2
      
      lfacs(1,jfpos(6)+2) = son(2)
      lfacs(2,jfpos(6)+2) = 3
      
      !Element 3 to 7
      lfacs(1,jfpos(3)+3) = son(7)
      lfacs(2,jfpos(3)+3) = 1
      
      lfacs(1,jfpos(7)+1) = son(3)
      lfacs(2,jfpos(7)+1) = 3
      
      !Element 4 to 8
      lfacs(1,jfpos(4)+3) = son(8)
      lfacs(2,jfpos(4)+3) = 2
      
      lfacs(1,jfpos(8)+2) = son(4)
      lfacs(2,jfpos(8)+2) = 3
      
      !Element 5 to 6
      lfacs(1,jfpos(5)+2) = son(6)
      lfacs(2,jfpos(5)+2) = 3
      
      lfacs(1,jfpos(6)+3) = son(5)
      lfacs(2,jfpos(6)+3) = 2
      
      !Element 5 to 8
      lfacs(1,jfpos(5)+4) = son(8)
      lfacs(2,jfpos(5)+4) = 4
      
      lfacs(1,jfpos(8)+4) = son(5)
      lfacs(2,jfpos(8)+4) = 4
      
      !Element 6 to 7
      lfacs(1,jfpos(6)+4) = son(7)
      lfacs(2,jfpos(6)+4) = 4
      
      lfacs(1,jfpos(7)+4) = son(6)
      lfacs(2,jfpos(7)+4) = 4
      
      !Element 7 to 8
      lfacs(1,jfpos(7)+2) = son(8)
      lfacs(2,jfpos(7)+2) = 3
      
      lfacs(1,jfpos(8)+3) = son(7)
      lfacs(2,jfpos(8)+3) = 2
      
      
      
      
   end subroutine

   subroutine LinearTetra_AdaptiveGetChildrenExternalFacesList2(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      integer(ip) :: i
      
      !there are 8 children with faces shared with the parent element
      inchild = 8
      !These elements are child 1 to 8
      do i = 1,8
         jjchild(i) = i
      enddo
      
      !First four children (external) have 3 external faces
      jjnface(1:4) = 3
      !Second (internal) children have 1 external face
      jjnface(5:8) = 1
      
      !Which are the external faces for each children
      jjface(1,1) = 1
      jjface(2,1) = 2
      jjface(3,1) = 4
      
      jjface(1,2) = 1
      jjface(2,2) = 2
      jjface(3,2) = 4
      
      jjface(1,3) = 1
      jjface(2,3) = 2
      jjface(3,3) = 4
      
      jjface(1,4) = 1
      jjface(2,4) = 2
      jjface(3,4) = 4
      
      jjface(1,5) = 3
      
      jjface(1,6) = 1
      
      jjface(1,7) = 3
      
      jjface(1,8) = 1
      
      
   end subroutine
   
   !-------------------------------------------------------------------
   !*******************************************************************
   !Variation 3
   
   subroutine LinearTetra_Add_node_to_child3(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
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
      elseif ((inode == 1 .and. jnode == 3) .or. (inode == 3 .and. jnode == 1)) then
         knode = 6
      elseif ((inode == 1 .and. jnode == 4) .or. (inode == 4 .and. jnode == 1)) then
         knode = 7    
      elseif ((inode == 2 .and. jnode == 3) .or. (inode == 3 .and. jnode == 2)) then
         knode = 8
      elseif ((inode == 2 .and. jnode == 4) .or. (inode == 4 .and. jnode == 2)) then
         knode = 9    
      elseif ((inode == 3 .and. jnode == 4) .or. (inode == 4 .and. jnode == 3)) then
         knode = 10        
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
         lnods(ispos+3) = ipoin
      
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
      case(6)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      
      case(7)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin

         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
      case(8)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
      case(9)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
      case(10)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      end select
   end subroutine
   
   subroutine LinearTetra_up_faces3(nface,ison,jface)
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
            elseif (jface == 4) then
               jface = 4
            else
               call runend('LinearTetra_up_faces')
            endif
         case (2)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then 
               jface = 3
            else
               call runend('LinearTetra_upfaces')
            endif
               
         case (3)
            if (jface == 1) then
               jface = 1   
            elseif (jface == 2) then
               jface = 3
            elseif (jface == 4) then 
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(4)
            if (jface == 1) then 
               jface = 3   
            elseif (jface == 2) then
               jface = 4
            elseif (jface == 4) then
               jface = 2
            else
               call runend('LinearTetra_upfaces')
            endif
         case(5)
            if (jface == 4) then 
               jface = 2   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(6)
            if (jface == 4) then 
               jface = 3   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(7)
            if (jface == 4) then
               jface = 1   
            else
               call runend('LinearTetra_upfaces')
            endif
         case(8)
            if (jface == 4) then
               jface = 4
            else
               call runend('LinearTetra_upfaces')
            endif
      end select
      
   end subroutine
   
   subroutine LinearTetra_down_faces3(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1
         facchild(2,1) = 1 !first son, first face
         facchild(1,2) = 2
         facchild(2,2) = 1 !second son, first face
         facchild(1,3) = 3
         facchild(2,3) = 1
         facchild(1,4) = 7
         facchild(2,4) = 4 
      case (2)
         facchild(1,1) = 1
         facchild(2,1) = 2 
         facchild(1,2) = 3
         facchild(2,2) = 4 
         facchild(1,3) = 4
         facchild(2,3) = 4 
         facchild(1,4) = 5
         facchild(2,4) = 4 
      case (3)
         facchild(1,1) = 2
         facchild(2,1) = 4 
         facchild(1,2) = 3
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 1 
         facchild(1,4) = 6
         facchild(2,4) = 4 
      case(4)
         facchild(1,1) = 1
         facchild(2,1) = 4 
         facchild(1,2) = 2
         facchild(2,2) = 2 
         facchild(1,3) = 4
         facchild(2,3) = 2 
         facchild(1,4) = 8
         facchild(2,4) = 4 
      end select

   end subroutine
  
   subroutine LinearTetra_ConnectInteriorFaces3(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      !Element 1 to 5
      lfacs(1,jfpos(1)+3) = son(5)
      lfacs(2,jfpos(1)+3) = 1
      
      lfacs(1,jfpos(5)+1) = son(1)
      lfacs(2,jfpos(5)+1) = 3
      
      !Element 2 to 6
      lfacs(1,jfpos(2)+3) = son(6)
      lfacs(2,jfpos(2)+3) = 2
      
      lfacs(1,jfpos(6)+2) = son(2)
      lfacs(2,jfpos(6)+2) = 3
      
      !Element 3 to 7
      lfacs(1,jfpos(3)+3) = son(7)
      lfacs(2,jfpos(3)+3) = 1
      
      lfacs(1,jfpos(7)+1) = son(3)
      lfacs(2,jfpos(7)+1) = 3
      
      !Element 4 to 8
      lfacs(1,jfpos(4)+3) = son(8)
      lfacs(2,jfpos(4)+3) = 2
      
      lfacs(1,jfpos(8)+2) = son(4)
      lfacs(2,jfpos(8)+2) = 3
      
      !Element 5 to 7
      lfacs(1,jfpos(5)+2) = son(7)
      lfacs(2,jfpos(5)+2) = 2
      
      lfacs(1,jfpos(7)+2) = son(5)
      lfacs(2,jfpos(7)+2) = 2
      
      !Element 5 to 8
      lfacs(1,jfpos(5)+3) = son(8)
      lfacs(2,jfpos(5)+3) = 1
      
      lfacs(1,jfpos(8)+1) = son(5)
      lfacs(2,jfpos(8)+1) = 3
      
      !Element 6 to 7
      lfacs(1,jfpos(6)+1) = son(7)
      lfacs(2,jfpos(6)+1) = 3
      
      lfacs(1,jfpos(7)+3) = son(6)
      lfacs(2,jfpos(7)+3) = 1
      
      !Element 6 to 8
      lfacs(1,jfpos(6)+3) = son(8)
      lfacs(2,jfpos(6)+3) = 3
      
      lfacs(1,jfpos(8)+3) = son(6)
      lfacs(2,jfpos(8)+3) = 3
      
      
      
      
   end subroutine

   subroutine LinearTetra_AdaptiveGetChildrenExternalFacesList3(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      integer(ip) :: i
      
      !there are 8 children with faces shared with the parent element
      inchild = 8
      !These elements are child 1 to 8
      do i = 1,8
         jjchild(i) = i
      enddo
      
      !First four children (external) have 3 external faces
      jjnface(1:4) = 3
      !Second (internal) children have 1 external face
      jjnface(5:8) = 1
      
      !Which are the external faces for each children
      jjface(1,1) = 1
      jjface(2,1) = 2
      jjface(3,1) = 4
      
      jjface(1,2) = 1
      jjface(2,2) = 2
      jjface(3,2) = 4
      
      jjface(1,3) = 1
      jjface(2,3) = 2
      jjface(3,3) = 4
      
      jjface(1,4) = 1
      jjface(2,4) = 2
      jjface(3,4) = 4
      
      jjface(1,5) = 4
      
      jjface(1,6) = 4
      
      jjface(1,7) = 4
      
      jjface(1,8) = 4
      
      
   end subroutine
   
   subroutine LinearTetra_GetInterpolationList(ivariation,InterpList)
      use typre
      implicit none
      integer(ip) :: ivariation
      type(InterpolationList) :: InterpList
      
      InterpList%n = 0
   end subroutine
   
   subroutine LinearTetra_GetEdgeNodesForInode(inode,nnode,nodelist)
      use typre
      implicit none
      integer(ip) :: inode, nnode, nodelist(3)
      
      nnode = 3
      
      if (inode == 1) then
         nodelist(1) = 2
         nodelist(2) = 3
         nodelist(3) = 4
      elseif (inode == 2) then
         nodelist(1) = 1
         nodelist(2) = 3
         nodelist(3) = 4
      elseif (inode == 3) then
         nodelist(1) = 1
         nodelist(2) = 2
         nodelist(3) = 4
      elseif (inode == 4) then
         nodelist(1) = 1
         nodelist(2) = 2
         nodelist(3) = 3   
      endif
   end subroutine
   
    subroutine LinearTetra_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
      use typre
      implicit none
      integer(ip) :: nInteriorNodes
      
      nInteriorNodes = 0
      
   end subroutine
   
   subroutine LinearTetra_GetInteriorPointFromParent(inode,jchild,jnode)
      use typre
      implicit none
      integer(ip) :: inode,jchild,jnode
      
      
      call runend('Wrong number of interior nodes')
   end subroutine   
   
   subroutine LinearTetra_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,iadd,newpo,PointerToChildren(*),pnods(*),lnods(*)
      
      call runend('Should not add interior node in tetra')
   end subroutine    
   
   subroutine  LinearTetra_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
      use typre
      implicit none
      integer(ip) :: nInTheFaceNodes
      
      nInTheFaceNodes = 0
   end subroutine
   
   subroutine LinearTetra_GetFaceForInTheFaceNode(inode,iface)
      use typre
      implicit none
      integer(ip) :: inode,iface
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine   
   
   subroutine LinearTetra_GetInTheFacePointFromParent(iface,kchild,knode)
      use typre
      implicit none
      integer(ip) :: iface,kchild,knode
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   subroutine LinearTetra_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,inode,kpoin,PointerToChildren(*),pnods(*),lnods(*)
      
      call runend('Adaptive: This subroutine should not be called')
   end subroutine
   
   subroutine LinearTetra_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,4),NewNodalArray(ndofn,4)
      integer(ip) :: ielty,ivariation,ichild
   
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
      elseif (ichild == 2) then
         NewNodalArray(:,1) = OldNodalArray(:,2)
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
      elseif (ichild == 3) then
         NewNodalArray(:,1) = OldNodalArray(:,3)
         NewNodalArray(:,2) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
      elseif (ichild == 4) then
         NewNodalArray(:,1) = OldNodalArray(:,4)
         NewNodalArray(:,2) = (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,4)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,4)+OldNodalArray(:,1))*0.5
      else
         select case (ivariation)
         case(1)
            if (ichild == 5) then
               NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
            elseif (ichild == 6) then
               NewNodalArray(:,1) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
            elseif (ichild == 7) then
               NewNodalArray(:,1) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
            elseif (ichild == 8) then
               NewNodalArray(:,1) = (OldNodalArray(:,4)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,4)+OldNodalArray(:,2))*0.5      
            endif   
         case(2)
            if (ichild == 5) then
               NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
            elseif (ichild == 6) then
               NewNodalArray(:,1) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
            elseif (ichild == 7) then
               NewNodalArray(:,1) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
            elseif (ichild == 8) then
               NewNodalArray(:,1) = (OldNodalArray(:,4)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,4)+OldNodalArray(:,2))*0.5      
            endif   
         case(3)
            if (ichild == 5) then
               NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
            elseif (ichild == 6) then
               NewNodalArray(:,1) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,4))*0.5
            elseif (ichild == 7) then
               NewNodalArray(:,1) = (OldNodalArray(:,3)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,3)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,2)+OldNodalArray(:,1))*0.5
            elseif (ichild == 8) then
               NewNodalArray(:,1) = (OldNodalArray(:,4)+OldNodalArray(:,1))*0.5
               NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
               NewNodalArray(:,3) = (OldNodalArray(:,4)+OldNodalArray(:,3))*0.5
               NewNodalArray(:,4) = (OldNodalArray(:,4)+OldNodalArray(:,2))*0.5      
            endif   
      
         end select
      endif
   
   
   end subroutine
   
   subroutine LinearTetra_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
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


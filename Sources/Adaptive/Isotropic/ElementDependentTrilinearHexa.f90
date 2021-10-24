module Mod_AdaptiveHexa
   use Mod_AdaptiveInterpolationList
   implicit none
   


contains

   subroutine TrilinearHexa_GetElementVariation(elcod,ivariation)
      use typre
      implicit none
      real(rp) :: elcod(:,:)
      integer(ip) :: ivariation
      
      ivariation = 1
      
   end subroutine

   subroutine TrilinearHexa_ElementGetDimensions(pnodb,nchild,nface,nchildface)
      use typre
      implicit none
      integer(ip), intent(out) :: pnodb,nchild,nface,nchildface
   
      pnodb = 4
      nchild = 8
      nface = 6
      nchildface = 4
   end subroutine
   
   subroutine TrilinearHexa_Add_node_to_child(nchild,inode,jnode,pnods,lnods,ipoin,PointerToChildren)
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
         knode = 9
      elseif ((inode == 1 .and. jnode == 4) .or. (inode == 4 .and. jnode == 1)) then
         knode = 10
      elseif ((inode == 1 .and. jnode == 5) .or. (inode == 5 .and. jnode == 1)) then
         knode = 11    
      elseif ((inode == 2 .and. jnode == 3) .or. (inode == 3 .and. jnode == 2)) then
         knode = 12
      elseif ((inode == 2 .and. jnode == 6) .or. (inode == 6 .and. jnode == 2)) then
         knode = 13    
      elseif ((inode == 3 .and. jnode == 4) .or. (inode == 4 .and. jnode == 3)) then
         knode = 14
      elseif ((inode == 3 .and. jnode == 7) .or. (inode == 7 .and. jnode == 3)) then
         knode = 15        
      elseif ((inode == 4 .and. jnode == 8) .or. (inode == 8 .and. jnode == 4)) then
         knode = 16           
      elseif ((inode == 5 .and. jnode == 6) .or. (inode == 6 .and. jnode == 5)) then
         knode = 17           
      elseif ((inode == 5 .and. jnode == 8) .or. (inode == 8 .and. jnode == 5)) then
         knode = 18           
      elseif ((inode == 6 .and. jnode == 7) .or. (inode == 7 .and. jnode == 6)) then
         knode = 19           
      elseif ((inode == 7 .and. jnode == 8) .or. (inode == 8 .and. jnode == 7)) then
         knode = 20           
      endif   
      
      select case (knode)
      case(1)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      case(2)
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
      case(3)
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      case(4)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
      case(5)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+5) = ipoin
      case(6)
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+6) = ipoin
      case(7)
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+7) = ipoin
      case(8)
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+8) = ipoin
      case(9)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
         
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      case(10)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
         
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      case(11)
         ielem = PointerToChildren(1)
         ispos = pnods(ielem)-1
         lnods(ispos+5) = ipoin
         
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+1) = ipoin
      case(12)
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
      case(13)
         ielem = PointerToChildren(5)
         ispos = pnods(ielem)-1
         lnods(ispos+6) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+2) = ipoin
      case(14)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
         
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
      case(15)
         ielem = PointerToChildren(7)
         ispos = pnods(ielem)-1
         lnods(ispos+7) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+3) = ipoin
      case(16)
         ielem = PointerToChildren(3)
         ispos = pnods(ielem)-1
         lnods(ispos+8) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+4) = ipoin
      case(17)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+6) = ipoin
         
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+5) = ipoin
      case(18)
         ielem = PointerToChildren(2)
         ispos = pnods(ielem)-1
         lnods(ispos+8) = ipoin
         
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+5) = ipoin
      case(19)
         ielem = PointerToChildren(6)
         ispos = pnods(ielem)-1
         lnods(ispos+7) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+6) = ipoin
      case(20)
         ielem = PointerToChildren(4)
         ispos = pnods(ielem)-1
         lnods(ispos+7) = ipoin
         
         ielem = PointerToChildren(8)
         ispos = pnods(ielem)-1
         lnods(ispos+8) = ipoin
      end select
      
   end subroutine
   
   subroutine TrilinearHexa_get_edges_from_child(pnods,lnods,nchild,PointerToChildren,nedge,ledger)
      use typre
      implicit none
      integer(ip), intent(in) :: pnods(*),lnods(*),nchild,PointerToChildren(nchild)
      integer(ip), intent(out) :: nedge,ledger(3,*)
      
      integer(ip) :: ispos(8)

      nedge = 12
      
      ispos(:) = pnods(PointerToChildren(:))-1
      
      !first edge (node 1 2 9)
      ledger(1,1) = lnods(ispos(1)+1)
      ledger(2,1) = lnods(ispos(5)+2)
      ledger(3,1) = lnods(ispos(1)+2)
      
      !second edge (node 1 4 10)
      ledger(1,2) = lnods(ispos(1)+1)
      ledger(2,2) = lnods(ispos(3)+4)
      ledger(3,2) = lnods(ispos(1)+4)
      
      !third edge (1 5 11)
      ledger(1,3) = lnods(ispos(1)+1)
      ledger(2,3) = lnods(ispos(2)+5)
      ledger(3,3) = lnods(ispos(1)+5)
      
      !fourth edge (2 3 12)
      ledger(1,4) = lnods(ispos(5)+2)
      ledger(2,4) = lnods(ispos(7)+3)
      ledger(3,4) = lnods(ispos(5)+3)

      !fifth edge (2 6 13)
      ledger(1,5) = lnods(ispos(5)+2)
      ledger(2,5) = lnods(ispos(6)+6)
      ledger(3,5) = lnods(ispos(5)+6)
      
      !6 edge (3 4 14
      ledger(1,6) = lnods(ispos(7)+3)
      ledger(2,6) = lnods(ispos(3)+4)
      ledger(3,6) = lnods(ispos(7)+4)
      
      !7 edge (3 7 15)
      ledger(1,7) = lnods(ispos(7)+3)
      ledger(2,7) = lnods(ispos(8)+7)
      ledger(3,7) = lnods(ispos(7)+7)
      
      !8 (4 8 16)
      ledger(1,8) = lnods(ispos(3)+4)
      ledger(2,8) = lnods(ispos(4)+8)
      ledger(3,8) = lnods(ispos(3)+8)
      
      !9 (5 6 17)
      ledger(1,9) = lnods(ispos(2)+5)
      ledger(2,9) = lnods(ispos(6)+6)
      ledger(3,9) = lnods(ispos(2)+6)
      
      !10 (5 8 18)
      ledger(1,10) = lnods(ispos(2)+5)
      ledger(2,10) = lnods(ispos(4)+8)
      ledger(3,10) = lnods(ispos(2)+8)
      
      !11 edge (6 7 19)
      ledger(1,11) = lnods(ispos(6)+6)
      ledger(2,11) = lnods(ispos(8)+7)
      ledger(3,11) = lnods(ispos(6)+7)
      
      !12 edge (7 8 20)
      ledger(1,12) = lnods(ispos(8)+7)
      ledger(2,12) = lnods(ispos(4)+8)
      ledger(3,12) = lnods(ispos(8)+8)
      
      

   end subroutine  

   subroutine TrilinearHexa_face_from_el(nface,pnodb,lnods,iface,sort,facnod)
      use typre
      implicit none

      integer(ip) :: iface,nface,pnodb
      integer(ip) :: facnod(pnodb,nface),lnods(*)
      integer(ip) :: sort
      
      integer(ip) :: faces(4,6), auxfaces(4),jface
      
      faces = reshape((/ 1,4,8,5, &
                         2,6,7,3, &
                         1,5,6,2, &
                         4,3,7,8, &
                         1,2,4,3, &
                         5,8,7,6  &
                     /),shape(faces))

      if (iface == 0) then
         do jface = 1,6
            auxfaces = faces(:,jface)
            facnod(:,jface) = lnods(auxfaces)
         enddo
         if (sort == 1) then
            do jface = 1,6
               call sortb(4,facnod(1,jface))
            enddo   
         endif
      else
         auxfaces = faces(:,iface)
         facnod(:,1) = lnods(auxfaces)
         if (sort == 1) call sortb(4,facnod(1,1))     
      endif   
   
   end subroutine

   subroutine TrilinearHexa_up_faces(nface,ison,jface)
      !this subroutine gives the new neighbour face in the parent when an element is gone
      use typre
      implicit none
      integer(ip):: ielty,nface,jface,ison
      select case (ison)
         case (1)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 3) then
               jface = 3
            elseif (jface == 5) then
               jface = 5   
            endif
         case (2)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 3) then
               jface = 3
            elseif (jface == 6) then
               jface = 6   
            endif
         case (3)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 4) then
               jface = 4
            elseif (jface == 5) then
               jface = 5   
            endif
         case (4)
            if (jface == 1) then
               jface = 1 
            elseif (jface == 4) then
               jface = 4
            elseif (jface == 6) then
               jface = 6   
            endif
         case (5)
            if (jface == 2) then
               jface = 2 
            elseif (jface == 3) then
               jface = 3
            elseif (jface == 5) then
               jface = 5   
            endif
         case (6)
            if (jface == 2) then
               jface = 2 
            elseif (jface == 3) then
               jface = 3
            elseif (jface == 6) then
               jface = 6   
            endif
         case (7)
            if (jface == 2) then
               jface = 2 
            elseif (jface == 4) then
               jface = 4
            elseif (jface == 5) then
               jface = 5   
            endif
         case (8)
            if (jface == 2) then
               jface = 2 
            elseif (jface == 4) then
               jface = 4
            elseif (jface == 6) then
               jface = 6   
            endif
         end select
      
   end subroutine

   subroutine TrilinearHexa_down_faces(jface,facchild)
   use typre
   implicit none
   integer(ip):: jface,facchild(2,*)
   !this subroutine gives the neighbour face in the child for a given face in the parent
      select case (jface)
      case (1)
         facchild(1,1) = 1    !first the son, then the face
         facchild(2,1) = 1 
         facchild(1,2) = 2
         facchild(2,2) = 1 
         facchild(1,3) = 3
         facchild(2,3) = 1 
         facchild(1,4) = 4
         facchild(2,4) = 1 
      case (2)
         facchild(1,1) = 5
         facchild(2,1) = 2 
         facchild(1,2) = 6
         facchild(2,2) = 2 
         facchild(1,3) = 7
         facchild(2,3) = 2 
         facchild(1,4) = 8
         facchild(2,4) = 2 
      case (3)
         facchild(1,1) = 1
         facchild(2,1) = 3 
         facchild(1,2) = 2
         facchild(2,2) = 3 
         facchild(1,3) = 5
         facchild(2,3) = 3 
         facchild(1,4) = 6
         facchild(2,4) = 3 
      case (4)
         facchild(1,1) = 3
         facchild(2,1) = 4
         facchild(1,2) = 4
         facchild(2,2) = 4
         facchild(1,3) = 7
         facchild(2,3) = 4
         facchild(1,4) = 8
         facchild(2,4) = 4
      case (5)
         facchild(1,1) = 1
         facchild(2,1) = 5 
         facchild(1,2) = 3
         facchild(2,2) = 5 
         facchild(1,3) = 5
         facchild(2,3) = 5 
         facchild(1,4) = 7
         facchild(2,4) = 5 
      case (6)
         facchild(1,1) = 2
         facchild(2,1) = 6 
         facchild(1,2) = 4
         facchild(2,2) = 6 
         facchild(1,3) = 6
         facchild(2,3) = 6 
         facchild(1,4) = 8
         facchild(2,4) = 6 
      end select

   end subroutine

   subroutine TrilinearHexa_nodes_to_face(pnodb,lnodb,iface)
      use typre
      implicit none
      integer(ip) :: ndime,pnode,pnodb,lnodb(pnodb),iface
      integer(ip) :: auxlnodb(4)

      auxlnodb = lnodb
      call sortb(4,auxlnodb)
      
      if (auxlnodb(1) == 1 .and. auxlnodb(2) == 4 .and. auxlnodb(3) == 5 .and. auxlnodb(4) == 8) then
         iface = 1
      elseif (auxlnodb(1) == 2 .and. auxlnodb(2) == 3 .and. auxlnodb(3) == 6 .and. auxlnodb(4) == 7) then
         iface = 2
      elseif (auxlnodb(1) == 1 .and. auxlnodb(2) == 2.and. auxlnodb(3) == 5 .and. auxlnodb(4) == 6) then
         iface = 3
      elseif (auxlnodb(1) == 3 .and. auxlnodb(2) == 4 .and. auxlnodb(3) == 7 .and. auxlnodb(4) == 8) then
         iface = 4
      elseif (auxlnodb(1) == 1 .and. auxlnodb(2) == 2 .and. auxlnodb(3) == 3 .and. auxlnodb(4) == 4) then
         iface = 5
      elseif (auxlnodb(1) == 5 .and. auxlnodb(2) == 6 .and. auxlnodb(3) == 7 .and. auxlnodb(4) == 8) then
         iface = 6   
      endif   

   end subroutine

   subroutine TrilinearHexa_ConnectInteriorFaces(jfpos,son,lfacs)
      use typre
      implicit none
      integer(ip) :: jfpos(*), son(*), lfacs(2,*)
      
      !interior faces
      !Elements 1 and 5
      lfacs(1,jfpos(1)+2) = son(5)
      lfacs(2,jfpos(1)+2) = 1
      
      lfacs(1,jfpos(5)+1) = son(1)
      lfacs(2,jfpos(5)+1) = 2
      
      !Elements 1 and 2           
      lfacs(1,jfpos(1)+6) = son(2)
      lfacs(2,jfpos(1)+6) = 5
      
      lfacs(1,jfpos(2)+5) = son(1)
      lfacs(2,jfpos(2)+5) = 6
      
      !Elements  1 and 3
      lfacs(1,jfpos(1)+4) = son(3)
      lfacs(2,jfpos(1)+4) = 3
      
      lfacs(1,jfpos(3)+3) = son(1)
      lfacs(2,jfpos(3)+3) = 4
      
      !Elements 2 and 4
      lfacs(1,jfpos(2)+4) = son(4)
      lfacs(2,jfpos(2)+4) = 3
      
      lfacs(1,jfpos(4)+3) = son(2)
      lfacs(2,jfpos(4)+3) = 4
      
      !Elements 2 and 6
      lfacs(1,jfpos(2)+2) = son(6)
      lfacs(2,jfpos(2)+2) = 1
      
      lfacs(1,jfpos(6)+1) = son(2)
      lfacs(2,jfpos(6)+1) = 2
      
      !Elements 3 and 4
      lfacs(1,jfpos(3)+6) = son(4)
      lfacs(2,jfpos(3)+6) = 5
      
      lfacs(1,jfpos(4)+5) = son(3)
      lfacs(2,jfpos(4)+5) = 6
      
      !Elements 3 and 7
      lfacs(1,jfpos(3)+2) = son(7)
      lfacs(2,jfpos(3)+2) = 1
      
      lfacs(1,jfpos(7)+1) = son(3)
      lfacs(2,jfpos(7)+1) = 2
      
      !Elements 4 and 8
      lfacs(1,jfpos(4)+2) = son(8)
      lfacs(2,jfpos(4)+2) = 1
      
      lfacs(1,jfpos(8)+1) = son(4)
      lfacs(2,jfpos(8)+1) = 2
      
      !Elements 5 and 6
      lfacs(1,jfpos(5)+6) = son(6)
      lfacs(2,jfpos(5)+6) = 5
      
      lfacs(1,jfpos(6)+5) = son(5)
      lfacs(2,jfpos(6)+5) = 6
      
      !Elements 5 and 7
      lfacs(1,jfpos(5)+4) = son(7)
      lfacs(2,jfpos(5)+4) = 3
      
      lfacs(1,jfpos(7)+3) = son(5)
      lfacs(2,jfpos(7)+3) = 4
      
      !Elements 6 and 8
      lfacs(1,jfpos(6)+4) = son(8)
      lfacs(2,jfpos(6)+4) = 3
      
      lfacs(1,jfpos(8)+3) = son(6)
      lfacs(2,jfpos(8)+3) = 4
      
      !Elements 7 and 8
      lfacs(1,jfpos(7)+6) = son(8)
      lfacs(2,jfpos(7)+6) = 5
      
      lfacs(1,jfpos(8)+5) = son(7)
      lfacs(2,jfpos(8)+5) = 6
      
   end subroutine


   subroutine TrilinearHexa_AdaptiveGetChildrenExternalFacesList(inchild,jjchild,jjnface,jjface)
      use typre
      implicit none
      integer(ip) :: inchild,jjchild(:), jjnface(:), jjface(:,:)
      
      !there are 3 children with faces shared with the parent element
      inchild = 8
      !These elements are child 1, 2 and 3
      jjchild(1) = 1
      jjchild(2) = 2
      jjchild(3) = 3
      jjchild(4) = 4
      jjchild(5) = 5
      jjchild(6) = 6
      jjchild(7) = 7
      jjchild(8) = 8
      !Each child has 2 external faces
      jjnface(1:8)  = 3
      !The external faces are the 2 first ones for all the children
      jjface(1,1) = 1
      jjface(2,1) = 3
      jjface(3,1) = 5
      
      jjface(1,2) = 1
      jjface(2,2) = 3
      jjface(3,2) = 6
      
      jjface(1,3) = 1
      jjface(2,3) = 4
      jjface(3,3) = 5
      
      jjface(1,4) = 1
      jjface(2,4) = 4
      jjface(3,4) = 6
      
      jjface(1,5) = 2
      jjface(2,5) = 3
      jjface(3,5) = 5
      
      jjface(1,6) = 2
      jjface(2,6) = 3
      jjface(3,6) = 6
      
      jjface(1,7) = 2
      jjface(2,7) = 4
      jjface(3,7) = 5
      
      jjface(1,8) = 2
      jjface(2,8) = 4
      jjface(3,8) = 6
   end subroutine
   
   subroutine TrilinearHexa_GetInterpolationList(ivariation,InterpList)
      use typre
      implicit none
      integer(ip) :: ivariation
      type(InterpolationList) :: InterpList
      
      integer(ip) :: auxlnods(8)    , facnods(4),ichild,iface,inode
      
      auxlnods = (/1_ip, 2_ip, 3_ip, 4_ip, 5_ip, 6_ip, 7_ip, 8_ip /)
      
      InterpList%n = 7
      
      !The in the faces ones
      iface = 1
      ichild = 1
      inode  = 8
      call AddFace
      
      iface = 2
      ichild = 5
      inode  = 7
      call AddFace
      
      iface = 3
      ichild = 1
      inode  = 6
      call AddFace
      
      iface = 4
      ichild = 7
      inode  = 8
      call AddFace
      
      iface = 5
      ichild = 1
      inode  = 3
      call AddFace
      
      iface = 6
      ichild = 6
      inode  = 8
      call AddFace
      
      
      
      !The interior one
      InterpList%childElementList(7) = 1
      InterpList%childObjectiveNodeList(7) = 7
      InterpList%nInterpolatingNodes(7) = 8
      InterpList%InterpolatingNodesList(1:8,7) = (/1_ip, 2_ip, 3_ip, 4_ip, 5_ip, 6_ip, 7_ip, 8_ip /)
      InterpList%InterpolatingCoefficientsList(1:8,7) = 0.125_rp
      
contains
      subroutine AddFace
         implicit none
         
         InterpList%childElementList(iface) =   ichild
         InterpList%childObjectiveNodeList(iface) = inode
         InterpList%nInterpolatingNodes(iface) = 4
         call TrilinearHexa_face_from_el(6_ip,4_ip,auxlnods,iface,0,InterpList%InterpolatingNodesList(1:4,iface))
         InterpList%InterpolatingCoefficientsList(1:4,iface) = 0.25_rp
      end subroutine
   
   end subroutine
   
   subroutine TrilinearHexa_GetEdgeNodesForInode(inode,nnode,nodelist)
      use typre
      implicit none
      integer(ip) :: inode, nnode,nodelist(3)
      
      nnode = 3
      
      if (inode == 1) then
         nodelist(1) = 2
         nodelist(2) = 4
         nodelist(3) = 5
      elseif (inode == 2) then
         nodelist(1) = 1
         nodelist(2) = 3
         nodelist(3) = 6
      elseif (inode == 3) then
         nodelist(1) = 2
         nodelist(2) = 4
         nodelist(3) = 7
      elseif (inode == 4) then
         nodelist(1) = 1
         nodelist(2) = 3
         nodelist(3) = 8   
      elseif (inode == 5) then
         nodelist(1) = 1
         nodelist(2) = 6
         nodelist(3) = 8
      elseif (inode == 6) then
         nodelist(1) = 2
         nodelist(2) = 5
         nodelist(3) = 7      
      elseif (inode == 7) then
         nodelist(1) = 3
         nodelist(2) = 6
         nodelist(3) = 8     
      elseif (inode == 8) then
         nodelist(1) = 4
         nodelist(2) = 5
         nodelist(3) = 7        
      endif
   end subroutine
   
   subroutine TrilinearHexa_GetNumberOfInteriorNodesToAdd(nInteriorNodes)
      use typre
      implicit none
      integer(ip) :: nInteriorNodes
      
      nInteriorNodes = 1
      
   end subroutine
   
   subroutine TrilinearHexa_GetInteriorPointFromParent(inode,jchild,jnode)
      use typre
      implicit none
      integer(ip) :: inode,jchild,jnode
      
      if (inode == 1) then
         jchild = 1
         jnode = 7
      else
         call runend('Wrong number of interior nodes')
      endif
   end subroutine   
   
   subroutine TrilinearHexa_AddNodeToChildInterior(ivariation,iadd,newpo,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,iadd,newpo,PointerToChildren(*),pnods(*),lnods(*)
      integer(ip) :: ispos(8), ichild,inode
      
      !Only one interior node needs to be added
      if (iadd == 1) then
         ichild = 1
         inode = 7
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 2
         inode = 3
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 3
         inode = 6
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 4
         inode = 2
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 5
         inode = 8
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 6
         inode = 4
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo

         ichild = 7
         inode = 5
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
         
         ichild = 8
         inode = 1
         ispos(ichild) = pnods(PointerToChildren(ichild))-1
         lnods(ispos(ichild)+inode) = newpo
      endif
   end subroutine    
   
   subroutine TrilinearHexa_GetNedge(nedge)
      use typre
      implicit none
      integer(ip) :: nedge
      
      nedge = 12
   end subroutine
   
   subroutine TrilinearHexa_GetEdges(iedge,lnodedge)
      use typre
      implicit none
      integer(ip) :: ielty,iedge,lnodedge(*)
      
      select case(iedge)
         case(1)
            lnodedge(1) = 1
            lnodedge(2) = 2
         case(2)
            lnodedge(1) = 1
            lnodedge(2) = 4
         case(3)
            lnodedge(1) = 1
            lnodedge(2) = 5
         case(4)
            lnodedge(1) = 2
            lnodedge(2) = 3
         case(5)
            lnodedge(1) = 2
            lnodedge(2) = 6
         case(6)
            lnodedge(1) = 3
            lnodedge(2) = 4
         case(7)
            lnodedge(1) = 3
            lnodedge(2) = 7
         case(8)
            lnodedge(1) = 4
            lnodedge(2) = 8   
         case(9)
            lnodedge(1) = 5
            lnodedge(2) = 6   
         case(10)
            lnodedge(1) = 5
            lnodedge(2) = 8   
         case(11)
            lnodedge(1) = 6
            lnodedge(2) = 7   
         case(12)
            lnodedge(1) = 7
            lnodedge(2) = 8   
            
      end select
      
   end subroutine
   
   subroutine TrilinearHexa_Up_Edges(ison,childedge,parentedge)
      use typre
      implicit none
      integer(ip) :: ielty,ison, childedge,parentedge

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
            case(3)
               parentedge = 3
            case(9)
               parentedge = 9
            case(10)
               parentedge = 10
            case default
               parentedge = -1
         end select
      elseif (ison == 3) then
         select case(childedge)
            case(2)
               parentedge = 2
            case(6)
               parentedge = 6
            case(8)
               parentedge = 8
            case default
               parentedge = -1
         end select
       elseif (ison == 4) then
         select case(childedge)
            case(10)
               parentedge = 10
            case(8)
               parentedge = 8
            case(12)
               parentedge = 12
            case default
               parentedge = -1
         end select   
      elseif (ison == 5) then
         select case(childedge)
            case(1)
               parentedge = 1
            case(4)
               parentedge = 4
            case(5)
               parentedge = 5
            case default
               parentedge = -1
         end select   
      elseif (ison == 6) then
         select case(childedge)
            case(9)
               parentedge = 9
            case(5)
               parentedge = 5
            case(11)
               parentedge = 11
            case default
               parentedge = -1
         end select   
      elseif (ison == 7) then
         select case(childedge)
            case(4)
               parentedge = 4
            case(7)
               parentedge = 7
            case(6)
               parentedge = 6
            case default
               parentedge = -1
         end select   
      elseif (ison == 8) then
         select case(childedge)
            case(7)
               parentedge = 7
            case(11)
               parentedge = 11
            case(12)
               parentedge = 12
            case default
               parentedge = -1
         end select      
      endif   
   end subroutine
   
  
   
   
   
   subroutine TrilinearHexa_GetNumberofInTheFaceNodesToAdd(nInTheFaceNodes)
      use typre
      implicit none
      integer(ip) :: nInTheFaceNodes
      
      nInTheFaceNodes = 6
      
   end subroutine
   
   subroutine TrilinearHexa_GetFaceForInTheFaceNode(inode,iface)
      use typre
      implicit none
      integer(ip) :: ielty,inode,iface
      
      iface = inode
   end subroutine
   
   subroutine TrilinearHexa_GetInTheFacePointFromParent(iface,kchild,knode)
      use typre
      implicit none
      integer(ip) :: iface,kchild,knode
      
      select case (iface)
         case(1)
            kchild = 1
            knode = 8
         case(2)
            kchild = 5
            knode = 7
         case(3)
            kchild = 1
            knode = 6
         case(4)
            kchild = 3
            knode = 7
         case(5)
            kchild = 1
            knode = 3
         case(6)
            kchild = 2
            knode = 7
      end select
      
   end subroutine
   
   subroutine TrilinearHexa_AddNodeToChildInTheFace(ivariation,inode,kpoin,PointerToChildren,pnods,lnods)
      use typre
      implicit none
      integer(ip) :: ivariation,inode,kpoin,PointerToChildren(*), pnods(*),lnods(*)
      
      integer(ip) :: ielem,ispos
      
      select case (inode)
         case(1)
            ielem = PointerToChildren(1)
            ispos = pnods(ielem)-1
            lnods(ispos+8) = kpoin
            
            ielem = PointerToChildren(2)
            ispos = pnods(ielem)-1
            lnods(ispos+4) = kpoin
            
            ielem = PointerToChildren(3)
            ispos = pnods(ielem)-1
            lnods(ispos+5) = kpoin
            
            ielem = PointerToChildren(4)
            ispos = pnods(ielem)-1
            lnods(ispos+1) = kpoin
         case(2)
            ielem = PointerToChildren(5)
            ispos = pnods(ielem)-1
            lnods(ispos+7) = kpoin
            
            ielem = PointerToChildren(6)
            ispos = pnods(ielem)-1
            lnods(ispos+3) = kpoin
            
            ielem = PointerToChildren(7)
            ispos = pnods(ielem)-1
            lnods(ispos+6) = kpoin
            
            ielem = PointerToChildren(8)
            ispos = pnods(ielem)-1
            lnods(ispos+2) = kpoin
         case(3)
            ielem = PointerToChildren(1)
            ispos = pnods(ielem)-1
            lnods(ispos+6) = kpoin
            
            ielem = PointerToChildren(2)
            ispos = pnods(ielem)-1
            lnods(ispos+2) = kpoin
            
            ielem = PointerToChildren(5)
            ispos = pnods(ielem)-1
            lnods(ispos+5) = kpoin
            
            ielem = PointerToChildren(6)
            ispos = pnods(ielem)-1
            lnods(ispos+1) = kpoin
         case(4)
            ielem = PointerToChildren(3)
            ispos = pnods(ielem)-1
            lnods(ispos+7) = kpoin
            
            ielem = PointerToChildren(4)
            ispos = pnods(ielem)-1
            lnods(ispos+3) = kpoin
            
            ielem = PointerToChildren(7)
            ispos = pnods(ielem)-1
            lnods(ispos+8) = kpoin
            
            ielem = PointerToChildren(8)
            ispos = pnods(ielem)-1
            lnods(ispos+4) = kpoin
         
         case(5)
            ielem = PointerToChildren(1)
            ispos = pnods(ielem)-1
            lnods(ispos+3) = kpoin
            
            ielem = PointerToChildren(3)
            ispos = pnods(ielem)-1
            lnods(ispos+2) = kpoin
            
            ielem = PointerToChildren(5)
            ispos = pnods(ielem)-1
            lnods(ispos+4) = kpoin
            
            ielem = PointerToChildren(7)
            ispos = pnods(ielem)-1
            lnods(ispos+1) = kpoin
         
         case(6)
            ielem = PointerToChildren(2)
            ispos = pnods(ielem)-1
            lnods(ispos+7) = kpoin
            
            ielem = PointerToChildren(4)
            ispos = pnods(ielem)-1
            lnods(ispos+6) = kpoin
            
            ielem = PointerToChildren(6)
            ispos = pnods(ielem)-1
            lnods(ispos+8) = kpoin
            
            ielem = PointerToChildren(8)
            ispos = pnods(ielem)-1
            lnods(ispos+5) = kpoin
         
      end select
   end subroutine
   
   subroutine TrilinearHexa_ParentToChildrenElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,8),NewNodalArray(ndofn,8)
      integer(ip) :: ielty,ivariation,ichild
   
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
         NewNodalArray(:,5) = (OldNodalArray(:,1)+OldNodalArray(:,5))*0.5
         NewNodalArray(:,6) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,5)+OldNodalArray(:,6))*0.25
         NewNodalArray(:,7) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,8) = (OldNodalArray(:,1)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,8))*0.25
         
      elseif (ichild == 2) then
         NewNodalArray(:,1) =  (OldNodalArray(:,1)+OldNodalArray(:,5))*0.5
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,5)+OldNodalArray(:,6))*0.25
         NewNodalArray(:,3) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,5) = OldNodalArray(:,5)
         NewNodalArray(:,6) =  (OldNodalArray(:,5)+OldNodalArray(:,6))*0.5
         NewNodalArray(:,7) = (OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,8) =  (OldNodalArray(:,5)+OldNodalArray(:,8))*0.5
         
      elseif (ichild == 3) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,4))*0.5
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
         NewNodalArray(:,4) = OldNodalArray(:,4)
         NewNodalArray(:,5) = (OldNodalArray(:,1)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,6) =  (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,7) = (OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,8) =  (OldNodalArray(:,4)+OldNodalArray(:,8))*0.5
      elseif (ichild == 4) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,2) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,4) = (OldNodalArray(:,4)+OldNodalArray(:,8))*0.5
         NewNodalArray(:,5) = (OldNodalArray(:,5)+OldNodalArray(:,8))*0.5
         NewNodalArray(:,6) = (OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,7) = (OldNodalArray(:,7)+OldNodalArray(:,8))*0.5
         NewNodalArray(:,8) = OldNodalArray(:,8)
      elseif (ichild == 5) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,2))*0.5
         NewNodalArray(:,2) = OldNodalArray(:,2)
         NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,5) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,5)+OldNodalArray(:,6))*0.25
         NewNodalArray(:,6) = (OldNodalArray(:,2)+OldNodalArray(:,6))*0.5
         NewNodalArray(:,7) = (OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,6)+OldNodalArray(:,7))*0.25
         NewNodalArray(:,8) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
      
      elseif (ichild == 6) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,5)+OldNodalArray(:,6))*0.25
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,6))*0.5
         NewNodalArray(:,3) = (OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,6)+OldNodalArray(:,7))*0.25
         NewNodalArray(:,4) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,5) = (OldNodalArray(:,5)+OldNodalArray(:,6))*0.5
         NewNodalArray(:,6) =  OldNodalArray(:,6)
         NewNodalArray(:,7) = (OldNodalArray(:,6)+OldNodalArray(:,7))*0.5
         NewNodalArray(:,8) = (OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         
      elseif (ichild == 7) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4))*0.25
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3))*0.5
         NewNodalArray(:,3) = OldNodalArray(:,3)
         NewNodalArray(:,4) = (OldNodalArray(:,3)+OldNodalArray(:,4))*0.5
         NewNodalArray(:,5) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,6) = (OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,6)+OldNodalArray(:,7))*0.25
         NewNodalArray(:,7) = (OldNodalArray(:,3)+OldNodalArray(:,7))*0.5
         NewNodalArray(:,8) = (OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
      elseif (ichild == 8) then
         NewNodalArray(:,1) = (OldNodalArray(:,1)+OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.125
         NewNodalArray(:,2) = (OldNodalArray(:,2)+OldNodalArray(:,3)+OldNodalArray(:,6)+OldNodalArray(:,7))*0.25
         NewNodalArray(:,3) = (OldNodalArray(:,3)+OldNodalArray(:,7))*0.5
         NewNodalArray(:,4) = (OldNodalArray(:,3)+OldNodalArray(:,4)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,5) = (OldNodalArray(:,5)+OldNodalArray(:,6)+OldNodalArray(:,7)+OldNodalArray(:,8))*0.25
         NewNodalArray(:,6) = (OldNodalArray(:,6)+OldNodalArray(:,7))*0.5
         NewNodalArray(:,7) = OldNodalArray(:,7)
         NewNodalArray(:,8) = (OldNodalArray(:,7)+OldNodalArray(:,8))*0.5
      endif
   
   end subroutine
   
   subroutine TrilinearHexa_ChildrenToParentElement(ndofn,OldNodalArray,NewNodalArray,ivariation,ichild)
      integer(ip) :: ndofn
      real(rp) :: OldNodalArray(ndofn,8),NewNodalArray(ndofn,8)
      integer(ip) :: ivariation,ichild
      
      if (ichild == 1) then
         NewNodalArray(:,1) = OldNodalArray(:,1)
      elseif (ichild == 2) then
         NewNodalArray(:,5) = OldNodalArray(:,5)
      elseif (ichild == 3) then
         NewNodalArray(:,4) = OldNodalArray(:,4)
      elseif (ichild == 4) then
         NewNodalArray(:,8) = OldNodalArray(:,8)
      elseif (ichild == 5) then
         NewNodalArray(:,2) = OldNodalArray(:,2)
      elseif (ichild == 6) then
         NewNodalArray(:,6) = OldNodalArray(:,6)
      elseif (ichild == 7) then
         NewNodalArray(:,3) = OldNodalArray(:,3)
      elseif (ichild == 8) then
         NewNodalArray(:,7) = OldNodalArray(:,7)   
      endif
   
   end subroutine

      
      
end module
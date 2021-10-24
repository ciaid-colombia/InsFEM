subroutine FindElement(a,pnode,lnode,ielem)
!-----------------------------------------------------------------------
!    This routine looks for ielem in in the list lnode using
!    a%pbopo and a%lbopo
!-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: pnode,lnode(pnode)
   integer(ip) :: ielem


   integer(ip) :: kelem,elemk,kpnode,inode,ipoin
   integer(ip), pointer :: klnode(:) => NULL()
   
   
   ielem=0
   ipoin = lnode(1)
   elementsLoop : do elemk = a%pelpo(ipoin),a%pelpo(ipoin+1)-1
      kelem = a%lelpo(elemk)
      
      call a%GetLnode(kelem,kpnode,klnode)
      if (pnode /= kpnode) cycle elementsLoop
      do inode = 1,pnode
         if (lnode(inode) /= klnode(inode)) cycle elementsLoop
      enddo
      ielem = kelem
      exit elementsLoop
   enddo elementsLoop


end subroutine FindElement

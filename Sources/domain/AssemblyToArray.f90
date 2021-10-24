subroutine AssemblyToArray(a,e,ndofn,elrhs,array)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   class(FiniteElement) :: e
   integer(ip) :: ndofn
   real(rp) :: elrhs(ndofn,*),array(ndofn,*)
   
   integer(ip) :: inode,ipoin,phang,jspos,jpoin,jhang
   real(rp) :: jweight
   
   !Usual
   if (a%kfl_HangingNodes .eqv. .false.) then
      !Do not vectorize since they can repeat!!
      do inode = 1,e%pnode
         array(:,e%lnods(inode)) = array(:,e%lnods(inode)) + elrhs(:,inode)
      enddo
      !Do not vectorize since they can repeat!
      !array(:,e%lnods(1:e%pnode)) = array(:,e%lnods(1:e%pnode)) + elrhs(:,1:e%pnode)
   
   !Hanging nodes
   else
      call a%AssemblyToArrayHanging(e,ndofn,elrhs,array)
   endif
end subroutine

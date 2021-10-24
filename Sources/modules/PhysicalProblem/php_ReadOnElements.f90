subroutine php_InitializeOnElementsConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
end subroutine   
   
! subroutine php_ReadOnElements(a,lnode)
!    use typre
!    use Mod_PhysicalProblem
!    implicit none
!    class(PhysicalProblem) :: a
!    integer(ip) :: lnode(:)
!    
!    integer(ip) :: ielem,pnode
!    
!    !Read element nodes
!    pnode=int(a%Listener%param(2))
!    lnode(1:pnode)=int(a%Listener%param(3:2+pnode))
!    
!    !To local numbering
!    call a%Mesh%Global2Local(pnode,lnode(1:pnode),lnode(1:pnode))
!    
!    !Find which is the boundary
!    call a%Mesh%FindElement(pnode,lnode(1:pnode),ielem)    
!    if(ielem==0) call runend('php_ReadOnElements: Element not found')
!    
!    !On Elements specific
!    a%gielem = ielem
!    a%gipsta=2+pnode
!    call a%SpecificReadOnElements
!    
!       
! end subroutine

subroutine php_ReadOnElements(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: ielem
   
   !Which is the element in local numbering
   ielem = a%Listener%param(1)
   call a%Mesh%ElementLocal2Initial%Global2Local(ielem,ielem)
   
   
   
   !On Elements specific
   a%gielem = ielem
   a%gipsta=1
   call a%SpecificReadOnElements
   
      
end subroutine


subroutine php_FinalizeOnElements(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
end subroutine

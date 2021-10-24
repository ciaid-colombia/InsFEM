module Mod_RpMeshAllocator
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   implicit none
   
contains

   subroutine AllocR1P(Memor,Mesh,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r1p), allocatable :: r1p_array(:)
      
      
      class(FiniteElement), pointer :: e => NULL() 
      integer(ip) :: ielem,nelem
      
      call Mesh%GetNelem(nelem)
      call Memor%alloc(nelem,r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         call Memor%palloc(e%pgaus,r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      call Mesh%ElementdeAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      
   end subroutine
   
   subroutine DeAllocR1P(Memor,Mesh,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r1p), allocatable :: r1p_array(:)
      
      
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ielem,nelem
      
      call Mesh%GetNelem(nelem)
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         call Memor%pdealloc(e%pgaus,r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      call Mesh%ElementdeAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      
      call Memor%dealloc(nelem,r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      
   end subroutine
   
   subroutine AllocR2P(Memor,Mesh,ndime1,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r2p), allocatable :: r1p_array(:)
      integer(ip) :: ndime1
      
      
      class(FiniteElement), pointer :: e => NULL() 
      integer(ip) :: ielem,nelem
      
      call Mesh%GetNelem(nelem)
      call Memor%alloc(nelem,r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         call Memor%palloc(ndime1,e%pgaus,r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      call Mesh%ElementdeAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      
   end subroutine
   
   subroutine DeAllocR2P(Memor,Mesh,ndime1,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r2p), allocatable :: r1p_array(:)
      integer(ip) :: ndime1
      
      
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ielem,nelem
      
      call Mesh%GetNelem(nelem)
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         call Memor%pdealloc(ndime1,e%pgaus,r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      call Mesh%ElementdeAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      
      call Memor%dealloc(nelem,r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      
   end subroutine
   
   subroutine AllocR3P(Memor,Mesh,ndime1,ndime2,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r3p), allocatable :: r1p_array(:)
      integer(ip) :: ndime1,ndime2
      
      
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ielem,nelem
      
      call Mesh%GetNelem(nelem)
      call Memor%alloc(nelem,r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         call Memor%palloc(ndime1,ndime2,e%pgaus,r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      call Mesh%ElementdeAlloc(e,Memor,'ForceClosedRule','plcd_memall')
      
   end subroutine
   
   subroutine DeAllocR3P(Memor,Mesh,ndime1,ndime2,r1p_array)
      type(MemoryMan) :: Memor
      class(FemMesh) :: Mesh
      type(r3p), allocatable :: r1p_array(:)
      integer(ip) :: ndime1,ndime2
      
      
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ielem,nelem
      
      elements : do ielem = 1,size(r1p_array,1)
         call Memor%pdealloc(size(r1p_array(ielem)%a,1),size(r1p_array(ielem)%a,2),size(r1p_array(ielem)%a,3),r1p_array(ielem)%a,'DPostprocess','plcd_memall')
      enddo elements
      
      call Memor%dealloc(size(r1p_array,1),r1p_array,'DPostprocess','MD_SetupPostprocess')
      
      
   end subroutine

   
end module Mod_RpMeshAllocator


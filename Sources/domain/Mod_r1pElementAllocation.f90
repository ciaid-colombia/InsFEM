module Mod_r1pElementAllocation
   use typre
   use Mod_Memor
   use Mod_Mesh
   implicit none

contains
   
   subroutine Allocater1pElement(Mesh,r1pArray,Memor)
      use typre
      implicit none
      class(FemMesh) :: Mesh
      type(r1p), allocatable     :: r1pArray(:)
      type(MemoryMan) :: Memor
      
      integer(ip) :: nelem,ielem,pnode
      integer(ip) :: counter, pgaus
      
      call Mesh%GetNelem(nelem)
      call Memor%alloc(nelem,r1pArray,'r1pArray','Allocater1pElement')
      counter = 0
      
      do ielem=1,nelem
         call Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(r1pArray(ielem)%a(pgaus))
         r1pArray(ielem)%a = 0.0_rp
         counter = counter + pgaus
      end do
      call Memor%allocObj(0,'r1pArray%a','Allocater1pElement',counter*rp)
      
   end subroutine
   
   
   subroutine DeAllocater1pElement(Mesh,r1pArray,Memor)
      use typre
      implicit none
      class(FemMesh) :: Mesh
      type(r1p), allocatable     :: r1pArray(:)
      type(MemoryMan) :: Memor
      
      integer(ip) :: nelem,ielem,pnode
      integer(ip) :: counter, pgaus
      
      counter = 0
      
      do ielem=1,size(r1pArray)
         pgaus = size(r1pArray(ielem)%a)
         deallocate(r1pArray(ielem)%a)
         counter = counter + pgaus
      end do
      call Memor%deallocObj(0,'r1pArray%a','Allocater1pElement',counter*rp)
      call Memor%dealloc(size(r1pArray),r1pArray,'r1pArray','Allocater1pElement')
   end subroutine
   
   
   subroutine Allocater2pElement(Mesh,ndime,r1pArray,Memor,string)
      use typre
      implicit none
      class(FemMesh) :: Mesh
      type(r2p), allocatable     :: r1pArray(:)
      type(MemoryMan) :: Memor
      character(*) :: string
      
      integer(ip) :: nelem,ielem,pnode,ndime
      integer(ip) :: counter, pgaus
      
      call Mesh%GetNelem(nelem)
      call Memor%alloc(nelem,r1pArray,string,'Allocater1pElement')
      counter = 0
      
      do ielem=1,nelem
         call Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(r1pArray(ielem)%a(ndime,pgaus))
         r1pArray(ielem)%a = 0.0_rp
         counter = counter + pgaus*ndime
      end do
      call Memor%allocObj(0,string,'Allocater1pElement',counter*rp)
      
   end subroutine
   
   
   subroutine DeAllocater2pElement(Mesh,ndime,r1pArray,Memor,string)
      use typre
      implicit none
      class(FemMesh) :: Mesh
      type(r2p), allocatable     :: r1pArray(:)
      type(MemoryMan) :: Memor
      character(*) :: string
      
      integer(ip) :: nelem,ielem,pnode,ndime
      integer(ip) :: counter, pgaus
      
      
      counter = 0
      
      do ielem=1,size(r1pArray)
         pgaus = size(r1pArray(ielem)%a,2)
         ndime = size(r1pArray(ielem)%a,1)
         deallocate(r1pArray(ielem)%a)
         counter = counter + pgaus*ndime
      end do
      call Memor%deallocObj(0,string,'Allocater1pElement',counter*rp)
      call Memor%dealloc(size(r1pArray),r1pArray,string,'Allocater1pElement')
   end subroutine
end module

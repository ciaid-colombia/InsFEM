module Mod_BroadCastBuffer
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_MPIObject
   implicit none
   private
   public BroadCastBuffer
   
   type :: i1pscalar
      integer(ip), pointer :: l => NULL()
   end type
   type :: r1pscalar
      real(rp), pointer :: a => NULL()
   end type
   type :: l1pscalar
      logical, pointer :: l => NULL()
   end type
   type :: c1p
      character, pointer :: l => NULL()
      integer(ip) :: nchar
   end type 
   type :: iap 
      integer(ip), pointer :: l(:) => NULL()
      integer(ip) :: ninte
   end type
   type :: iap2 
      integer(ip), pointer :: l(:,:) => NULL()
      integer(ip) :: ninte
   end type
   type :: rap 
      real(rp),    pointer :: a(:) => NULL()
      integer(ip) :: nreal
   end type
   type :: rap2
      real(rp),    pointer :: a(:,:) => NULL()
      integer(ip) :: nreal
   end type

   type, extends(MPIObject) :: BroadCastBuffer
   
      integer(ip) :: nz = 0
      integer(ip), allocatable :: ia(:)
      character(1), allocatable :: ja(:)
      integer(ip), allocatable :: SendType(:)
   
      type(iap), allocatable :: ExternalIPointers(:)         
      type(iap2), allocatable :: ExternalI2Pointers(:) 
      type(rap), allocatable :: ExternalRPointers(:)
      type(rap2), allocatable :: ExternalR2Pointers(:)
      type(l1p), allocatable :: ExternalLPointers(:)
      type(c1p), allocatable :: ExternalCPointers(:)
      
      type(i1pscalar), allocatable :: ExternalIPointersScalar(:)
      type(r1pscalar), allocatable :: ExternalRPointersScalar(:)
      type(l1pscalar), allocatable :: ExternalLPointersScalar(:)
   
contains
      procedure :: Initialize
      procedure :: AddInteger
      procedure :: AddLogical
      procedure :: AddReal
      procedure :: BroadCast
      procedure :: CheckSizes
      procedure :: Dealloc
      
      procedure :: AddInteger1
      procedure :: AddLogical1
      procedure :: AddReal1
      procedure :: AddCharacter
      procedure :: AddInteger2
      procedure :: AddReal2
      
      generic :: Add => AddInteger,AddInteger2,AddReal,AddReal2,AddLogical,AddInteger1,AddReal1,AddLogical1,AddCharacter
   end type

contains   
   subroutine Initialize(a,nBufferEntries,BufferSize)
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip) :: NBufferEntries,BufferSize
      
      a%nz = 0
      call a%Memor%alloc(nBufferEntries+1,a%ia,'ia','BroadCastBuffer')
      call a%Memor%alloc(rp*BufferSize,a%ja,'ja','BroadCastBuffer')
      allocate(a%ExternalRPointers(nBufferEntries))
      call a%Memor%allocObj(0,'ExternalRPointers','BroadCastBuffer',nBufferEntries)
      allocate(a%ExternalLPointers(nBufferEntries))
      call a%Memor%allocObj(0,'ExternalLPointers','BroadCastBuffer',nBufferEntries)
      call a%Memor%alloc(nBufferEntries,a%SendType,'SendType','BroadCastBuffer')
      
      allocate(a%ExternalIPointers(nBufferEntries))
      allocate(a%ExternalI2Pointers(nBufferEntries))
      allocate(a%ExternalCPointers(nBufferEntries))
      allocate(a%ExternalIPointersScalar(nBufferEntries))
      allocate(a%ExternalRPointersScalar(nBufferEntries))
      allocate(a%ExternalR2Pointers(nBufferEntries))
      allocate(a%ExternalLPointersScalar(nBufferEntries))
      call a%Memor%allocObj(0,'ExternalPointersScalar','BroadCastBuffer',5*nBufferEntries)  
      
      
      
      a%ia(1) = 1_ip
      
   end subroutine

   subroutine AddInteger(a,intarray)
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip), target :: intarray(:)
      
      integer(ip) :: nchar,ispos            
      
      nchar = size(intarray)*ip  

      call a%CheckSizes(nchar)     
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,intarray,a%ja(ispos+1:ispos+nchar))   
      
      !External Pointer
      a%ExternalIPointers(a%nz)%l => intarray
      a%ExternalIPointers(a%nz)%ninte = nchar
      a%SendType(a%nz) = 0
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine
   
   subroutine AddInteger2(a,intarray)
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip), target :: intarray(:,:)
      
      integer(ip) :: nchar,ispos            
      
      nchar = size(intarray)*ip  

      call a%CheckSizes(nchar)     
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,intarray,a%ja(ispos+1:ispos+nchar))   
      
      !External Pointer
      a%ExternalI2Pointers(a%nz)%l => intarray
      a%ExternalI2Pointers(a%nz)%ninte = nchar
      a%SendType(a%nz) = 7
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine   
   
   subroutine AddReal(a,realarray)
      implicit none
      class(BroadCastBuffer) :: a
      real(rp), target :: realarray(:)
      
      integer(ip) :: nchar,ispos            
      
      nchar = size(realarray)*rp  

      call a%CheckSizes(nchar)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,realarray,a%ja(ispos+1:ispos+nchar))   
      
      !External Pointer
      a%ExternalRPointers(a%nz)%a => realarray
      a%ExternalRPointers(a%nz)%nreal = nchar
      a%SendType(a%nz) = 1
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine
   
   subroutine AddReal2(a,realarray)
      implicit none
      class(BroadCastBuffer) :: a
      real(rp), target :: realarray(:,:)
      
      integer(ip) :: nchar,ispos            
      
      nchar = size(realarray)*rp  

      call a%CheckSizes(nchar)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,realarray,a%ja(ispos+1:ispos+nchar))   
      
      !External Pointer
      a%ExternalR2Pointers(a%nz)%a => realarray
      a%ExternalR2Pointers(a%nz)%nreal = nchar
      a%SendType(a%nz) = 8
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine   
   
   subroutine AddLogical(a,logicalarray)
      implicit none
      class(BroadCastBuffer) :: a
      logical, target :: logicalarray(:)
      
      integer(ip) :: i,ispos,n
      
      n = size(logicalarray)
      call a%CheckSizes(n)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      do i = 1,n
         if (logicalarray(i) .eqv. .true.) then
            a%ja(ispos+i) = "1"
         else
            a%ja(ispos+i) = "0"
         endif
      enddo
      
      !External Pointer
      a%ExternalLPointers(a%nz)%l => LogicalArray
      a%SendType(a%nz) = 2
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+n
   end subroutine  
   
   subroutine AddInteger1(a,integerarray)
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip), target :: integerarray
      
      integer(ip) :: nchar,ispos            
      
      nchar = ip  
      
      call a%CheckSizes(nchar)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,integerarray,a%ja(ispos+1:ispos+nchar))   
      
      !External Pointer
      a%ExternalIPointersScalar(a%nz)%l => integerarray
      a%SendType(a%nz) = 3
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine
   
   subroutine AddReal1(a,realarray)
      implicit none 
      class(BroadCastBuffer) :: a
      real(rp), target :: realarray
      
      integer(ip) :: nchar,ispos            
      
      nchar = rp  
      
      call a%CheckSizes(nchar)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(nchar,realarray,a%ja(ispos+1:ispos+nchar))
      
      !External Pointer
      a%ExternalRPointersScalar(a%nz)%a => realarray
      a%SendType(a%nz) = 4
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+nchar
   end subroutine
   
   subroutine AddLogical1(a,logicalarray)
      use typre
      implicit none
      class(BroadCastBuffer) :: a
      logical, target :: logicalarray
      
      integer(ip) :: ispos
      integer(ip), parameter :: n = 1
      
      call a%CheckSizes(n)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      if (logicalarray .eqv. .true.) then
         a%ja(ispos+1) = "1"
      else
         a%ja(ispos+1) = "0"
      endif
      
      !External Pointer
      a%ExternalLPointersScalar(a%nz)%l => logicalarray
      a%SendType(a%nz) = 5
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+n
   end subroutine
   
   subroutine AddCharacter(a,n,chararray)
      use typre
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip) :: n
      character, target :: chararray
      
      integer(ip) :: ispos
      
      
      call a%CheckSizes(n)
      
      a%nz = a%nz+1
      ispos = a%ia(a%nz)-1
      call CharacterCopy(n,chararray,a%ja(ispos+1))
      
      !External Pointer
      a%ExternalCPointers(a%nz)%l => chararray
      a%ExternalCPointers(a%nz)%nchar = n
      a%SendType(a%nz) = 6
      
      !Update next ia
      a%ia(a%nz+1) = a%ia(a%nz)+n
   end subroutine   
    
  
   
   subroutine BroadCast(a)
      use typre
      use MPI
      implicit none
      class(BroadCastBuffer) :: a
      
      integer(ip) :: i, ispos,pnode,inode
      
      integer(ip) :: ierr
      
      call MPI_BCAST(a%ja, a%ia(a%nz+1), MPI_CHARACTER, a%MPIroot, a%MPIcomm, ierr)
      
      !Split to each array
      do i = 1,a%nz
         ispos = a%ia(i)-1
         pnode = a%ia(i+1)-a%ia(i)
         if (a%SendType(i) == 0) then
            call CharacterCopy(a%ExternalIPointers(i)%ninte,a%ja(ispos+1:ispos+pnode),a%ExternalIPointers(i)%l)
         elseif (a%SendType(i) == 7) then
            call CharacterCopy(a%ExternalI2Pointers(i)%ninte,a%ja(ispos+1:ispos+pnode),a%ExternalI2Pointers(i)%l)   
         elseif(a%SendType(i) == 1) then
            call CharacterCopy(a%ExternalRPointers(i)%nreal,a%ja(ispos+1:ispos+pnode),a%ExternalRPointers(i)%a)
         elseif(a%SendType(i) == 8) then
            call CharacterCopy(a%ExternalR2Pointers(i)%nreal,a%ja(ispos+1:ispos+pnode),a%ExternalR2Pointers(i)%a)   
         elseif(a%SendType(i) == 2) then
            do inode = 1,pnode
               if (a%ja(ispos+inode) == "1") then
                  a%ExternalLPointers(i)%l(inode) = .true.
               else
                  a%ExternalLPointers(i)%l(inode) = .false.
               endif
            enddo
         elseif (a%SendType(i) == 3) then
            call CharacterCopy(ip,a%ja(ispos+1),a%ExternalIPointersScalar(i)%l)
         elseif (a%SendType(i) == 4) then
            call CharacterCopy(rp,a%ja(ispos+1),a%ExternalRPointersScalar(i)%a)
         elseif (a%SendType(i) == 5) then   
            if (a%ja(ispos+1) == "1") then
               a%ExternalLPointersScalar(i)%l = .true.
            else
               a%ExternalLPointersScalar(i)%l = .false.
            endif
         elseif (a%SendType(i) == 6) then   
            call CharacterCopy(a%ExternalCPointers(i)%nchar,a%ja(ispos+1),a%ExternalCPointers(i)%l)
         endif
      enddo
   
   
   end subroutine
      
   subroutine CheckSizes(a,n)
      use typre
      implicit none
      class(BroadCastBuffer) :: a
      integer(ip) :: n
      
      if (a%nz+1 > size(a%ia)-1) call runend('BroadCastBuffer, too many entries')
      if (a%ia(a%nz+1)+n > size(a%ja)) call runend('BroadCastBuffer, exceeded buffer size')
   end subroutine
   
   subroutine Dealloc(a)
      use typre
      implicit none
      class(BroadCastBuffer) :: a
      
      integer(ip) :: BufferSize,nBufferEntries
      
      BufferSize = size(a%ja)
      nBufferEntries = size(a%ia)-1
      
      call a%Memor%dealloc(nBufferEntries+1,a%ia,'ia','BroadCastBuffer')
      call a%Memor%dealloc(BufferSize,a%ja,'ja','BroadCastBuffer')
      deallocate(a%ExternalRPointers)
      call a%Memor%deallocObj(0,'ExternalRPointers','BroadCastBuffer',nBufferEntries)
      deallocate(a%ExternalLPointers)
      call a%Memor%deallocObj(0,'ExternalLPointers','BroadCastBuffer',nBufferEntries)
      call a%Memor%dealloc(nBufferEntries,a%SendType,'SendType','BroadCastBuffer')
      
      deallocate(a%ExternalIPointers)
      deallocate(a%ExternalI2Pointers)
      deallocate(a%ExternalR2Pointers)
      deallocate(a%ExternalCPointers)
      deallocate(a%ExternalIPointersScalar)
      deallocate(a%ExternalRPointersScalar)
      deallocate(a%ExternalLPointersScalar)
      call a%Memor%deallocObj(0,'ExternalPointersScalar','BroadCastBuffer',5*nBufferEntries) 
      
      
   end subroutine
   
end module

subroutine CharacterCopy(n,character1,character2)
   use typre
   implicit none
   integer(ip) :: n
   character :: character1(n), character2(n)
   
   character2(1:n) = character1(1:n)
end subroutine

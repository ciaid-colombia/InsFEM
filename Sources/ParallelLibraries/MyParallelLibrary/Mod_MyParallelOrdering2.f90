module Mod_MyParallelOrdering2
   use Mod_LinkedList
   use Mod_HeapSortExternal
   use Mod_ParallelOrderingInterface 
   implicit none
   private 
   public MyParallelOrdering2
   
        
   !In this parallel ordering, a sorted binary search is used for the inverse mapping
   type, extends(ParallelOrderingInterface) :: MyParallelOrdering2
      private

      integer(ip), allocatable :: isLocal2Global(:)
      integer(ip), allocatable :: ProcessorList(:)
      
      integer(ip), allocatable :: isGlobal2Local(:,:)
      
contains   
   
      procedure :: Init => Init
      procedure :: Local2GlobalScalar => Local2GlobalScalar
      procedure :: Local2GlobalVector => Local2GlobalVector
      procedure :: Global2LocalScalar => Global2LocalScalar
      procedure :: Global2LocalVector => Global2LocalVector
      procedure :: GetProcScalar => GetProcScalar
      procedure :: GetProcVector => GetProcVector
      procedure :: Deallocate
      
   end type
   
contains

   subroutine Init(a,npoin,isLocal2Global,ProcessorList,Memor)
      use typre
      use Mod_Memor
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: npoin
      integer(ip) :: isLocal2Global(npoin),ProcessorList(npoin)
      type(MemoryMan) :: Memor
      type(MemoryMan), pointer :: Memor2 => NULL()
      
      integer(ip) :: gpoin,ipoin,iprime,modulus,n
      integer(ip),allocatable :: auxOrdering(:)


   
      call Memor%alloc(npoin,a%isLocal2Global,'isLocal2Global','init')
      call Memor%alloc(npoin,a%ProcessorList,'ProcessorList','init')
      a%isLocal2Global = isLocal2Global
      a%ProcessorList = ProcessorList
      
      call Memor%alloc(npoin,2,a%isGlobal2Local,'isGlobal2Local','init')
      
      
      call Memor%alloc(npoin,auxOrdering,'auxOrdering','init')
      
      
      call heapsortAuxArray(npoin,a%isLocal2Global,auxOrdering)
      
      a%isGlobal2Local(:,1) = a%isLocal2Global(auxOrdering(:))
      a%isGlobal2Local(:,2) = auxOrdering(:)
      
      
      call Memor%dealloc(npoin,auxOrdering,'auxOrdering','init')

   end subroutine
      
   subroutine Local2GlobalVector(a,npoin,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
   
      GlobalArray = a%isLocal2Global(LocalArray)
   
   end subroutine
   
   subroutine Local2GlobalScalar(a,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      GlobalArray = a%isLocal2Global(LocalArray)
   
   end subroutine
   
   subroutine Global2LocalVector(a,npoin,GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
      
      integer(ip) :: ipoin
   
      do ipoin = 1,npoin
         call Global2LocalScalar(a,GlobalArray(ipoin),LocalArray(ipoin))
      enddo
   end subroutine
   
   subroutine Global2LocalScalar(a,GlobalArray, LocalArray)
      use typre
      implicit none
      class(MyParallelOrdering2) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      integer(ip) :: modulus,ipos
      
      
      ipos = binarySearch_I (a%isGlobal2Local(:,1), GlobalArray)
      LocalArray = a%isGlobal2Local(ipos,2)
   end subroutine
   
   subroutine GetProcVector(a,npoin,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), ProcNumber(npoin)
   
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine
   
   subroutine GetProcScalar(a,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(MyParallelOrdering2) :: a
      integer(ip) :: LocalArray, ProcNumber
   
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine

   subroutine Deallocate(a,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(MyParallelOrdering2) :: a
      type(MemoryMan) :: Memor
   
      call Memor%dealloc(size(a%isLocal2Global),a%isLocal2Global,'isLocal2Global','init')
      call Memor%dealloc(size(a%ProcessorList),a%ProcessorList,'ProcessorList','init')
      
      call Memor%dealloc(size(a%isGlobal2Local,1),2,a%isGlobal2Local,'isGlobal2Local','init')
      
   end subroutine
   
end module

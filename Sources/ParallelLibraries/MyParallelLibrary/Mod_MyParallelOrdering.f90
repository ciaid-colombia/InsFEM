module Mod_MyParallelOrdering
   use Mod_LinkedList
   use Mod_HeapSortExternal
   use Mod_ParallelOrderingInterface 
   implicit none
   private 
   public MyParallelOrdering
   
   integer(ip), parameter :: primeNumberList(20) = (/23, 61, 101, 199, 401, 827, 1619, 4273, 10781, 20903, 41957, 97883, 200003, 404941, 949001, 2020693, 4250957, 9624157, 21168611, 49979687 /)
   
   type InverseMap
      integer(ip) :: n = 0
      type(LinkedListHeader) :: Header
      integer(ip), pointer :: Sorted(:,:) => NULL()

   end type
      
   !In this parallel ordering, a hashtable approach is followed for the inverse mapping (global2local)
   !the hash function is defined by computing the modulus against prime numbers
   !when there are collisions, a sorted binary search algorithm is used
   type, extends(ParallelOrderingInterface) :: MyParallelOrdering
      private

      integer(ip), allocatable :: isLocal2Global(:)
      integer(ip), allocatable :: ProcessorList(:)
      
      integer(ip) :: primeNumber
      
      type(InverseMap), allocatable :: imap(:)
      type(LinkedListStorage) :: imapStorage
      integer(ip), allocatable :: auxStorage(:) 
      
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

   subroutine Init(a,MPIcomm,npoin,isLocal2Global,ProcessorList,Memor)
      use typre
      use Mod_Memor
      implicit none
      
      class(MyParallelOrdering) :: a
      integer(ip) :: npoin,MPIcomm
      integer(ip) :: isLocal2Global(npoin),ProcessorList(npoin)
      type(MemoryMan), target :: Memor
      
      integer(ip) :: gpoin,ipoin,iprime,modulus,n
      integer(ip),allocatable :: auxForSort(:)

      
   
      call Memor%alloc(npoin,a%isLocal2Global,'isLocal2Global','init')
      call Memor%alloc(npoin,a%ProcessorList,'ProcessorList','init')
      a%isLocal2Global = isLocal2Global
      a%ProcessorList = ProcessorList
      
      !Now we need to build the inverse mapping
      !Choose the most convenient prime number
      a%primeNumber = 0
      do iprime = 1,size(primeNumberList)
         if (primeNumberList(iprime) > npoin/30) then
            a%primeNumber = primeNumberList(iprime)
            exit
         endif
      enddo
      if (a%primeNumber == 0) a%primeNumber = primeNumberList(size(primeNumberList))
      
      !Allocate imap and header for lists
      allocate(a%imap(a%primeNumber))
      call Memor%allocObj(0,'imap','init',a%primeNumber)
      
      !Allocate storage for lists
      call a%imapStorage%Initialize(2*npoin,'NoRemove','Repeat',Memor)
      
      !add to lists
      do ipoin = 1,npoin
         gpoin = a%isLocal2Global(ipoin)
         modulus = mod(gpoin,a%primeNumber)+1
         
         call a%imapStorage%Add(a%imap(modulus)%Header,gpoin)
         call a%imapStorage%Add(a%imap(modulus)%Header,ipoin)
         a%imap(modulus)%n = a%imap(modulus)%n+1
      enddo
      
      !Extract from lists to arrays and sort
      !Allocate memory
      call Memor%alloc(2*npoin,a%auxStorage,'auxStorage','init')
      call Memor%alloc(npoin,auxForSort,'auxForSort','init')
      n = 0
      do iprime = 1,a%primeNumber
         if (a%imap(iprime)%n /= 0) then
            call localsub1(a%imap(iprime),a%imapStorage,n,a%auxStorage,auxForSort)
         endif
      enddo
      
      call a%imapStorage%Dealloc
      
      call Memor%dealloc(npoin,auxForSort,'auxForSort','init')

   end subroutine
      
   subroutine localsub1(imap2,imapStorage,n,auxStorage, auxForSort)
      type(InverseMap), target :: imap2
      type(LinkedListStorage) :: imapStorage
      integer(ip) :: n
      integer(ip), target :: auxStorage(*)
      integer(ip) :: auxForSort(*)
      
      integer(ip) :: ipoin,ipoin2
      integer(ip) :: n2,i2
      
      integer(ip), target :: superauxarray(2,imap2%n)

      !Extract
      imap2%sorted(1:imap2%n,1:2) => auxStorage(2*n+1:2*n+2*imap2%n)
      
      call imapStorage%ListToArray(imap2%Header,superauxarray)
      
      imap2%sorted(:,1) = superauxarray(1,:)
      imap2%sorted(:,2) = superauxarray(2,:)
      
      !Update general counter
      n = n + imap2%n
      
      !If necessary, sort
      if (imap2%n > 1) then

         call HeapSortAuxArray(imap2%n,imap2%Sorted(:,1),auxForSort(1:imap2%n))
         
         imap2%Sorted(1:imap2%n,1) = imap2%Sorted(auxForSort(1:imap2%n),1)
         imap2%Sorted(1:imap2%n,2) = imap2%Sorted(auxForSort(1:imap2%n),2)
         
         
         !Avoid repeated entries
         n2 = 1
         auxForSort(1) = 1
         do i2 = 2,imap2%n
            if (imap2%Sorted(i2,1) /= imap2%Sorted(i2-1,1)) then
               n2 = n2+1
               auxForSort(n2) = i2
            else
               call runend('repeated entries')
            
            endif
         enddo
         imap2%n = n2
         imap2%Sorted(1:imap2%n,1) = imap2%Sorted(auxForSort(1:imap2%n),1)
         imap2%Sorted(1:imap2%n,2) = imap2%Sorted(auxForSort(1:imap2%n),2)
      endif
      
      
      
   end subroutine
  
  
   
   subroutine Local2GlobalVector(a,npoin,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
   
      GlobalArray = a%isLocal2Global(LocalArray)
   
   end subroutine
   
   subroutine Local2GlobalScalar(a,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      GlobalArray = a%isLocal2Global(LocalArray)
   
   end subroutine
   
   subroutine Global2LocalVector(a,npoin,GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(MyParallelOrdering) :: a
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
      class(MyParallelOrdering) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      integer(ip) :: modulus,ipos
      
      modulus = mod(GlobalArray,a%primeNumber)+1
      
      if (a%imap(modulus)%n == 1) then
         if (a%imap(modulus)%Sorted(1,1) == GlobalArray) then
            LocalArray = a%imap(modulus)%Sorted(1,2)
         else
            LocalArray = 0
         endif 
      else
         ipos = binarySearch_I (a%imap(modulus)%Sorted(1:a%imap(modulus)%n,1), GlobalArray)
         if (ipos /= 0) then
            LocalArray = a%imap(modulus)%Sorted(ipos,2)      
         else
            LocalArray = 0
         endif
         
      endif
      
      
   end subroutine
   
   subroutine GetProcVector(a,npoin,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(MyParallelOrdering) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), ProcNumber(npoin)
   
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine
   
   subroutine GetProcScalar(a,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(MyParallelOrdering) :: a
      integer(ip) :: LocalArray, ProcNumber
   
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine

   subroutine Deallocate(a,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(MyParallelOrdering) :: a
      type(MemoryMan) :: Memor
   
      call Memor%dealloc(size(a%isLocal2Global),a%isLocal2Global,'isLocal2Global','init')
      call Memor%dealloc(size(a%ProcessorList),a%ProcessorList,'ProcessorList','init')
      deallocate(a%imap)
      call Memor%deallocObj(0,'imap','init',a%primeNumber)
      call Memor%dealloc(size(a%auxStorage),a%auxStorage,'auxStorage','init')
      
   end subroutine
   
end module

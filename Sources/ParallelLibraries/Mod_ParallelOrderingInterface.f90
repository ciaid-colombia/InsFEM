module Mod_ParallelOrderingInterface
   implicit none
   private
   public ParallelOrderingInterface
   
   type ParallelOrderingInterface
      private
   
   
contains   
   
      procedure :: Init
      procedure :: Local2GlobalScalar
      procedure :: Local2GlobalVector
      procedure :: Global2LocalScalar
      procedure :: Global2LocalVector
      procedure :: GetProcScalar
      procedure :: GetProcVector
      procedure :: Deallocate
      
      generic :: Local2Global => Local2GlobalScalar, Local2GlobalVector
      generic :: Global2Local => Global2LocalScalar, Global2LocalVector
      generic :: GetProc => GetProcScalar, GetProcVector
   
   end type

contains

   subroutine Init(a,MPIcomm,npoin,isLocal2Global,ProcessorList,Memor)
      use typre
      use Mod_Memor
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: npoin,MPIcomm
      integer(ip) :: isLocal2Global(npoin),ProcessorList(npoin)
      type(MemoryMan), target :: Memor
   
      call runend('ParallelOrdering, InitNoProcList procedure not implemented')
   
   end subroutine
   
   subroutine Local2GlobalVector(a,npoin,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
   
      call runend('ParallelOrdering, GetGlobalFromLocal procedure not implemented')
   
   end subroutine
   
   subroutine Local2GlobalScalar(a,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      call runend('ParallelOrdering, GetGlobalFromLocal procedure not implemented')
   
   end subroutine
   
   subroutine Global2LocalVector(a,npoin,GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
   
      call runend('ParallelOrdering, GetLocalFromGlobal procedure not implemented')
   end subroutine
   
   subroutine Global2LocalScalar(a,GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      call runend('ParallelOrdering, GetLocalFromGlobal procedure not implemented')
   end subroutine
   
   subroutine GetProcVector(a,npoin,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), ProcNumber(npoin)
   
      call runend('ParallelOrdering, GetProcNumber procedure not implemented')
   end subroutine
   
   subroutine GetProcScalar(a,LocalArray, ProcNumber)
      use typre
      implicit none
      
      class(ParallelOrderingInterface) :: a
      integer(ip) :: LocalArray, ProcNumber
   
      call runend('ParallelOrdering, GetProcNumber procedure not implemented')
   end subroutine

   subroutine Deallocate(a,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(ParallelOrderingInterface) :: a
      type(MemoryMan) :: Memor
   
      call runend('ParallelOrdering, Deallocate procedure not implemented')
   
   end subroutine
   
end module

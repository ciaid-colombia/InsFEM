module Mod_MyParallelLibrary
   implicit none
   
   
contains

   subroutine CreateMyCommunicator(ParallelCommunicator, Memor)  !FACTORY
      use Mod_ParallelCommunicatorInterface
      use Mod_MyCommunicator
      use Mod_Memor
      use typre
      implicit none
      class(ParallelCommunicatorInterface), pointer :: ParallelCommunicator
      type(MemoryMan) :: Memor
      integer(ip) :: istat
      
      class(MyCommunicator), pointer :: iParallelCommunicator
      
      allocate(iParallelCommunicator,stat=istat)   
      call Memor%allocObj(istat,'ParallelCommunicator','CreateCommunicatorPETSc',1)
      
      ParallelCommunicator => iParallelCommunicator
   
   end subroutine
   
   subroutine CreateMyOrdering(ParallelOrdering,Memor)
      use Mod_ParallelOrderingInterface
      use Mod_MyParallelOrdering
      use Mod_Memor
      use typre
      implicit none
      class(ParallelOrderingInterface), pointer :: ParallelOrdering
      type(MemoryMan) :: Memor
      
      integer(ip) :: istat
      class(MyParallelOrdering), pointer :: iParallelOrdering => NULL()
      
      allocate(iParallelOrdering,stat=istat)   
      call Memor%allocObj(istat,'ParallelOrdering','CreateOrderingPETSc',1)
      
      ParallelOrdering => iParallelOrdering
   end subroutine
   
   
end module

module Mod_ParallelLibraryInterface
   use typre
   use Mod_ParallelOrderingInterface
   use Mod_ParallelCommunicatorInterface
   use Mod_ParallelSystemInterface
   implicit none
   private
   public ParallelLibraryInterface

   type :: ParallelLibraryInterface
      private

   contains

      procedure :: Initialize
      procedure :: Finalize
      procedure :: CreatePartitioner
      procedure :: CreateRebalancePartitioner
      procedure :: DeallocatePartitioner
      procedure :: CreateOrdering
      procedure :: CreateCommunicator
      procedure :: CreateSystem
      procedure :: DeallocateOrdering
      procedure :: DeallocateCommunicator
      procedure :: DeallocateSystem
      procedure :: GetMaximumMemory

   end type

contains

   subroutine Initialize(a,namda)
      implicit none

      class(ParallelLibraryInterface) :: a
      character(*) :: namda

      call runend('ParallelLibrary, Initialize procedure not implemented')
   end subroutine

   subroutine Finalize(a)
      implicit none

      class(ParallelLibraryInterface) :: a

      call runend('ParallelLibrary, Finalize procedure not implemented')
   end subroutine

   subroutine CreatePartitioner(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(ParallelLibraryInterface) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, CreatePartitioner not implemented')
   end subroutine

   subroutine CreateRebalancePartitioner(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(ParallelLibraryInterface) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, CreatePartitioner not implemented')
   end subroutine

   subroutine DeallocatePartitioner(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(ParallelLibraryInterface) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, DeallocatePartitioner not implemented')
   end subroutine

   subroutine CreateOrdering(a,ParallelOrdering, Memor)   !FACTORY 
      use typre
      use Mod_ParallelOrderingInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelOrderingInterface), pointer  :: ParallelOrdering
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, CreateOrdering procedure not implemented')

   end subroutine

   subroutine DeallocateOrdering(a,ParallelOrdering, Memor)   !FACTORY 
      use typre
      use Mod_ParallelOrderingInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelOrderingInterface), pointer  :: ParallelOrdering
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, DeallocateOrdering procedure not implemented')

   end subroutine

   subroutine CreateCommunicator(a,ParallelCommunicator, Memor)   !FACTORY 
      use typre
      use Mod_ParallelCommunicatorInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelCommunicatorInterface), pointer  :: ParallelCommunicator
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, CreateCommunicator procedure not implemented')

   end subroutine

   subroutine DeallocateCommunicator(a,ParallelCommunicator, Memor)   !FACTORY 
      use typre
      use Mod_ParallelCommunicatorInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelCommunicatorInterface), pointer  :: ParallelCommunicator
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, DeallocateCommunicator procedure not implemented')

   end subroutine

   subroutine CreateSystem(a,ParallelSystem,Memor)   !FACTORY 
      use typre
      use Mod_ParallelSystemInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelSystemInterface), pointer  :: ParallelSystem
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, CreateSystem procedure not implemented')

   end subroutine

   subroutine DeallocateSystem(a,ParallelSystem, Memor)   !FACTORY 
      use typre
      use Mod_ParallelSystemInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      class(ParallelSystemInterface), pointer  :: ParallelSystem
      type(MemoryMan) :: Memor

      call runend('ParallelLibrary, DeallocateSystem procedure not implemented')

   end subroutine

   subroutine GetMaximumMemory(a,maxmemo)   !FACTORY 
      use typre
      use Mod_ParallelSystemInterface
      use Mod_Memor
      implicit none

      class(ParallelLibraryInterface) :: a
      integer(8) :: maxmemo

      call runend('ParallelLibrary, maxmemo not implemented')

   end subroutine



end module

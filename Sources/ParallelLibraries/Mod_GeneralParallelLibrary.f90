module Mod_GeneralParallelLibrary
   use Mod_ParallelLibraryInterface
   use Mod_MyCommunicator
   use Mod_PetscLibrary
   use Mod_Memor
   use Mod_MyParallelLibrary
   implicit none
   private
   public GeneralParallelLibrary, GeneralParallelLibrary_Const
   
   type, extends(ParallelLibraryInterface) :: GeneralParallelLibrary
      character(6) :: PartitionerType = 'Petsc'
      character(6) :: RebalancePartitionerType = 'Petsc'
      character(6) :: CommunicatorType = 'Petsc'
      !character(6) :: OrderingType = 'My    '
      character(6) :: OrderingType = 'Petsc'
   
contains   
   
   procedure :: SetPartitionerType !Choose between petsc and zoltan as partitioners
   procedure :: SetRebalancePartitionerType !Choose between petsc and zoltan as partitioners
   procedure :: SetCommunicatorType
   procedure :: SetOrderingType
   
   procedure :: Initialize => InitializeGeneral
   procedure :: Finalize   => FinalizeGeneral
   procedure :: CreatePartitioner => CreatePartitionerGeneral
   procedure :: CreateRebalancePartitioner => CreateRebalancePartitionerGeneral
   procedure :: DeallocatePartitioner => DeallocatePartitionerGeneral
   procedure :: CreateOrdering => CreateOrderingGeneral
   procedure :: CreateCommunicator => CreateCommunicatorGeneral
   procedure :: CreateSystem => CreateSystemGeneral
   procedure :: DeallocateOrdering => DeallocateOrderingGeneral
   procedure :: DeallocateCommunicator => DeallocateCommunicatorGeneral
   procedure :: DeallocateSystem => DeallocateSystemGeneral
   procedure :: GetMaximumMemory => GetMaximumMemoryGeneral
   
   
   

   end type
   
   interface GeneralParallelLibrary_Const
      procedure constructor
   
   end interface

contains
   function constructor()
      type(GeneralParallelLibrary), pointer :: constructor
      
      allocate(constructor)
   end function

   subroutine SetPartitionerType(a,PartitionerLib)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      character(6) :: PartitionerLib

      a%PartitionerType = PartitionerLib
   end subroutine

   subroutine SetRebalancePartitionerType(a,PartitionerLib)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      character(6) :: PartitionerLib

      a%RebalancePartitionerType = PartitionerLib
   end subroutine

   subroutine SetCommunicatorType(a,CommunicatorType)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      character(6) :: CommunicatorType

      a%CommunicatorType = CommunicatorType
   end subroutine


   subroutine SetOrderingType(a,OrderingType)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      character(6) :: OrderingType
      
      a%OrderingType = OrderingType
   end subroutine
   
   !--------------------------------------------------------


   subroutine InitializeGeneral(a,namda)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      character(*) :: namda

      !Petsc
      call InitializePETSc

   end subroutine

   subroutine FinalizeGeneral(a)
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a

      !Petsc
      call FinalizePETSc

   end subroutine

   subroutine CreateOrderingGeneral(a, ParallelOrdering, Memor)  !FACTORY
      use Mod_ParallelOrderingInterface
      use Mod_Memor
      use typre
      implicit none

      class(GeneralParallelLibrary) :: a
      class(ParallelOrderingInterface), pointer :: ParallelOrdering
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      if (a%OrderingType .eq. 'Petsc') then
         call CreateOrderingPETSc(ParallelOrdering, Memor)
      elseif (a%OrderingType .eq. 'My') then
         call CreateMyOrdering(ParallelOrdering,Memor)
      endif
      

   end subroutine

   subroutine CreateCommunicatorGeneral(a, ParallelCommunicator, Memor)  !FACTORY
      use Mod_ParallelCommunicatorInterface
      use Mod_Memor
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelCommunicatorInterface), pointer :: ParallelCommunicator
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      if (a%CommunicatorType .eq. 'Petsc') then
         call CreateCommunicatorPETSc(ParallelCommunicator, Memor)
      elseif (a%CommunicatorType .eq. 'My') then
         call CreateMyCommunicator(ParallelCommunicator, Memor)
      endif

   end subroutine

   subroutine CreatePartitionerGeneral(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      call CreatePartitioner(ParallelPartitioner,Memor,a%PartitionerType)
   end subroutine

   subroutine CreateRebalancePartitionerGeneral(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface

      use Mod_Memor
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      call CreatePartitioner(ParallelPartitioner,Memor,a%RebalancePartitionerType)
   end subroutine

   subroutine CreatePartitioner(ParallelPartitioner,Memor,PartitionerType)
      use typre
      use Mod_ParallelPartitionerInterface
      implicit none
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor
      character(6) :: PartitionerType

      if (PartitionerType .eq. 'Petsc') then
         call CreatePartitionerPetsc(ParallelPartitioner,Memor)
      else
         call runend('CreatePartitionerGeneral: Wrong type of parallel partitioner')
      endif

   end subroutine

   subroutine DeallocatePartitionerGeneral(a,ParallelPartitioner,Memor)
      use typre
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor

      integer(ip) :: istat

      deallocate(ParallelPartitioner,stat = istat)
      call Memor%DeallocObj(istat,'ParallelPartitioner','DeallocatePartitionerGeneral',1)
   end subroutine

   subroutine CreateSystemGeneral(a,ParallelSystem,Memor)  !FACTORY
      use Mod_ParallelSystemInterface
      use Mod_Memor
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelSystemInterface), pointer :: ParallelSystem
      type(MemoryMan) :: Memor

         call CreateSystemPETSc(ParallelSystem, Memor)

   end subroutine

   subroutine DeallocateOrderingGeneral(a, ParallelOrdering, Memor)  !FACTORY
      use Mod_ParallelOrderingInterface
      use Mod_Memor
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelOrderingInterface), pointer :: ParallelOrdering
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      deallocate(ParallelOrdering,stat=istat)   
      call Memor%deallocObj(istat,'ParallelOrdering','DeallocateOrderingGeneral',1)
   end subroutine

   subroutine DeallocateCommunicatorGeneral(a, ParallelCommunicator, Memor)  !FACTORY
      use Mod_ParallelCommunicatorInterface
      use Mod_Memor
      use typre
      implicit none

      class(GeneralParallelLibrary) :: a
      class(ParallelCommunicatorInterface), pointer :: ParallelCommunicator
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      deallocate(ParallelCommunicator,stat=istat)   
      call Memor%deallocObj(istat,'ParallelCommunicator','DeallocateCommunicatorGeneral',1)
   end subroutine

   subroutine DeallocateSystemGeneral(a, ParallelSystem, Memor)  !FACTORY
      use Mod_ParallelSystemInterface
      use Mod_Memor
      use typre
      implicit none
      class(GeneralParallelLibrary) :: a
      class(ParallelSystemInterface), pointer :: ParallelSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      deallocate(ParallelSystem,stat=istat)   
      call Memor%deallocObj(istat,'ParallelSystem','DeallocateSystemGeneral',1)

   end subroutine

   subroutine GetMaximumMemoryGeneral(a,maxmemo)   !FACTORY 
      use typre
      use Mod_Memor
      implicit none
      class(GeneralParallelLibrary) :: a
      integer(8) :: maxmemo

      call GetMaximumMemoryPETSC(maxmemo)  
   end subroutine

end module

module Mod_PETScLibrary
   use typre
#include <petsc/finclude/petscsys.h>
   use petscsys
   implicit none
   integer(ip) :: IsInitializedPetsc = 0
   
contains
   
   subroutine InitializePETSc
      implicit none
      integer :: ierr
      
      !we ensure that it is initialized only once
      if (IsInitializedPetsc == 0) then
         call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
         call PetscMemorySetGetMaximumUsage(ierr)
      endif
      isInitializedPetsc = isInitializedPetsc+1   
   end subroutine
   
   subroutine FinalizePETSc
      implicit none
      integer :: ierr
      
      !We ensure that it is finalized only once
      isInitializedPetsc = isInitializedPetsc-1
      if (isInitializedPetsc == 0) then
         call PetscFinalize(ierr)
      endif   
   end subroutine
   
   subroutine CreateOrderingPETSc(ParallelOrdering, Memor)  !FACTORY
      use Mod_PETScOrdering
      use Mod_ParallelOrderingInterface
      use Mod_Memor
      use typre
      implicit none
      class(ParallelOrderingInterface), pointer :: ParallelOrdering
      type(MemoryMan) :: Memor
      integer(ip) :: istat
      
      class(PETScOrdering), pointer :: iParallelOrdering
      
      allocate(iParallelOrdering,stat=istat)   
      call Memor%allocObj(istat,'ParallelOrdering','CreateOrderingPETSc',1)
      
      ParallelOrdering => iParallelOrdering
      
   end subroutine
   
   subroutine CreateCommunicatorPETSc(ParallelCommunicator, Memor)  !FACTORY
      use Mod_PETScCommunicator
      use Mod_ParallelCommunicatorInterface
      use Mod_Memor
      use typre
      implicit none
      class(ParallelCommunicatorInterface), pointer :: ParallelCommunicator
      type(MemoryMan) :: Memor
      integer(ip) :: istat
      
      class(PETScCommunicator), pointer :: iParallelCommunicator => NULL()
      
      allocate(iParallelCommunicator,stat=istat)   
      call Memor%allocObj(istat,'ParallelCommunicator','CreateCommunicatorPETSc',1)
      
      ParallelCommunicator => iParallelCommunicator
   
   end subroutine
   
   subroutine CreatePartitionerPETSc(ParallelPartitioner,Memor)
      use typre
      use Mod_PetscPartitioner
      use Mod_ParallelPartitionerInterface
      use Mod_Memor
      implicit none
      class(ParallelPartitionerInterface), pointer :: ParallelPartitioner
      type(MemoryMan) :: Memor
      character(6) :: PartitionerType
      
      integer(ip) :: istat
      class(PetscPartitioner1), pointer :: iParallelPartitioner => NULL()
      
      allocate(iParallelPartitioner,stat = istat)
      ParallelPartitioner => iParallelPartitioner
      call Memor%allocObj(istat,'ParallelPartitioner','CreatePartitionerPetsc',1)
      
   end subroutine
   
   subroutine CreateSystemPETSc(ParallelSystem,Memor)  !FACTORY
      use Mod_PETScSystem
      use Mod_ParallelSystemInterface
      use Mod_Memor
      use typre
      implicit none
      class(ParallelSystemInterface), pointer :: ParallelSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat    
      
      class(PETScSystem), pointer :: iParallelSystem => NULL()
      
      allocate(iParallelSystem,stat=istat)   
      call Memor%allocObj(istat,'ParallelSystem','CreateSystemPETSc',1)
      
      ParallelSystem => iParallelSystem
   
   end subroutine
   
   subroutine GetMaximumMemoryPETSC(maxmemo)   !FACTORY 
      use typre
      use Mod_Memor
#include <petscversion.h>
#include <petsc/finclude/petscsys.h>
   use petscsys
      implicit none
      integer(8) :: maxmemo
      integer(ip) :: ierr
      real(rp) :: mem
      
      maxmemo = 0
      !call PetscMallocGetMaximumUsage(maxmemo,ierr)
      call PetscMemoryGetMaximumUsage(mem,ierr)
      maxmemo = mem
   end subroutine  
   
end module









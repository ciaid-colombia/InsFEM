module Mod_PETScOrdering
   use typre
   use Mod_ParallelOrderingInterface
#include <petsc/finclude/petscis.h>
   use petscis
   implicit none
   private
   public PETScOrdering
   
   type , extends(ParallelOrderingInterface) :: PETScOrdering
  
      integer(ip), allocatable ::  IsLocal2Global(:),ProcessorList(:)
      ISLocalToGlobalMapping   ::  PIsLocal2Global
      integer(ip) :: npoin, MPIcomm
      
contains   
      procedure :: Init                => InitPETSc
      procedure :: InitNdofn
      procedure :: Local2GlobalVector  => Local2GlobalVectorPETSc
      procedure :: Local2GlobalScalar  => Local2GlobalScalarPETSc
      procedure :: Global2LocalVector  => Global2LocalVectorPETSc
      procedure :: Global2LocalScalar  => Global2LocalScalarPETSc
      procedure :: GetProcScalar       => GetProcScalarPETSc
      procedure :: GetProcVector       => GetProcVectorPETSc
      procedure :: Deallocate          => DeallocatePETSc
   
   end type

contains

   subroutine InitPETSc(a,MPIcomm,npoin,isLocal2Global,ProcessorList,Memor)
      use typre
      use Mod_Memor
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: npoin,IsLocal2Global(npoin),ProcessorList(npoin),MPIcomm
      type(MemoryMan), target :: Memor
      integer :: ierr
   
      a%MPIcomm = MPIcomm
      a%npoin = npoin
      call Memor%alloc(npoin,a%IsLocal2Global,'IsLocal2Global','PETScOrdering_INIT')
      a%IsLocal2Global = IsLocal2Global-1

      call ISLocalToGlobalMappingCreate(a%MPIcomm,1_ip,npoin,a%IsLocal2Global,PETSC_COPY_VALUES,a%PIsLocal2Global,ierr)

      call Memor%alloc(a%npoin,a%ProcessorList,'ProcessorList','PETScOrdering_INIT')
      a%ProcessorList(1:a%npoin) = ProcessorList(1:a%npoin)
   end subroutine
   
   subroutine InitNdofn(a,MPIComm,npoin,isLocal2Global,ProcessorList,Memor,ndofn)
      use typre
      use Mod_Memor
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: npoin,IsLocal2Global(npoin),ProcessorList(npoin),ndofn,MPIComm
      type(MemoryMan), target :: Memor
      integer :: ierr
   
      a%MPIcomm = MPIcomm
      a%npoin = npoin
      call Memor%alloc(npoin,a%IsLocal2Global,'IsLocal2Global','PETScOrdering_INIT')
      a%IsLocal2Global = IsLocal2Global-1
      
      call Memor%alloc(a%npoin,a%ProcessorList,'ProcessorList','PETScOrdering_INIT')
      a%ProcessorList(1:a%npoin) = ProcessorList(1:a%npoin)

      call ISLocalToGlobalMappingCreate(a%MPIcomm,ndofn,npoin,a%IsLocal2Global,PETSC_COPY_VALUES,a%PIsLocal2Global,ierr)
   end subroutine
   
   subroutine Local2GlobalVectorPETSc(a,npoin,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
   
      GlobalArray = a%IsLocal2Global(LocalArray) + 1
   
   end subroutine
   
    subroutine Local2GlobalScalarPETSc(a,LocalArray, GlobalArray)
      use typre
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: LocalArray, GlobalArray
   
      GlobalArray = a%IsLocal2Global(LocalArray) + 1
   end subroutine
   
   subroutine Global2LocalVectorPETSc(a,npoin,GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), GlobalArray(npoin)
      integer(ip) :: idummy,ierr
   
      !Transform lnods to local numbering
      LocalArray = GlobalArray- 1
      call ISGlobalToLocalMappingApply(a%PIsLocal2Global,IS_GTOLM_MASK,npoin,LocalArray,idummy,LocalArray,ierr)
      LocalArray = LocalArray + 1 
   end subroutine
   
   subroutine Global2LocalScalarPETSc(a, GlobalArray, LocalArray)
      use typre
      implicit none
      
      class(PETScOrdering) :: a
      integer(ip) :: LocalArray, GlobalArray, auxLocalArray(1)
      integer(ip) :: idummy,ierr
   
      !Transform lnods to local numbering
      auxLocalArray(1) = GlobalArray- 1
      call ISGlobalToLocalMappingApply(a%PIsLocal2Global,IS_GTOLM_MASK,1,auxLocalArray,idummy,auxLocalArray,ierr)
      auxLocalArray(1) = auxLocalArray(1) + 1
      LocalArray = auxLocalArray(1)
   end subroutine
   
   subroutine GetProcVectorPETSc(a,npoin,LocalArray, ProcNumber)
      use typre
      implicit none
      class(PETScOrdering) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalArray(npoin), ProcNumber(npoin)
      
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine
   
   subroutine GetProcScalarPETSc(a,LocalArray, ProcNumber)
      use typre
      implicit none
      class(PETScOrdering) :: a
      integer(ip) :: LocalArray, ProcNumber
      
      ProcNumber = a%ProcessorList(LocalArray)
   end subroutine
   
   subroutine DeallocatePETSc(a,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(PETScOrdering) :: a
      type(MemoryMan) :: Memor
      
      integer(ip) :: iaux
      integer :: ierr
      
      call Memor%dealloc(a%npoin,a%IsLocal2Global,'IsLocal2Global','PETScOrdering_DEALLOCATEIS')
      call ISLocalToGlobalMappingDestroy(a%PIslocal2Global,ierr)
      if (allocated(a%ProcessorList)) then
         call Memor%dealloc(a%npoin,a%ProcessorList,'ProcessorList','PETScOrdering_DEALLOCATEIS')
      endif
   
   end subroutine

end module

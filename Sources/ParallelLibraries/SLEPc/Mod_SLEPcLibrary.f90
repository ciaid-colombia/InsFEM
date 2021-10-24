module Mod_SLEPcLibrary
   use typre
   use Mod_Memor
   use Mod_MyCommunicator
   use Mod_ParallelSystemInterface
   use Mod_EigenLibraryInterface
   use Mod_EigenSystemInterface
   use Mod_SLEpcSystem
#ifdef SLEPC
#include <petscversion.h>
#include <slepc/finclude/slepc.h>
#include <petsc/finclude/petscsys.h>
   use petscsys
   private
   public SLEPcLibrary

   type,abstract, extends(EigenLibraryInterface) :: SLEPcLibrary

      character(3) :: solvertype

   contains   

      procedure :: Initialize               => InitializeSLEPc
      procedure :: Finalize                 => FinalizeSLEPc
      procedure :: DeallocateSystem         => DeallocateSystemSLEPc
      procedure :: CreateSystemProj         => CreateSystemPETSc 

   end type

contains

   subroutine InitializeSLEPc(a,namda)
      implicit none

      class(SLEPcLibrary) :: a
      integer ::  ierr
      character(*) :: namda

      a%solvertype = namda
      call SLEPcFinalize(ierr)
      call SLEPcInitialize(PETSC_NULL_CHARACTER,ierr)
      call PETScMemorySetGetMaximumUsage(ierr)

   end subroutine

   subroutine FinalizeSLEPc(a)
      implicit none

      class(SLEPcLibrary) :: a
      integer :: ierr

      call SLEPcFinalize(ierr)

   end subroutine


   subroutine DeallocateSystemSLEPc(a, EigenSystem, Memor)  !FACTORY
      implicit none
      class(SLEPcLibrary) :: a
      class(EigenSystemInterface), pointer :: EigenSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      deallocate(EigenSystem,stat=istat)   
      call Memor%deallocObj(istat,'EigenSystem','DeallocateSystemSLEPc',1)

   end subroutine
   
   subroutine CreateSystemPETSc(a,ParallelSystem,Memor)  !FACTORY
      implicit none
      class(SLEPcLibrary) :: a
      class(ParallelSystemInterface), pointer :: ParallelSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat    
      
      class(BasisPETScSystem), pointer :: iParallelSystem => NULL()
      
      allocate(iParallelSystem,stat=istat)   
      call Memor%allocObj(istat,'ParallelSystem','CreateSystemPETSc',1)
      
      ParallelSystem => iParallelSystem
   
   end subroutine
   
#endif  

end module


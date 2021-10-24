module Mod_EPSSystem
   use typre
   use Mod_Memor
   use Mod_PETScSystem
   use Mod_SLEPcSystem
   use Mod_ParallelSystemInterface
   use Mod_EigenSystemInterface
#ifdef SLEPC
#include <slepc/finclude/slepceps.h>  
#include <slepc/finclude/slepcbv.h>  
   use slepceps
   use slepcbv
   implicit none
   private
   public EPSSystem
   
   type, extends(SLEPcSystem) :: EPSSystem

      Mat            :: A_Matrix,B_Matrix
      Vec            :: xr,xi     !Re and Im part of EigVector
      EPS            :: eps
      ST             :: st
      EPSType        :: tname
      PetscReal      :: eigr,eigi !Re and Im part of EigVector
      PetscReal      :: tol
      PetscInt       :: nev,maxit,its

   contains

      procedure :: InitBuild              => EPSInitSystem
      procedure :: SetMatrixPointersBuild => EPSsetMatrixPointersBuild
      procedure :: SolveSystem            => EPSSolveSys
      procedure :: SpecificWriter         => EPSWrite
      procedure :: Deallocate             => EPSSystemDeallocate
   end type

contains

   subroutine EPSInitSystem(a,isnap,SnapInterval,nsnap)
      implicit none
      class(EPSSystem) :: a
      class(ParallelSystemInterface),pointer:: LS => NULL()
      integer(ip)      :: iaux,nsnap,SnapInterval
      integer(ip), target :: isnap
      PetscErrorCode   :: ierr
      
      a%SnapInterval = SnapInterval
      a%isnap => isnap
      a%nsnap = 0

      call VecCreate(a%MPIcomm,a%xr,ierr)
      call VecSetFromOptions(a%xr,ierr)
      call VecSetSizes(a%xr,PETSC_DECIDE,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(a%xr,ierr)

      call VecCreate(a%MPIcomm,a%xi,ierr)
      call VecSetFromOptions(a%xi,ierr)
      call VecSetSizes(a%xi,PETSC_DECIDE,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(a%xi,ierr)

      !Set up solver
      call EPSCreate(a%MPIcomm,a%eps,ierr)
      call a%Memor%allocObj(ierr,'EPSSolver','EPSSystemallocate',1)

      call EPSSetOptionsPrefix(a%eps,adjustl(trim(a%exmod)),ierr)
      call EPSSetFromOptions(a%eps,ierr)
  
  end subroutine

  subroutine EPSsetMatrixPointersBuild(a)
      implicit none
      class(EPSSystem) :: a
      class(ParallelSystemInterface),pointer:: LS => NULL()
      PetscErrorCode :: ierr

      LS => a%LinearSystem
      select type(LS)
      class is(PETScSystem)

          call MatConvert(LS%K_Matrix   ,MATMPIAIJ,MAT_INITIAL_MATRIX,a%B_Matrix,ierr)
          call MatConvert(LS%Mass_Matrix,MATMPIAIJ,MAT_INITIAL_MATRIX,a%A_Matrix,ierr)

      end select

  end subroutine

  subroutine EPSSolveSys(a)
      implicit none
      class(EPSSystem) :: a
      integer(ip)      :: i,iaux
      PetscReal        :: err
      Vec              :: u_aux1
      PetscErrorCode   :: ierr

      call SLEPcInitialize(PETSc_NULL_CHARACTER,ierr)

      call EPSSetOperators(a%eps,a%A_Matrix,a%B_Matrix,ierr)

      !Solve the system
      call EPSSolve(a%eps,ierr) 

      call EPSGetIterationNumber(a%eps,a%its,ierr)
      call EPSGetST(a%eps,a%st,ierr)
      call EPSGetType(a%eps,a%tname,ierr)
      call EPSGetDimensions(a%eps,a%nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
      call EPSGetTolerances(a%eps,a%tol,a%maxit,ierr)

      !Display solution and clean up
      call EPSGetConverged(a%eps,a%nconv,ierr)
      a%ndofr  = a%nconv

      call BVCreate(a%MPIcomm,a%Basis,ierr)
      call BVSetFromOptions(a%Basis,ierr)
      call BVSetSizes(a%Basis,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,a%nconv,ierr)
      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%allocObj(ierr,'Basis','EPSSystem',iaux)

      call a%Memor%alloc(a%nconv,a%sigma,'sigma','EPSSystem')

      call VecCreate(a%MPIcomm,u_aux1,ierr)
      call VecSetFromOptions(u_aux1,ierr)
      call VecSetSizes(u_aux1,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(u_aux1,ierr)

      do i=0,a%nconv-1

          call EPSGetEigenpair(a%eps,i,a%sigma(i+1),a%eigi,u_aux1,a%xi,ierr)
          call BVInsertVec(a%Basis,i,u_aux1,ierr)

      end do

  end subroutine

  subroutine EPSWrite(a)
     implicit none
     class(EPSSystem) :: a
     PetscErrorCode :: ierr

     call EPSView(a%eps,a%outviewer,ierr)

  end subroutine

  subroutine EPSSystemDeallocate(a)
      implicit none
      class(EPSSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call EPSReset(a%eps,ierr)
      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%deallocObj(ierr,'Basis','EPSSystem',iaux)

      call EPSDestroy(a%eps,ierr)
      call a%Memor%deallocObj(ierr,'EPSSolver','EPSSystemDeallocate',1)
      
      call a%Memor%dealloc(a%nconv,a%sigma,'sigma','EPSSystem')

      !Clean up after solve
      call VecDestroy(a%xr,ierr)
      call VecDestroy(a%xi,ierr)

   end subroutine
#endif

end module

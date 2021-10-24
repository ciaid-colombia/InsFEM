module Mod_SVDSystem
   use typre
   use Mod_PETScSystem
   use Mod_SLEPcSystem
   use Mod_ParallelSystemInterface
   use Mod_EigenSystemInterface
#ifdef SLEPC
#include <slepc/finclude/slepcsvd.h>  
#include <slepc/finclude/slepcbv.h>  
   use slepcsvd
   use slepcbv
   implicit none
   private
   public SVDSystem

   type, extends(SLEPcSystem) :: SVDSystem

      SVD          :: svd
      KSP          :: SOLVER1
      Mat          :: auxMass,LMass,invLMass

   contains

      procedure :: InitBuild           => SVDInitSystem
      
      procedure :: SnapshotsMean       => SVDSnapshotsMean
      procedure :: SnapshotsMeanToZero => SVDSnapshotsMeanToZero
      procedure :: SetSnapshots        => SVDSetSnapshots
      
      procedure :: SetLumpedMass         => SVDSetLumpedMass
      procedure :: AssemblySnapshotsMass => SVDAssemblySnapshotsMass
      procedure :: MassMultBasis         => SVDMassMultBasis
      
      procedure :: AssemblySnapshots   => SVDAssemblySnapshots

      procedure :: SolveSystem         => SVDSolveSys
      procedure :: SpecificWriter      => SVDWrite
      procedure :: LMassDeallocate     => SVDLMassDeallocate
      procedure :: Deallocate          => SVDSystemDeallocate
   end type

contains

   subroutine SVDInitSystem(a,isnap,SnapInterval,nsnap)
      implicit none
      class(SVDSystem) :: a
      integer(ip), target       :: isnap
      integer(ip)      :: iaux,nsnap,SnapInterval
      PetscErrorCode   :: ierr

      a%SnapInterval = SnapInterval
      a%isnap => isnap
      a%nsnap = nsnap
      
      call MatCreate(a%MPIcomm,a%Snapshots,ierr)
      call MatSetFromOptions(a%Snapshots,ierr)
      call MatSetSizes(a%Snapshots,a%ndofn*a%npoinLocal,PETSC_DECIDE,a%ndofn*a%gnpoin,a%nsnap,ierr)
      call MatSetUp(a%Snapshots,ierr)
      iaux = a%ndofn*a%gnpoin*a%nsnap
      call a%Memor%allocObj(ierr,'Snapshots','SVDSystem',iaux)
      
      call a%Memor%alloc(a%ndofn,a%npoinLocal,a%SnapshotMean,'SnapshotsMean','SVDSystem')
      a%SnapshotMean = 0.0_rp

      call SVDCreate(a%MPIcomm,a%svd,ierr)
      call a%Memor%allocObj(ierr,'SVDSolver','SVDSystem',1)

      call SVDSetOptionsPrefix(a%svd,adjustl(trim(a%exmod)),ierr)
      call SVDSetFromOptions(a%svd,ierr)

   end subroutine

   subroutine SVDSetLumpedMass(a,vmass)
      implicit none
      class(SVDSystem) :: a
      integer(ip)      :: iaux,ipoin,ispos,idofn
      real(rp)         :: vmass(*)
      PC               :: PETSC_PC
      PetscErrorCode   :: ierr

      call MatCreate(a%MPIcomm,a%LMass,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofn*a%npoinLocal
      call a%Memor%allocObj(ierr,'LMass','SLEPcSystem',iaux)
      call MatSetFromOptions(a%LMass,ierr)
      call MatSetSizes(a%LMass,a%ndofn*a%npoinLocal,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,a%ndofn*a%gnpoin,ierr)
      call MatSetUp(a%LMass,ierr)
      call MatZeroEntries(a%LMass,ierr)

      call MatCreate(a%MPIcomm,a%invLMass,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofn*a%npoinLocal
      call a%Memor%allocObj(ierr,'invLMass','SLEPcSystem',iaux)
      call MatSetFromOptions(a%invLMass,ierr)
      call MatSetSizes(a%invLMass,a%ndofn*a%npoinLocal,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,a%ndofn*a%gnpoin,ierr)
      call MatSetUp(a%invLMass,ierr)
      call MatZeroEntries(a%invLMass,ierr)

      do ipoin = 1, a%npoinLocal
         do idofn = 1, a%ndofn
            ispos = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(ipoin))+idofn
            call MatSetValues(a%LMass,1,ispos-1,1,ispos-1,sqrt(1.0_rp/vmass(ipoin)),INSERT_VALUES,ierr)
         end do
      end do
      
      do ipoin = 1, a%npoinLocal
         do idofn = 1, a%ndofn
            ispos = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(ipoin))+idofn
            call MatSetValues(a%invLMass,1,ispos-1,1,ispos-1,sqrt(vmass(ipoin)),INSERT_VALUES,ierr)
         end do
      end do
      
      call MatAssemblyBegin(a%LMass,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%LMass,MAT_FINAL_ASSEMBLY,ierr)
      
      call MatAssemblyBegin(a%invLMass,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%invLMass,MAT_FINAL_ASSEMBLY,ierr)

   end subroutine
   
   subroutine SVDSetSnapshots(a,auxvec)
      implicit none
      class(SVDSystem) :: a
      real(rp)         :: auxvec(a%ndofn,a%npoinLocal)
      integer(ip)      :: ipoin,ispos(a%npoinLocal),idofn
      PetscErrorCode   :: ierr

      a%SnapshotMean = a%SnapshotMean + auxvec 

      ispos(1:a%npoinLocal) = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(1:a%npoinLocal))
      do idofn = 1, a%ndofn
         call MatSetValues(a%Snapshots,a%npoinLocal,ispos+idofn-1,1,a%isnap-1,auxvec(idofn,:),INSERT_VALUES,ierr)
      enddo

   end subroutine

   subroutine SVDSnapshotsMean(a,nsnap)
      implicit none
      class(SVDSystem) :: a
      integer(ip)      :: nsnap,isnap,ipoin,idofn,ispos(a%npoinLocal)
      PetscErrorCode   :: ierr

      a%SnapshotMean = a%SnapshotMean/nsnap
     
      ispos(1:a%npoinLocal) = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(1:a%npoinLocal))
      do isnap =0,a%nsnap-1
         do idofn = 1, a%ndofn
            call MatSetValues(a%Snapshots,a%npoinLocal,ispos+idofn-1,1,isnap,-1.0_rp*a%SnapshotMean(idofn,:),ADD_VALUES,ierr)
         enddo
      enddo
      
   end subroutine

   subroutine SVDSnapshotsMeanToZero(a)
      implicit none
      class(SVDSystem) :: a

      a%SnapshotMean = 0.0_rp
      
   end subroutine

   subroutine SVDWrite(a)
      implicit none
      class(SVDSystem) :: a
      PetscErrorCode :: ierr

      call SVDView(a%svd,a%outviewer,ierr)

   end subroutine

   subroutine SVDAssemblySnapshots(a)
      implicit none
      class(SVDSystem)     :: a
      PetscErrorCode       :: ierr
      
      call MatAssemblyBegin(a%Snapshots,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%Snapshots,MAT_FINAL_ASSEMBLY,ierr)
      
   end subroutine

   subroutine SVDAssemblySnapshotsMass(a)
      implicit none
      class(SVDSystem) :: a
      Mat              :: auxmat
      PetscErrorCode   :: ierr

      call MatMatMult(a%LMass,a%Snapshots,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,auxmat,ierr)
      call MatZeroEntries(a%Snapshots,ierr)
      call MatCopy(auxmat,a%Snapshots,DIFFERENT_NONZERO_PATTERN,ierr)
      call MatAssemblyBegin(a%Snapshots,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%Snapshots,MAT_FINAL_ASSEMBLY,ierr)
      
   end subroutine

   subroutine SVDSolveSys(a)
      implicit none
      class(SVDSystem) :: a
      integer(ip)      :: i,iaux,niter,nreas
      real(rp)         :: cputim2,cputim1
      Vec              :: u_aux1,v_aux1
      PetscErrorCode   :: ierr

      call cpu_time(cputim1)

      call SLEPcInitialize(PETSc_NULL_CHARACTER,ierr)

      call SVDSetOperators(a%svd,a%Snapshots,PETSC_NULL_MAT,ierr)
      call SVDSolve(a%svd,ierr)
      call SVDGetConverged(a%svd,a%nconv,ierr)
      a%ndofr  = a%nconv

      if (a%nconv .le. 0) then
          call runend('SVD found 0 basis vectors')
      endif

      call BVCreate(a%MPIcomm,a%Basis,ierr)
      call BVSetFromOptions(a%Basis,ierr)
      call BVSetSizes(a%Basis,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,a%nconv+1,ierr)
      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%allocObj(ierr,'Basis','SVDSystem',iaux)

      call BVCreate(a%MPIcomm,a%RBasis,ierr)
      call BVSetFromOptions(a%RBasis,ierr)
      call BVSetSizes(a%RBasis,PETSC_DECIDE,a%nconv,a%nconv+1,ierr)
      iaux = a%nconv*a%nconv
      call a%Memor%allocObj(ierr,'RBasis','SVDSystem',iaux)

      call a%Memor%alloc(a%nconv,a%sigma,'sigma','SVDSystem')

      call SVDGetBV(a%svd,a%RBasis,a%Basis,ierr)
      call BVSetActiveColumns(a%Basis,0,a%nconv,ierr)
      call BVSetActiveColumns(a%RBasis,0,a%nconv,ierr)

      do i =0, a%ndofr -1
         call SVDGetSingularTriplet(a%svd,i,a%sigma(i+1),PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
      end do
      call SVDGetIterationNumber(a%svd,niter,ierr);
      call SVDGetConvergedReason(a%svd,nreas,ierr);

      call cpu_time(cputim2)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_solve,500)
         write(a%lun_solve,100) niter
         write(a%lun_solve,101) nreas
         write(a%lun_solve,102) cputim2-cputim1
         write(a%lun_solve,501)
         if (a%kfl_flush == 1) call flush(a%lun_solve)
      endif

      500 format('SVD solver INFO')
      100 format('   Number of iterations: ',i12)
      101 format('   Converged Reason:     ',i12)
      102 format('   Solve time:           ',e12.6)
      501 format('------------------------------')
     
   end subroutine
   
   subroutine SVDMassMultBasis(a)
      implicit none
      class(SVDSystem) :: a
      BV               :: auxBasis
      PetscErrorCode   :: ierr

      call BVDuplicate(a%Basis,auxBasis,ierr)
      call BVCopy(a%Basis,auxBasis,ierr)
      call BVSetActiveColumns(auxBasis,0,a%nconv,ierr)
      call BVMatMult(auxBasis,a%invLMass,a%Basis,ierr)
      call BVDestroy(auxBasis,ierr)
      
   end subroutine

   subroutine SVDLMassDeallocate(a)
      implicit none
      class(SVDSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call MatDestroy(a%LMass,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofn*a%npoinLocal
      call a%Memor%deallocObj(ierr,'LMass','SLEPcSystem',iaux)
      call MatDestroy(a%invLMass,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofn*a%npoinLocal
      call a%Memor%deallocObj(ierr,'invLMass','SLEPcSystem',iaux)

   end subroutine

   subroutine SVDSystemDeallocate(a)
      implicit none
      class(SVDSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call SVDReset(a%svd,ierr)

      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%deallocObj(ierr,'Basis','SVDSystem',iaux)
      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%deallocObj(ierr,'RBasis','SVDSystem',iaux)

      call SVDDestroy(a%svd,ierr)
      call a%Memor%deallocObj(ierr,'SVDSolver','SVDSystem',1)
      
      call a%Memor%dealloc(a%ndofn,a%npoinLocal,a%SnapshotMean,'SnapshotsMean','SVDSystem')
      
      call a%Memor%dealloc(a%nconv,a%sigma,'sigma','SVDSystem')
      
      call MatDestroy(a%Snapshots,ierr)
      iaux = a%ndofn*a%gnpoin*a%nsnap
      call a%Memor%deallocObj(ierr,'Snapshots','SVDSystem',iaux)

   end subroutine

#endif

end module

module Mod_SLEPcSystem
   use typre
   use Mod_Memor
   use Mod_PETScSystem
   use Mod_PETScOrdering
   use Mod_ParallelSystemInterface
   use Mod_ParallelOrderingInterface
   use Mod_ParallelCommunicatorInterface
   use Mod_EigenSystemInterface
#ifdef SLEPC
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>
#include <slepc/finclude/slepcsvd.h>  
#include <slepc/finclude/slepcbv.h>  
#include <petsc/finclude/petscksp.h>
   use petscksp
   use petscmat
   use petscvec
   use slepcsvd
   use slepcbv
   implicit none
   private
   public SLEPcSystem, BasisPETScSystem

   type, extends(EigenSystemInterface) :: SLEPcSystem

      class(ParallelSystemInterface), pointer :: LinearSystem => NULL()
      class(PETScOrdering), pointer           :: LocalOrdering => NULL()

      type(MemoryMan), pointer :: Memor => NULL()

      type(tMat), pointer   :: AMatrix => NULL()
      type(tVec), pointer   :: RHS => NULL(), RHS_pointer => NULL(), Basis_pointer(:) => NULL(), Prec_pointer(:) => NULL()
      type(tVec), pointer   :: BasisU(:),ModBasis(:)
      type(tMat)            :: BasisV,APhi,Snapshots
      type(tVec)            :: SnapMean,ModRHS,RHSPhi,UNKNO
      type(tBV)             :: Basis,RBasis
      type(tKSP)            :: SOLVER
      PetscViewer    :: outviewer

      character(150) :: exmod,optionsfile
      integer(ip)    :: lun_solve,sizecut,numbasis,nsnap,SnapInterval

      integer(ip), pointer  :: isnap => NULL()
      real(rp), allocatable :: SnapshotMean(:,:),sigma(:)

   contains

      procedure :: SetLinearSystem     => SLEPcSetLinearSystem
      procedure :: SetMatrixPointers   => SLEPcSetMatrixPointers
      procedure :: SetRHSPointers      => SLEPcSetRHSPointers
      procedure :: SetOrderingPointers => SLEPcSetOrderingPointers
     
      procedure :: Init                => SLEPcInit
      procedure :: InitSystem          => SLEPcInitSolve
      procedure :: InitOperators       => SLEPcInitOperators
      procedure :: SetBasisSize        => SLEPcSetBasisSize
      procedure :: ResetBasisPointer   => SLEPcResetBasisPointer
      procedure :: BasisPointer
      procedure :: PrecPointer
      procedure :: RHSPointer
      
      procedure :: PetrovGalerkinProj  => SLEPcPetrovGalerkinProjection
      procedure :: GalerkinProj        => SLEPcGalerkinProjection

      procedure :: MeanToRHS           => SLEPcMeanToRHS

      procedure :: ToZero              => SLEPcToZero
      procedure :: SetupOperators      => SLEPcSetupOperators
      procedure :: MatProjection       => SLEPcMatProjection
      procedure :: VecProjection       => SLEPcVecProjection
      procedure :: SolveProjSystem     => SLEPcSolveProjSystem

      procedure :: WriteInfo           => SLEPcWriteInfo

      procedure :: SetNConvergence     => SLEPcSetConv
      procedure :: SetBasis            => SLEPcSetBasis
      procedure :: SetMean             => SLEPcSetMean
      procedure :: SetSigma            => SLEPcSetSigma
      procedure :: OrthogonalizeBasis  => SLEPcOrthogonalizeBasis
      procedure :: OrthoDeallocate     => SLEPcOrthogonalizeDeallocate

      procedure :: GetSigma            => SLEPcGetSigma
      procedure :: GetMean             => SLEPcGetMean
      procedure :: GetBasis            => SLEPcGetBasis
      procedure :: GetBasisR           => SLEPcGetBasisR
      procedure :: GetSolution         => SLEPcGetSolution
      procedure :: GetNConvergence     => SLEPcGetConv

      procedure :: DeallocOperators    => SLEPcDeallocateOperators
      procedure :: DeallocateRun       => SLEPcSystemDeallocate
      
      procedure :: SetBasisToLS        => SLEPcSetBasisToLS
   end type

   type, extends(PETScSystem) :: BasisPETScSystem
      type(tMat)               :: PMATRIX
      type(tVec)               :: PUNKNO,PRHS,PMATRIXD
      type(tVec), pointer      :: BasisU(:) => NULL(), auxBasis(:) => NULL(), ModBasis(:) => NULL()
      real(rp), pointer :: Basis(:,:,:) => NULL()
      integer           :: ndofr
      logical           :: kfl_massMatrix           !Mass matrix projection in basis
   contains
      procedure :: Init2       => Proj_Init
      procedure :: Solve       => Proj_Solve
      procedure :: GetSolution => Proj_GetSolution
      procedure :: Deallocate  => Proj_Deallocate
   end type

contains

   subroutine SLEPcInit(a,gnpoin,npoin,npoinLocal,ndofn,optionsfile,exmod,lun_solve,Memor)
      implicit none

      class(SLEPcSystem) :: a
      type(MemoryMan), target :: Memor
      integer(ip)    :: gnpoin,npoin,npoinLocal,ndofn,lun_solve,optionslen
      character(150) :: exmod,optionsfile
      character(120) :: ErrorString
      character(10)  :: string
      PetscErrorCode :: ierr
      
      a%Memor => Memor
      a%npoinLocal = npoinLocal
      a%npoin = npoin
      a%gnpoin = gnpoin
      a%ndofn = ndofn

      a%optionsfile = adjustl(trim(optionsfile))
      string = adjustl(trim(exmod))
      a%exmod = adjustl(trim(string))//'_'

      a%lun_solve = lun_solve

      optionslen = len(adjustl(trim(optionsfile)))
      optionsfile = adjustl(trim(optionsfile))
      call PetscOptionsInsertFile(a%MPIcomm,PETSC_NULL_OPTIONS,optionsfile(1:optionslen),PETSC_TRUE,ierr) 

      if (ierr /= 0) then
         ErrorString = 'SLEPc: '//adjustl(trim(optionsfile))//' not found'
         call runend(ErrorString)
      endif

   end subroutine

   subroutine SLEPcSetLinearSystem(a,LinearSystem)
      implicit none
      class(SLEPcSystem) :: a
      class(ParallelSystemInterface), pointer :: LinearSystem

      a%LinearSystem => LinearSystem

   end subroutine

   subroutine SLEPcSetMatrixPointers(a)
      implicit none
      class(SLEPcSystem) :: a
      class(ParallelSystemInterface),pointer:: LS => NULL()

      LS => a%LinearSystem
      select type(LS)
      class is(PETScSystem)
         a%AMatrix => LS%MATRIX
      end select

   end subroutine

   subroutine SLEPcSetRHSPointers(a)
      implicit none
      class(SLEPcSystem) :: a
      class(ParallelSystemInterface),pointer:: LS => NULL()

      LS => a%LinearSystem
      select type(LS)
      class is(PETScSystem)
         a%RHS => LS%RHS
      end select

      call a%RHSPointer(a%RHS)

   end subroutine

   subroutine SLEPcSetOrderingPointers(a)
      implicit none
      class(SLEPcSystem) :: a
      class(ParallelSystemInterface),pointer:: LS => NULL()

      LS => a%LinearSystem
      select type(LS)
      class is(PETScSystem)
         a%LocalOrdering => LS%LocalOrdering
      end select

   end subroutine

   subroutine SLEPcWriteInfo(a,fil_outpu)
      implicit none
      class(SLEPcSystem) :: a 
      character(150)   :: fil_outpu
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,a%outviewer,ierr)
      call PetscViewerFileSetMode(a%outviewer,FILE_MODE_WRITE,ierr)
      call PetscViewerASCIIOpen(a%MPIcomm,fil_outpu,a%outviewer,ierr)

      call a%SpecificWriter

      call PetscViewerDestroy(a%outviewer,ierr)

   end subroutine

   subroutine SLEPcSetConv(a,nconv)
      implicit none
      class(SLEPcSystem) :: a 
      integer(ip)        :: nconv

      a%ndofr = nconv
      a%nconv = nconv

   end subroutine

   subroutine SLEPcSetMean(a,SnapMean)
      implicit none
      class(SLEPcSystem) :: a 
      integer(ip)        :: ipoin,ispos,idofn
      real(rp)           :: SnapMean(a%ndofn,a%npoinLocal)
      PetscErrorCode     :: ierr

      do ipoin = 1, a%npoinLocal
         do idofn = 1, a%ndofn
            ispos = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(ipoin))+idofn
            call VecSetValues(a%SnapMean,1,ispos-1,SnapMean(idofn,ipoin),INSERT_VALUES,ierr)
         enddo
      enddo
      call VecAssemblyBegin(a%SnapMean,ierr)
      call VecAssemblyEnd(a%SnapMean,ierr)
    
   end subroutine

   subroutine SLEPcSetBasis(a,idofr,Basis)
      implicit none
      class(SLEPcSystem) :: a 
      integer(ip)        :: idofr,ipoin,ispos,idofn
      real(rp), target   :: Basis(a%ndofn,a%npoinLocal)
      PetscErrorCode     :: ierr

      call VecZeroEntries(a%BasisU(idofr),ierr)
      do ipoin = 1, a%npoinLocal
         do idofn = 1, a%ndofn
            ispos = (a%ndofn)*(a%LocalOrdering%IsLocal2Global(ipoin))+idofn
            call VecSetValues(a%BasisU(idofr),1,ispos-1,Basis(idofn,ipoin),INSERT_VALUES,ierr)
         enddo
      enddo
      call VecAssemblyBegin(a%BasisU(idofr),ierr)
      call VecAssemblyEnd(a%BasisU(idofr),ierr)

      call a%BasisPointer(a%BasisU)
      call a%PrecPointer(a%BasisU)

   end subroutine

   subroutine SLEPcSetSigma(a,sigma)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: sigma(a%nconv) 

      call a%Memor%alloc(a%nconv,a%sigma,'sigma','SLEPcSystem')
      a%sigma = sigma

   end subroutine
   
   subroutine SLEPcOrthogonalizeBasis(a)
      implicit none
      class(SLEPcSystem) :: a 
      integer(ip)        :: idofr,iaux
      real(rp)           :: aux(a%ndofn*a%npoinLocal)
      PetscErrorCode     :: ierr

      call BVCreate(a%MPIcomm,a%Basis,ierr)
      call BVSetFromOptions(a%Basis,ierr)
      call BVSetSizes(a%Basis,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,a%nconv,ierr)
      iaux = a%ndofn*a%gnpoin*a%nconv
      call a%Memor%allocObj(ierr,'BasisBV','SLEPcSystem',iaux)
      do idofr =1,a%ndofr
         call BVInsertVec(a%Basis,idofr-1,a%BasisU(idofr),ierr)
      enddo

      call BVOrthogonalize(a%Basis,PETSC_NULL_MAT,ierr)

      do idofr =1,a%ndofr
         call BVCopyVec(a%Basis,idofr-1,a%BasisU(idofr),ierr)
      enddo

   end subroutine

   subroutine SLEPcGetSigma(a,sigma)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: sigma(a%nconv) 

      sigma = a%sigma

   end subroutine

   subroutine SLEPcGetMean(a,Mean)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: Mean(a%ndofn,a%npoinLocal)

      Mean = a%SnapshotMean

   end subroutine

   subroutine SLEPcGetBasis(a,U)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: i
      real(rp)           :: U(a%ndofn*a%npoinLocal,a%ndofr)
      Vec                :: u_aux
      PetscErrorCode     :: ierr

      call VecCreate(a%MPIcomm,u_aux,ierr)
      call VecSetFromOptions(u_aux,ierr)
      call VecSetSizes(u_aux,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(u_aux,ierr)
      do i=0,a%ndofr-1
         call VecPlaceArray(u_aux,U(:,i+1),ierr)
         call BVCopyVec(a%Basis,i,u_aux,ierr)  !GetColumn may be faster
         call VecResetArray(u_aux,ierr)
      end do
      call VecDestroy(u_aux,ierr)

   end subroutine

   subroutine SLEPcGetBasisR(a,V)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: i
      real(rp)           :: V(a%ndofr,a%ndofr)
      Vec                :: v_aux,v_auxs
      VecScatter         :: ToAllScatter
      PetscErrorCode     :: ierr

      call VecCreate(a%MPIcomm,v_aux,ierr)
      call VecSetFromOptions(v_aux,ierr)
      call VecSetSizes(v_aux,PETSC_DECIDE,a%ndofr,ierr)
      call VecSetUp(v_aux,ierr)

      call VecScatterCreateToAll(v_aux,ToAllScatter,v_auxs,ierr)

      do i=0,a%ndofr-1
         call BVCopyVec(a%RBasis,i,v_aux,ierr)  !GetColumn may be faster
         call VecPlaceArray(v_auxs,V(:,i+1),ierr)
         call VecScatterBegin(ToAllScatter,v_aux,v_auxs,INSERT_VALUES,SCATTER_FORWARD,ierr)
         call VecScatterEnd(ToAllScatter,v_aux,v_auxs,INSERT_VALUES,SCATTER_FORWARD,ierr)
         call VecResetArray(v_auxs,ierr)
      end do
      call VecDestroy(v_aux,ierr)
      call VecDestroy(v_auxs,ierr)
      call VecScatterDestroy(ToAllScatter,ierr)

   end subroutine

   subroutine SLEPcToZero(a)
      implicit none
      class(SLEPcSystem) :: a
      PetscErrorCode     :: ierr

      call MatZeroEntries(a%APhi,ierr)
      call VecZeroEntries(a%RHSPhi,ierr)
   end subroutine

   subroutine SLEPcSetupOperators(a)
      implicit none
      class(SLEPcSystem) :: a
      PetscErrorCode     :: ierr

      call MatAssemblyBegin(a%AMatrix,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%AMatrix,MAT_FINAL_ASSEMBLY,ierr)
      call VecAssemblyBegin(a%RHS,ierr)
      call VecAssemblyEnd(a%RHS,ierr)
      call VecAssemblyBegin(a%RHS_pointer,ierr)
      call VecAssemblyEnd(a%RHS_pointer,ierr)

   end subroutine

   subroutine SLEPcInitOperators(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: iaux
      PetscErrorCode     :: ierr

      call VecCreate(a%MPIcomm,a%RHSPhi,ierr)
      call VecSetFromOptions(a%RHSPhi,ierr)
      call VecSetSizes(a%RHSPhi,PETSC_DECIDE,a%ndofr,ierr)
      call VecSetUp(a%RHSPhi,ierr)
      call VecAssemblyBegin(a%RHSPhi,ierr)
      call VecAssemblyEnd(a%RHSPhi,ierr)
      iaux = a%ndofr
      call a%Memor%allocObj(ierr,'RHSPHI','SLEPcSystem',iaux)
      
      call VecCreate(a%MPIcomm,a%UNKNO,ierr)
      call VecSetFromOptions(a%UNKNO,ierr)
      call VecSetSizes(a%UNKNO,PETSC_DECIDE,a%ndofr,ierr)
      call VecSetUp(a%UNKNO,ierr)
      call VecAssemblyBegin(a%UNKNO,ierr)
      call VecAssemblyEnd(a%UNKNO,ierr)
      iaux = a%ndofr
      call a%Memor%allocObj(ierr,'UNKNO','SLEPcSystem',iaux)

      call MatCreate(a%MPIcomm,a%APhi,ierr)
      call MatSetFromOptions(a%APhi,ierr)
      call MatSetSizes(a%APhi,PETSC_DECIDE,PETSC_DECIDE,a%ndofr,a%ndofr,ierr)
      call MatSetUp(a%APhi,ierr)
      call MatZeroEntries(a%APhi,ierr)
      call MatAssemblyBegin(a%APhi,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%APhi,MAT_FINAL_ASSEMBLY,ierr)
      iaux = a%ndofr*a%ndofr
      call a%Memor%allocObj(ierr,'APHI','SLEPcSystem',iaux)
   end subroutine

   subroutine SLEPcSetBasisSize(a,ndofr)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: ndofr
     
      a%ndofr = ndofr
   end subroutine

   subroutine BasisPointer(a,BasisU)
      implicit none
      class(SLEPcSystem) :: a
      Vec, target        :: BasisU(:)

      a%Basis_pointer => BasisU

   end subroutine

   subroutine PrecPointer(a,BasisU)
      implicit none
      class(SLEPcSystem) :: a
      Vec, target        :: BasisU(:)

      a%Prec_pointer => BasisU

   end subroutine

   subroutine SLEPcResetBasisPointer(a)
      implicit none
      class(SLEPcSystem) :: a

      call a%BasisPointer(a%BasisU)

   end subroutine

   subroutine RHSPointer(a,RHSp)
      implicit none
      class(SLEPcSystem) :: a
      Vec, target        :: RHSp

      a%RHS_pointer => RHSp

   end subroutine

   subroutine SLEPcPetrovGalerkinProjection(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: idofr
      PetscErrorCode     :: ierr

      do idofr = 1, a%ndofr
         call MatMult(a%AMatrix,a%BasisU(idofr),a%ModBasis(idofr),ierr)
      end do
      call a%BasisPointer(a%ModBasis)
      call a%PrecPointer(a%ModBasis)

   end subroutine

   subroutine SLEPcGalerkinProjection(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: idofr
      PetscErrorCode     :: ierr

      do idofr = 1, a%ndofr
         call MatMult(a%AMatrix,a%BasisU(idofr),a%ModBasis(idofr),ierr)
      end do
      call a%BasisPointer(a%ModBasis)

   end subroutine

   subroutine SLEPcMeanToRHS(a)
      implicit none
      class(SLEPcSystem) :: a
      Vec                :: auxvec,u_aux
      PetscErrorCode     :: ierr

      call VecDuplicate(a%RHS_pointer,auxvec,ierr)
      call VecCopy(a%RHS_pointer,auxvec,ierr)

      call MatMult(a%AMatrix,a%SnapMean,a%ModRHS,ierr)
      call VecAYPX(a%ModRHS,-1.0_rp,auxvec,ierr)
      
      call VecDestroy(auxvec,ierr)
      
      call a%RHSPointer(a%ModRHS)

   end subroutine

   subroutine SLEPcMatProjection(a)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: auxmat(a%ndofr)
      integer(ip)        :: idofr,j
      PetscErrorCode     :: ierr

      do idofr = 1,a%ndofr
         call VecMDot(a%Basis_pointer(idofr),a%ndofr,a%Prec_pointer,auxmat,ierr)
         call MatSetValues(a%APhi,a%ndofr,(/(j,j=0,a%ndofr-1)/),1,idofr-1,auxmat,INSERT_VALUES,ierr)
      end do
      call MatAssemblyBegin(a%APhi,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%APhi,MAT_FINAL_ASSEMBLY,ierr)

   end subroutine

   subroutine SLEPcVecProjection(a)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: RHSPhi(a%ndofr)
      integer(ip)        :: j
      PetscErrorCode     :: ierr
       
      call VecMDot(a%RHS_pointer,a%ndofr,a%Prec_pointer,RHSPhi,ierr)
      call VecSetValues(a%RHSPhi,a%ndofr,(/(j,j=0,a%ndofr-1)/),RHSPhi,INSERT_VALUES,ierr)
      call VecAssemblyBegin(a%RHSPhi,ierr)
      call VecAssemblyEnd(a%RHSPhi,ierr)

   end subroutine

   subroutine SLEPcInitSolve(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: iaux
      PC                 :: PETSC_PC
      PetscErrorCode     :: ierr

      call VecCreate(a%MPIcomm,a%ModRHS,ierr)
      call VecSetFromOptions(a%ModRHS,ierr)
      call VecSetSizes(a%ModRHS,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(a%ModRHS,ierr)
      iaux = a%ndofn*a%npoinLocal
      call a%Memor%allocObj(ierr,'ModRHS','SLEPcSystem',iaux)
      call VecZeroEntries(a%ModRHS,ierr)

      call VecDuplicateVecsF90(a%ModRHS,a%ndofr,a%BasisU,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%allocObj(ierr,'Basis','SLEPcSystem',iaux)

      call VecDuplicateVecsF90(a%ModRHS,a%ndofr,a%ModBasis,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%allocObj(ierr,'ModBasis','SLEPcSystem',iaux)

      call VecCreate(a%MPIcomm,a%SnapMean,ierr)
      call VecSetFromOptions(a%SnapMean,ierr)
      call VecSetSizes(a%SnapMean,a%ndofn*a%npoinLocal,a%ndofn*a%gnpoin,ierr)
      call VecSetUp(a%SnapMean,ierr)
      iaux = a%ndofn*a%npoinLocal
      call a%Memor%allocObj(ierr,'SnapMean','SLEPcSystem',iaux)
      call VecZeroEntries(a%SnapMean,ierr)
      
      call KSPCreate(a%MPIcomm,a%SOLVER,ierr)
      call KSPSetOptionsPrefix(a%SOLVER,adjustl(trim(a%exmod)),ierr)
      call KSPSetFromOptions(a%SOLVER,ierr)
      call a%Memor%allocObj(ierr,'Solver','SLEPcSystem',1)
      call KSPGetPC(a%SOLVER,PETSC_PC,ierr)
      call PCSetOptionsPrefix(PETSC_PC,adjustl(trim(a%exmod)),ierr)
      call PCSetFromOptions(PETSC_PC,ierr)

   end subroutine

   subroutine SLEPcSolveProjSystem(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip)        :: niter,ConvergedReason,i,j
      real(rp)           :: rnorm
      real(rp)           :: cputim2,cputim1
      real(rp), pointer  :: matpointer(:,:) => NULL(),vecpointer(:) => NULL()
      PetscErrorCode     :: ierr

      call cpu_time(cputim1)
      call KSPSetOperators(a%SOLVER,a%APhi,a%APhi,ierr)
      call KSPSolve(a%SOLVER,a%RHSPhi,a%UNKNO,ierr)

      call KSPGetConvergedReason(a%SOLVER,ConvergedReason,ierr)
      if (ierr < 0) then
         write(*,*) 'Warning, ROM solver did not converge'
      endif
      call KSPGetIterationNumber(a%SOLVER,niter,ierr)
      call KSPGetResidualNorm(a%SOLVER,rnorm,ierr)

      call cpu_time(cputim2)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_solve,500)
         write(a%lun_solve,100) niter
         write(a%lun_solve,101) rnorm
         write(a%lun_solve,102) cputim2-cputim1
         write(a%lun_solve,501)
         if (a%kfl_flush == 1) call flush(a%lun_solve)
      endif

      500 format('ROM solver INFO')
      100 format('   Number of iterations: ',i12)
      101 format('   Residual norm:        ',e12.6)
      102 format('   Solve time:           ',e12.6)
      501 format('------------------------------')
   end subroutine

   subroutine SLEPcGetSolution(a,Solution)
      implicit none
      class(SLEPcSystem) :: a
      real(rp)           :: Solution(a%ndofr)
      VecScatter         :: ToAllScatter
      Vec                :: UNKAll
      PetscErrorCode     :: ierr

      Solution = 0.0_rp
      call VecScatterCreateToAll(a%UNKNO,ToAllScatter,UNKAll,ierr)
      call VecPlaceArray(UNKAll,Solution,ierr)
      call VecScatterBegin(ToAllScatter,a%UNKNO,UNKAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(ToAllScatter,a%UNKNO,UNKAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecResetArray(UNKAll,ierr)
      call VecDestroy(UNKAll,ierr)
      call VecScatterDestroy(ToAllScatter,ierr)

   end subroutine

   subroutine SLEPcGetConv(a,ndofr)
      implicit none
      class(SLEPcSystem) :: a 
      integer(ip)      :: ndofr

      ndofr = a%nconv
  
   end subroutine

   subroutine SLEPcDeallocateOperators(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call MatDestroy(a%APhi,ierr)
      call VecDestroy(a%RHSPhi,ierr)
      iaux = a%ndofr*a%ndofr
      call a%Memor%deallocObj(ierr,'APHI','SLEPcSystem',iaux)
      iaux = a%ndofr
      call a%Memor%deallocObj(ierr,'RHSPHI','SLEPcSystem',iaux)
      call VecDestroy(a%UNKNO,ierr)
      call a%Memor%deallocObj(ierr,'UNKNO','SLEPcSystem',iaux)

      call KSPReset(a%SOLVER,ierr)

   end subroutine

   subroutine SLEPcSystemDeallocate(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call VecDestroy(a%ModRHS,ierr)
      iaux = a%ndofn*a%npoinLocal
      call a%Memor%deallocObj(ierr,'ModRHS','SLEPcSystem',iaux)
      
      call VecDestroy(a%SnapMean,ierr)
      iaux = a%ndofn*a%npoinLocal
      call a%Memor%deallocObj(ierr,'SnapMean','SLEPcSystem',iaux)
      
      call KSPDestroy(a%SOLVER,ierr)
      call a%Memor%deallocObj(ierr,'Solver','SLEPcSystem',1)

      call VecDestroyVecsF90(a%ndofr,a%BasisU,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%deallocObj(ierr,'Basis','SLEPcSystem',iaux)
      call VecDestroyVecsF90(a%ndofr,a%ModBasis,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%deallocObj(ierr,'ModBasis','SLEPcSystem',iaux)
      
      call a%Memor%dealloc(a%nconv,a%sigma,'sigma','SLEPcSystem')

   end subroutine
   
   subroutine SLEPcOrthogonalizeDeallocate(a)
      implicit none
      class(SLEPcSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      call BVDestroy(a%Basis,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%deallocObj(ierr,'BasisBV','SLEPcSystem',iaux)

   end subroutine


   subroutine SLEPcSetBasisToLS(a,LS,Basis)
      implicit none
      class(SLEPcSystem) :: a
      class(ParallelSystemInterface),pointer:: LS
      real(rp), target   :: Basis(a%ndofn,a%npoinLocal,a%ndofr)
      PetscErrorCode :: ierr

      select type(LS)
      class is(BasisPETScSystem)
         LS%BasisU => a%BasisU
         LS%Basis => Basis
         LS%ndofr = a%ndofr
         LS%ndofr  = a%ndofr
         LS%ndofn  = a%ndofn
         LS%npoinLocal = a%npoinLocal
      end select

   end subroutine

   subroutine Proj_Init(a,kfl_massMatrix,kfl_Precondition)
      implicit none
      class(BasisPETScSystem) :: a
      real(rp)    :: auxmat(a%ndofr)
      integer(ip) :: idofr,iaux,j
      integer(ip) :: kfl_Precondition    !Use a Precondition projection
      logical     :: kfl_massMatrix      !Mass matrix projection in basis
      PetscErrorCode :: ierr

      a%kfl_massMatrix = kfl_massMatrix
      call VecCreate(a%MPIcomm,a%PRHS,ierr)
      call VecSetFromOptions(a%PRHS,ierr)
      call VecSetSizes(a%PRHS,PETSC_DECIDE,a%ndofr,ierr)
      call VecSetUp(a%PRHS,ierr)
      call VecZeroEntries(a%PRHS,ierr)
      call VecAssemblyBegin(a%PRHS,ierr)
      call VecAssemblyEnd(a%PRHS,ierr)
      iaux = a%ndofr
      call a%Memor%allocObj(ierr,'PRHS','SLEPcSystem',iaux)
      
      call VecCreate(a%MPIcomm,a%PUNKNO,ierr)
      call VecSetFromOptions(a%PUNKNO,ierr)
      call VecSetSizes(a%PUNKNO,PETSC_DECIDE,a%ndofr,ierr)
      call VecSetUp(a%PUNKNO,ierr)
      call VecZeroEntries(a%PUNKNO,ierr)
      call VecAssemblyBegin(a%PUNKNO,ierr)
      call VecAssemblyEnd(a%PUNKNO,ierr)
      iaux = a%ndofr
      call a%Memor%allocObj(ierr,'PUNKNO','SLEPcSystem',iaux)

      call VecDuplicateVecsF90(a%RHS,a%ndofr,a%auxBasis,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%allocObj(ierr,'auxBasis','SLEPcSystem',iaux)

      call MatAssemblyBegin(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)
      !Projection of the matrix to the rom space
      do idofr = 1, a%ndofr
         call MatMult(a%MATRIX,a%BasisU(idofr),a%auxBasis(idofr),ierr)
      end do
      if (a%kfl_massMatrix) then
         call VecCreate(a%MPIcomm,a%PMATRIXD,ierr)
         call VecSetFromOptions(a%PMATRIXD,ierr)
         call VecSetSizes(a%PMATRIXD,PETSC_DECIDE,a%ndofr,ierr)
         call VecSetUp(a%PMATRIXD,ierr)
         call VecZeroEntries(a%PMATRIXD,ierr)
         call VecAssemblyBegin(a%PMATRIXD,ierr)
         call VecAssemblyEnd(a%PMATRIXD,ierr)
         iaux = a%ndofr
         call a%Memor%allocObj(ierr,'PMATRIXD','SLEPcSystem',iaux)
         do idofr = 1,a%ndofr
            !Diagonal form
            call VecDot(a%auxBasis(idofr),a%BasisU(idofr),auxmat(idofr),ierr)
            !call VecSetValue(a%PMATRIXD,idofr-1,1.0_rp,INSERT_VALUES,ierr)   !Identity matrix form
            call VecSetValue(a%PMATRIXD,idofr-1,1.0_rp/auxmat(idofr),INSERT_VALUES,ierr)  !Full multiplication
         end do
         a%ModBasis => a%BasisU
      elseif (kfl_Precondition .eq. 1) then
         call MatCreate(a%MPIcomm,a%PMATRIX,ierr)
         call MatSetFromOptions(a%PMATRIX,ierr)
         call MatSetSizes(a%PMATRIX,PETSC_DECIDE,PETSC_DECIDE,a%ndofr,a%ndofr,ierr)
         call MatSetUp(a%PMATRIX,ierr)
         call MatZeroEntries(a%PMATRIX,ierr)
         call MatAssemblyBegin(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         iaux = a%ndofr*a%ndofr
         call a%Memor%allocObj(ierr,'PMATRIX','SLEPcSystem',iaux)
         do idofr = 1,a%ndofr
            !Non-diagonal form
            call VecMDot(a%auxBasis(idofr),a%ndofr,a%auxBasis,auxmat,ierr)
            call MatSetValues(a%PMATRIX,a%ndofr,(/(j,j=0,a%ndofr-1)/),1,idofr-1,auxmat,INSERT_VALUES,ierr)
         end do
         a%ModBasis => a%auxBasis
         call MatAssemblyBegin(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
      else
         call MatCreate(a%MPIcomm,a%PMATRIX,ierr)
         call MatSetFromOptions(a%PMATRIX,ierr)
         call MatSetSizes(a%PMATRIX,PETSC_DECIDE,PETSC_DECIDE,a%ndofr,a%ndofr,ierr)
         call MatSetUp(a%PMATRIX,ierr)
         call MatZeroEntries(a%PMATRIX,ierr)
         call MatAssemblyBegin(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         iaux = a%ndofr*a%ndofr
         call a%Memor%allocObj(ierr,'PMATRIX','SLEPcSystem',iaux)
         do idofr = 1,a%ndofr
            !Non-diagonal form
            call VecMDot(a%auxBasis(idofr),a%ndofr,a%BasisU,auxmat,ierr)
            call MatSetValues(a%PMATRIX,a%ndofr,(/(j,j=0,a%ndofr-1)/),1,idofr-1,auxmat,INSERT_VALUES,ierr)
         end do
         a%ModBasis => a%BasisU
         call MatAssemblyBegin(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(a%PMATRIX,MAT_FINAL_ASSEMBLY,ierr)
      end if

   end subroutine
   
   subroutine Proj_Solve(a,unkno)
      implicit none
      class(BasisPETScSystem) :: a
      real(rp)    :: unkno(a%ndofn,a%npoinLocal),auxunkno(a%ndofr),rnorm,RHSVec(a%ndofr)
      integer(ip) :: niter,idofn,idofr,ipoin,ConvergedReason,j
      real(rp)    :: cputim2, cputim1
      PetscErrorCode :: ierr

      call cpu_time(cputim1)
      
      call VecAssemblyBegin(a%RHS,ierr)
      call VecAssemblyEnd(a%RHS,ierr)
      !Projection of the RHS to the rom space
      call VecMDot(a%RHS,a%ndofr,a%ModBasis,RHSVec,ierr)
      call VecSetValues(a%PRHS,a%ndofr,(/(j,j=0,a%ndofr-1)/),RHSVec,INSERT_VALUES,ierr)
      call VecAssemblyBegin(a%PRHS,ierr)
      call VecAssemblyEnd(a%PRHS,ierr)

      if (a%kfl_massMatrix) then
         call VecPointwiseMult(a%PUNKNO,a%PMATRIXD,a%PRHS,ierr)
      else
         call KSPSetOperators(a%SOLVER,a%PMATRIX,a%PMATRIX,ierr)
         call KSPSolve(a%SOLVER,a%PRHS,a%PUNKNO,ierr)
         call KSPGetConvergedReason(a%SOLVER,ConvergedReason,ierr)
         if (ierr < 0) then
            write(*,*) 'Warning, ROM mass matrix solver did not converge'
         endif
         call KSPGetIterationNumber(a%SOLVER,niter,ierr);
         call KSPGetResidualNorm(a%SOLVER,rnorm,ierr)
      end if

      unkno = 0.0_rp
      call VecPlaceArray(a%UNKNO,unkno,ierr)
      call a%GetSolution(auxunkno)
      do idofr = 1, a%ndofr
         call VecAXPY(a%UNKNO,auxunkno(idofr),a%BasisU(idofr),ierr)
      end do
      call VecResetArray(a%UNKNO,ierr)
      call VecResetArray(a%UNKNO,ierr)

      call cpu_time(cputim2)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_solve,500)
         write(a%lun_solve,100) niter
         write(a%lun_solve,101) rnorm
         write(a%lun_solve,102) cputim2-cputim1
         write(a%lun_solve,501)
         if (a%kfl_flush == 1) call flush(a%lun_solve)
      endif

      500 format('ROM mass matrix solver INFO')
      100 format('   Number of iterations: ',i12)
      101 format('   Residual norm:        ',e12.6)
      102 format('   Solve time:           ',e12.6)
      501 format('------------------------------')

   end subroutine
   
   subroutine Proj_GetSolution(a,Solution)
      implicit none
      class(BasisPETScSystem) :: a
      real(rp)           :: Solution(a%ndofr)
      VecScatter         :: ToAllScatter
      Vec                :: UNKAll
      PetscErrorCode     :: ierr

      Solution = 0.0_rp
      call VecScatterCreateToAll(a%PUNKNO,ToAllScatter,UNKAll,ierr)
      call VecPlaceArray(UNKAll,Solution,ierr)
      call VecScatterBegin(ToAllScatter,a%PUNKNO,UNKAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecScatterEnd(ToAllScatter,a%PUNKNO,UNKAll,INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecResetArray(UNKAll,ierr)
      call VecDestroy(UNKAll,ierr)
      call VecScatterDestroy(ToAllScatter,ierr)

   end subroutine

   subroutine Proj_Deallocate(a)
      implicit none
      class(BasisPETScSystem) :: a
      integer(ip) :: iaux
      PetscErrorCode :: ierr

      if (a%kfl_massMatrix) then
         call VecDestroy(a%PMATRIXD,ierr)
         iaux = a%ndofr
         call a%Memor%deallocObj(ierr,'PMATRIXD','SLEPcSystem',iaux)
      else
         call MatDestroy(a%PMATRIX,ierr)
         iaux = a%ndofr*a%ndofr
         call a%Memor%deallocObj(ierr,'PMATRIX','SLEPcSystem',iaux)
      end if
      call VecDestroy(a%PRHS,ierr)
      iaux = a%ndofr
      call a%Memor%deallocObj(ierr,'PRHS','SLEPcSystem',iaux)
      call VecDestroy(a%PUNKNO,ierr)
      call a%Memor%deallocObj(ierr,'PUNKNO','SLEPcSystem',iaux)

      call VecDestroyVecsF90(a%ndofr,a%auxBasis,ierr)
      iaux = a%ndofn*a%npoinLocal*a%ndofr
      call a%Memor%deallocObj(ierr,'auxBasis','SLEPcSystem',iaux)

      call MatDestroy(a%MATRIX,ierr);
      iaux = a%MatrixSize
      call a%Memor%deallocObj(ierr,'MATRIX','PETScSystem_Deallocate',iaux)

      call KSPDestroy(a%SOLVER,ierr);
      call a%Memor%deallocObj(ierr,'Solver','PETScSystem_Deallocate',1)

      call VecDestroy(a%UNKNO,ierr);
      call a%Memor%deallocObj(ierr,'UNKNO','PETScSystem_Deallocate',1)

      call VecDestroy(a%RHS,ierr);
      iaux = a%RhsSize
      call a%Memor%deallocObj(ierr,'RHS','PETScSystem_Deallocate',iaux)

      call a%LocalOrdering%Deallocate(a%Memor)
   end subroutine
#endif

end module

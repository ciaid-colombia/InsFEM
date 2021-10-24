module Mod_PETScSystem
   use typre
   use Mod_ParallelSystemInterface
   use Mod_ParallelOrderingInterface
   use Mod_Memor
   use Mod_Debugging
   use Mod_PETScOrdering
   use Mod_Element
#include <petsc/finclude/petscksp.h>
   use petscksp
   implicit none
   private
   public PETScSystem

   type, extends(ParallelSystemInterface) :: PETScSystem
      Mat :: MATRIX, auxMATRIX,Mass_Matrix,K_Matrix
      Vec :: RHS, auxRHS
      Vec :: UNKNO
      type(PETScOrdering) :: LocalOrdering
      KSP :: SOLVER

      integer(ip)    :: MatrixSize, RhsSize
      character(150) :: exmod,optionsfile
      integer(ip)    :: lun_solve
      logical       :: eigen = .false.

   contains

      procedure :: Init                 => PETScInitSystem
      procedure :: ToZero               => PETScSystemToZero
      procedure :: ToZeroMatrix         => PETScSystemToZeroMatrix
      procedure :: ToZeroRhs            => PETScSystemToZeroRhs
      procedure :: Assembly             => PETScAssembly
      procedure :: AssemblyRowBlock     => PETScAssemblyRowBlock
      procedure :: AssemblyElmat        => PETSCAssemblyElmat
      procedure :: AssemblyEigenMat     => PETScAssemblyEigenMatrix
      procedure :: AssemblyElrhs        => PETScAssemblyElrhs
      procedure :: AssemblyPointRhs     => PETScAssemblyPointRhs
      procedure :: AssemblyForceGhosts  => PETScAssemblyForceGhosts
      procedure :: Solve                => PETScSolve
      procedure :: Deallocate           => PETScSystemDeallocate
      procedure :: MatGetDiagonal       => PETScMatGetDiagonal
      procedure :: CopyToRHS            => PETScCopyToRHS
      procedure :: Matvec               => PETScMatvec
      procedure :: RhsDot               => PETScRhsDot
      procedure :: createMatrix
      procedure :: setMatrix
      procedure :: setfiles
      procedure :: LocalDimensionsAndOrdering
      procedure :: createVector
      procedure :: createSolver
      procedure :: createPreConditioner
   end type

contains

   subroutine PETScInitSystem(a,npoinLocal,npoinGhost,gnpoin,ndofn,issym,ia,ja,&
         &LocalOrdering,Memor,optionsfile,exmod,lun_solve,opti)
      implicit none
      class(PETScSystem) :: a
      type(MemoryMan), target   :: Memor
      class(ParallelOrderingInterface) :: LocalOrdering
      integer(ip) :: gnpoin,npoinLocal,npoinGhost,ndofn,issym,ia(*),ja(*)
      character(150) :: exmod,optionsfile
      integer(ip)   :: lun_solve
      integer(ip) :: ipoin,jpoin,poinj,iaux,npoin
      integer(ip), allocatable :: d_nnz(:), o_nnz(:),iwa(:),iwa2(:)
      PetscErrorCode     :: ierr
      integer(ip) :: MPIsize
      integer(ip) :: ufields(3),pfields
      character(5), optional :: opti

      a%Memor => Memor
      call a%setfiles(exmod,optionsfile,lun_solve)

      call a%LocalDimensionsAndOrdering(gnpoin,ndofn,npoinLocal,npoinGhost,LocalOrdering)
      
      !Non symmetric matrix
      if (issym == 0) then

         call a%createMatrix(a%MATRIX,'MATRIX',npoinLocal,gnpoin,ndofn,ia,ja)

         if(present(opti)) then
            if(opti == 'EIGEN') then
               !Special for eigenvalue problems
               a%eigen = .true.
               call a%createMatrix(a%Mass_Matrix,'MMATRI',npoinLocal,gnpoin,ndofn,ia,ja)
               call a%createMatrix(a%K_Matrix   ,'KMATRI',npoinLocal,gnpoin,ndofn,ia,ja)

               !call MatSetUp(a%Mass_Matrix,ierr)
               call a%setMatrix(a%Mass_Matrix,npoinLocal,gnpoin,ndofn,ia,ja)
               !call MatSetUp(a%K_Matrix,ierr)
               call a%setMatrix(a%K_Matrix   ,npoinLocal,gnpoin,ndofn,ia,ja)
            end if 
         end if 

         call a%setMatrix(a%MATRIX,npoinLocal,gnpoin,ndofn,ia,ja)

         !Unkno PETSC vector
         call a%createVector(npoinLocal,gnpoin,ndofn)

         !Solver 
         call a%createSolver

         !Preconditioner
         call a%createPreConditioner()

      else
         call runend('PETScSystem, Initialize: Not ready for symmetric matrices')

      endif
   end subroutine
   
   subroutine setfiles(a,exmod,optionsfile,lun_solve)
      use MPI
      implicit none
      class(PetscSystem) :: a
      character(150) :: exmod,optionsfile
      integer(ip)   :: lun_solve
      integer(ip) :: ierr,MPIsize,optionslen
      character(120) :: ErrorString
      
      a%exmod = adjustl(trim(exmod))//"_"
      a%optionsfile = adjustl(trim(optionsfile))
      a%lun_solve = lun_solve
      
      !Set OptionsFile for PETSC
      optionslen = len(adjustl(trim(optionsfile)))
      optionsfile = adjustl(trim(optionsfile))
      call PetscOptionsInsertFile(a%MPIcomm,PETSC_NULL_OPTIONS,optionsfile(1:optionslen),PETSC_TRUE,ierr)

      if (ierr /= 0) then
         ErrorString = 'PETScInitSystem: '//adjustl(trim(optionsfile))//' not found'
         call runend(ErrorString)
      endif
   end subroutine
   
   subroutine LocalDimensionsAndOrdering(a,gnpoin,ndofn,npoinLocal,npoinGhost,LocalOrdering)
      implicit none
      class(PetscSystem) :: a
      integer(ip) :: gnpoin,ndofn,npoinLocal,npoinGhost
      class(ParallelOrderingInterface) :: LocalOrdering
      PetscErrorCode     :: ierr
      integer(ip), allocatable :: d_nnz(:), o_nnz(:),iwa(:),iwa2(:)
      integer(ip) :: MPIsize,ipoin,npoin
      
      a%gnpoin = gnpoin
      a%ndofn = ndofn
      a%npoinLocal = npoinLocal
      a%npoinGhost = npoinGhost

      !Initialize local numbering
      npoin = npoinLocal+npoinGhost
      call a%Memor%alloc(npoin,iwa,'iwa','PETScSystem_Init')
      call a%Memor%alloc(npoin,iwa2,'iwa2','PETScSystem_Init')

      do ipoin = 1,npoin
         iwa(ipoin) = ipoin
         iwa2(ipoin) = ipoin
      enddo

      call LocalOrdering%Local2Global(npoin,iwa,iwa)
      call LocalOrdering%GetProc(npoin,iwa2,iwa2)
      call a%LocalOrdering%InitNdofn(a%MPIComm,npoin,iwa,iwa2,a%Memor,ndofn)
      call a%Memor%dealloc(npoin,iwa,'iwa','PETScSystem_Init')
      call a%Memor%dealloc(npoin,iwa2,'iwa2','PETScSystem_Init')
   end subroutine
      
   subroutine createVector(a,npoinLocal,gnpoin,ndofn)
      implicit none
      class(PETScSystem) :: a
      integer(ip) :: gnpoin,npoinLocal,ndofn
      PetscErrorCode     :: ierr
      integer(ip) :: iaux

      call VecCreateMPIWithArray(a%MPIcomm,ndofn,npoinLocal*ndofn,gnpoin*ndofn,PETSC_NULL_SCALAR,a%UNKNO,ierr)
      iaux = 1
      call a%Memor%allocObj(ierr,'UNKNO','PETScSystem_Init',iaux)
 
      !Right Hand Side Vector
      call VecDuplicate(a%UNKNO,a%RHS,ierr)
      call VecSetBlockSize(a%RHS,ndofn,ierr);
   
      call VecSetLocalToGlobalMapping(a%RHS,a%LocalOrdering%PIsLocal2Global,ierr)
      a%RhsSize = npoinLocal*ndofn*rp
      call a%Memor%allocObj(ierr,'RHS','PETScSystem_Init',a%RhsSize)

   end subroutine
 
   subroutine createSolver(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr

      call KSPCreate(a%MPIcomm,a%SOLVER,ierr)
      call KSPSetOptionsPrefix(a%SOLVER,adjustl(trim(a%exmod)),ierr);
      call KSPSetFromOptions(a%SOLVER,ierr)
      call a%Memor%allocObj(ierr,'Solver','PETScSystem_Init',1)

   end subroutine

   subroutine createPreConditioner(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr
      PC                 :: PETSC_PC
      call KSPGetPC(a%SOLVER,PETSC_PC,ierr)
      call PCSetOptionsPrefix(PETSC_PC,adjustl(trim(a%exmod)),ierr);
      call PCSetFromOptions(PETSC_PC,ierr)

   end subroutine

   subroutine createMatrix(a,auxMat,matName,npoinLocal,gnpoin,ndofn,ia,ja)
      implicit none
      class(PETScSystem) :: a
      Mat :: auxMat
      character(6)    :: matName
      integer(ip)     :: gnpoin,npoinLocal,ndofn,ia(*),ja(*)
      PetscErrorCode     :: ierr

      call MatCreate(a%MPIcomm,auxMat,ierr)
      call MatSetsizes(auxMat,npoinLocal*ndofn,npoinLocal*ndofn,gnpoin*ndofn,gnpoin*ndofn,ierr)
      call MatSetBlockSize(auxMat,ndofn,ierr);

      a%MatrixSize = npoinLocal*ip + (ia(npoinLocal+1)-1)*(ndofn*ndofn+1)*rp
      call a%Memor%allocObj(ierr,matName,'PETScSystem_Init',a%MatrixSize)

   end subroutine

   subroutine setMatrix(a,auxMat,npoinLocal,gnpoin,ndofn,ia,ja)
      implicit none
      class(PETScSystem) :: a
      Mat :: auxMat
      integer(ip) :: gnpoin,npoinLocal,ndofn,ia(*),ja(*)
      double  precision :: info(MAT_INFO_SIZE)
      double  precision :: mal, nz_a
      integer(ip) :: ipoin,jpoin,poinj,iaux,npoin
      integer(ip), allocatable :: d_nnz(:), o_nnz(:),iwa(:),iwa2(:)
      PetscErrorCode     :: ierr

      !Left Hand Side MATRIX
      call a%Memor%alloc(npoinLocal,d_nnz,'d_nnz','PETScSystem_Init')
      call a%Memor%alloc(npoinLocal,o_nnz,'o_nnz','PETScSystem_Init')

      d_nnz = 1
      do ipoin = 1,npoinLocal
         do poinj = ia(ipoin),ia(ipoin+1)-1
            jpoin = ja(poinj)
            if (jpoin <= npoinLocal) then
               d_nnz(ipoin) = d_nnz(ipoin)+1
            else
               o_nnz(ipoin) = o_nnz(ipoin)+1
            endif
         enddo
      enddo
      call MatSetOptionsPrefix(auxMat,adjustl(trim(a%exmod)),ierr);
      call MatSetFromOptions(auxMat,ierr);

      call MatXAIJSetPreallocation(auxMat,ndofn,d_nnz,o_nnz,d_nnz,d_nnz,ierr)
      call MatSetLocalToGlobalMapping(auxMat,a%LocalOrdering%PISLocal2Global,a%LocalOrdering%PIsLocal2Global,ierr)
      call MatSetOption(auxMat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,ierr)

      call a%Memor%dealloc(npoinLocal,d_nnz,'d_nnz','PETScSystem_Init')
      call a%Memor%dealloc(npoinLocal,o_nnz,'o_nnz','PETScSystem_Init')

      call MatGetInfo(auxMat,MAT_LOCAL,info,ierr)
      mal = info(MAT_INFO_MALLOCS)
      nz_a = info(MAT_INFO_NZ_ALLOCATED)

   end subroutine  

   subroutine PETScSystemToZero(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr

      call PETScSystemToZeroRhs(a)
      call MatZeroEntries(a%MATRIX,ierr)
   end subroutine
   
   subroutine PETScSystemToZeroMatrix(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr

      call MatZeroEntries(a%MATRIX,ierr)
   end subroutine

   subroutine PETScSystemToZeroRhs(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr

      call VecSet(a%RHS,0.0_rp,ierr)
   end subroutine

   subroutine PETScAssemblyForceGhosts(a,e,elmat,elrhs)
      implicit none
      class(PETScSystem)  :: a
      class(FiniteElement) :: e
      real(rp), target    :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)

      real(rp)    :: welmat(a%ndofn,e%pnode,a%ndofn,e%pnode),welrhs(a%ndofn,e%pnode)
      integer(ip) :: ilnods(e%pnode),jlnods(e%pnode),rlnods(e%pnode),inode,jdime,jnode,idime
      PetscErrorCode     :: ierr

      integer(ip) :: iaux(20)
     
      ilnods = e%lnods(1:e%pnode)-1
      
      !RHS
      call VecSetValuesBlockedLocal(a%RHS,e%pnode,ilnods,elrhs,ADD_VALUES,ierr)

      !ROW_ORIENTED
      forall (jdime = 1:a%ndofn,jnode = 1:e%pnode,idime=1:a%ndofn,inode = 1:e%pnode)
         welmat(jdime,jnode,idime,inode) = elmat(idime,inode,jdime,jnode)
      end forall
      call MatSetValuesBlockedLocal0(a%MATRIX,e%pnode,ilnods,e%pnode,ilnods,welmat,ADD_VALUES,ierr)

   end subroutine

   subroutine PETScAssembly(a,e,elmat,elrhs)
      implicit none
      class(PETScSystem)  :: a
      class(FiniteElement) :: e
      real(rp), target    :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)

      real(rp)    :: welmat(a%ndofn,e%pnode,a%ndofn,e%pnode)
      integer(ip) :: ilnods(e%pnode),jlnods(e%pnode),inode,jdime,jnode,idime
      PetscErrorCode     :: ierr

      call PETScAssemblyElmat(a,e,elmat)

      call PETScAssemblyElrhs(a,e,elrhs)
   end subroutine
   
   subroutine PETScAssemblyRowBlock(a,e,irowBlock,elmat,elrhs)
      implicit none
      class(PETScSystem)  :: a
      class(FiniteElement) :: e
      integer(ip) :: irowblock
      real(rp), target    :: elmat(a%ndofn,1,a%ndofn,e%mnode), elrhs(a%ndofn,1)

      call PETScAssemblyElmatRowBlock(a,e,irowBlock,elmat)

      call PETScAssemblyElrhsRowBlock(a,e,irowBlock,elrhs)
   end subroutine
   
   subroutine PETScAssemblyElMat(a,e,elmat)
      implicit none
      class(PETScSystem)  :: a
      class(FiniteElement) :: e
      real(rp), target    :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)

      real(rp)    :: welmat(a%ndofn,e%pnode,a%ndofn,e%pnode)
      integer(ip) :: ilnods(e%pnode),jlnods(e%pnode),inode,jdime,jnode,idime
      PetscErrorCode     :: ierr

      do inode = 1,e%pnode
         if (e%lnods(inode) <= a%npoinLocal) then
            ilnods(inode) = e%lnods(inode)-1
         else
            ilnods(inode) = -1
         endif
      enddo
      jlnods = e%lnods(1:e%pnode)-1

      !ROW_ORIENTED
      forall (jdime = 1:a%ndofn,jnode = 1:e%pnode,idime=1:a%ndofn,inode = 1:e%pnode)
         welmat(jdime,jnode,idime,inode) = elmat(idime,inode,jdime,jnode)
      end forall
      call MatSetValuesBlockedLocal0(a%MATRIX,e%pnode,ilnods,e%pnode,jlnods,welmat,ADD_VALUES,ierr)
   
   end subroutine
   
   subroutine PETScAssemblyElMatRowBlock(a,e,irowBlock,elmat)
      implicit none
      class(PETScSystem)  :: a
      class(FiniteElement) :: e
      real(rp), target    :: elmat(a%ndofn,1,a%ndofn,e%mnode)
      integer(ip) :: irowBlock

      real(rp)    :: welmat(a%ndofn,e%pnode,a%ndofn,1)
      integer(ip) :: ilnods(1),jlnods(e%pnode),inode,jdime,jnode,idime
      PetscErrorCode     :: ierr

      if (irowBlock <= a%npoinLocal) then
         ilnods = irowBlock -1
      else
         return
      endif
      jlnods = e%lnods(1:e%pnode)-1

      !ROW_ORIENTED
      forall (jdime = 1:a%ndofn,jnode = 1:e%pnode,idime=1:a%ndofn,inode = 1:1)
         welmat(jdime,jnode,idime,inode) = elmat(idime,inode,jdime,jnode)
      end forall

      call MatSetValuesBlockedLocal0(a%MATRIX,1,ilnods,e%pnode,jlnods,welmat,ADD_VALUES,ierr)

   end subroutine
   
   

   subroutine PETScAssemblyEigenMatrix(a,e,massmat,kmat)
      implicit none
      class(PETScSystem) :: a
      type(FiniteElement)   :: e
      real(rp), target :: massmat(a%ndofn,e%mnode,a%ndofn,e%mnode)
      real(rp), target :: kmat(a%ndofn,e%mnode,a%ndofn,e%mnode)
      real(rp)    :: wmassmat(a%ndofn,e%pnode,a%ndofn,e%pnode)
      real(rp)    :: wkmat(a%ndofn,e%pnode,a%ndofn,e%pnode)
      integer(ip) :: ilnods(e%pnode),jlnods(e%pnode),inode,jdime,jnode,idime
      integer     :: ierr
      integer(ip) :: iaux(20)

      do inode = 1,e%pnode
         if (e%lnods(inode) <= a%npoinLocal) then
            ilnods(inode) = e%lnods(inode)-1
         else
            ilnods(inode) = -1
         endif
      enddo
      jlnods = e%lnods(1:e%pnode)-1

      !ROW_ORIENTED
      forall (jdime = 1:a%ndofn,jnode = 1:e%pnode,idime=1:a%ndofn,inode = 1:e%pnode)
         wmassmat(jdime,jnode,idime,inode) = massmat(idime,inode,jdime,jnode)
            wkmat(jdime,jnode,idime,inode) =    kmat(idime,inode,jdime,jnode)
      end forall

      call MatSetValuesBlockedLocal0(a%Mass_Matrix,e%pnode,ilnods,e%pnode,jlnods,wmassmat,ADD_VALUES,ierr)
      call MatSetValuesBlockedLocal0(a%K_Matrix   ,e%pnode,ilnods,e%pnode,jlnods,   wkmat,ADD_VALUES,ierr)

   end subroutine

   subroutine PETScAssemblyPointRhs(a,rhSize,pointId,pointRhs)
      implicit none
      class(PETScSystem) :: a
      integer(ip) :: rhSize,pointId(rhSize)
      real(rp)    :: pointRhs(a%ndofn,rhSize)
      integer     :: nod,count,nd
      integer(ip) :: localIndex(rhSize)
      real(rp)    :: localPointRhs(a%ndofn,rhSize)
      PetscErrorCode     :: ierr

      call VecSetValuesBlockedLocal(a%RHS,rhSize,pointId-1,pointRhs,ADD_VALUES,ierr)

   end subroutine

   subroutine PETScAssemblyElrhs(a,e,elrhs)
      implicit none
      class(PETScSystem) :: a
      class(FiniteElement)   :: e
      real(rp), target :: elrhs(a%ndofn,e%mnode)
      real(rp)    :: welrhs(a%ndofn,e%pnode)
      integer(ip) :: ilnods(e%pnode),rlnods(e%pnode),inode,nodecount
      PetscErrorCode     :: ierr

      nodecount = 0
      do inode = 1,e%pnode
         if (e%lnods(inode) <= a%npoinLocal) then
            ilnods(inode) = e%lnods(inode)-1
            nodecount = nodecount+1
            welrhs(:,nodecount) = elrhs(:,inode)
            rlnods(nodecount) = ilnods(inode)
         else
            ilnods(inode) = -1
         endif
      enddo
      
      !RHS
      call VecSetValuesBlockedLocal(a%RHS,nodecount,rlnods,welrhs,ADD_VALUES,ierr)
      
   end subroutine
   
   subroutine PETScAssemblyElrhsRowBlock(a,e,irowBlock,elrhs)
      implicit none
      class(PETScSystem) :: a
      class(FiniteElement)   :: e
      real(rp), target :: elrhs(a%ndofn,1)
      integer(ip) :: irowBlock
      
      integer(ip) :: ilnods(1),rlnods(e%pnode),inode,nodecount
      PetscErrorCode     :: ierr

      
      if (irowBlock >a%npoinLocal) return
      ilnods(1) = irowBlock-1
      
      !RHS
      call VecSetValuesBlockedLocal(a%RHS,1,ilnods,elrhs,ADD_VALUES,ierr)
      
   end subroutine

   subroutine PETScCopyToRHS(a,rhsvec)
      class(PETScSystem) :: a
      real(rp) :: rhsvec(a%ndofn,a%npoinLocal)
      integer(ip) :: j
      PetscErrorCode     :: ierr

      call VecAssemblyBegin(a%RHS,ierr)
      call VecAssemblyEnd(a%RHS,ierr)
      call VecSetValuesBlockedLocal(a%RHS,a%npoinLocal,(/(j,j=0,a%npoinLocal-1)/),rhsvec,INSERT_VALUES,ierr)
   
   end subroutine

   subroutine PETScSolve(a,unkno)
      implicit none
      class(PETScSystem) :: a
      real(rp)    :: unkno(a%ndofn,a%npoinLocal)
      PetscViewer :: viewer, rhsviewer, unknoviewer
      integer(ip) :: niter
      real(rp)    :: rnorm
      PetscErrorCode     :: ierr
      integer(ip), save :: ipass = 0
      real(rp) :: cputim2, cputim1, cputim3, cputim4

      integer(ip) :: ConvergedReason

      call cpu_time(cputim1)

      call cpu_time(cputim3)
      call MatAssemblyBegin(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)

      if(a%eigen) then
         call MatAssemblyBegin(a%K_Matrix,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(  a%K_Matrix,MAT_FINAL_ASSEMBLY,ierr)

         call MatAssemblyBegin(a%Mass_Matrix,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(  a%Mass_Matrix,MAT_FINAL_ASSEMBLY,ierr)
      end if

      call VecAssemblyBegin(a%RHS,ierr)
      call VecAssemblyEnd(a%RHS,ierr)
      call cpu_time(cputim4)

      !Set ipass in order to debug matrix   
      ipass = ipass+1
      if (deb_PostProcessMatrix == 1) then
         call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.m",viewer,ierr);
         call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
         call MatView(a%MATRIX,viewer,ierr);

         call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vec.m",rhsviewer,ierr);
         call PetscViewerSetFormat(rhsviewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
         call VecView(a%RHS,rhsviewer,ierr);

         call MPI_BARRIER(a%MPIcomm,ierr)

         stop
      endif

      call cpu_time(cputim3)
      !Continue
      
      call VecPlaceArray(a%UNKNO,unkno,ierr)

      call KSPSetOperators(a%SOLVER,a%MATRIX,a%MATRIX,ierr)

      call cpu_time(cputim4)

      call cpu_time(cputim3)

      call KSPSolve(a%SOLVER,a%RHS,a%UNKNO,ierr)

      call cpu_time(cputim4)

      call cpu_time(cputim3)
      call KSPGetConvergedReason(a%SOLVER,ConvergedReason,ierr)
      if (ierr < 0) then
         write(*,*) 'Warning, solver did not converge'
      endif

      call cpu_time(cputim4)

      call cpu_time(cputim3)
      call KSPGetIterationNumber(a%SOLVER, niter,ierr);
      call KSPGetResidualNorm(a%SOLVER,rnorm,ierr)
      call VecResetArray(a%UNKNO,ierr)
      call cpu_time(cputim4)
      call cpu_time(cputim2)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_solve,100) niter
         write(a%lun_solve,101) rnorm
         write(a%lun_solve,102) cputim2-cputim1
         if (a%kfl_flush == 1) call flush(a%lun_solve)
      endif

      100 format('      Number of iterations: ',i12)
      101 format('      Residual norm:        ',e12.6)
      102 format('      Solve time:           ',e12.6)

   end subroutine

   subroutine PETScMatGetDiagonal(a,ipoin,diag,iwa,rwa)
      implicit none
      class(PETScSystem) :: a
      integer(ip) :: ipoin
      real(rp) :: diag(a%ndofn)
      integer(ip) :: iwa(*)
      real(rp)    :: rwa(*)
      PetscErrorCode     :: ierr
      integer(ip) :: gipoin,ipos,icol,idofn,irow,ncols

      do idofn = 1,a%ndofn
         call a%LocalOrdering%Local2Global(ipoin,gipoin)
         irow = a%ndofn*(gipoin-1)+idofn-1
         call MatGetRow(a%MATRIX,irow,ncols,iwa,rwa,ierr)
         ipos = 0
         icol = -1
         do while (icol /= irow)
            ipos = ipos+1
            icol = iwa(ipos)
         enddo
         diag(idofn) = rwa(ipos)
         call MatRestoreRow(a%MATRIX,irow,ncols,iwa,rwa,ierr)
      enddo

   end subroutine   

   subroutine PETScSystemDeallocate(a)
      implicit none
      class(PETScSystem) :: a
      PetscErrorCode     :: ierr
      integer(ip) :: iaux

      call MatDestroy(a%MATRIX,ierr);
      iaux = a%MatrixSize
      call a%Memor%deallocObj(ierr,'MATRIX','PETScSystem_Deallocate',iaux)

      if(a%eigen) then
         call MatDestroy(a%Mass_Matrix,ierr);
         iaux = a%MatrixSize
         call a%Memor%deallocObj(ierr,'MMATRI','PETScSystem_Deallocate',iaux)

         call MatDestroy(a%K_Matrix,ierr);
         iaux = a%MatrixSize
         call a%Memor%deallocObj(ierr,'KMATRI','PETScSystem_Deallocate',iaux)
      end if

      call KSPDestroy(a%SOLVER,ierr);
      call a%Memor%deallocObj(ierr,'Solver','PETScSystem_Deallocate',1)

      call VecDestroy(a%UNKNO,ierr);
      call a%Memor%deallocObj(ierr,'UNKNO','PETScSystem_Deallocate',1)

      call VecDestroy(a%RHS,ierr);
      iaux = a%RhsSize
      call a%Memor%deallocObj(ierr,'RHS','PETScSystem_Deallocate',iaux)

      call a%LocalOrdering%Deallocate(a%Memor)
   end subroutine

   subroutine PETScMatvec(a,vec,vecout)
      implicit none
      class(PetscSystem) :: a
      real(rp) :: vec(a%ndofn,a%npoinLocal), vecout(a%ndofn,a%npoinLocal)
      PetscErrorCode     :: ierr

      call MatAssemblyBegin(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(a%MATRIX,MAT_FINAL_ASSEMBLY,ierr)

      !We use the RHS and UNKNO headers
      call VecPlaceArray(a%RHS,vec,ierr)
      call VecPlaceArray(a%UNKNO,vecout,ierr)

      call MatMult(a%MATRIX,a%RHS,a%UNKNO,ierr)

      call VecResetArray(a%RHS,ierr)
      call VecResetArray(a%UNKNO,ierr)
   end subroutine

   subroutine PETScRhsDot(a,vec,vdot)
      implicit none
      class(PetscSystem) :: a
      real(rp) :: vec(a%ndofn,a%npoinLocal), vdot
      PetscErrorCode     :: ierr
      
      !VecAssembly, shouldn't do anything
      call VecAssemblyBegin(a%RHS,ierr)
      call VecAssemblyEnd(a%RHS,ierr)

      !We use unkno header
      call VecPlaceArray(a%UNKNO,vec,ierr)

      !Compute dot against RHS
      call VecDot(a%RHS,a%UNKNO,vdot,ierr);

      !We give unkno back
      call VecResetArray(a%UNKNO,ierr)

   end subroutine
end module

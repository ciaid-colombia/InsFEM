module Mod_EigenSystemInterface
   use typre
   use Mod_Memor
   use Mod_ParallelOrderingInterface
   use Mod_ParallelSystemInterface
   use Mod_ReadWrite
   implicit none
   private
   public EigenSystemInterface
   
   type :: EigenSystemInterface
      !For svd solver
      integer(ip) :: ndofn,gndofn,nconv,ndofr,npoin,npoinLocal,gnpoin
      integer(ip) :: kfl_flush = 0   !Default is do not flush
      integer(ip) :: MPIsize, MPIrank, MPIroot,MPIcomm

   contains

      procedure :: SetMPI
      procedure :: SetFlush

      procedure :: SetBasisToLS

      procedure :: SetLinearSystem
      procedure :: SetMatrixPointers
      procedure :: SetMatrixPointersBuild
      procedure :: SetRHSPointers
      procedure :: SetOrderingPointers
      
      procedure :: Init
      procedure :: InitSystem
      procedure :: InitOperators
      procedure :: SetBasisSize
      procedure :: ResetBasisPointer

      procedure :: MeanToRHS
      procedure :: PetrovGalerkinProj
      procedure :: GalerkinProj

      procedure :: ToZero
      procedure :: SetupOperators
      procedure :: MatProjection
      procedure :: VecProjection
      procedure :: SolveProjSystem

      procedure :: WriteInfo

      procedure :: SetNConvergence
      procedure :: SetBasis
      procedure :: SetMean
      procedure :: SetSigma
      procedure :: OrthogonalizeBasis
      procedure :: OrthoDeallocate

      procedure :: GetSigma
      procedure :: GetMean
      procedure :: GetBasis
      procedure :: GetBasisR
      procedure :: GetSolution
      procedure :: GetNConvergence

      procedure :: DeallocOperators
      procedure :: DeallocateRun
      
      !Specific solver
      procedure :: InitBuild
      procedure :: SolveSystem
      procedure :: SpecificWriter
      procedure :: LMassDeallocate
      procedure :: Deallocate

      procedure :: SnapshotsMean
      procedure :: SnapshotsMeanToZero
      procedure :: SetSnapshots
      
      procedure :: AssemblySnapshots
      procedure :: SetLumpedMass
      procedure :: AssemblySnapshotsMass
      procedure :: MassMultBasis

   end type

contains

   subroutine SetMPI(a,comm,size,root,irank)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: irank,size,root,comm
      
      a%MPIcomm = comm
      a%MPIrank = irank
      a%MPIsize = size
      a%MPIroot = root
   end subroutine
   
   subroutine SetFlush(a,kfl_flush)
      use typre
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: kfl_flush
      
      a%kfl_flush = kfl_flush
   end subroutine
      
   subroutine InitBuild(a,isnap,SnapInterval,nsnap)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip)      :: nsnap,SnapInterval
      integer(ip), target :: isnap

      call runend('EigenSystem, InitBuild procedure not implemented')

   end subroutine  

   subroutine SolveSystem(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SolveSystem procedure not implemented')

   end subroutine  

   subroutine SetLumpedMass(a,vmass)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)         :: vmass(*)

      call runend('EigenSystem, SetLumpedMass procedure not implemented')

   end subroutine

   subroutine SnapshotsMean(a,nsnap)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: nsnap

      call runend('EigenSystem, SnapshotsMean procedure not implemented')

   end subroutine

   subroutine SnapshotsMeanToZero(a)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: nsnap

      call runend('EigenSystem, SnapshotsMeanToZero procedure not implemented')

   end subroutine

   subroutine SetSnapshots(a,auxvec)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp) :: auxvec(a%ndofn,a%npoinLocal)

      call runend('EigenSystem, SetSnapshots procedure not implemented')

   end subroutine

   subroutine AssemblySnapshots(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, AssemblySnapshots procedure not implemented')

   end subroutine

   subroutine AssemblySnapshotsMass(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, AssemblySnapshotsMass procedure not implemented')

   end subroutine

   subroutine MassMultBasis(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, MassMultBasis procedure not implemented')

   end subroutine

   subroutine Init(a,gnpoin,npoin,npoinLocal,ndofn,optionsfile,exmod,lun_solve,Memor)
      implicit none
      class(EigenSystemInterface) :: a
      type(MemoryMan), target :: Memor
      integer(ip)    :: gnpoin,npoin,npoinLocal,ndofn,lun_solve
      character(150) :: exmod,optionsfile

      call runend('EigenSystem, Init procedure not implemented')

   end subroutine  

   subroutine SetLinearSystem(a,LinearSystem)
      implicit none
      class(EigenSystemInterface) :: a
      class(ParallelSystemInterface), pointer :: LinearSystem

      call runend('EigenSystem, setLinearSystem procedure not implemented')

   end subroutine

   subroutine SetBasisToLS(a,LS,Basis)
      implicit none
      class(EigenSystemInterface) :: a
      class(ParallelSystemInterface), pointer :: LS
      real(rp), target   :: Basis(a%ndofn,a%npoinLocal,a%ndofr)

      call runend('EigenSystem, SetBasisToLS procedure not implemented')

   end subroutine

   subroutine SetMatrixPointersBuild(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetMatrixPointersBuild procedure not implemented')

   end subroutine

   subroutine SetMatrixPointers(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetMatrixPointers procedure not implemented')

   end subroutine

   subroutine SetRHSPointers(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetRHSPointers procedure not implemented')

   end subroutine

   subroutine SetOrderingPointers(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetOrdringPointers procedure not implemented')

   end subroutine

   subroutine InitSystem(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, InitSystem procedure not implemented')

   end subroutine

   subroutine InitOperators(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, InitOperators procedure not implemented')

   end subroutine

   subroutine ToZero(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, ToZero procedure not implemented')

   end subroutine

   subroutine SetupOperators(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetupOperators procedure not implemented')

   end subroutine

   subroutine MatProjection(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, MatProjection procedure not implemented')

   end subroutine

   subroutine VecProjection(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, VecProjection procedure not implemented')

   end subroutine

   subroutine SolveProjSystem(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SolveProjSystem procedure not implemented')

   end subroutine

   subroutine SetBasisSize(a,ndofr)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip)                 :: ndofr

      call runend('EigenSystem, SetBasisSize procedure not implemented')

   end subroutine

   subroutine ResetBasisPointer(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, ResetBasisPointer procedure not implemented')

   end subroutine

   subroutine MeanToRHS(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, MeanToRHS procedure not implemented')

   end subroutine

   subroutine PetrovGalerkinProj(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, PetrovGalerkinProj procedure not implemented')

   end subroutine

   subroutine GalerkinProj(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, GalerkinProj procedure not implemented')

   end subroutine

   subroutine SetWriter(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SetWriter procedure not implemented')

   end subroutine

   subroutine DeallocateWriter(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, DeallocateWriter procedure not implemented')

   end subroutine

   subroutine WriteInfo(a,fil_outpu)
      implicit none
      class(EigenSystemInterface) :: a
      character(150)   :: fil_outpu

      call runend('EigenSystem, WriteInfo procedure not implemented')

   end subroutine
   
   subroutine SpecificWriter(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, SpecificWriter procedure not implemented')

   end subroutine

   subroutine SetNConvergence(a,nconv)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: nconv

      call runend('EigenSystem, SetNConvergence procedure not implemented')

   end subroutine

   subroutine SetBasis(a,idofr,Basis)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip)      :: idofr
      real(rp), target :: Basis(a%ndofn,a%npoinLocal)

      call runend('EigenSystem, SetBasis procedure not implemented')

   end subroutine

   subroutine SetMean(a,SnapMean)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)       :: SnapMean(a%ndofn,a%npoinLocal)

      call runend('EigenSystem, SetSnapMean procedure not implemented')

   end subroutine

   subroutine SetSigma(a,sigma)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)           :: sigma(a%nconv) 

      call runend('EigenSystem, SetSigma procedure not implemented')

   end subroutine
   
   subroutine OrthogonalizeBasis(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, OrthogonalizeBasis procedure not implemented')

   end subroutine

   subroutine OrthoDeallocate(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, OrthoDeallocate procedure not implemented')

   end subroutine

   subroutine GetSigma(a,sigma)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)           :: sigma(a%nconv) 

      call runend('EigenSystem, GetSigma procedure not implemented')

   end subroutine

   subroutine GetMean(a,Mean)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)         :: Mean(a%ndofn,a%npoinLocal)

      call runend('EigenSystem, GetMean procedure not implemented')

   end subroutine

   subroutine GetBasis(a,U)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)         :: U(a%ndofn*a%npoinLocal,a%ndofr)

      call runend('EigenSystem, GetBasis procedure not implemented')

   end subroutine

   subroutine GetBasisR(a,V)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)         :: V(a%ndofr,a%ndofr)

      call runend('EigenSystem, GetBasisR procedure not implemented')

   end subroutine

   subroutine GetSolution(a,Solution)
      implicit none
      class(EigenSystemInterface) :: a
      real(rp)           :: Solution(a%ndofr)

      call runend('EigenSystem, GetSolution procedure not implemented')

   end subroutine

   subroutine GetNConvergence(a,ndofr)
      implicit none
      class(EigenSystemInterface) :: a
      integer(ip) :: ndofr

      call runend('EigenSystem, GetNConvergence procedure not implemented')

   end subroutine

   subroutine LMassDeallocate(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, LMassDeallocate procedure not implemented')

   end subroutine  

   subroutine Deallocate(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, Deallocate procedure not implemented')

   end subroutine  

   subroutine DeallocOperators(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, DeallocOperators procedure not implemented')

   end subroutine  
   
   subroutine DeallocateRun(a)
      implicit none
      class(EigenSystemInterface) :: a

      call runend('EigenSystem, DeallocateRun procedure not implemented')

   end subroutine  
   
end module

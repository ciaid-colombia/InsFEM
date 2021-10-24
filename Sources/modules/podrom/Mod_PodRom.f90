module Mod_PodRom
   use typre
   use Mod_Timer
   use Mod_ReadWrite
   use Mod_MPIObject
   use Mod_ParallelLibraryInterface
   use Mod_ParallelSystemInterface
   use Mod_EigenLibraryInterface
   use Mod_EigenSystemInterface
   use Mod_Mesh
   use Mod_AdaptiveInterface
   use Mod_MeshInterpolator
   use Mod_PhysicalProblem
   implicit none
   private
   public PodRomProblem, PodRomProblem_Const
   
   type ROMTime
      type(Timer) :: Input
      type(Timer) :: Assembly
      type(Timer) :: Solve
      type(Timer) :: Output
      type(Timer) :: ReadSnapshots
      type(Timer) :: SetSnapshots
   end type

   type, extends(MPIObject) :: PodRomProblem

      type(ROMTime)                :: Timer
   
      integer(ip) :: lun_outro, lun_solro          !Output unit

      class(PhysicalProblem), pointer :: Problem => NULL()   !PhysicalProblem over which I am running the ROM
      character(5)                    :: ProblemType
      class(FemMesh), pointer         :: Mesh => NULL(), OldMesh => NULL()  !Mesh over which the PhysicalProblem runs
      class(Reader),pointer           :: Readerpr => NULL()
      class(Writer),pointer           :: Writerpr => NULL()
      character(150)                  :: exmod
      character (len=150),dimension(40):: nameBasis
      
      !Folder for storing the rom basis
      character(150) :: RestartFolder, OldRestartFolder
      character(150) :: fil_basis, fil_snap, fil_mesh
      integer(ip)    :: lun_basis, lun_snap, lun_mesh
      character(150) :: oldnamda      
      
      !Current time
      real(rp) :: ctime
      
      class(EigenLibraryInterface), pointer   :: EigenLibrary => NULL()
      class(EigenSystemInterface) , pointer   :: EigenSystem  => NULL()
      
      character(3) :: kfl_itask             !FOM, ROM, SNP, BAS
      character(3) :: kfl_eigentype         !SVD, EPS
      
      integer(ip)  :: nsnap                 !Total Number of Snapshots
      integer(ip)  :: isnap                 !Current Snapshot
      integer(ip)  :: wsnap                 !Saved number of snapshots
      integer(ip)  :: SnapshotsInterval     !Take snapshots every # of steps
      
      integer(ip)  :: numBasisFiles         !#Multiple basis files
      integer(ip)  :: inumBasis             !Current snapshots from multiple basis files
      
      integer(ip)  :: ndofn                 !Number of degrees of freedom
      integer(ip)  :: ndofr                 !Number of degrees of freedom for the reduced order model
      integer(ip)  :: nconv                 !Max number of degrees of freedom for the reduced order model
      integer(ip)  :: nstep                 !Number of time steps to be performed by the reduced-order model
      integer(ip)  :: nNonLinearIterations  !Number of non-linear iterations to be performed
      integer(ip)  :: basis_number          !number of basis vectors
      real(rp)     :: cotol                 !Convergence tolerance for ROM
      real(rp)     :: basis_energy          !Relative energy retained by the basis
      real(rp)     :: bEnergyT              !Total energy
      real(rp), allocatable :: bEnergy(:)   !Vector energy
      integer(ip)  :: npara                 !Number of parameters
      real(rp)     :: param                 !Parameter for the run
      real(rp)     :: parfun(2,2)           !Snapshot mean function
      
      real(rp), allocatable :: Basis(:,:,:)
      real(rp), allocatable :: SnapMean(:,:)
      real(rp), allocatable :: Snapshots(:,:)
      real(rp), allocatable :: FOMSolution(:,:)
      real(rp), allocatable :: sigma(:)
      real(rp), allocatable :: OldMeshBasis(:,:,:)
      real(rp), allocatable :: OldMeshSnapMean(:,:)

      !Options
      integer(ip) :: kfl_Precondition   !Use a Precondition projection

      logical :: kfl_SubstractMean  !Substract the mean before building the basis
      logical :: kfl_massMatrix     !Mass matrix projection in basis
      logical :: kfl_BkpOpen        !Snapshots old
      logical :: kfl_outBasis       !Basis output
      logical :: kfl_SaveSnap       !Save data
      logical :: kfl_basis_number   !by vector number instead of energy

      logical :: kfl_snapshotSpecific = .false. !Should I call snapshots specific


      !For Adaptive Mesh refiner
      class(AdaptiveRefinerInterface), pointer :: Refiner => NULL()
      character(6) :: RefinerErrorEstimator 
      character(9) :: RefinerErrorCriteria
      real(rp)     :: RefinerErrorLimits(2)

      !Interpolator
      type(Interpolator), pointer :: Int_Restart => NULL()
      integer(ip)  :: kfl_inter=0             ! Restart with interpolation
      logical      :: kfl_updateBasis=.false. ! Interpolate from OldMesh Basis during run
      integer(ip)  :: UpdateFreq              ! Frequency of interpolation
  
       
   
contains
      procedure :: SetRestartFolder
      procedure :: SetOldRestartFile
      procedure :: SetOldInputFile
      procedure :: EigenSystemMemall  => rom_EigenSystemMemall
      procedure :: GetEigenSystem     => rom_GetEigenSystem
      procedure :: SetPhysicalProblem => rom_SetPhysicalProblem
      procedure :: SetIOPointers      => rom_SetIOPointers 
      procedure :: SetEigenSystem     => rom_SetEigenSystem
      
      procedure :: SnapshotsToSystem         => rom_SnapshotsToSystem
      procedure :: SolutionToPhysicalProblem => rom_SolutionToPhysicalProblem

      procedure :: SetPointers => rom_SetPointers

      procedure :: Turnon    => rom_Turnon
      procedure :: Endste    => rom_Endste
      procedure :: Turnof    => rom_Turnof
      procedure :: GetItask  => rom_GetItask
      procedure :: Setste    => rom_Setste
      procedure :: Begste    => rom_Begste
      procedure :: getste    => rom_getste
      procedure :: Doiter    => rom_Doiter
      procedure :: Outerr    => rom_Outerr
      procedure :: Openfi    => rom_Openfi
      procedure :: Closefi   => rom_closefi
      
      procedure :: Memall      => rom_Memall
      procedure :: DeallocAll  => rom_Dealloc
      
      procedure :: BuildSystem => rom_BuildSystem
      
      procedure :: Readat      => rom_Readat
      procedure :: ReadatMPI   => rom_ReadatMPI
    
      procedure :: SVD         => rom_SVD
      procedure :: BasisFilter => rom_basisFilter
      procedure :: Interpolate => rom_Interpolate

      procedure :: OutputTimes => rom_OutputTimes
      procedure :: GetTimes

      procedure :: InputData   => rom_InputData
      procedure :: Input       => rom_Input
      procedure :: Output      => rom_Output
      procedure :: InputBasis  => rom_InputBasis
      procedure :: Postpr      => rom_postpr

      procedure :: SetAdaptiveRefiner
      procedure :: SetUpdateFreq
      procedure :: SetInterpolation
      procedure :: SetInterp
      procedure :: SetOldMesh
      procedure :: Refine                => rom_Refine
      procedure :: GetRefinementCriteria => rom_GetRefinementCriteria
      
      procedure :: ComputeMassMatrix       => rom_ComputeMassMatrix
      procedure :: ComputeLumpedMassMatrix => rom_ComputeLumpedMassMatrix
      
   end type
   
  
    interface 

       subroutine rom_SetPointers(a)
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
       end subroutine

       subroutine rom_Turnon(a)
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
       end subroutine

       subroutine rom_Endste(a,kfl_gotim)
          use typre
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
          integer(ip) :: kfl_gotim
       end subroutine

       subroutine rom_Turnof(a)
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
       end subroutine

       subroutine rom_Setste(a,dtinv,ctime)
          use typre
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
          real(rp) :: dtinv, ctime
       end subroutine   

       subroutine rom_Begste(a)
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
       end subroutine

       subroutine rom_getste(a,dtinv)
          use typre
          import PodRomProblem
          implicit none
          class(PodRomProblem) :: a
          real(rp) :: dtinv
       end subroutine

       subroutine rom_Doiter(a)
          import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Openfi(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_closefi(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Outerr(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Memall(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Dealloc(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_SVD(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_basisFilter(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Interpolate(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_postpr(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Output(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_InputData(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_Input(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_InputBasis(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_BuildSystem(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem), target :: a
      end subroutine
      
      subroutine rom_SnapshotsToSystem(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem), target :: a
      end subroutine
      
      subroutine rom_SolutionToPhysicalProblem(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem), target :: a
      end subroutine
      
      subroutine rom_Readat(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine
      
      subroutine rom_ReadatMPI(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine

      subroutine rom_EigenSystemMemall(a)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine

      subroutine rom_OutputTimes(a,task)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
         character(6) :: task
      end subroutine

      subroutine rom_Refine(a,itask)
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine rom_GetRefinementCriteria(a,markel)
         use typre
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine rom_ComputeMassMatrix(a)
         use typre
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine

      subroutine rom_ComputeLumpedMassMatrix(a)
         use typre
         import PodRomProblem
         implicit none
         class(PodRomProblem) :: a
      end subroutine

   end interface

  interface PodRomProblem_Const
      procedure constructor
  end interface PodRomProblem_Const

  contains

      function constructor()
          class(PodRomProblem), pointer :: constructor

          allocate(constructor)

      end function constructor
   !-----------------------------------------------------------------------------------
   !ROM subroutines
   subroutine rom_SetPhysicalProblem(a,Problem,ProblemType,ProblemAlgor)
      class(PodRomProblem) :: a
      class(PhysicalProblem), target :: Problem
      character(5)                   :: ProblemType
      integer(ip) :: ProblemAlgor,tn,ndime

      !Pointer to the PhysicalProblem
      a%Problem => Problem
      a%ProblemType = ProblemType
      
      call a%Problem%GetName(a%exmod)
      !Pointer to its Mesh
      call a%Problem%GetMesh(a%Mesh)

      call a%Mesh%GetNdime(ndime)
      tn=(ndime*(ndime+1))/2
      
      call a%Problem%GetNdofn(a%ndofn)
      if (a%ProblemType == 'BOUSS') then
          a%ndofn = a%ndofn + 1
      endif
      if (a%ProblemType == 'IRRSP') then 
          a%ndofn = a%ndofn + tn + 1
          a%kfl_snapshotSpecific = .true.
      endif
      call rom_SetVariableName(a,ProblemAlgor)

   end subroutine

   subroutine rom_SetVariableName(a,ProblemAlgor)
      class(PodRomProblem) :: a
      integer(ip)             :: idime,ndime,tn,iten,ProblemAlgor
      character (len=150),dimension(3):: vel,disp
      character (len=150),dimension(6):: stress
     
      call a%Mesh%GetNdime(ndime)
      tn=(ndime*(ndime+1))/2

      vel(1)    = 'Velocity X'
      vel(2)    = 'Velocity Y'
      vel(3)    = 'Velocity Z'
      disp(1)   = 'Displacement X'
      disp(2)   = 'Displacement Y'
      disp(3)   = 'Displacement Z'
      if(ndime == 2) then
          stress(1) = 'Stress XX'
          stress(2) = 'Stress YY'
          stress(3) = 'Stress XY'
      else
          stress(1) = 'Stress XX'
          stress(2) = 'Stress YY'
          stress(3) = 'Stress ZZ'
          stress(4) = 'Stress YZ'
          stress(5) = 'Stress XZ'
          stress(6) = 'Stress XY'
      endif

      if (a%ProblemType == 'BOUSS') then
         do idime=1,ndime
            a%nameBasis(idime) = vel(idime)
         end do
         a%nameBasis(ndime+1) = 'Pressure'
         a%nameBasis(ndime+2) = 'Temperature'
      elseif (a%ProblemType == 'LMACH') then
         do idime=1,ndime
            a%nameBasis(idime) = vel(idime)
         end do
         a%nameBasis(ndime+1) = 'Temperature'
         a%nameBasis(ndime+2) = 'Pressure'
      elseif (a%ProblemType == 'NSTIN') then
          if(ProblemAlgor==1) then                          !Standard NS
              do idime=1,ndime
                  a%nameBasis(idime) = vel(idime)
              end do
              a%nameBasis(ndime+1) = 'Pressure'
          elseif (ProblemAlgor==3) then                     !SUP NS
              do iten  = 1,tn
                  a%nameBasis(iten)  = stress(iten)
              end do
              do idime = 1+tn,tn+ndime
                  a%nameBasis(idime) = vel(idime-tn)
              end do

              a%nameBasis(ndime+tn+1) = 'Pressure'
          endif
      elseif (a%ProblemType == 'NSCOM') then
         a%nameBasis(1) = 'Pressure'
         do idime=1,ndime
            a%nameBasis(idime+1) = vel(idime)
         end do
         a%nameBasis(ndime+2) = 'Temperature'
      elseif (a%ProblemType == 'TEMPE') then
         a%nameBasis(1) = 'Temperature'
     elseif (a%ProblemType == 'SOLID') then
         if(ProblemAlgor==1) then                           !Standard Solid
             do idime=1,ndime
                 a%nameBasis(idime) = disp(idime)
             end do
         elseif (ProblemAlgor==2 .or. ProblemAlgor==3) then !SUP linear/NeoHookean Solid
             do idime = 1,ndime
                 a%nameBasis(idime) = disp(idime)
             end do
             do iten  = ndime+1,tn+ndime
                 a%nameBasis(iten) = stress(iten-ndime)
             end do

             a%nameBasis(ndime+tn+1) = 'Pressure'
         elseif (ProblemAlgor==4) then                      !UP NeoHookean Solid
             do idime = 1,ndime
                 a%nameBasis(idime) = disp(idime)
             end do
             a%nameBasis(ndime+1) = 'Pressure'
         endif

     elseif (a%ProblemType == 'IRRSP') then !When the irreducible solid pretends to be SUP
             do idime = 1,ndime
                 a%nameBasis(idime) = disp(idime)
             end do
             do iten  = ndime+1,tn+ndime
                 a%nameBasis(iten) = stress(iten-ndime)
             end do

             a%nameBasis(ndime+tn+1) = 'Pressure'
     end if

   end subroutine

   subroutine rom_SetIOPointers(a,Writerpr,Readerpr)
      use typre
      implicit none
      
      class(PodRomProblem) :: a
      class(Writer),target   :: Writerpr
      class(Reader),target   :: Readerpr
      
      a%Readerpr => Readerpr
      a%Writerpr => Writerpr
   end subroutine

   subroutine rom_GetEigenSystem(a,EigenSystem)
      use typre
      implicit none
      class(PodRomProblem) :: a
      class(EigenSystemInterface), pointer :: EigenSystem
      
      EigenSystem => a%EigenSystem 
   end subroutine
   

   subroutine rom_SetEigenSystem(a)
      class(PodRomProblem) :: a
      
      !Pointer to the PhysicalProblem eigensystem
      call a%Problem%SetEigenSystem(a%EigenSystem)
      
   end subroutine

   subroutine rom_GetItask(a,itask)
      use typre
      implicit none
      class(PodRomProblem) :: a
      integer(ip) :: itask

      if (a%kfl_itask == 'FOM' .or. a%kfl_itask == 'SNP') then
         itask = 0
      else
         itask = 1
      end if
   end subroutine
   
   subroutine rom_NULLSUB(a)
      use typre
      implicit none
      class (PodRomProblem) :: a
      
   end subroutine
   
   subroutine SetRestartFolder(a,RestartFolder) 
      use Mod_Int2str
      implicit none
      class(PodRomProblem) :: a
      character(150) :: RestartFolder
      if (a%kfl_ReadType == -1) then
         call runend('PhysicalProblem: setting restart folder before setting read type')
      endif
      
      if (a%kfl_ReadType == 0) then
         a%RestartFolder = RestartFolder
      elseif (a%kfl_ReadType == 1) then
         a%RestartFolder = trim(RestartFolder)//'/rst'//int2str(a%MPIrank)
      else 
      endif
   end subroutine
   
   subroutine SetOldRestartFile(a,OldRestartF)
      use Mod_Int2str
      implicit none
      class(PodRomProblem) :: a
      character(150) :: OldRestartF
      
      a%OldRestartFolder = OldRestartF
   end subroutine    

   subroutine SetOldInputFile(a,namda)
      implicit none
      class(PodRomProblem) :: a
      character(150) :: namda
      
      a%oldnamda = namda
   end subroutine

   subroutine SetAdaptiveRefiner(a,Refiner)
      use typre
      implicit none
      class(PodRomProblem) :: a
      class(AdaptiveRefinerInterface), target :: Refiner
      
      a%Refiner => Refiner
   end subroutine

   subroutine SetInterpolation(a)
      use typre
      implicit none
      class(PodRomProblem) :: a

      a%kfl_inter = 1
   end subroutine

   subroutine SetUpdateFreq(a,UpdateFreq)
      use typre
      implicit none
      class(PodRomProblem) :: a
      integer(ip) :: UpdateFreq

      a%kfl_updateBasis = .true.
      a%UpdateFreq = UpdateFreq
   end subroutine

   subroutine SetInterp(a,Interp)
      use typre
      implicit none
      class(PodRomProblem) :: a
      type(Interpolator), target     :: Interp
      
      a%Int_Restart => Interp
   end subroutine
   
   subroutine SetOldMesh(a)
      use typre
      implicit none
      class(PodRomProblem) :: a
      
      a%OldMesh => a%Problem%OldMesh
   end subroutine
   
   subroutine GetTimes(a,cpu_modul)
      use typre
      implicit none
      class(PodRomProblem) :: a
      real(rp) :: cpu_modul(4)
      
      call a%Timer%Input%GetValue(cpu_modul(1))
      call a%Timer%Assembly%GetValue(cpu_modul(2))
      call a%Timer%Solve%GetValue(cpu_modul(3))
      call a%Timer%Output%GetValue(cpu_modul(4))
   end subroutine

end module

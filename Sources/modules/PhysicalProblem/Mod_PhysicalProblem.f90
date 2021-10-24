module Mod_PhysicalProblem
   use typre
   use Mod_Timer
   use Mod_MpiObject
   use Mod_ParallelSystemInterface
   use Mod_ParallelLibraryInterface
   use Mod_EigenSystemInterface
   use Mod_postpr
   use Mod_ReadWrite
   use Mod_Mesh
   use Mod_AdaptiveInterface
   use Mod_MeshInterpolator
   use Mod_CutMesh
   implicit none
   private
   public PhysicalProblem, php_TimesToCpuModulArray

   type PhysicalTime
      type(Timer) :: Total
      type(Timer) :: Turnon
      type(Timer) :: Getste
      type(Timer) :: Begste
      type(Timer) :: Doiter
      type(Timer) :: BuildLinearSystem
      type(Timer) :: SolveLinearSystem
      type(Timer) :: Endite
      type(Timer) :: Endste
      type(Timer) :: Output
      type(Timer) :: Turnof
      type(Timer) :: Restar
      type(Timer) :: Refine
      type(Timer) :: ProjectArraysOUM
   end type

   type, abstract, extends(MpiObject) ::  PhysicalProblem

   type(Interpolator), pointer :: Int_Restart => NULL()
   type(Interpolator)          :: TrackingInterpolator

   !Miscellaneous
   class(FemMesh), pointer           :: Mesh => NULL(), OldMesh => NULL()
   class(PostprFile),pointer         :: FilePostpr => NULL()
   class(Reader),pointer             :: Readerpr => NULL()
   class(Writer),pointer             :: Writerpr => NULL()
   class(ParallelLibraryInterface), pointer   :: ParallelLibrary => NULL()

   type(PhysicalTime)                :: Timer
   character(66)                     :: exmod, namod
   character(150)                    :: oldnamda
   character(150)                    :: RestartFolder
   character(150)                    :: OldRestartFolder
   character(5)                      :: extramatrix

   !Linear System
   class(ParallelSystemInterface), pointer   :: LinearSystem => NULL()

   class(EigenSystemInterface), pointer      :: EigenSystem => NULL()
   real(rp), allocatable                     :: unkno(:,:)
   !For using in the specific routines
   integer(ip) :: &
            gipoin,&
            giboun,gbouni,&
            gipsta,gaux_nboun,gmnodb,gielem
   ! Logical units
   integer(ip) ::&
            lun_outph, &             ! Output for the physical problem
            lun_conve, &             ! Convergence
            lun_cpconve, &           ! Coupled Convergence
            lun_solve, &             ! Solver
            lun_rstar, &             ! Restart file
            lun_rsta2, &             ! Restart file
            lun_nolin, &             ! Non-linearity
            lun_adapt, &             ! Adaptivity
            lun_error, &             ! Error exact
            lun_ersgs                ! Error SGS

   ! General
   integer(ip) ::    &
            istep,      &            ! Number of time steps run
            ndofn,      &            ! Number of d.o.f
            ndofr,      &            ! Number of reduced d.o.f
            ndofbc,     &            ! Number of degrees of freedom for reading boundary conditions (standard)
            ndofbcstart=0, &         ! First component of the elemental matrix were the BCs are applied 
            nsmax

   real(rp)    :: &
            ctime = 0.0_rp, &        ! Current time
            bctime = 0.0_rp,&        ! Time for boundary conditions (important for Crank-Nicolson)
            timef = 0.0_rp

   !Physical Problem - PHYSICAL_PROBLEM ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !Physical Problem Definition - PROBLEM_DEFINTION -------------------------------------------
      integer(ip)  :: kfl_timei            ! Existence of du/dt or any temporal derivative
      integer(ip)  :: kfl_SwitchOff        ! Disconnection of elements
   !Exact solution
      integer(ip)  :: kfl_exacs = 0        ! Exact solution
      integer(ip)  :: exact                !Exact solution
      real(rp)     :: expar(10)            !Exact solution parameters
   !Physical Problem Properties - PROPERTIES --------------------------------------------------
      !nothing here, defined in the specific physical problem
   !Boundary conditions and initial conditions - BOUNDARY_CONDITIONS +++++++++++++++++++++++++++++

   !Boundary Conditions Read Strategy
   integer(ip) :: kfl_BoundaryConditionsReadStrategy = 0

   integer(ip), allocatable :: kfl_fixno(:,:) !Nodal fixity (1/ndime/ndofn,npoin)
   real(rp), allocatable :: bvess(:,:,:)      !Essential bc values (1/ndime/ndofn,npoin,)

   integer(ip), allocatable :: kfl_fixbo(:)     !Element boundary fixity (nboun)
   type(r1p)  , allocatable :: bvnat(:)         !Natural bc values (nboun)
   integer(ip) :: bvnat_coun  !Counts the dimension of the bvnat array

   integer(ip),  allocatable :: &
            kfl_funno(:),&              !Functions for node bc (npoin)
            kfl_funbo(:),&              !Functions for onboundaries bc (nboun)
            kfl_funty(:,:)              !Function type and number of parameters (ifunc,1:2)
   type(r1p), allocatable :: funpa(:)   !Function parameters (ifunc,1:a%kfl_funty(ifunc,2))

   integer(ip) ::&
            kfl_incnd=0,&               !Initial condition type
            kfl_conbc                   !Constant b.c.

   real(rp) :: &
            picnd(10),&                 !Initial condition parameters
            updbcn                      !Parameter to update the natural bc

   integer(ip) :: ManufacturedBoundaryCondition

!Numerical treatment - NUMERICAL_TREATMENT +++++++++++++++++++++++++++++++++++++++++++++++++++
   !outer loop iterative coupling
   integer(ip)  :: kfl_normc          ! Norm of convergence outer iterations (iterative coupling)
   integer(ip)  :: kfl_initRelax=2    ! Norm of convergence outer iterations (iterative coupling)
   real(rp)     :: resid = 0.0_rp     ! Residual for outer iterations (iterative coupling)
   real(rp)     :: cpres = 0.0_rp     ! Residual for coupling iterations (iterative coupling)
   !time marching loop
   integer(ip)  :: kfl_stead  = 0     ! Steady-state has been reached
   real(rp)     :: &
            dtinv ,&                  ! 1/dt
            dtinv2,&                  ! 1/dt^2
            dtcri ,&                  ! Critical time step
            dtime                     ! Time step
   character(5) :: &
         kfl_tsche_1st_datafile,&     ! Wanted Temporal Scheme: 1st time derivative
         kfl_tsche_1st_current ,&     ! Actual Temporal Scheme, 1st time derivative
         kfl_tsche_2nd_datafile,&     ! Wanted Temporal Scheme: 2nd time derivative
         kfl_tsche_2nd_current        ! Actual Temporal Scheme, 2nd time derivative
   integer(ip)  :: neule                     ! Number of Euler time steps
   integer(ip)  :: neule_2nd                 ! Number of Euler time steps 2nd order
   integer(ip)  :: kfl_tsche_change   ! detects changes in time scheme
   integer(ip)  :: kfl_elmat_datafile ! compute/re-compute element matrix (computed at first time step and on time scheme change)
   integer(ip)  :: kfl_elmat_current
   integer(ip)  :: kfl_tiacc          ! Temporal accuracy OBSOLETE
   integer(ip)  :: kfl_twost          ! Two steps BDF2 scheme for time derivative OBSOLETE
   !inner loop non-linearity
   integer(ip)  :: maxit              ! Max # of iterations of non-linearity loop
   integer(ip)  :: kfl_goite          ! Keep iterating non-linearity loop
   integer(ip)  :: itera  = 0         ! Internal iteration counter (non-linearity loop)
   integer(ip)  :: cpiter = 0         ! Case coupling iteration counter (coupling loop)
   real(rp)     :: cotol              ! Convergence tolerance (non-linearity loop)
   real(rp)     :: cptol              ! Convergence tolerance (coupled case loop)
   integer(ip)  :: kfl_linea          ! Linearization (RHS=0, Picard=1, Newton=2)
   integer(ip)  :: npica              ! Number of Picard iteration (Newton's lin.)
   !outer Loop non-linearity
   integer(ip) :: OutIiter

   logical :: kfl_SkipSystemSolve = .false.
   logical :: kfl_error = .false.  ! print error SGS

   integer(ip) :: ncomp,oldncomp                     ! Number of components, unknown vectors

   integer(ip) :: kfl_ProjectionType        ! Projection Type: 0: lumped 1: L2

   character(15) :: &
         EndLoopQuadrature       ! Quadrature rule for endste and endite loops: 'DefaultRule'. 'ForceClosedRule'

   real(rp) ::&
            sstol   ,&                  ! Steady state tolerance (time marching loop)
            safet   ,&                  ! Safety factor for time step
            err01(2),&                  ! L1 error u
            err02(2),&                  ! L2 error u
            err0i(2),&                  ! Linf error u
            err11(2),&                  ! L1 error grad(u)
            err12(2),&                  ! L2 error grad(u)
            err1i(2)                    ! Linf error grad(u)

!Output & postprocess - OUTPUT_&_POST-PROCESS +++++++++++++++++++++++++++++++++++++++++++++++
   !General Post-processing
   integer(ip)  :: npp_inits            ! Postprocess initial step
   real(rp)     :: pos_tinit            ! Postprocess initial time
   integer(ip)  :: npp_stepi(49)        ! Postprocess step intervals
   integer(ip)  :: pos_alrea(49)
   real(rp)     :: pos_times(10,49)     ! Postprocess times for u,p, etc.
   !Tracking of points
   integer(ip)  :: nptra = 0            ! Number of points to be tracked
   real(rp), allocatable :: cptra(:,:)  ! Coordinates of points to be tracked


   !For Adaptive Mesh refiner
   class(AdaptiveRefinerInterface), pointer :: Refiner => NULL()
   character(6) :: RefinerErrorEstimator
   character(9) :: RefinerErrorCriteria
   real(rp)     :: RefinerErrorLimits(2)
   integer(ip)  :: kfl_SkipLinearSystemRefinement = 0
   logical      :: kfl_adap = .false.

   integer(ip) :: kfl_MPIComType = 0

   !For FixedMeshes
   !Cut Mesh type
   type(CutMesh), pointer      :: CutMesh => NULL()

   !Eliminate or not the columns when applying the boundary conditions
   logical :: kfl_DeleteDirichletColumns = .true.

   !Restrictions
   integer(ip) :: kfl_restrictions = 0
   integer(ip) :: nrest
   real(rp), pointer :: RHSin(:) => NULL()         ! Restriction values coming from another mesh
   real(rp)          :: RHSout(5)                  ! Restriction being sent to another mesh
   character (len = 5), dimension(5) :: keyrest

   !Restart
   character(150) :: fil_rstar
   integer(ip) :: kfl_preli,&                     ! Preliminary for this module
                  kfl_rstar=0,&                   ! Restart within the running module
                  kfl_inter=0                     ! Restart with interpolation
   real(rp), allocatable :: interp_veloc(:,:,:)   ! ALE velocity interpolated to the new mesh
   real(rp), allocatable :: restar_veloc(:,:,:),& ! ALE velocity coming from the old mesh
                            restar_disp(:,:,:)    ! ALE displacement coming from the old mesh

   logical :: kfl_coupconv = .false.    ! Did coupled case converge?
   logical :: kfl_docoupconv = .false.  ! Should I check coupled case conv?
   logical :: kfl_printBC= .true.       ! Should I print BC at the beginning of the run?
   
   !Multiple communicators
   integer(ip) :: kfl_multicomm = 1
   integer(ip) :: MulticommColor = 0

   !For ROM-SVD
   logical :: kfl_eigensolve
   logical :: kfl_isSlave = .false. ! Is this physical problem slave to another one?


contains

      procedure :: SetMesh
      procedure :: SetPostprFile
      procedure :: SetIOPointers
      procedure :: SetMPICommunicationsType
      procedure :: SetGlobalIteration
      procedure :: SetCouplingIteration
      procedure :: SetMulticommData

      !Restart
      procedure :: SetRestart
      procedure :: SetInterpolation
      procedure :: SetRestartFolder
      procedure :: SetOldInputFile
      procedure :: SetOldRestartFile
      procedure :: SetOldMesh
      procedure :: SetInterp

      procedure :: SetParallelLibrary
      procedure :: OpenFiles    => php_openfi
      procedure :: CloseDataFile=> php_closedatafi
      procedure :: CloseFiles   => php_closefi
      procedure :: Turnon       => php_turnon
      procedure :: Reaous       => php_reaous
      procedure :: Reaphy       => php_reaphy
      procedure :: Reanut       => php_reanut
      procedure :: ReaMPI       => php_reaMPI
      procedure :: Reabcs       => php_reabcs
      procedure :: Bounod       => php_bounod
      procedure :: Setste       => php_setste
      procedure :: SetCtime     => php_SetCtime
      procedure :: Doiter       => php_doiter
      procedure :: LinearSystemMemall => php_LinearSystemMemall
      procedure :: LinearSystemTurnof => php_LinearSystemTurnof
      procedure :: BuildLinearSystem  => php_BuildLinearSystem
      procedure :: Solite       => php_solite
      procedure :: Turnof       => php_turnof
      procedure :: Begite       => php_begite
      procedure :: Begste       => php_begste
      procedure :: Updbcs       => php_updbcs
      procedure :: Iniunk       => php_iniunk
      procedure :: Endste       => php_Endste
      procedure :: lev1Output   => php_Output
      procedure :: Restart      => php_restar
      procedure :: OutputTimes  => php_OutputTimes
      procedure :: Endite       => php_Endite
      procedure :: Outerr       => php_Outerr
      procedure :: Project      => php_Project
      procedure :: GetMesh      => php_GetMesh
      procedure :: GetName      => php_GetName
      procedure :: GetNdofn     => php_GetNdofn
      procedure :: GetNdofr     => php_GetNdofr
      procedure :: GetIstep     => php_GetIstep
      procedure :: GetItera     => php_GetItera

      procedure :: GetUnkno   => php_GetUnkno
      procedure :: SetUnknoA1 => php_SetUnknoA1
      procedure :: SetUnknoA2 => php_SetUnknoA2
      generic   :: SetUnkno   => SetUnknoA1,SetUnknoA2


      procedure :: SpecificGetUnkno     => php_SpecificGetUnkno
      procedure :: SetInitialConditions => php_SetInitialConditions
      procedure :: GetInitialConditions => php_GetInitialConditions

      procedure :: Initializations => php_Initializations
      procedure :: SetLinearSystem => php_SetLinearSystem
      procedure :: GetLinearSystem => php_GetLinearSystem
      procedure :: SetEigenSystem  => php_SetEigenSystem
      procedure :: GetTimeScheme1stOrder_Datafile  => php_GetTimeScheme1stOrder_Datafile
      procedure :: SetCutMesh => php_SetCutMesh
      procedure :: GetCutMesh => php_GetCutMesh
      procedure :: ProjectArraysOntoUndeformedMesh => php_ProjectArraysOntoUndeformedMesh
      procedure :: AdvectArraysOntoUndeformedMesh  => php_AdvectArraysOntoUndeformedMesh
      procedure :: SpecificProjectArraysOUM
      procedure :: SpecificAdvectArraysOUM
      procedure :: Project_disc

      procedure :: OpenFilesRestart
      procedure :: CloseFilesRestart
      procedure :: GeneralRestart
      procedure :: GaussArrayRestart

      procedure :: MakeSteady => php_MakeSteady
      procedure :: SkipSystemSolve => php_SkipSystemSolve

      !boundary conditions
      procedure :: OTAReadOnNodesDo => php_OTAReadOnNodesDo
      procedure :: OTAReadOnBoundariesDo => php_OTAReadOnBoundariesDo
      procedure :: OTAReadOnElementsDo => php_OTAReadOnElementsDo
      procedure :: OTAReadSourceDo => php_OTAReadSourceDo


      procedure :: ReadOnNodes => php_ReadOnNodes
      procedure :: ReadOnBoundaries => php_ReadOnBoundaries
      procedure :: ReadOnElements => php_ReadOnElements
      procedure :: ReadSource => php_ReadSource

      procedure :: SpecificReabcs => php_SpecificReabcsRUNEND
      procedure :: SpecificReadOnNodes => php_RUNENDSUB
      procedure :: SpecificReadSource  => php_RUNENDSUB
      procedure :: SpecificReadOnBoundaries => php_RUNENDSUB
      procedure :: SpecificReadOnElements => php_RUNENDSUB
      procedure :: SpecificReadOnFunctions => php_NULLSUB_ifunc
      procedure :: SpecificScatterOnFunctions => php_NULLSUB
      procedure :: SpecificManufacturedBoundaryConditions => php_RUNENDSUB

      !Parameters
      procedure :: SetParameters => php_SetParameters
      procedure :: GetParameters => php_GetParameters

      procedure :: GetTime
      procedure :: WriteTimes
      procedure :: GetNamod
      procedure :: GetResidual
      procedure :: GetCouplingResidual
      procedure :: GetCoupledConvFlag
      procedure :: SetCoupledConvFlag
      procedure :: SetAdaptiveFlag

      !Adaptive Mesh Refinement
      procedure :: SetAdaptiveRefiner
      procedure :: Refine => php_Refine
      procedure :: PreRefine => php_NULLSUB
      procedure :: SpecificRefine => php_RUNENDSUBitaskchar
      procedure :: GetRefinementCriteria => php_RUNEND_RefCri
      procedure :: SpecificSubscalesRefCriteria  => php_RUNEND_SubRef

      !Restrictions
      procedure :: GetRHSarray
      procedure :: SetRHSarray

      !Restrictions
      procedure :: relaxVector=>php_relaxVector

      !Deferred Procedures
      procedure(SetExmod),       deferred       :: SetExmod
      procedure(SetNdofn),       deferred       :: SetNdofn
      procedure(SetNdofbc),      deferred       :: SetNdofbc
      procedure(SpecificReaous), deferred       :: SpecificReaous
      procedure(SpecificReaphy), deferred       :: SpecificReaphy
      procedure(SpecificReanut), deferred       :: SpecificReanut
      procedure(SpecificReaMPI), deferred       :: SpecificReaMPI

      procedure(SpecificTurnon), deferred       :: SpecificTurnon
      procedure(SpecificBounod), deferred       :: SpecificBounod
      procedure(SpecificSolite), deferred       :: SpecificSolite
      procedure(SpecificTurnof), deferred       :: SpecificTurnof
      procedure(SpecificBegite), deferred       :: SpecificBegite
      procedure(SpecificBegste), deferred       :: SpecificBegste
      procedure(SpecificUpdbcs), deferred       :: SpecificUpdbcs
      procedure(SpecificIniunk), deferred       :: SpecificIniunk
      procedure(SpecificExaerr), deferred       :: SpecificExaerr
      procedure(SpecificEndste), deferred       :: SpecificEndste
      procedure(SpecificEndite), deferred       :: SpecificEndite
      procedure(SpecificCrankNicolsonEndste), deferred       :: SpecificCrankNicolsonEndste

      procedure(SpecificOuterr), deferred :: SpecificOuterr
      procedure(Memall), deferred :: Memall
      procedure(Getste), deferred :: Getste
      procedure(CvgUnk), deferred :: Cvgunk
      procedure(Output), deferred :: Output
      procedure(Elmope), deferred :: Elmope
      procedure(Bouope), deferred :: Bouope
      procedure(SpecificRestart), deferred :: SpecificRestart

      procedure :: PoinOpe => php_NULLSUB

   end type

   interface

      subroutine php_openfi(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_closefi(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_closedatafi(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_turnon(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_reaous(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_reaphy(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_reanut(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_reaMPI(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_reabcs(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_bounod(a)
          use typre
          import PhysicalProblem
          implicit none
          class(PhysicalProblem) :: a
       end subroutine

      subroutine php_doiter(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Operations to be done at the beginning of a time step
      end subroutine

      subroutine php_LinearSystemMemall(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Allocation of unknown and LinearSystem
      end subroutine

      subroutine php_LinearSystemTurnof(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Allocation of unknown and LinearSystem
      end subroutine

      subroutine php_solite(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Allocation of unknown and LinearSystem
      end subroutine

      subroutine php_BuildLinearSystem(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Construct the LinearSystem
      end subroutine

      subroutine php_turnof(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a


      end subroutine

      subroutine php_begite(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Operations to be done at the beginning of a iteration
      end subroutine

      subroutine php_begste(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Operations to be done at the beginning of a time step
      end subroutine

      subroutine php_updbcs(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Update of the boundary conditions if non-constant
      end subroutine

      subroutine php_iniunk(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Initialization of the unknown vectors
      end subroutine

      subroutine php_Endste(a,kfl_gotim)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(out) :: kfl_gotim

         !Operations to be done at the end of a time step
      end subroutine

      subroutine php_Output(a)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_restar(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem)  :: a
         integer(ip), intent(in) :: itask

         !Operations to be done at the end of a time step
      end subroutine

      subroutine php_Outerr(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

         !Output of errors
      end subroutine


      subroutine php_OutputTimes(a)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Output of times at the end of the step
      end subroutine

      subroutine php_Endite(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

         !Operations to be done at the end of a iteration
      end subroutine

      subroutine php_Initializations(a)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

      end subroutine

      subroutine php_Refine(a,itask)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         character(6) :: itask

         !Adaptive Mesh Refinement
      end subroutine

      subroutine php_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
         use typre
         use Mod_MeshInterpolator
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         type(Interpolator) :: Interp
         integer(ip) :: itask
      end subroutine
      
      subroutine php_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
         use typre
         use Mod_Advector
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         type(Advector) :: Advect
         integer(ip) :: itask
      end subroutine

      subroutine OpenFilesRestart(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine

      subroutine CloseFilesRestart(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine

      subroutine GeneralRestart(a,array,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in)    :: itask
         real(rp) :: array(:,:)
      end subroutine

      subroutine GaussArrayRestart(a,array,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in)    :: itask
         type(r3p) :: array(:)
      end subroutine

      subroutine php_ReadOnNodes(a)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

       subroutine php_ReadOnBoundaries(a,knodb)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: knodb(:)

      end subroutine

       subroutine php_ReadSource(a)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

       subroutine php_ReadOnElements(a)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_OTAReadOnNodesDo(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_OTAReadOnBoundariesDo(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_OTAReadOnElementsDo(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_OTAReadSourceDo(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      end subroutine

      subroutine php_relaxVector(a,itask,w_r,w_rmax,r_i1,r_i21)
         import
         implicit none

         class(PhysicalProblem) :: a
         character(5) :: itask
         real(rp) :: w_r,w_rmax
         real(rp) :: r_i1(:,:),r_i21(:,:)

      end subroutine

   end interface

   abstract interface

      subroutine SetExmod(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

         !This subroutine sets the Physical Problem extension (.nsi)
         !a%exmod = 'nsi'

      end subroutine

      subroutine SetNdofn(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

         !This subroutine sets the Physical Problem Number of degrees of freedom
         !a%ndofn = ndime+1

      end subroutine

      subroutine SetNdofbc(a)
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

         !This subroutine sets the Physical Problem Number of degrees of freedom for setting boundary conditions
         !a%ndofn = ndime

         !It needs not coincide with ndofn (NavierStokes: ndofn = ndime+1, ndofbc = ndime)

      end subroutine

      subroutine SpecificReaous(a,itask)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificReaphy(a,itask)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificReanut(a,itask)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificReaMPI(a)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a

      end subroutine

      subroutine SpecificReabcs(a,itask)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificReabcsMPI(a,itask)
         use typre
         import PhysicalProblem
         implicit none

         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificBounod(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificSolite(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine SpecificTurnon(a)
         use typre
         use Mod_Element
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Turnon operations for the specific physical problem
      end subroutine

      subroutine SpecificTurnof(a)
         use typre
         use Mod_Element
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Turnof operations for the specific physical problem
      end subroutine

      subroutine SpecificBegste(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Operations to be done at the beginning of a time step
      end subroutine

      subroutine SpecificBegite(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Operations to be done at the beginning of a iteration
      end subroutine

       subroutine SpecificUpdbcs(a)
          import PhysicalProblem
          implicit none
          class(PhysicalProblem) :: a

          !Overwrite boundary conditions using exact values
       end subroutine

       subroutine SpecificExaerr(a)
          import PhysicalProblem
          implicit none
          class(PhysicalProblem) :: a

       end subroutine

      subroutine SpecificIniunk(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Initialization of the unknown vectors
      end subroutine

      subroutine SpecificEndste(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

         !Operations to be done at the end of a time step
      end subroutine

      subroutine SpecificEndite(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

         !Operations to be done at the end of a iteration
      end subroutine

      subroutine SpecificCrankNicolsonEndste(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Crank-Nicolson Update for the end of the step
      end subroutine


      subroutine Elmope(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         !Computation of the LinearSystem, volume contribution
      end subroutine

      subroutine Bouope(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Computation of the LinearSystem, boundary contribution
      end subroutine

      subroutine SpecificOuterr(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Error output

      end subroutine

      subroutine Memall(a)
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a

         !Memory allocation
      end subroutine

      subroutine Getste(a,dtinv)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         real(rp) :: dtinv

         !Compute the critical time step, return dtinv
      end subroutine

      subroutine Cvgunk(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in) :: itask

         !Convergence computations
      end subroutine

      subroutine Output(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip) :: itask

         !Output results
      end subroutine

      subroutine SpecificRestart(a,itask)
         use typre
         import PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in) :: itask

         !Restart calculations
      end subroutine

   end interface

contains

   subroutine SetOldMesh(a,Mesh)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(FemMesh),target      :: Mesh

      a%OldMesh => Mesh
   end subroutine

   subroutine SetMesh(a,Mesh)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(FemMesh),target      :: Mesh

      a%Mesh => Mesh
   end subroutine

   subroutine SetPostprFile(a,FilePostpr)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(PostprFile),target   :: FilePostpr

      a%FilePostpr => FilePostpr
   end subroutine

   subroutine SetIOPointers(a,Writerpr,Readerpr)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(Writer),target   :: Writerpr
      class(Reader),target   :: Readerpr

      a%Readerpr => Readerpr
      a%Writerpr => Writerpr
   end subroutine

   subroutine SetMPICommunicationsType(a,kfl_MPIComType)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: kfl_MPIComType

      a%kfl_MPIComType = kfl_MPIComType
   end subroutine

   subroutine SetInterp(a,Interp)
      use typre
      use Mod_MeshInterpolator
      implicit none
      class(PhysicalProblem) :: a
      type(Interpolator), target     :: Interp

      a%Int_Restart => Interp
   end subroutine

   subroutine SetRestart(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a

      a%kfl_rstar = 1
   end subroutine

   subroutine SetInterpolation(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a

      a%kfl_inter = 1
   end subroutine

   subroutine SetParallelLibrary(a,ParallelLibrary)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(ParalleLLibraryInterface), target :: ParallelLibrary

      a%ParallelLibrary => ParallelLibrary
   end subroutine

   subroutine GetNdofbc(a,ndofbc)
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: ndofbc

      ndofbc = a%ndofbc
   end subroutine

   subroutine php_setste(a,dtinv,ctime,timef,nsmax)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: dtinv,ctime,timef
      integer(ip) :: nsmax

      a%dtinv = dtinv
      a%dtinv2= dtinv*dtinv
      a%ctime = ctime
      a%bctime= ctime
      a%dtime = 1.0/dtinv
      a%kfl_coupconv = .false.
      a%timef = timef
      a%nsmax = nsmax

   end subroutine

   subroutine php_SetCtime(a,ctime)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: ctime

      a%ctime = ctime
   end subroutine

   subroutine GetTime(a,cpu_modul)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: cpu_modul

      call a%Timer%Total%GetValue(cpu_modul)

   end subroutine

   subroutine php_TimesToCpuModulArray(a,cpu_modul)
      class(PhysicalProblem) :: a
      real(rp) :: cpu_modul(14)

      call a%Timer%Total%GetValue(cpu_modul(1))
      call a%Timer%Turnon%GetValue(cpu_modul(2))
      call a%Timer%Getste%GetValue(cpu_modul(3))
      call a%Timer%Begste%GetValue(cpu_modul(4))
      call a%Timer%Doiter%GetValue(cpu_modul(5))
      call a%Timer%BuildLinearSystem%GetValue(cpu_modul(6))
      call a%Timer%SolveLinearSystem%GetValue(cpu_modul(7))
      call a%Timer%Endite%GetValue(cpu_modul(8))
      call a%Timer%Endste%GetValue(cpu_modul(9))
      call a%Timer%Output%GetValue(cpu_modul(10))
      call a%Timer%Turnof%GetValue(cpu_modul(11))
      call a%Timer%Restar%GetValue(cpu_modul(12))
      call a%Timer%Refine%GetValue(cpu_modul(13))
      call a%Timer%ProjectArraysOUM%GetValue(cpu_modul(14))

   end subroutine


   subroutine WriteTimes(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: cpu_modul(14)
      real(rp) :: cpuminim,cpudenom

      integer(ip) :: ifield

      call php_TimesToCpuModulArray(a,cpu_modul)

      cpuminim = 1.0d-6
      cpudenom = max(cpu_modul(1),cpuminim)
      write(a%lun_outpu,101) (cpu_modul(ifield),&
            100.0_rp*(cpu_modul(ifield))/cpudenom,ifield = 2,14)

      101 format(&
            10x,'     Turnon                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Getste                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Begste                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Doiter                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'       Build Linear System  :       ',F12.2,' (',F6.2,' % )',/,&
            10x,'       Solve Linear System  :       ',F12.2,' (',F6.2,' % )',/,&
            10x,'       Endite               :       ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Endste                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Output                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Turnof                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Restar                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Refine                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Project between meshes :   ',F12.2,' (',F6.2,' % )')
   end subroutine


   subroutine GetNamod(a,namod)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: namod

      namod(1:6) = a%namod(1:6)

   end subroutine

   subroutine GetResidual(a,residual)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: residual

      residual = a%resid
   end subroutine

   subroutine GetCouplingResidual(a,residual)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: residual

      residual = a%cpres
   end subroutine

   subroutine GetCoupledConvFlag(a,kfl_coupconv)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      logical :: kfl_coupconv

      kfl_coupconv= a%kfl_coupconv
   end subroutine

   subroutine SetCoupledConvFlag(a,kfl_docoupconv)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      logical :: kfl_docoupconv

      a%kfl_docoupconv = kfl_docoupconv
   end subroutine

   subroutine SetAdaptiveFlag(a,kfl_adap)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: kfl_adap

      if (kfl_adap == 1) a%kfl_adap = .true.
   end subroutine

   subroutine php_RUNENDSUB(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a

      call runend('A procedure in PhysicalProblem is pointing to the RUNENDSUB')

   end subroutine

   subroutine php_RUNENDSUBitaskchar(a,itask)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask

      call runend('A procedure in PhysicalProblem is pointing to the RUNENDSUB')

   end subroutine

   subroutine php_RUNENDSUBitask(a,itask)
      use typre
      implicit none
      class(PhysicalProblem)  :: a
      integer(ip), intent(in) :: itask

      call runend('A procedure in PhysicalProblem is pointing to the RUNENDSUB')

   end subroutine

   subroutine php_SpecificReabcsRUNEND(a,itask,kflag)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag

      call runend('A procedure in PhysicalProblem is pointing to the php_SpecificReabcsRUNEND')

   end subroutine

   subroutine php_Project(a,ndofn,array)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: ndofn
      real(rp)    :: array(ndofn,*)

      if (a%kfl_ProjectionType == 0) then
         call a%Mesh%Smooth(ndofn,array)
      elseif (a%kfl_ProjectionType == 1) then
         call a%Mesh%Project(ndofn,array)
      endif
   end subroutine

   subroutine php_GetMesh(a,Mesh)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(FemMesh), pointer :: Mesh

      Mesh => a%Mesh
   end subroutine

   subroutine php_GetName(a,exmod)
      use typre
      implicit none
      class(PhysicalProblem)   :: a
      character(66)  :: exmod

      exmod = a%exmod
   end subroutine

   subroutine php_GetNdofn(a,ndofn)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: ndofn

      ndofn = a%ndofn
   end subroutine

   subroutine php_GetNdofr(a,ndofr)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: ndofr

      ndofr = a%ndofr
   end subroutine

   subroutine php_GetUnkno(a,unkno)
      use typre
      implicit none
      class(PhysicalProblem), target :: a
      real(rp)    :: unkno(:,:)
      integer(ip) :: idofn

      do idofn = 1,a%ndofn
         unkno(idofn,:) = a%unkno(idofn,:)
      end do
   end subroutine

   subroutine php_SetUnknoA1(a,unkno,idofn)
      use typre
      implicit none
      class(PhysicalProblem), target :: a
      real(rp)    :: unkno(:)
      integer(ip) :: idofn

      a%unkno(idofn,:) = unkno
      call a%Mesh%ArrayCommunicator%GhostCommunicate(idofn,a%unkno(idofn,:))
   end subroutine

   subroutine php_SetUnknoA2(a,unkno)
      use typre
      implicit none
      class(PhysicalProblem), target :: a
      real(rp) :: unkno(:,:)

      a%unkno = unkno
      call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%unkno)
   end subroutine

   subroutine php_SpecificGetUnkno(a,unkno)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: unkno(:,:)

      call runend('SpecificGetUnkno not ready for this physical problem')
   end subroutine

   subroutine php_SetInitialConditions(a,unkno)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: unkno(:,:,:)

      call runend('SetInitialConditions not ready for this physical problem')
   end subroutine

   subroutine php_GetInitialConditions(a,unkno)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: unkno(:,:,:)

      call runend('GetInitialConditions not ready for this physical problem')
   end subroutine

   subroutine php_GetIstep(a,istep)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: istep

      istep = a%istep
   end subroutine

   subroutine php_GetItera(a,itera)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: itera

      itera = a%itera
   end subroutine

   subroutine php_GetTimeScheme1stOrder_Datafile(a,kfl_tsche_1st_datafile)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      character(5) :: kfl_tsche_1st_datafile

      kfl_tsche_1st_datafile = a%kfl_tsche_1st_datafile
   end subroutine

   subroutine php_SetLinearSystem(a,LinearSystem)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(ParallelSystemInterface), target :: LinearSystem

      a%LinearSystem => LinearSystem
   end subroutine

   subroutine php_GetLinearSystem(a,LinearSystem)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(ParallelSystemInterface), pointer :: LinearSystem

      LinearSystem => a%LinearSystem
   end subroutine

   subroutine php_SetEigenSystem(a,EigenSystem)
      class(PhysicalProblem) :: a
      class(EigenSystemInterface), target :: EigenSystem

      !Pointer to the PodRomProblem eigensystem
      a%EigenSystem=> EigenSystem

   end subroutine

   !> This subroutine is used to set parameters from the exterior <br>
   !> It needs to be specifically implemented from each physical problem
   subroutine php_SetParameters(a,iparam,param)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: iparam
      real(rp)    :: param

      call runend('SetParameters not defined for this physical problem')
   end subroutine

   !> This subroutine is used to get parameters from the exterior <br>
   !> It needs to be specifically implemented from each physical problem
   subroutine php_GetParameters(a,iparam,param)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: iparam
      real(rp)    :: param

      call runend('GetParameters not defined for this physical problem')
   end subroutine

   subroutine SetAdaptiveRefiner(a,Refiner)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      class(AdaptiveRefinerInterface), target :: Refiner

      a%Refiner => Refiner
   end subroutine

   subroutine SetRestartFolder(a,RestartFolder)
      use Mod_Int2str
      implicit none
      class(PhysicalProblem) :: a
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

   subroutine php_RUNEND_RefCri(a,markel)
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: markel(*)

      call runend('GetRefinementCriteria not defined for current module')
   end subroutine

   subroutine php_RUNEND_SubRef(a,error,TotalEstimatedError)
      implicit none
      class(PhysicalProblem) :: a
      real(rp) :: error(:), TotalEstimatedError

      call runend('SubscalesRefCriteria not defined for current module')
   end subroutine

   subroutine php_SetCutMesh(a,CMesh)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      type(CutMesh), target :: CMesh

      a%CutMesh => CMesh
   end subroutine

   subroutine php_GetCutMesh(a,CMesh)
      use typre
      implicit none
      class(PhysicalProblem), target :: a
      type(CutMesh), pointer :: CMesh

      CMesh => a%CutMesh
   end subroutine

   subroutine SpecificProjectArraysOUM(a,Interp,itask)
      use typre
      use Mod_MeshInterpolator
      implicit none
      class(PhysicalProblem) :: a
      type(Interpolator) :: Interp
      integer(ip) :: itask

      call runend('SpecificProjectArraysOUM not implemented')
   end subroutine
   
   subroutine SpecificAdvectArraysOUM(a,Advect,itask)
      use typre
      use Mod_Advector
      implicit none
      class(PhysicalProblem) :: a
      type(Advector) :: Advect
      integer(ip) :: itask

      call runend('SpecificAdvectArraysOUM not implemented')
   end subroutine

   subroutine GetRHSarray(a,RHS)
      use typre
      implicit none
      class(PhysicalProblem), target :: a
      real(rp), pointer :: RHS(:)

      RHS => a%RHSout
   end subroutine

   subroutine SetRHSarray(a,RHS)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      real(rp), target :: RHS(:)

      a%RHSin => RHS
   end subroutine

   subroutine SetOldInputFile(a,namda)
      implicit none
      class(PhysicalProblem) :: a
      character(150) :: namda

      a%oldnamda = namda
   end subroutine

   subroutine SetOldRestartFile(a,OldRestartF)
      use Mod_Int2str
      implicit none
      class(PhysicalProblem) :: a
      character(150) :: OldRestartF

      a%OldRestartFolder = OldRestartF
   end subroutine

   subroutine php_NULLSUB_ifunc(a,ifunc)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: ifunc
   end subroutine

   subroutine php_MakeSteady(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a

      a%kfl_stead = 1
   end subroutine

   subroutine php_SkipSystemSolve(a,doskip)
      use typre
      implicit none
      class(PhysicalProblem) :: a
      logical :: doskip

      a%kfl_SkipSystemSolve = doskip
   end subroutine

   subroutine php_NULLSUB(a)
      use typre
      implicit none
      class(PhysicalProblem) :: a

   end subroutine

   subroutine Project_disc(a,ndofn,array,kfl_fixno)
      use typre
      use Mod_Element
      implicit none
      class(PhysicalProblem) :: a
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ndofn,idofn,ipoin,npoin,ielem,nelem,inode,igaus,disc
      real(rp)    :: array(:,:),dvol
      real(rp), allocatable :: vmass_disc(:)
      integer(ip) :: kfl_perio,kfl_HangingNodes,kfl_fixno(:)

      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%Alloc(npoin,vmass_disc,'vmass_disc','Project_disc')
      call a%Mesh%ElementAlloc(e,a%Memor,'ForceClosedRule','Project_disc')
      vmass_disc = 0.0_rp

      elements : do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         call e%elmdel
         disc = 0_ip
         do inode=1,e%pnode
            ipoin = e%lnods(inode)
            if (kfl_fixno(ipoin)==-3_ip) then
               disc = 1_ip
            endif
         enddo

         if (disc==0) then
            gauss_points: do igaus=1,e%pgaus
               e%igaus = igaus
               call e%elmder
               dvol = e%weigp(e%igaus)*e%detjm
               do inode=1,e%pnode
                  ipoin = e%lnods(inode)
                  vmass_disc(ipoin) = vmass_disc(ipoin) + dvol*e%shape(inode,e%igaus)
               enddo
            enddo gauss_points
         endif
      enddo elements

      call a%Mesh%ElementDealloc(e,a%Memor,'ForceClosedRule','Project_disc')

      !Communicate between subdomains
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,vmass_disc)

      !Periodic Boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(1_ip,vmass_disc)

      !Store the inverse of vmass
      do ipoin=1,npoin
         if (vmass_disc(ipoin) /= 0.0_rp) then
            vmass_disc(ipoin)=1.0_rp/vmass_disc(ipoin)
         endif
      end do
      do ipoin=1,npoin
         array(:,ipoin) = array(:,ipoin)*vmass_disc(ipoin)
      enddo

      !Hanging nodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(ndofn,array)

      !Communicate Ghosts
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndofn,array)

      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(ndofn,array)

   end subroutine

   subroutine SetGlobalIteration(a,iiter)
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: iiter

      a%OutIiter = iiter
   end subroutine

   subroutine SetCouplingIteration(a,cpiter)
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: cpiter

      a%cpiter = cpiter
   end subroutine

   subroutine SetMulticommData(a,kfl_multicomm,MulticommColor)
      implicit none
      class(PhysicalProblem) :: a
      integer(ip) :: kfl_multicomm,MulticommColor
      
      a%kfl_multicomm = kfl_multicomm
      a%MulticommColor = MulticommColor
   end subroutine
   

end module

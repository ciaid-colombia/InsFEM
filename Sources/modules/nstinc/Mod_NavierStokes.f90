module Mod_NavierStokes
   use typre
   use Mod_PhysicalProblem
   use Mod_CutMesh
   use Mod_MaterialProperties
   use Mod_TurbulentInletBoundaryConditions
   use Mod_TurbulentBodyForces
   use Mod_ExitBodyForces
   use Mod_PorousMedia
   use Mod_PlasmaActuator
   implicit none
   private
   public NavierStokesProblem,NavierStokesProblem_Const, zensi
   
   real(rp), parameter   :: &
     zensi = epsilon(0.0_rp)        ! zero     
  
   type, extends(PhysicalProblem) ::  NavierStokesProblem

      !Material Properties type 
      type(MatProperties) :: MatProp(5)         ! Maximum number of materials 
      integer(ip)  :: nmat                      ! Number of diferents material
      
      ! Logical units
      integer(ip) ::&
            lun_thpre, lun_wallo, lun_anatu,&
            lun_kespe, lun_mval,&
            lun_trap
      integer(ip), allocatable   :: lun_force(:)
      
      ! Physical problem
      real(rp), allocatable :: &
            veloc(:,:,:), &                     !Velocity field
            veloc_cp(:,:), &                    !Velocity field FSI coupling 
            press(:,:),&                        !Pressure field
            press_cp(:),&                       !Pressure field FSI coupling
            vorti(:,:),&                        !Vorticity
            qfac(:),&                           !Q-factor
            force(:,:),&                        !Body force
            momen(:,:),&                        !Body moment
            AnalyticalVelocity(:,:,:),&         !Exact velocity
            AnalyticalPressure(:,:)             !Exact pressure
            
      real(rp), pointer ::   &                  !For setting Dirichlet Boundary conditions on embedded methods
            DirichletVelocity(:,:) => NULL()
            
      integer(ip) :: lun_press,lun_avg
      character(150) :: fil_press,fil_avg
 
      type(r1p), allocatable :: &
            prsgs(:), &                         ! Pressure subgrid scales
            residualP(:)                        ! Pressure residual at gp, for postprocess
      type(r2p), allocatable :: &
            residualU(:),&                      ! Residual at gauss point, for postprocess
            residualGraP2(:),&                  ! Residual at gauss point, for postprocess, Split-OSS
            bvesgs(:)                           ! Velocity subgrid scales in boundaries
      type(r3p), allocatable :: &
            vesgs(:), &                         ! Velocity subgrid scales
            vesgs2(:)                           ! Velocity subgrid scales needed in dyn Split-Oss 
      integer(ip) :: ResidualSize
      
      integer(ip) ::&
            kfl_visco,&                         ! Viscous term
            kfl_advec,&                         ! Existence of (u.grad)u
            kfl_advco,&                         ! Oseen problem
            kfl_cotem,&                         ! Coupling with TEMPER
            kfl_confi,&                         ! Confined flow
            kfl_cotur,&                         ! Turbulence model
            kfl_stabm,&                         ! Stabilization method: 1: VMS, 0 SUPG, -1 GLS
            kfl_bousg,&                         ! Stabilization in the subscales
            kfl_colev,&                         ! Level set Coupling
            kfl_fsurf,&                         ! Free surface flag 0 don't use 1 on -1 off
            kfl_fsurfLapla, &                   ! Free surface with laplacian problem outside
            kfl_fsurfDirichlet = 0,&            ! Enforce Dirichlet Boundary Conditions on the free surface
            kfl_fsurfDirichletMethod = 1, &     ! Dirichlet method: 1: Nitshce, 2: Strong, 3: LLM
            kfl_SurfaceTension,&                ! Surface tension flag
            kfl_EnrichElem,&                    ! Enriched element used in level set
            kfl_EnrichPressure,&                ! =0 off  =1 The enrichment is in pressure
            kfl_EnrichVelocity, &               ! =0 off  =1 The enrichment is in velocity            
            kfl_PressTempSubscale               ! =0 off =1 Activate dynamic pressure term

      real(rp) ::&
            grnor    ,&                         ! Gravity norm
            gravi(3) ,&                         ! Gravity vector
            turbu(2) ,&                         ! Turbulence parameters
            advco(3) ,&                         ! Constant Convection term
            fcons    ,&                         ! Convection term
            fvins    ,&                         ! Viscous term
            boube    ,&                         ! Boussinesq volume expansion
            boutr    ,&                         ! Boussinesq reference temperature
            bougr(3)                            ! Boussinesq gravity field   

      type(r1p), allocatable :: &
            viscarray(:)                        ! Gauss point Viscosity 
      
      ! Output & postprocess
      integer(ip) ::&
            kfl_outfm         =-999,&           ! Output of forces and moments
            kfl_postBtract    =-999             ! Output boundary traction
     logical :: kfl_openforce = .true.          ! Open file for forces

      real(rp) ::&
            adimf( 3),&                         ! Adimensionalization factor for forces
            adimm( 3),&                         ! Adimensionalization factor for moments
            origm( 3)                           ! Moment origin (application) point
      
      !Statistics data
      real(rp) :: StatisticsStartTime = 0.0_rp
      real(rp) :: StatisticsCumulativeTime = 0.0_rp
      real(rp),allocatable  ::   &
         timav(:,:,:),&                         ! velocity time average
         timap(:,:),&                           ! pressure  time average
         rmsv(:,:),&                            ! R.M.S. (deviation) of velocity
         rmsp(:),&                              ! R.M.S. (deviation) of pressure
         turbi(:),&
         timaTractions(:,:,:)
         
      ! Numerical treatment
      integer(ip) ::&
            mtrit    ,&                         ! Maximum Number of iterations of tracking
            kfl_penal,&                         ! Penalization
            kfl_wtemp,&                         ! Weighting of du/dt
            kfl_wlapl,&                         ! Weighting of nu Lap(u)
            kfl_repro,&                         ! Stabilization based on residual projection
            kfl_repro_SkipFE,&                  ! Skip Finite Element Part of the Solution (from both LHS and RHS)
            kfl_shock,&                         ! Shock capturing type
            kfl_trasg,&                         ! Tracking of subgrid scale
            kfl_tacsg,&                         ! Time accuracy of subscales
            kfl_nolsg,&                         ! Non-linearity of the subscales
            kfl_nolsgNewtonRaphson,&            ! Newton Raphson iterative scheme for non-linear subgrid scales
            kfl_inist,&                         ! Initial Stokes
            kfl_local=0,&                       ! Local system of reference
            kfl_dispa,&                         ! Compute dissipation: 0 NO, 1:Only if postprocess, 2: Always
            kfl_tausm,&                         ! Use a smoothed tau from the previous iteration
            kfl_ExchangeLocationWallLaw,&       ! ExchangeLocationWallLaw
            kfl_FixWallLawDeltaForAll,&
            kfl_hdifumin,&                      ! Use minimum h for the diffusive stabilization term
            nodpr,&                             ! Node on which pressure is prescribed
            kfl_adapsgs                         ! 1: scaled L_2 norm, 2: entropy function
           
!           kfl_local_nsi,&                     ! Local system of reference
!           maxit_nsf,&                         ! Max # of iterations for the fractional step iterations
             
      real(rp), allocatable :: &
            Tausmo(:,:)                         !Smoothed Tau array (timom and tidiv)
            
      real(rp) ::&
            staco(4) ,&                         ! Stability constants
            shock(2) ,&                         ! Shock capturing parameters
            prepr=0.0_rp, &                     ! Prescribed pressure
            tosgs,    &                         ! Subgrid velocity internal tolerance
            relsg,    &                         ! Subgrid velocity relaxation parameter
            resip,    &                         ! Residual for outer iterations (p)
            weigh,    &                         ! Weight of dU/dt in the residual
            penal, &                            ! Penalty term if penalization is used
            WallLawDeltaValue,&
            erph1t(2), &                        ! Error p L2(H1)
            eruh1t(2), &
            erulit(2), &                        ! Error p Linf(L2)
            erplit(2)
            
      !Statistics
      integer(ip) ::  &
            nmean, &                            ! Number of statistics entries, for taking the mean 
            nmean_walllaw, &
            AVGnsteps
      
      real(rp) :: &      
            tamin,&                             ! Minimum tau
            tamax,&                             ! Maximum tau
            tamea,&                             ! Mean tau
            remin,&                             ! Minimum Re
            remax,&                             ! Maximum Re
            remea,&                             ! Mean Re
            ypmin,&                             ! Minimum y+
            ypmax,&                             ! Maximum y+
            ypmean                              ! Mean y+

      !Boundary conditions
      integer(ip), allocatable :: &
         kfl_fixrs(:),&                         ! Reference system for the BV
         kfl_bours(:)
      
      !Output and Postprocess
      real(rp), allocatable :: &
         veltp(:,:),&                           ! Velocty at tracking points
         pretp(:)                               ! Pressure at tracking points

      character(150) :: fil_trap

      !Dissipation and residual projection
      real(rp), allocatable :: &
         repro(:,:),&                           ! Residual Projection
         grprj(:,:,:),&                         ! Gradient finite element projection
         Dissipation(:)                         ! Sum of Numerical, turbulent and viscous dissipation

      type(r2p), allocatable :: NodalArrayDataForRefinement(:)

      !GaussPointDissipations
      type(r1p), allocatable :: GPDissipation(:)
      type(r1p), allocatable :: divergence(:)
         
      !Coupling with temperature
      integer(ip) :: kfl_ExternalTemperature = 0    
      real(rp), pointer :: &
         tempe(:) => NULL()                           ! Coupling with temperature problem
      !Non-Linear Subscales for Temperature   
      integer(ip) :: kfl_ExternalTemperatureSGS = 0   ! External temperature subgrid scales      
      type(r2p), pointer :: &
            tesgs(:)=> NULL()                         ! Temperature Subgrid scales
         
!       real(rp), allocatable :: &
!             tausm(:,:)                              ! Smoothed Tau for fractional step   

      !For reading boundary conditions, to be used by MPIroot
      integer(ip), allocatable :: &
         gkfl_bours(:),&                        ! Boundary reference system
         gkfl_fixrs(:),&                        ! Reference system for the BV
         aux_kfl_bours(:)

      !For spatially averaging velocity, along 1 dimension
      integer(ip) :: Avg1DIdime 
      
      !Coupling with LevelSet
      real(rp), pointer :: &
!          level(:)         => NULL() ,&        ! Coupling with levelset problem      
         CutGradient(:,:) => NULL()
         
      !Stabilizing free surface
      integer(ip) :: kfl_StabilizeFreeSurface = 0
      real(rp)    :: StabilizeFreeSurface_Param = 0.0_rp
      real(rp), allocatable :: StabilizeFreeSurface_VelocityGradient(:,:,:), StabilizeFreeSurface_PressureGradientAndForce(:,:)

      !Subrelaxation
      real(rp) :: subrelax
      
      !Boundary coupling
      integer(ip) :: kfl_computeTractions = 0
      logical :: doRobin = .false.
      real(rp) :: bstco(3)                      ! Stability constants for DG-DD
      real(rp), allocatable ::&
         PointwiseSource(:,:),&                 ! For FSI coupling
         btraction(:,:),&                       ! Outward traction
         btau(:)                                ! Boundary stabilization parameter (out)
      real(rp), pointer ::&
         eveloc(:,:) => NULL(),&                ! Boundary velocity for coupling
         etraction(:,:) => NULL(),&             ! Inward traction
         etau(:) => NULL()                      ! Boundary stabilization parameter (in)
      !Relaxation boundary coupling
      real(rp), allocatable :: velocDD(:,:,:)
      real(rp), allocatable :: tractDD(:,:,:)
      real(rp), allocatable :: bres(:,:,:)
      real(rp) :: relax = 1.0_rp                ! Relaxation param for mesh disp
      real(rp) :: relax_max = 1.0_rp            ! Relaxation param for mesh disp
      logical  :: kfl_doAitken  = .false.       ! Used to know Aitken relaxation should be used
      real(rp)    ::  alfa_robin = 0.0_rp
      
      !FrameOfReference (FOR) acceleration
      integer(ip) :: kfl_FORAcceleration = 0
      integer(ip) :: kfl_FORFunty
      real(rp)    :: FORParam(6)
      real(rp)    :: FORAcceleration(3) = 0.0_rp
      integer(ip) :: kfl_FORAxesRotation = 0         ! Non-inertial rotating frame
      real(rp)    :: FORAxesAngularVeloc(3) = 0.0_rp ! Angular velocity of rotating axis
      
      !For suction boundary conditions
      integer(ip) :: kfl_SuctionDirichletBC = 0

      !For spring type boundary conditions
      !It accounts for the presence of an elastic shell in the boundary
      integer(ip) :: kfl_ElasticBoundary
      real(rp), allocatable :: EB_Displacement(:,:)
      real(rp) :: EB_Stiffness,EB_density,EB_Damping,EB_ShearStiffness
      integer(ip) :: EB_NDelaySteps = 0

      integer(ip) :: sizeOfPointWise
      integer(ip) :: dimOfPointWise
      
      !Turbulent Inlet Boundary Conditions
      integer(ip) :: kfl_TurbulentInletBoundaryConditions = 0
      type(TurbulentInletBoundaryConditions) :: TIBC
      
      !Turbulent body forces
      integer(ip) :: kfl_TurbulentBodyForces = 0
      type(TurbulentInletForce) :: TBF
      
      !Exit Body Forces
      integer(ip) :: kfl_ExitBodyForces = 0
      type(ExitBodyForce) :: EBF
      
      !Extra Initial Viscosity
      integer(ip) :: kfl_ExtraInitialViscosity = 0, EIV_nsteps = 100
      real(rp) :: EIV_visco
      
      !Porous Media
      integer(ip) :: kfl_PorousMedia = 0
      type(PorousMedia) :: PME
      
      !Plasma Actuator
      integer(ip) :: kfl_Plasma = 0
      type(PlasmaAct) :: PAC
      
      !RVE periodic BC
      integer(ip) :: kfl_RVE
      real(rp) :: RVETensor(3,3)

      !Postprocessing External Forces
      real(rp), allocatable :: ExternalForcesArray(:,:)
      
      !Coriolis Force
      integer(ip) :: kfl_CoriolisForce = 0
      real(rp) :: CoriolisW(3) = 0.0_rp

contains

      procedure :: SetExmod             => nsi_SetExmod
      procedure :: SetNdofn             => nsi_SetNdofn
      procedure :: SetNdofbc            => nsi_SetNdofbc
      procedure :: SpecificReaphy       => nsi_reaphy
      procedure :: SpecificReanut       => nsi_reanut
      procedure :: SpecificReaous       => nsi_reaous
      procedure :: SpecificReaMPI       => nsi_reampi
      procedure :: NstincSpecificReampi => nsi_NULLSUB
      procedure :: SpecificBounod    => nsi_bounod
      procedure :: SpecificSolite    => nsm_solite
      procedure :: SpecificTurnon    => nsi_turnon
      procedure :: SpecificTurnof    => nsi_turnof
      procedure :: SpecificBegite    => nsi_begite
      procedure :: SpecificBegste    => nsi_begste
      procedure :: SpecificIniunk    => nsi_iniunk
      procedure :: SpecificEndste    => nsi_endste
      procedure :: SpecificEndite    => nsi_endite
      procedure :: SpecificCrankNicolsonEndste => nsi_CrankNicolsonEndste      
      procedure :: SpecificUpdbcs    => nsi_updbcs 
      procedure :: SpecificExaerr    => nsi_exaerr
      
      procedure :: SpecificReabcs => nsi_Reabcs
      procedure :: SpecificReadOnNodes => nsi_ReadOnNodes
      procedure :: SpecificReadSource =>  nsi_ReadSource
      procedure :: SpecificReadOnBoundaries => nsi_ReadOnBoundaries
      procedure :: SpecificReadOnFunctions => nsi_ReadOnFunctions
      
      procedure :: SpecificOuterr    => nsi_outerr
      procedure :: Memall            => nsi_memall
      procedure :: Getste            => nsi_getste
      procedure :: Cvgunk            => nsi_cvgunk
      procedure :: InnerResiduals    => nsi_InnerResiduals
      procedure :: Output            => nsi_output
      procedure :: PointTracking     => nsi_outtpo
      procedure :: SpecificRestart   => nsi_restar
      
      procedure :: Ifconf            => nsi_ifconf
      procedure :: Elmope            => nsm_Elmope
      procedure :: Bouope            => nsm_Bouope
      procedure :: EndElmope         => nsm_EndElmope
      procedure :: EndBouope         => nsm_EndBouope
      
      procedure :: InitStats         => nsm_InitStats
      procedure :: InGaussStats      => nsm_InGaussStats
      procedure :: FinalizeStats     => nsm_FinalizeStats
      
      procedure :: GetPhysicalParameters => nsi_GetPhysicalParameters
      procedure :: GetPhysicalParametersPointers
      
      procedure :: SetVelocityArray
      procedure :: GetVelocityArray

      procedure :: GetPressureArray
      procedure :: GetSigmaArray
      procedure :: GetVelocitySGS
      
      procedure :: GetDissipationArray
      
      procedure :: SetTemperatureArray
      procedure :: SetTemperatureSGS

      procedure :: SetTractionFlag
      procedure :: GetTractionArray
      procedure :: SetExternalTractionArray
      procedure :: GetBoundaryTau
      procedure :: SetBoundaryTau
     
      !ROM 
      procedure :: SpecificGetUnkno     => nsi_GetUnkno
      
      procedure :: SetParameters => nsi_SetParameters
      procedure :: GetParameters => nsi_GetParameters
     
      !LevelSet
      procedure :: GetEnrichedFlag
      procedure :: SetCutGradient
      procedure :: SetDirichletVelocity
      
      !Adaptivity
      procedure :: GetRefinementCriteria         => nsi_GetRefinementCriteria
      procedure :: PreRefine                     => nsi_Prerefine
      procedure :: SpecificRefine                => nsi_Refine
      procedure :: SpecificSubscalesRefCriteria  => nsi_SubscalesRefCriteria
      
      procedure :: SpecificProjectArraysOUM      => nsi_ProjectArraysOntoUndeformedMesh
      procedure :: SpecificAdvectArraysOUM       => nsi_AdvectArraysOntoUndeformedMesh

      procedure :: NstincSpecificOutput          => nsi_NULLSUBitask
      procedure :: NstincSpecificRestart         => nsi_NULLSUBrestart
      procedure :: GaussRestart                  => nsi_gaussRestart
      !RVE
      procedure :: RVE => nsi_RVE

   end type

   interface

      subroutine nsi_ReadOnFunctions(a,ifunc)
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: ifunc
      end subroutine
   
      subroutine nsi_ReadSource(a)
         import 
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine

      subroutine nsi_SetExmod(a)
         import
         implicit none         
         class(NavierStokesProblem) :: a      
      end subroutine
      
      subroutine nsi_SetNdofn(a)
         import
         implicit none         
         class(NavierStokesProblem) :: a      
      end subroutine
      
      subroutine nsi_SetNdofbc(a)
         import
         implicit none         
         class(NavierStokesProblem) :: a      
      end subroutine
     
      subroutine nsi_reaphy(a,itask)
         use typre
         import
         implicit none         
         integer(ip) :: itask
         class(NavierStokesProblem) :: a      
      end subroutine
   
      subroutine nsi_reanut(a,itask)
         use typre
         import
         implicit none         
         integer(ip) :: itask
         class(NavierStokesProblem) :: a      
      end subroutine
   
      subroutine nsi_reaous(a,itask)
         use typre
         import
         implicit none         
         integer(ip) :: itask
         class(NavierStokesProblem) :: a      
      end subroutine
   
      subroutine nsi_reampi(a)
         use typre
         import
         implicit none         
         class(NavierStokesProblem) :: a      
      end subroutine

      subroutine nsi_bounod(a,itask)
         use typre
         import
         implicit none
         integer(ip) :: itask
         class(NavierStokesProblem) :: a
      end subroutine
   
      subroutine nsi_outerr(a)
         import
         implicit none
         class(NavierStokesProblem) :: a    
      end subroutine

      subroutine nsi_exaerr(a)
         import
         implicit none
         class(NavierStokesProblem) :: a    
      end subroutine      
      
      subroutine nsi_memall(a)
         import
         implicit none
         class(NavierStokesProblem) :: a      
      end subroutine
      
      subroutine nsi_getste(a,dtinv)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp) :: dtinv
      end subroutine
      
      subroutine nsi_begste(a)
         import
         implicit none
         class(NavierStokesProblem) :: a      
      end subroutine
      
      subroutine nsi_begite(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_endite(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
      end subroutine
      
      subroutine nsi_cvgunk(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine
      
      subroutine nsi_InnerResiduals(a,rinsi,rprnsi)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp) :: rinsi,rprnsi
      end subroutine

      subroutine nsi_outtpo(a)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_output(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
      end subroutine
      
      subroutine nsi_endste(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
      end subroutine
            
      subroutine nsi_ifconf(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_iniunk(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_updbcs(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine      
      
      subroutine nsm_solite(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
      end subroutine
      
      subroutine nsm_Elmope(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsm_Bouope(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine

      subroutine nsi_turnon(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine

      subroutine nsi_turnof(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_CrankNicolsonEndste(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_restar(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine

      subroutine nsi_gaussRestart(a,itask)
         use typre
         import 
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip), intent(in) :: itask
         
      end subroutine
      
      subroutine nsm_InitStats(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsm_InGaussStats(a,acden,acvis,gpvno,chale,timom)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp) :: acden,acvis,gpvno,chale(2),timom
      end subroutine
      
      subroutine nsm_FinalizeStats(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine GetPhysicalParametersPointers(a,imat,acden,acvis)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp), pointer    :: acden, acvis
         integer(ip) :: imat
      end subroutine

      subroutine nsi_GetPhysicalParameters(a,imat,acden,acvis)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp)    :: acden, acvis
         integer(ip) :: imat 
      end subroutine

      subroutine nsi_Reabcs(a,itask,kflag)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine
      
      subroutine nsi_ReadOnNodes(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_ReadOnBoundaries(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsm_EndElmope(NSProblem,task)
         import
         implicit none
         class(NavierStokesProblem), target :: NSProblem
         character(6) :: task
      end subroutine
      
      subroutine nsm_EndBouope(NSProblem,task)
         import
         implicit none
         class(NavierStokesProblem), target :: NSProblem
         character(6) :: task
      end subroutine
      
      subroutine nsi_GetRefinementCriteria(a,markel)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: markel(*)
      end subroutine
      
      subroutine nsi_SubscalesRefCriteria(a,error,TotalEstimatedError)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         real(rp) :: error(:), TotalEstimatedError
      end subroutine

      subroutine nsi_Prerefine(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_Refine(a,itask)
         import
         implicit none
         class(NavierStokesProblem) :: a
         character(6) :: itask
      end subroutine
      
      subroutine nsi_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
         use typre
         use Mod_MeshInterpolator
         import
         implicit none
         class(NavierStokesProblem) :: a
         type(Interpolator) :: Interp
         integer(ip) :: itask
      end subroutine

      subroutine nsi_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
         use typre
         use Mod_Advector
         import
         implicit none
         class(NavierStokesProblem) :: a
         type(Advector) :: Advect
         integer(ip) :: itask
      end subroutine
  
      subroutine nsi_ManufacturedBoundaryConditions(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsi_RVE(a,itask)
         use typre
         import
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
      end subroutine
      
   end interface
   
  interface NavierStokesProblem_Const
      procedure constructor
  end interface NavierStokesProblem_Const

contains
  
   function constructor()
      class(NavierStokesProblem), pointer :: constructor

      allocate(constructor)
   end function constructor
   
   subroutine nsi_NULLSUB(a)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
   
   end subroutine

   subroutine nsi_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      integer(ip) :: itask
      
   end subroutine

   subroutine nsi_NULLSUBrestart(a,itask)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine

   subroutine GetVelocityArray(a,veloc)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp), pointer :: veloc(:,:)
      
      veloc => a%veloc(:,:,1)
   end subroutine

   subroutine GetVelocitySGS(a,vesgs)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      type(r3p), pointer :: vesgs(:)
      
      vesgs => a%vesgs
   end subroutine
   
   subroutine GetEnrichedFlag(a,kfl_EnrichElem)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      integer(ip), pointer :: kfl_EnrichElem
      
      kfl_EnrichElem => a%kfl_EnrichElem
   end subroutine   
   
   subroutine GetPressureArray(a,press)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp),pointer :: press(:)
      
      press => a%press(:,1)
   end subroutine
   
    subroutine GetSigmaArray(a,sigma)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp),pointer :: sigma(:,:)
      integer(ip) :: kfl_cosig
      sigma => NULL()
   end subroutine
   
   subroutine GetDissipationArray(a,dissipation)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp), pointer :: dissipation(:)
      
      if (a%kfl_dispa == 2) then
         dissipation => a%dissipation
      else
         dissipation => null()
      endif   
   end subroutine
   
   subroutine SetTemperatureArray(a,tempe)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      real(rp), target :: tempe(:)
      
      a%tempe => tempe
      a%kfl_ExternalTemperature = 1
   end subroutine
   
   subroutine SetTemperatureSGS(a,tesgs)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      type(r2p), target :: tesgs(:)
      
      a%tesgs => tesgs
      a%kfl_ExternalTemperatureSGS = 1
   end subroutine
   
   subroutine SetTractionFlag(a)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      
      a%kfl_computeTractions = 1
   end subroutine
   
   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      real(rp), target :: veloc(:,:)
      
      a%eveloc => veloc
   end subroutine
   
   subroutine nsi_GetUnkno(a,unkno)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      real(rp) :: unkno(:,:)
           
      if (a%kfl_cotem == 1) unkno(a%ndofn+1,:) = a%tempe
   end subroutine

   !Implementation of php_SetParameters for NavierStokes
   subroutine nsi_SetParameters(a,iparam,param)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      integer(ip) :: iparam
      real(rp)    :: param
      
      !For setting the viscosity
      if (iparam == 1) then
         a%MatProp(1)%visco = param
      endif
   end subroutine
   
   ! Implementation of php_GetParameters for NavierStokes
   subroutine nsi_GetParameters(a,iparam,param)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      integer(ip) :: iparam
      real(rp)    :: param
      
      !For getting the viscosity
      if (iparam == 1) then
         param = a%MatProp(1)%visco
      endif
   end subroutine

   subroutine GetTractionArray(a,trac)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp), pointer :: trac(:,:)

      trac => a%btraction
   end subroutine

   subroutine SetExternalTractionArray(a,trac)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      real(rp), target :: trac(:,:)

      a%etraction => trac
   end subroutine

   subroutine GetBoundaryTau(a,tau)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp), pointer :: tau(:)
      
      tau => a%btau
   end subroutine

   subroutine SetBoundaryTau(a,tau)
      use typre
      implicit none
      class(NavierStokesProblem), target :: a
      real(rp), target :: tau(:)
      
      a%etau => tau
   end subroutine

   subroutine SetCutGradient(a,CutGradient)
      use typre
      implicit none
      class(NavierStokesProblem) :: a
      real(rp), target :: CutGradient(:,:)
   
      a%CutGradient => CutGradient
   end subroutine
   
   subroutine SetDirichletVelocity(a,DirichletVelocity)
      implicit none
      class(NavierStokesProblem) :: a
      real(rp), target :: DirichletVelocity(:,:)
      
      a%DirichletVelocity => DirichletVelocity
   end subroutine

end module Mod_NavierStokes


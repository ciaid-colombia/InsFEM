module Mod_Temperature
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_tem_TemperatureMaterial
   implicit none
   private
   public TemperatureProblem,TemperatureProblem_Const,zetem
   
   real(rp), parameter :: &
       zetem = epsilon(1.0_rp)
   
   type, extends(PhysicalProblem) :: TemperatureProblem
   
      !Physical Problem
      real(rp), allocatable :: tempe(:,:),heatf(:)
      
      ! Logical units
      integer(ip), allocatable   :: lun_force(:)
      
      !TemperatureMaterials
      integer(ip) :: NumberOfMaterials
      type(TemperatureMaterial), allocatable :: Materials(:)
      integer(ip), allocatable :: ElementMaterials(:)
      
      !Problem data and flags
      integer(ip) ::&
         kfl_advec,&              ! Existence of (u.grad)T
         kfl_joule,&              ! Existence of Joule effect
         kfl_sourc,&              ! Existence and type of source term
         kfl_cotur,&              ! Coupling with a turbulence model
         kfl_wtemp,&              ! Weighting of du/dt
         kfl_wlapl,&              ! Weighting of nu Lap(u)
         kfl_repro,&              ! Stabilization based on residual projection
         kfl_shock,&              ! Shock capturing type 
         kfl_trasg,&              ! Tracking of subgrid scale
         kfl_tacsg,&              ! Time accuracy of subscales
         kfl_nolsg,&              ! Non-linearity of the subgrid scales
         kfl_splot,&              ! Output for 3D gnuplot flag
         kfl_dispa,&              ! Calculate dissipation
         kfl_stabm,&              ! Stabilization Method: 1 VMS, 0 SUPG, -1 GLS
         ipsou(2),&               ! Source integer parameters
         kfl_CouplingThreeField,& ! Coupling with Sigmaup model
         kfl_adapsgs              ! 1: scaled L_2 norm, 2: entropy function
      
      real(rp) ::&
         sourc(2),&               ! Source term value
         rpsou(10),&              ! Source real parameters
         staco( 3),&              ! Stability constants
         shock    ,&              ! Shock capturing parameter
         relsg,    &              ! Subgrid velocity relaxation parameter
         turbu,    &              ! Turbulence parameters
         relax                    ! Relaxation factor
      real(rp),     allocatable ::&
         PointwiseSource(:)       ! Power of sources (Q)
        
      !Numerical
      real(rp), allocatable :: &
         repro(:), &                  ! Residual Projection
         dissipation(:),&             ! Dissipation
         gradient(:,:),&
         grprj(:,:),&                 ! Gradient finite element projection
         ViscousDissipation(:)  
      
      type(r2p), allocatable :: NodalArrayDataForRefinement(:)
         
      type(r1p), allocatable ::  &
         ShockCapturingViscosity(:), &
         sigmatermarray(:)
         
      type(r2p), allocatable :: &
         GradientGaussPoints(:)
         
         
      type(r2p),allocatable ::  &
       tesgs(:)              ! subgrid scale temperature TEMPER   
      
      !External forces 
      real(rp), allocatable :: &
            tfsou(:)                 ! Source time function parameters
      real(rp) ::&
            react                    ! Reaction term

      ! Boundary conditions and turbulence
      real(rp) :: & 
            visco,&                  ! Viscosity for the Joule effect (mu) and wall law
            prtur,&                  ! Turbulent Prandtl number  
            delta                    ! Distance to the wall    
            
      ! Output & postprocess
      integer(ip) :: kfl_outfm                  ! Output of forces and moments
      logical :: kfl_openforce = .true.          ! Open file for forces

      integer(ip) :: &
         kfl_ExternalAdvection=0  ! External advection (Nstinc, for instance)   
      !Velocity pointer, obtained from another object (NSI, for instance)
      real(rp), pointer     :: veloc(:,:) => NULL()
      real(rp), pointer     :: sigmaNS(:,:) => NULL()   !Coupling with three field problem
      integer(ip) :: kfl_ExternalAdvectionSGS=0
      type(r3p), pointer    :: vesgs(:) => NULL()    
      
      !For spatially averaging temperature, along 1 dimension
      integer(ip) :: Avg1DIdime
      
      !For Dust Transport
      integer(ip) :: DT_kfl_DustTransport = 0
      real(rp)    :: DT_ParticleSize = 1e-6
      real(rp)    :: DT_MediumDensity = 1.18
      real(rp)    :: DT_MediumDynamicViscosity =1.8e-5
      real(rp)    :: DT_GravityForce(3) = (/ 0.0_rp, -9.81_rp, 0.0_rp /)
      !real(rp)    :: DT_SurfaceWetness
      !To be computed
      real(rp)    :: DT_DustDepositionVelocityNorm
      real(rp)    :: DT_DustDepositionVelocityVector(3) = 0.0_rp
      real(rp)    :: DT_ThresholdFrictionVelocity
      
      real(rp), allocatable :: SmoothedVelocityGradient(:,:,:)
      
      
             
contains
      
      procedure :: SetExmod          => tem_SetExmod
      procedure :: SetNdofn          => tem_SetNdofn
      procedure :: SetNdofbc         => tem_SetNdofbc
      procedure :: SpecificReaphy    => tem_reaphy
      procedure :: SpecificReanut    => tem_reanut
      procedure :: SpecificReaous    => tem_reaous
      procedure :: SpecificReaMPI    => tem_reampi
      procedure :: SpecificBounod    => tem_bounod
      procedure :: SpecificReabcs => tem_Reabcs
      procedure :: SpecificManufacturedBoundaryConditions => tem_ManufacturedBoundaryConditions
      procedure :: SpecificReadOnNodes => tem_ReadOnNodes
      procedure :: SpecificReadOnBoundaries => tem_ReadOnBoundaries
      procedure :: SpecificReadOnElements =>   tem_ReadOnElements
      
      procedure :: SpecificOuterr            => tem_outerr
      procedure :: Memall            => tem_memall
      procedure :: SpecificIniunk    => tem_iniunk
      
      procedure :: Getste            => tem_getste
      procedure :: SpecificBegste    => tem_begste
      procedure :: SpecificBegite    => tem_begite
      procedure :: SpecificSolite    => tem_solite
      procedure :: SpecificEndite    => tem_endite
      procedure :: EnditeElmope      => tem_EnditeElmope
      procedure :: SpecificEndste    => tem_endste
      procedure :: SpecificCrankNicolsonEndste => tem_CrankNicolsonEndste
      Procedure :: SpecificUpdbcs    => tem_NULLSUB 
      Procedure :: SpecificExaerr    => tem_NULLSUB
      
      procedure :: Cvgunk            => tem_cvgunk
      procedure :: Output            => tem_output
      procedure :: Elmope            => tem_elmope
      procedure :: Bouope            => tem_bouope
      
      procedure :: SpecificTurnon    => tem_turnon
      procedure :: SpecificTurnof    => tem_turnof
     
      procedure :: SpecificRestart   => tem_restar
      
      procedure :: GetPhysicalParameters => tem_GetPhysicalParameters
      
      !Boundary operations 
      procedure :: EndBouope => tem_EndBouope
    
      !ROM
      procedure :: SpecificGetUnkno     => tem_GetUnkno
     
      !Coupling
      procedure :: SetVelocityArray
      procedure :: SetSigmaNSArray
      
      procedure :: GetTemperatureArray
      procedure :: GetDissipationArray
      
      procedure :: GetTemperatureSGS
      procedure :: SetVelocitySGS
      
      !For Adaptive refinement
      procedure :: PreRefine                     => tem_Prerefine
      procedure :: SpecificRefine                => tem_Refine
      procedure :: GetRefinementCriteria         => tem_GetRefinementCriteria
      procedure :: SpecificSubscalesRefCriteria  => tem_SubscalesRefCriteria
      

      
   end type
   
   interface

      subroutine tem_SetExmod(a) 
         import TemperatureProblem
         implicit none
         
         
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_SetNdofn(a)
         import TemperatureProblem
         implicit none
         
         
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_SetNdofbc(a)
         import TemperatureProblem
         implicit none
         
         
         class(TemperatureProblem) :: a
      
      end subroutine
   
      subroutine tem_reaphy(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         
         integer(ip) :: itask
         class(TemperatureProblem) :: a
      
      end subroutine
   
      subroutine tem_reanut(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         
         integer(ip) :: itask
         class(TemperatureProblem) :: a
      
      end subroutine
   
      subroutine tem_reaous(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         
         integer(ip) :: itask
         class(TemperatureProblem) :: a
      
      end subroutine
   
      subroutine tem_reampi(a)
         use typre
         import TemperatureProblem
         implicit none
         
         class(TemperatureProblem) :: a
      
      end subroutine
         
      subroutine tem_bounod(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         
         integer(ip) :: itask
         class(TemperatureProblem) :: a
      
      end subroutine
   
      subroutine tem_outerr(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      
      end subroutine
      
      subroutine tem_memall(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      
      end subroutine
      
      subroutine tem_getste(a,dtinv)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         real(rp) :: dtinv
      
      
      end subroutine
      
      subroutine tem_begste(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      
      end subroutine
      
      subroutine tem_begite(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      
      end subroutine
      
      subroutine tem_endite(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: itask
      
      
      end subroutine
      
      subroutine tem_cvgunk(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip), intent(in) :: itask
      
      
      end subroutine
      
      subroutine tem_output(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: itask
      
      
      end subroutine
      
      subroutine tem_endste(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: itask
      
      end subroutine
            
      subroutine tem_iniunk(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
 
      subroutine tem_updbcs(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine  
      
      
      subroutine tem_solite(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine tem_elmope(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_EnditeElmope(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_bouope(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_turnon(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_turnof(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_CrankNicolsonEndste(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         
      end subroutine
      
      subroutine tem_restar(a,itask)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine
            
      subroutine tem_GetPhysicalParameters(a,ielem,acden,acsph,actco,acrea,acsou)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: ielem
         real(rp) :: acden, acsph,actco,acrea,acsou
         
      end subroutine
      
      subroutine tem_reabcs(a,itask,kflag)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      
      end subroutine
      
      subroutine tem_ReadOnNodes(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_ReadOnBoundaries(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine
      
      subroutine tem_ReadOnElements(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine

      subroutine tem_Prerefine(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      end subroutine
      
      subroutine tem_Refine(a,itask)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         character(6) :: itask
      end subroutine
      
      subroutine tem_GetRefinementCriteria(a,markel)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine tem_SubscalesRefCriteria(a,error,TotalEstimatedError)
         use typre
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
         real(rp) :: error(:),TotalEstimatedError
      end subroutine

      subroutine tem_ManufacturedBoundaryConditions(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a
      
      end subroutine

      subroutine tem_EndBouope(a)
         import TemperatureProblem
         implicit none
         class(TemperatureProblem) :: a     
      end subroutine

   end interface

  interface TemperatureProblem_Const
      procedure constructor
  end interface TemperatureProblem_Const
   
contains

   function constructor()
      class(TemperatureProblem), pointer :: constructor

      allocate(constructor)
   end function constructor

   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(TemperatureProblem) :: a
      real(rp), target :: veloc(:,:)
      
      a%veloc => veloc
      a%kfl_ExternalAdvection = 1
   end subroutine
   
   subroutine SetSigmaNSArray(a,sigmaNS)
      use typre
      implicit none
      class(TemperatureProblem) :: a
      real(rp), target :: sigmaNS(:,:)

       a%sigmaNS => sigmaNS
   end subroutine
   
   subroutine GetTemperatureArray(a,tempe)
      use typre
      implicit none
      class(TemperatureProblem), target :: a
      real(rp), pointer :: tempe(:)
      
      if(allocated(a%tempe)) then

          tempe => a%tempe(:,1)
      end if
   end subroutine
   
   subroutine GetDissipationArray(a,dissipation)
      use typre
      implicit none
      class(TemperatureProblem), target :: a
      real(rp), pointer :: dissipation(:)
      
      if (a%kfl_dispa == 2) then
         dissipation => a%dissipation
      else
         dissipation => null()
      endif   
   end subroutine
   
   subroutine GetTemperatureSGS(a,tesgs)
      use typre
      implicit none
      class(TemperatureProblem), target :: a
      integer(ip) :: kfl_trasg
      type(r2p), pointer :: tesgs(:)
      
      tesgs => a%tesgs

   end subroutine
   
   subroutine SetVelocitySGS(a,vesgs)  
      use typre
      implicit none
      class(TemperatureProblem) :: a
      type(r3p), target :: vesgs(:)
      
      a%vesgs => vesgs
      a%kfl_ExternalAdvectionSGS = 1
   end subroutine
   
   subroutine tem_NULLSUB(a)
      implicit none
      class(TemperatureProblem) :: a
      
   end subroutine
   
   subroutine tem_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(TemperatureProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine
   
   subroutine tem_GetUnkno(a,unkno)
      use typre
      implicit none
      class(TemperatureProblem) :: a
      real(rp) :: unkno(:,:)
            
   end subroutine


end module

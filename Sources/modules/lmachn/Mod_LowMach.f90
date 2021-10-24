module Mod_LowMach
   use typre
   use Mod_PhysicalProblem
   implicit none
   private
   public LowMachProblem,LowMachProblem_Const, zelmn

   real(rp), parameter   :: &
     zelmn = epsilon(0.0_rp)        ! zero
    
   type, extends(PhysicalProblem) ::  LowMachProblem
      
      ! Logical units
      integer(ip), allocatable   :: lun_force(:)
      integer(ip)                :: lun_trap         ! Tracking of points
      character(150)             :: fil_trap         ! Tracking of points
      
      ! Physical problem
      real(rp), allocatable :: &
         press(:,:),&                        ! Pressure field
         veloc(:,:,:),&                      ! Velocity field
         tempe(:,:),&                        ! Temperature field
         itemp(:),&                          ! Initial temperature
         pther(:),&                          ! Thermodynamic Pressure field
         densf(:),&                          ! Density
         force(:,:),&                        ! Body force
         momen(:,:),&                        ! Body moment
         heatf(:)                            ! Heat flux

      logical :: kfl_pther                   ! Pther not constant
      type(r1p), allocatable :: &
         prsgs(:)                            ! Pressure subgrid scales
      type(r2p), allocatable :: &
         tesgs(:)                            ! Temperature subgrid scales
      type(r3p), allocatable :: &
         vesgs(:)                            ! Velocity subgrid scales
           
      type(r2p), allocatable :: &
         residualU(:)                        ! Velocity residual at gp for postprocess
      type(r1p), allocatable :: &
         residualP(:),&                      ! Pressure residual at gp for postprocess
         residualT(:)                        ! Temperature residual at gp for postprocess

      integer(ip) ::&
            kfl_eqnst,&                      ! Equation of state
            kfl_visco,&                      ! Viscous term
            kfl_advec,&                      ! Existence of (u.grad)u
            kfl_confi,&                      ! Confined flow
            kfl_stabm,&                      ! Stabilization method: 1: VMS, 0 SUPG, -1 GLS   
            kfl_sourc                        ! heat sources are off

      real(rp) ::&
            grnor    ,&                      ! Gravity norm
            gravi(3) ,&                      ! Gravity vector
            fvins                            ! Viscous term

            
      ! Physical properties
      real(rp) ::&
            itpre,&                          ! Initial thermodynamic pressure
            cphea,&                          ! Specific heat constant pressure (Cp)
            tcond,&                          ! Thermal conduction (k)
            visco,&                          ! Viscosity (mu)
            densi,&                          ! Initial density (rho)
            sgasc,&                          ! Ideal equation specific gas constant(R)
            texpc,&                          ! Thermal expansion coeficcient (alpha)
            react,&                          ! Reaction term
            sourc                            ! source term

      real(rp),     allocatable ::&
         PointwiseSource(:)                  ! heat sources pointwise (Q)

      real(rp), allocatable :: &
            tfsou(:)                         ! Source time function paramters

      type(r1p), allocatable :: &
         viscarray(:)                        ! Gauss point Viscosity

      ! Boundary conditions
      real(rp) :: &
         delta,&                             ! Distance to the wall
         epspe,&                             ! Penalization for confined flow
         nlstep                              ! step coefficient nl

      ! Volume (to be safe by each processor)
      real(rp) :: &
         dvolt                               ! Total volume

      ! Output & postprocess
      integer(ip) :: &
         kfl_outfm                           ! Output of forces and moments
     logical :: kfl_openforce = .true.       ! Open file for forces

      real(rp) :: &
         adimf(3),&                          ! Adimensionalization factor for forces
         adimm(3),&                          ! Adimensionalization factor for moments
         adimh,   &                          ! Adimensionalization factor for heat flux
         origm(3)                            ! Moment origin (application) point
            
      ! Numerical treatment
      integer(ip) ::&
         mtrit    ,&                         ! Maximum Number of iterations of tracking
         kfl_repro,&                         ! Stabilization based on residual projection
         kfl_trasg,&                         ! Tracking of subgrid scale
         kfl_tacsg,&                         ! Time accuracy of subscales
         kfl_nolsgScheme,&                   ! Linearization model for the subscale
         kfl_nolsg,&                         ! Non-linearity of the subscales
         kfl_linop,&                         ! Linearization method used
         kfl_dispa,&                         ! Compute dissipation: 0 NO, 1:Only if postprocess, 2: Always
         nodpr,&                             ! Node on which pressure is prescribed
         kfl_acite,&                         ! # picard iterations in N-R
         kfl_adapsgs                         ! 1: scaled L_2 norm, 2: entropy function
           
      real(rp) ::&
         staco(4) ,&                         ! Stability constants
         prepr    ,&                         ! Prescribed pressure
         tosgs,    &                         ! Subgrid velocity internal tolerance
         relsg,    &                         ! Subgrid velocity relaxation parameter
         rilmn,    &                         ! Linearization residual fe
         resiv,    &                         ! Residual for outer iterations (v)
         resip,    &                         ! Residual for outer iterations (p)
         resit,    &                         ! Residual for outer iterations (t)
         rhsnorm                             ! Residual norm Newton method

      !Statistics
      integer(ip) ::  &
         nmean                               ! Number of statistics entries, for taking the mean 
      
      real(rp) :: &      
         tamin,&                             ! Minimum tau
         tamax,&                             ! Maximum tau
         tamea,&                             ! Mean tau
         remin,&                             ! Minimum Re
         remax,&                             ! Maximum Re
         remea                               ! Mean Re

      !Output and Postprocess
      real(rp), allocatable :: &
         veltp(:,:),&                        ! Velocty at tracking points
         pretp(:),&                          ! Pressure at tracking points
         temtp(:,:),&                        ! Temperature at tracking points
         tprtp(:,:)                          ! Thermodynamic pressure at tracking points

      !Auxiliar working arrays
      real(rp), allocatable :: &
         repro(:,:), &                       ! Residual Projection
         grprj(:,:,:)                        ! Gradient finite element projection
      type(r2p), allocatable :: NodalArrayDataForRefinement(:)
       

contains

      !specific deferred procedures from PhysicalProblem
      procedure :: SetExmod          => lmn_SetExmod
      procedure :: SetNdofn          => lmn_SetNdofn
      procedure :: SetNdofbc         => lmn_SetNdofbc

      procedure :: SpecificReaphy    => lmn_reaphy
      procedure :: SpecificReanut    => lmn_reanut
      procedure :: SpecificReaous    => lmn_reaous
      procedure :: SpecificReaMPI    => lmn_reampi

      procedure :: SpecificReabcs =>  lmn_Reabcs
      procedure :: SpecificReadOnNodes =>  lmn_ReadOnNodes
      procedure :: SpecificReadOnBoundaries =>  lmn_ReadOnBoundaries

      procedure :: SpecificBounod    => lmn_NULLSUBitask
      procedure :: SpecificSolite    => lmn_NULLSUBitask!lmn_solite
      procedure :: SpecificTurnon    => lmn_NULLSUB!lmn_turnon 
      procedure :: SpecificTurnof    => lmn_turnof
      procedure :: SpecificBegite    => lmn_begite
      procedure :: SpecificBegste    => lmn_begste
      procedure :: SpecificUpdbcs    => lmn_NULLSUB!lmn_updbcs 
      procedure :: SpecificIniunk    => lmn_iniunk
      procedure :: SpecificEndste    => lmn_endste
      procedure :: SpecificEndite    => lmn_endite
      procedure :: SpecificCrankNicolsonEndste => lmn_CrankNicolsonEndste
      procedure :: SpecificExaerr    => lmn_NULLSUB
       
      
      procedure :: SpecificOuterr    => lmn_outerr
      procedure :: Memall            => lmn_memall
      procedure :: Getste            => lmn_getste
      procedure :: Cvgunk            => lmn_cvgunk
      procedure :: Output            => lmn_output
      procedure :: Elmope            => lmn_elmope
      procedure :: Bouope            => lmn_bouop0
      procedure :: SpecificRestart   => lmn_restar
      
           !own procedures of the class
      procedure :: InnerResiduals    => lmn_InnerResiduals
      procedure :: Ifconf            => lmn_ifconf
      
      procedure :: InitStats         => lmn_InitStats
      procedure :: InGaussStats      => lmn_InGaussStats
      procedure :: FinalizeStats     => lmn_FinalizeStats
     
      procedure :: GetPhysicalParameters => lmn_GetPhysicalParameters
      procedure :: ComputeThermPressure  => lmn_ComputeThermPressure
      procedure :: ComputeDensity        => lmn_ComputeDensity
      
      !Boundary operations 
      procedure :: EndBouope => lmn_EndBouope
      
      ! ROM
      procedure :: SpecificGetUnkno     => lmn_GetUnkno
      
      !For Adaptive refinement
      procedure :: PreRefine                     => lmn_Prerefine
      procedure :: SpecificRefine                => lmn_Refine
      procedure :: GetRefinementCriteria         => lmn_GetRefinementCriteria
      procedure :: SpecificSubscalesRefCriteria  => lmn_SubscalesRefCriteria
    
   end type

   interface
     
      subroutine lmn_SetExmod(a)
         import LowMachProblem
         implicit none         
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_SetNdofn(a)
         import LowMachProblem
         implicit none         
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_SetNdofbc(a)
         import LowMachProblem
         implicit none         
         class(LowMachProblem) :: a      
      end subroutine
     
      subroutine lmn_reaphy(a,itask)
         use typre
         import LowMachProblem
         implicit none         
         integer(ip) :: itask
         class(LowMachProblem) :: a      
      end subroutine
   
      subroutine lmn_reanut(a,itask)
         use typre
         import LowMachProblem
         implicit none         
         integer(ip) :: itask
         class(LowMachProblem) :: a      
      end subroutine
   
      subroutine lmn_reaous(a,itask)
         use typre
         import LowMachProblem
         implicit none         
         integer(ip) :: itask
         class(LowMachProblem) :: a      
      end subroutine
   
      subroutine lmn_reaMPI(a)
         use typre
         import LowMachProblem
         implicit none         
         class(LowMachProblem) :: a      
      end subroutine

      subroutine lmn_Reabcs(a,itask,kflag)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine
 
      subroutine lmn_ReadOnNodes(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
      end subroutine
      
      subroutine lmn_ReadOnBoundaries(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
      end subroutine

      subroutine lmn_outerr(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine

      subroutine lmn_memall(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_getste(a,dtinv)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         real(rp) :: dtinv      
      end subroutine
      
      subroutine lmn_begste(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_begite(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_endite(a,itask)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip) :: itask     
      end subroutine
      
      subroutine lmn_cvgunk(a,itask)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip), intent(in) :: itask   
      end subroutine
      
      subroutine lmn_InnerResiduals(a,rvlmn,rplmn,rtlmn)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         real(rp) :: rvlmn,rplmn,rtlmn
      end subroutine
      
      subroutine lmn_output(a,itask)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip) :: itask    
      end subroutine
      
      subroutine lmn_endste(a,itask)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip) :: itask    
      end subroutine
            
      subroutine lmn_ifconf(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a 
      end subroutine
      
      subroutine lmn_iniunk(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
      end subroutine 
      
!      subroutine lmn_updbcs(a)
!         import LowMachProblem
!         implicit none
!         class(LowMachProblem) :: a
!      end subroutine      

!      subroutine lmn_solite(a,itask)
!         use typre
!         import LowMachProblem
!         implicit none
!         class(LowMachProblem) :: a
!         integer(ip) :: itask      
!      end subroutine
      
      subroutine lmn_elmope(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_bouop0(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a     
      end subroutine
      
      subroutine lmn_turnof(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine

      subroutine lmn_CrankNicolsonEndste(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a       
      end subroutine
      
      subroutine lmn_restar(a,itask)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip), intent(in) :: itask       
      end subroutine
      
      subroutine lmn_InitStats(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_InGaussStats(a,acden,acvis,gpvno,chale,timom)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         real(rp) :: acden,acvis,gpvno,chale(2),timom
      end subroutine
      
      subroutine lmn_FinalizeStats(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a      
      end subroutine
      
      subroutine lmn_GetPhysicalParameters(a,acden,acvis,acsph,actco,actex,acgas,acpth,acrea,acsou)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         real(rp), optional    :: acden,acvis,acsph,actco,actex,acgas,acrea,acsou
         real(rp), optional    :: acpth(*) 
      end subroutine
      
      subroutine lmn_ComputeDensity(a)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem)  :: a
      end subroutine
      
      subroutine lmn_ComputeThermPressure(a)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem)  :: a  
      end subroutine

      subroutine lmn_EndBouope(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a     
      end subroutine
      
      subroutine lmn_Prerefine(a)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
      end subroutine
      
      subroutine lmn_Refine(a,itask)
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         character(6) :: itask
      end subroutine
      
      subroutine lmn_GetRefinementCriteria(a,markel)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         integer(ip) :: markel(*)
      end subroutine
 
      subroutine lmn_SubscalesRefCriteria(a,error,TotalEstimatedError)
         use typre
         import LowMachProblem
         implicit none
         class(LowMachProblem) :: a
         real(rp) :: error(:), TotalEstimatedError
      end subroutine

   end interface

    
  interface LowMachProblem_Const
      procedure constructor
  end interface LowMachProblem_Const
   
contains

function constructor()
    class(LowMachProblem), pointer :: constructor

    allocate(constructor)

end function constructor

   subroutine lmn_NULLSUB(a)
      use typre
      implicit none
      class(LowMachProblem) :: a
   
   end subroutine

   subroutine lmn_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(LowMachProblem) :: a
      integer(ip) :: itask
      
   end subroutine

   subroutine lmn_GetUnkno(a,unkno)
      use typre
      implicit none
      class(LowMachProblem) :: a
      real(rp) :: unkno(:,:)
            
   end subroutine

end module Mod_LowMach

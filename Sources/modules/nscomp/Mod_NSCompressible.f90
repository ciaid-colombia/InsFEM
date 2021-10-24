module Mod_NSCompressible
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_PhysicalProblem
   implicit none
   private
   public NSCompressibleProblem, zensi, ncsgs
   
   real(rp), parameter   :: &
     zensi = epsilon(0.0_rp)        ! zero
   
    integer(ip), parameter :: ncsgs = 3
  
   type, extends(PhysicalProblem) ::  NSCompressibleProblem
      
      ! Logical units
      integer(ip), allocatable   :: lun_force(:)
      integer(ip) ::&
             lun_trap                 ! Tracking of points

      character(150) ::&
             fil_trap               ! Tracking of points
      
      ! Physical problem
      real(rp), allocatable :: &
         densf(:,:),&               !Density field
         momen(:,:,:),&             !Momentum field
         energ(:,:),&               !Total energy field
         veloc(:,:,:),&             !Velocity field
         press(:,:),&               !Pressure field
         tempe(:,:)                 !Teperature field

      type(r2p), allocatable :: &
            cosgs(:),&                          !Continuity subgrid scales
            ensgs(:)                            !Energy subgrid scales

      type(r3p), allocatable :: &
            mosgs(:)                            !Momentum subgrid scales

      type(r1p), allocatable :: &
            coResGP(:),&                        !Continuity residual at Gauss point
            enResGP(:)                          !Energy residual at Gauss point

      type(r2p), allocatable :: &
            moResGP(:)                          !Momentum residual at Gauss point

      integer(ip) ::&
            kfl_visco,&                         ! Viscous term
            kfl_advec,&                         ! Existence of A_j d_j(U)
            kfl_confi,&                         ! Confined flow
            kfl_sourc,&                         ! Existence of momentum or energy sources
            kfl_react,&                         ! Existence of Reaction terms
            kfl_stabm,&                         ! Stabilization method: 1: VMS, 0 SUPG, -1 GLS
            kfl_jacgr                         ! VMS: 0: A·grad V, 1: div(AV)

      real(rp) ::&
            grnor    ,&                         ! Gravity norm
            gravi(3) ,&                         ! Gravity vector
            srce                                ! Heat sources
            
      ! Physical properties
      real(rp) ::&
             cphea,&                             !Specific heat constant pressure (Cp)
             cvhea,&                             !Specific heat constant volume(Cv)
             tcond,&                             !Thermal conduction (k)
             visco,&                             !Viscosity (mu)
             lawdep(10)!,&			 !Ideal gas eqn Parameters

      integer(ip) :: &
            lawde                               ! Eq of state 0=undefined 1=Ideal >2=other

      !Relative Thermodynamic Properties
      real(rp) ::&
             relpre,&                            ! Relative pressure     
             reltem                              ! Relative temperature


     type(r1p), allocatable :: &
            scdiffarray(:,:),&                     ! Gauss point shock capturing added diffusivities
            sgdiffarray(:,:)                       ! Gauss point subgrid added diffusivities

      ! Damping
     integer(ip) :: &
            ndamp                               ! Number of damping zones 0=undefined >0=exists

     real(rp)::&
            dampco(3,10),&                         ! Coeficients of damping along (X,Y,Z)
            dampxo(3,10),&                         ! Coordinates of initial vertex rectangular region
            dampxf(3,10),&                         ! Coordinates of final vextex rectangular region
            dampff(5,10)                        ! Unperturbed values of fields
      ! Radial Damping
     integer(ip) :: &
            rdamp,&                              ! Radial damping zone 0=undefined 1=exists
            rdampex                              ! Exponential of damping

     real(rp)::&
            rdampco,&                           ! Coeficient of damping along radius
            rdampxo(3),&                        ! Coordinates of origin
            rdampro,&                           ! Initial radius of damping region
            rdampff(5)                          ! Unperturbed values of fields

      ! Boundary conditions

      real(rp) :: &
            delta,&                              ! Distance to the wall
            epspe                                ! Penalization for confined flow

       integer(ip) ::&
             kfl_outfm,&                         ! Output of forces and moments
             kfl_nscnbc                          ! NSC Non-Reflecting Boundary Conditions
       logical :: kfl_openforce = .true.          ! Open file for forces
 
       real(rp), allocatable :: &
             forfl(:,:),&                         !Body forces field
             momfl(:,:),&                         !Body moments field
             hfxfl(:),&                           !Body heat flux field
             force(:,:),&                         !Body force
             formo(:,:)                           !Body moment

       real(rp) ::&
             adimf( 4),&                         ! Adimensionalization factor for forces
             adimm( 3),&                         ! Adimensionalization factor for moments
             origm( 3)                           ! Moment origin (application) point

       real(rp),allocatable  ::   &
             timad(:,:),&                          ! density  time average
             timav(:,:,:),&                        ! velocity time average
             timap(:,:),&                          ! pressure  time average
             timat(:,:),&                          ! temperature  time average
             rmsd(:),&                           ! R.M.S. (deviation) of density
             rmsv(:,:),&                         ! R.M.S. (deviation) of velocity 
             rmsp(:),&                           ! R.M.S. (deviation) of pressure 
             rmst(:),&                           ! R.M.S. (deviation) of temperature 
             turbi(:)                            ! Turbulence intensity = R.M.S(U)/|U| 
      
      ! Numerical treatment
      integer(ip) ::&
            kfl_trasg,&                         ! Tracking of subgrid scale
            kfl_tacsg,&                         ! Time accuracy of subscales
            kfl_nolsg,&                         ! Non-linearity of the subscales
            kfl_repro,&                         ! Stabilization based on residual projection
            kfl_shock,&                         ! Shock capturing off(0), residual(1), gradient orthogonal projection (2)
            nodpr,&                             ! Node on which pressure is prescribed
            kfl_sctyp,&                         ! Shock capturing type: isotropic(0), anisotropic(1)
            ErrorEstimatorTypeOfSubscales = 0   !1: scaled L_2 norm, 2: entropy function
           
             
      real(rp) ::&
            staco(4) ,&                         ! Stability constants
            prepr    ,&                         ! Prescribed pressure
            shock                           ! Shock capturing parameter
      real(rp), allocatable :: &
           grprj(:,:,:)                         !Gradient finite element projection

      real(rp), allocatable :: &
         repro(:,:)                             ! Residual Projection

      !FSI coupling
      real(rp),     allocatable ::&
         btraction(:,:)             ! For FSI coupling

      !Boundary VMass for NRBCS
      real(rp),     allocatable ::&
         bvmass(:)             ! For boundary non reflecting problem

      !Restart from incompressible flow
      integer(ip) :: kfl_nsirstar
      real(rp)    :: speed,molmass,densi
 
      !Adaptivity
      type(r2p), allocatable :: NodalArrayDataForRefinement(:)
contains

      !specific deferred procedures from PhysicalProblem
      procedure :: SetExmod          => nsc_SetExmod
      procedure :: SetNdofn          => nsc_SetNdofn
      procedure :: SetNdofbc         => nsc_SetNdofbc
      procedure :: SpecificReaphy    => nsc_reaphy
      procedure :: SpecificReanut    => nsc_reanut
      procedure :: SpecificReaous    => nsc_reaous
      procedure :: SpecificReaMPI    => nsc_reampi
      procedure :: SpecificReabcs => nsc_Reabcs
      procedure :: SpecificReadOnNodes => nsc_ReadOnNodes
      procedure :: SpecificReadOnBoundaries => nsc_ReadOnBoundaries
      procedure :: SpecificTurnof    => nsc_turnof
      procedure :: SpecificBegite    => nsc_begite
      procedure :: SpecificBegste    => nsc_begste
      procedure :: SpecificUpdbcs    => nsc_updbcs 
      procedure :: SpecificIniunk    => nsc_iniunk
      procedure :: SpecificEndite    => nsc_endite  
      procedure :: SpecificEndste    => nsc_endste
      procedure :: SpecificExaerr    => nsc_exaerr
      procedure :: SpecificOuterr    => nsc_outerr
      procedure :: SpecificSolite    => nsc_NULLSUBitask
      procedure :: SpecificTurnon    => nsc_NULLSUB
      procedure :: SpecificBounod    => nsc_NULLSUBitask
      procedure :: SpecificCrankNicolsonEndste => nsc_NULLSUB
      procedure :: SpecificRestart   => nsc_restar
      procedure :: Cvgunk            => nsc_NULLSUBitaskin
      procedure :: Elmope            => nsc_NULLSUB
      procedure :: Bouope            => nsc_NULLSUB
      procedure :: Memall            => nsc_memall
      procedure :: Getste            => nsc_getste  
      procedure :: Output            => nsc_output
      procedure :: InnerResiduals    => nsc_InnerResiduals
      procedure :: Ifconf            => nsc_ifconf
      procedure :: EndBouope      => nsc_EndBouope 
      procedure :: SpecificManufacturedBoundaryConditions => nsc_ManufacturedBoundaryConditions
      procedure :: CalculatePrimitivesLogic   => nsc_CalculatePrimitivesLogic
      procedure :: CalculatePrimitives   => nsc_CalculatePrimitives
      procedure :: CalculateConservatives   => nsc_CalculateConservatives
      procedure :: CalculateStatistics   => nsc_CalculateStatistics
      procedure :: FinalizeStats         => nsc_FinalizeStats
      procedure :: GetPhysicalParameters => nsc_GetPhysicalParameters
      procedure :: GetVelocityArray
      procedure :: GetTractionArray
      procedure :: SetVelocityArray
      procedure :: SetPressureArray
      procedure :: SetTemperatureArray
      procedure :: CalculateBoundaryVmass => nsc_CalculateBoundaryVmass

      procedure :: SpecificNSCompReaphy => nsc_SpecificNSCompReaphyRUNEND
      procedure :: SpecificNSCompReanut => nsc_SpecificNSCompReanutRUNEND
      procedure :: SpecificNSCompReaMPI => nsc_SpecificNSCompReaMPIRUNEND
      procedure :: SpecificNSCompReabcs => nsc_SpecificNSCompReabcsRUNEND
      procedure :: SpecificNSCompMemall => nsc_SpecificNSCompMemallRUNEND
      procedure :: SpecificNSCompTurnof => nsc_SpecificNSCompTurnofRUNEND
      procedure :: SpecificNSCompBegste => nsc_SpecificNSCompBegsteRUNEND
      procedure :: SpecificNSCompBegite => nsc_SpecificNSCompBegiteRUNEND
      procedure :: SpecificNSCompEndite => nsc_SpecificNSCompEnditeRUNEND
      procedure :: SpecificNSCompEndste => nsc_SpecificNSCompEndsteRUNEND
      procedure :: SpecificNSCompExaerr => nsc_SpecificNSCompExaerrRUNEND
      procedure :: SpecificNSCompRefine => nsc_SpecificNSCompRefineRUNEND
      procedure :: SpecificNSCompRefinementCriteria  => nsc_SpecificNSCompRefinementCriteriaRUNEND
      procedure :: SpecificSubscalesRefCriteria  => nsc_SpecificRefCriteriaRUNEND2
      procedure :: SpecificExactSolRefCriteria  => nsc_SpecificRefCriteriaRUNEND
      procedure :: SpecificNSCompRestar => nsc_SpecificNSCompRestarRUNEND

      !For Adaptive refinement
      procedure :: SpecificRefine => nsc_Refine
      procedure :: GetRefinementCriteria => nsc_GetRefinementCriteria
      procedure :: PreRefine  => nsc_Prerefine

   end type

   interface
     
      subroutine nsc_SetExmod(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a      

      end subroutine
      
      subroutine nsc_SetNdofn(a)
         import NSCompressibleProblem
         implicit none         
         class(NSCompressibleProblem) :: a      

      end subroutine
      
      subroutine nsc_SetNdofbc(a)
         import NSCompressibleProblem
         implicit none         
         class(NSCompressibleProblem) :: a      

      end subroutine
     
      subroutine nsc_reaphy(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none         
         integer(ip) :: itask
         class(NSCompressibleProblem) :: a      

      end subroutine
   
      subroutine nsc_reanut(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none         
         integer(ip) :: itask
         class(NSCompressibleProblem) :: a      

      end subroutine
   
      subroutine nsc_reaous(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none         
         integer(ip) :: itask
         class(NSCompressibleProblem) :: a      
      end subroutine
   
      subroutine nsc_reaMPI(a)
         use typre
         import NSCompressibleProblem
         implicit none         
         class(NSCompressibleProblem) :: a      

      end subroutine

      subroutine nsc_Reabcs(a,itask,kflag)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      
      end subroutine
      
      subroutine nsc_ReadOnNodes(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

      subroutine nsc_ReadOnBoundaries(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

        subroutine nsc_exaerr(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
     
      end subroutine 

      subroutine nsc_outerr(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
     
      end subroutine
      
      subroutine nsc_memall(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a      

      end subroutine
      
      subroutine nsc_getste(a,dtinv)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         real(rp) :: dtinv
    
      end subroutine
      
      subroutine nsc_begste(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a      
      end subroutine
      
      subroutine nsc_begite(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine
      
      subroutine nsc_endite(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip) :: itask
      
      end subroutine

      subroutine nsc_endste(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine nsc_output(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine nsc_iniunk(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine
      
      subroutine nsc_ifconf(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a 
      end subroutine
      
      subroutine nsc_updbcs(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      end subroutine      
      
      subroutine nsc_EndBouope(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine
      
      subroutine nsc_turnof(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

      subroutine nsc_restar(a,itask)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip), intent(in) :: itask
         
      end subroutine
      
      subroutine nsc_CalculatePrimitivesLogic(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

      subroutine nsc_CalculatePrimitives(a,ncomp)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip), intent(in) :: ncomp
      
      end subroutine

      subroutine nsc_CalculateStatistics(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

      subroutine nsc_FinalizeStats(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine
      
      subroutine nsc_GetPhysicalParameters(a,acvis,actco,accph,accvh)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         real(rp) :: acvis, actco, accph, accvh
         
      end subroutine

      subroutine nsc_ManufacturedBoundaryConditions(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      end subroutine

      subroutine nsc_GetRefinementCriteria(a,markel)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip) :: markel(*)
      end subroutine
      
      subroutine nsc_Refine(a,itask)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine nsc_CalculateConservatives(a,ncomp)
         use typre
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
         integer(ip), intent(in) :: ncomp
      
      end subroutine

      subroutine nsc_CalculateBoundaryVmass(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      
      end subroutine

      subroutine nsc_Prerefine(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      
      
      end subroutine

   end interface

contains

   subroutine nsc_InnerResiduals(a,drnsc,mrnsc,ernsc)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp) :: drnsc,mrnsc,ernsc
      real(rp) :: zensi = epsilon(0.0_rp)        ! zero
   
      integer(ip) :: ndime, npoinLocal
      
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPIHeterogeneous(2_ip,a%ndofn,1,npoinLocal,1,1,1,a%unkno,a%densf(1,1),drnsc,zensi,a%MPIcomm,a%MPIroot)     
      call vecresMPIHeterogeneous(2_ip,a%ndofn,ndime,npoinLocal,2,1,ndime,a%unkno,a%momen(1,1,1),mrnsc,zensi,a%MPIcomm,a%MPIroot)
      call vecresMPIHeterogeneous(2_ip,a%ndofn,1,npoinLocal,a%ndofn,1,1,a%unkno,a%energ(1,1),ernsc,zensi,a%MPIcomm,a%MPIroot)     
   
   end subroutine

   subroutine GetVelocityArray(a,veloc)
      use typre
      implicit none
      class(NSCompressibleProblem), target :: a
      real(rp), pointer :: veloc(:,:)
      
      veloc => a%veloc(:,:,1)
   end subroutine
      
   subroutine GetTractionArray(a,trac)
      use typre
      implicit none
      class(NSCompressibleProblem), target :: a
      real(rp), pointer :: trac(:,:)
      
      trac => a%btraction(:,:)
   end subroutine

   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp), target :: veloc(:,:)
      
      a%veloc(:,:,1) = veloc

   end subroutine

   subroutine SetPressureArray(a,press)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp), target :: press(:)
      
      a%press(:,1) = press

   end subroutine

   subroutine SetTemperatureArray(a,tempe)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp), target :: tempe(:)
     
      a%tempe(:,1) = tempe

   end subroutine

   subroutine nsc_SpecificNSCompReaphyRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Reaphy is pointing to the RUNEND')
   end subroutine

   subroutine nsc_SpecificNSCompReanutRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Reanut is pointing to the RUNEND')
   end subroutine

   subroutine nsc_SpecificNSCompReaMPIRUNEND(a)
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP ReaMPI is pointing to the RUNEND')
   end subroutine

   subroutine nsc_SpecificNSCompReabcsRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Reabcs is pointing to the RUNEND')
   end subroutine

   subroutine nsc_SpecificNSCompMemallRUNEND(a)
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP Memall is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompTurnofRUNEND(a)
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP Turnof is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompBegsteRUNEND(a)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP Begste is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompBegiteRUNEND(a)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP Begite is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompEnditeRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Endite is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompEndsteRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Endste is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompExaerrRUNEND(a)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
   
      call runend('NSCOMP Exaerr is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompRefineRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      character(6) :: itask
   
      call runend('NSCOMP Refine is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompRefinementCriteriaRUNEND(a,markel)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: markel(*)
   
      call runend('NSCOMP RefinementCriteria is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificRefCriteriaRUNEND2(a,error,TotalEstimatedError)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp) :: error(:), TotalEstimatedError
   
      call runend('NSCOMP RefCriteria is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificRefCriteriaRUNEND(a,error)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      real(rp) :: error(:)
   
      call runend('NSCOMP RefCriteria is pointing to the RUNEND')

   end subroutine

   subroutine nsc_SpecificNSCompRestarRUNEND(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
   
      call runend('NSCOMP Restar is pointing to the RUNEND')

   end subroutine

   subroutine nsc_NULLSUB(a)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
   
   end subroutine

   subroutine nsc_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip) :: itask
      
   end subroutine
   
   subroutine nsc_NULLSUBitaskin(a,itask)
      use typre
      implicit none
      class(NSCompressibleProblem) :: a
      integer(ip),intent(in) :: itask
      
   end subroutine

end module Mod_NSCompressible


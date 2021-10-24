module Mod_Solids
   use typre
   use Mod_PhysicalProblem
   implicit none
   private
   public SolidsProblem,SolidsProblem_Const, zesld

   real(rp), parameter   :: &
     zesld = epsilon(0.0_rp)        ! zero
    
   type, extends(PhysicalProblem) ::  SolidsProblem

       integer(ip) :: udofn = 0              !Displacement degrees of freedom
       ! Physical problem
      ! Physical problem
      real(rp), allocatable :: &
               disp(:,:,:)    ,&
              sigma(:,:,:)    ,&
              press(:,:)

       real(rp), allocatable ::&
               disp_cp(:,:),&                !Past values for coupling iteration
               veloc(:,:,:),&                !For use with newmark integrator
               accel(:,:,:),&
               btraction_nodal(:,:),&        !Traction in the solid boundary, used to check FSI convergence
               c_elastic(:,:)                !Elastic tensor

       real(rp), allocatable :: &
               stress(:,:),&                 !Stress calculated in endste
               strain(:,:)                   !Strain calculated in endste

       type(r1p), allocatable :: &
               pStrain(:)                    !Principal strains 

       type(r1p), allocatable :: &
               press_g(:)                    !Gauss press
       type(r2p), allocatable :: &
               sigma_g(:)                    !Gauss sigma

       type(r2p), allocatable :: &
               extForce(:)      ,&           !External force, calculated in begite
               stress_g(:)      ,&           !Gauss stress calculated in endste
               strain_g(:)      ,&           !Gauss strain calculated in endste
               printPStrain(:)  ,&           !Print principal strains calculated in endste for ALE coupling
               btraction(:)                  !Boundary traction

       real(rp) ::&
               force_factor = 0.0_rp,&       !Actual force percentage
               delta_force  = 0.0_rp,&       !Delta to be added to factor
               grnor        = 0.0_rp,&       !Gravity norm
               gravi(3)     = 0.0_rp,&       !Gravity vector
               traction(3)  = 0.0_rp         !Source term body traction

       !Material type (linear, nonlinear)
       character(5) :: sld_type   = 'NONE '
       !constitutive equation 
       character(5) :: sld_model  = 'NONE '
       character(5) :: sld_strain = 'NONE '

       !Linear elasticity
       real(rp) ::&
               densi   = 0.0_rp,&            !Density (rho)
               young   = 0.0_rp,&            !Young's modulus (E)
               poisson = 0.0_rp              !Poisson modulus (nu)

       !Non-Linear elasticity
       real(rp) ::&
               bulkMod  = 0.0_rp,&           !Bulk modulus K 
               lambda   = 0.0_rp,&           !lambda (lame param) 
               mu       = 0.0_rp             !mu     (Poisson's ratio)

       !Newmark constants
       real(rp) ::&
               omega = 0.0_rp   ,&           !Omega(usually 1/2)
               beta  = 0.0_rp   ,&           !Beta (usually 1/4)
               gamma = 0.0_rp                !Gamma (usually 1/2)

       ! Logical units
       integer(ip)    :: lun_trap = 0           !Used for point tracking
       character(150) :: fil_trap

       ! Output & postprocess
       integer(ip) ::&
               nbody     = 0    ,&                !# of bodies
               kfl_outfm = 0                 !Output of forces and moments

       ! Numerical treatment
       integer(ip) ::&
               kfl_local = 0    ,&           !Local system of reference
               kfl_traction = 0              !Distributed tractions exist ?

       !For reading boundary conditions, to be used by MPIroot
       integer(ip), allocatable :: &
               kfl_bours(:),&                !Boundary reference system
               kfl_fixrs(:)                  !Reference system for the BV

       !For point wise sources, e.g applied nodal force
       logical    :: kfl_sourcesPresent= .false.  !Used to know if there are PointSources
       integer(ip):: kfl_constPointForce = 9      !1 = constant, 0 = impulse        
       real(rp)   :: kfl_PointForceTime  = 0.0_rp !Duration of force

       !For solution with EPS eigen value solver
       logical::kfl_eigenter = .true.

       real(rp),    allocatable :: pwSource(:,:)
       integer(ip), allocatable :: pwSourceId(:)
       integer(ip) :: pwSourceSize = 0, pwSourceDim = 0

       !FSI used variables
       integer(ip):: kfl_imposedTraction = 0_ip    !Output of imposed tractions
       logical    :: kfl_foldedElements  = .false. !Check if folded elements

       !Relaxation factors 
       real(rp) :: sld_dispRelax         =1.0_rp   !Relaxation for displacements

       !Tractions used with FSI
       real(rp), pointer     :: trac(:,:) => NULL()

       !Used when ALE owns the solid module
       logical :: kfl_modProperties = .false. !Used when coupled with ALE to deform mesh smoothly
       logical :: kfl_saveStrainsALE= .false. !Used when coupled with ALE to deform mesh smoothly
       logical :: kfl_printPrincipalStresses= .false. !Used when coupled with ALE to deform mesh smoothly

       !Other flags
       logical :: kfl_FSI_restart = .false.   !Used in sld_begste for fsi restart

       logical :: kfl_printStrain  = .false.   !Used to know how to print stress and strain
       logical :: kfl_printStress  = .false.   !Used to know how to print stress and strain

       logical :: kfl_NodalStress  = .false.   !Used to know how to print stress and strain
       logical :: kfl_NodalStrain  = .false.   !Used to know how to print stress and strain

       logical :: kfl_GaussStress  = .false.   !Used to know how to print stress and strain
       logical :: kfl_GaussStrain  = .false.   !Used to know how to print stress and strain

       logical :: kfl_printSigma      = .false.   !Used to know how to print stress and strain
       logical :: kfl_printNodalSigma = .false.   !Used to know how to print stress and strain
       logical :: kfl_printGaussSigma = .false.   !Used to know how to print stress and strain

       logical :: kfl_printPress      = .false.   !Used to know how to print stress and strain
       logical :: kfl_printNodalPress = .false.   !Used to know how to print stress and strain
       logical :: kfl_printGaussPress = .false.   !Used to know how to print stress and strain




   contains

      !specific deferred procedures from PhysicalProblem
      !Degrees of freedom
      procedure :: SetExmod                 => sld_SetExmod
      procedure :: SetNdofn                 => sld_SetNdofn
      procedure :: SetNdofbc                => sld_SetNdofbc
      !Readers
      procedure :: SpecificReaphy           => sld_reaphy
      procedure :: SpecificReanut           => sld_reanut
      procedure :: SolidSpecificReanut      => sld_NULLSUBitask
      procedure :: SpecificReaous           => sld_reaous
      procedure :: SolidSpecificReaous      => sld_NULLSUBitask
      procedure :: SpecificReaMPI           => sld_reampi
      procedure :: SolidSpecificReampi      => sld_NULLSUB
      procedure :: SpecificReabcs           => sld_Reabcs
      procedure :: SpecificReadOnNodes      => sld_ReadOnNodes
      procedure :: SpecificReadSource       => sld_ReadSource
      procedure :: SpecificReadOnBoundaries => sld_ReadOnBoundaries
      !Preprocess and postprocess of linear system
      procedure :: SpecificRestart          => sld_restar
      procedure :: SolidSpecificRestart     => sld_NULLSUBrestart
      procedure :: SpecificBounod           => sld_bounod
      procedure :: SpecificSolite           => sld_NULLSUBitask
      procedure :: SpecificTurnon           => sld_turnonGeneral
      procedure :: SolidSpecificTurnon      => sld_turnonSpecific
      procedure :: SpecificTurnof           => sld_turnofGeneral
      procedure :: SolidSpecificTurnof      => sld_turnofSpecific 
      procedure :: SpecificBegite           => sld_begite
      procedure :: SpecificBegste           => sld_begste
      procedure :: SpecificUpdbcs           => sld_NULLSUB
      procedure :: SpecificIniunk           => sld_iniunk
      procedure :: SpecificEndste           => sld_endste
      procedure :: SpecificEndite           => sld_endite
      procedure :: SpecificRefine           => sld_refine
      procedure :: SolidSpecificRefine      => sld_NULLSUBitaskc
      procedure :: SolidSpecificRefineGauss => sld_refineGauss
      procedure :: SpecificCrankNicolsonEndste => sld_NULLSUB
      procedure :: UpdateDynamicComponents  => sld_updateDynamicComp
      !Elemental procedures and solution
      procedure :: Memall                   => sld_memallGeneral
      procedure :: SolidSpecificMemall      => sld_memallSpecific
      procedure :: InnerResiduals           => sld_InnerResiduals
      procedure :: Getste                   => sld_NULLSUBdtinv
      procedure :: Cvgunk                   => sld_cvgunk
      procedure :: preElmope                => sld_preElmope
      procedure :: Elmope                   => sld_Elmope
      procedure :: EndElmope                => sld_EndElmope
      procedure :: Bouope                   => sld_bouop0
      procedure :: ModifyPointForces        => sld_NullsubPointForcesIn
      procedure :: ModifyBouopeRHS          => sld_NullsubRhsIn
      procedure :: EndBouope                => sld_EndBouope
      !Output and error output
      procedure :: Output                   => sld_output
      procedure :: SolidSpecificOutput      => sld_NULLSUBitask
      procedure :: PointTracking            => sld_outtpo
      procedure :: SpecificOuterr           => sld_outerr
      procedure :: SpecificExaerr           => sld_NULLSUB
      !own procedures of the class
      procedure :: GetPhysicalParameters    => sld_GetPhysicalParameters
      procedure :: sld_calcVelandAccel
      !Boundary operations 
      !Setters
      procedure :: SetExternalTractionArray
      procedure :: SetAitkenParam
      !Getters
      procedure :: GetDisplacementArray
      procedure :: GetVelocityArray
      procedure :: GetInternalTractionArray
      procedure :: GetFoldedElements     => sld_GetFoldedElements
      !ROM
      procedure :: SpecificGetUnkno      => sld_GetUnkno
      procedure :: GetRefinementCriteria => sld_GetRefinementCriteria
      !Assembly
      procedure :: GetMatrixOrganization => sld_GetMatrixOrganization
      
   end type

   interface
     
      subroutine sld_SetExmod(a)
         import 
         implicit none         
         class(SolidsProblem) :: a      
      end subroutine
      
      subroutine sld_SetNdofn(a)
         import 
         implicit none         
         class(SolidsProblem) :: a      
      end subroutine
      
      subroutine sld_SetNdofbc(a)
         import 
         implicit none         
         class(SolidsProblem) :: a      
      end subroutine
     
      subroutine sld_reaphy(a,itask)
         use typre
         import 
         implicit none         
         integer(ip) :: itask
         class(SolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reanut(a,itask)
         use typre
         import 
         implicit none         
         integer(ip) :: itask
         class(SolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reaous(a,itask)
         use typre
         import SolidsProblem
         implicit none         
         integer(ip) :: itask
         class(SolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reaMPI(a)
         use typre
         import 
         implicit none         
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_Reabcs(a,itask,kflag)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine
 
      subroutine sld_ReadOnNodes(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine

      subroutine sld_ReadSource(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine
      
      subroutine sld_ReadOnBoundaries(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine

      subroutine sld_bounod(a,itask)
         use typre
         import 
         implicit none      
         integer(ip) :: itask
         class(SolidsProblem) :: a      
      end subroutine
   
      subroutine sld_outerr(a)
         import
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_memallGeneral(a)
         import
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_memallSpecific(a)
         import
         implicit none
         class(SolidsProblem) :: a      
      end subroutine
      
      subroutine sld_begste(a)
         import
         implicit none
         class(SolidsProblem) :: a      
      end subroutine
      
      subroutine sld_begite(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_refine(a,itask)
         import 
         implicit none
         class(SolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sld_refineGauss(a,itask)
         import 
         implicit none
         class(SolidsProblem) :: a
         character(6) :: itask
      end subroutine
      
      subroutine sld_endite(a,itask)
         use typre
         import
         implicit none
         class(SolidsProblem) :: a
         integer(ip) :: itask     
      end subroutine
      
      subroutine sld_cvgunk(a,itask)
         use typre
         import
         implicit none
         class(SolidsProblem) :: a
         integer(ip), intent(in) :: itask   
      end subroutine
      
      subroutine sld_InnerResiduals(a,resld,srsld,prsld)
          use typre
          import
          implicit none
          class(SolidsProblem) :: a
          real(rp) :: resld
          real(rp), optional :: srsld,prsld
      end subroutine

      subroutine sld_outtpo(a)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine
      
      subroutine sld_output(a,itask)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine
      
      subroutine sld_endste(a,itask)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine
            
      subroutine sld_iniunk(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine 

      subroutine sld_preElmope(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_Elmope(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine
      
      subroutine sld_bouop0(a)
         import 
         implicit none
         class(SolidsProblem) :: a     
      end subroutine
      
      subroutine sld_turnonGeneral(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_turnonSpecific(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_turnofGeneral(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_turnofSpecific(a)
         import 
         implicit none
         class(SolidsProblem) :: a      
      end subroutine

      subroutine sld_GetPhysicalParameters(a,J,sz,densi,c_elas,ielem)
         use typre
         import 
         implicit none
         class(SolidsProblem), target :: a
         integer(ip) :: ndime,sz
         real(rp) :: J,densi
         real(rp)    :: c_elas(sz,sz)
         integer(ip), optional :: ielem
      end subroutine
      
      subroutine sld_calcVelandAccel(a,dtinv)
         use typre
         import 
         implicit none
         class(SolidsProblem), target :: a
         real(rp) :: dtinv
      end subroutine

      subroutine sld_EndElmope(a,task)
         import 
         implicit none
         class(SolidsProblem) :: a
         character(6) :: task
      end subroutine
      
      subroutine sld_EndBouope(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine

      subroutine sld_updateDynamicComp(a)
         import 
         implicit none
         class(SolidsProblem) :: a
      end subroutine

      subroutine sld_restar(a,itask)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
         integer(ip), intent(in) :: itask       
         !Restart calculations
      end subroutine

      subroutine sld_GetRefinementCriteria(a,markel)
         use typre
         import 
         implicit none
         class(SolidsProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

   end interface

    
  interface SolidsProblem_Const
      procedure constructor
  end interface SolidsProblem_Const

  contains

      function constructor()
          class(SolidsProblem), pointer :: constructor

          allocate(constructor)

      end function constructor

   subroutine sld_NULLSUB(a)
      use typre
      implicit none
      class(SolidsProblem) :: a
   
   end subroutine

   subroutine sld_NullsubPointForcesIn(a,ndof,pfsize,pf)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip), intent(in) :: ndof,pfsize
      real(rp), intent(inout) :: pf(ndof,pfsize)
   
   end subroutine

   subroutine sld_NullsubRhsIn(a,ndof,mnode,rhs)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip), intent(in) :: ndof,mnode
      real(rp), intent(inout) :: rhs(ndof,mnode)
   
   end subroutine

   subroutine sld_NULLSUBrestart(a,itask)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine

   subroutine sld_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip) :: itask
      
   end subroutine
   
   subroutine sld_NULLSUBitaskc(a,itask)
      use typre
      implicit none
      class(SolidsProblem) :: a
      character(6) :: itask
      
   end subroutine

   subroutine sld_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine
   
    subroutine sld_NULLSUBitaskkflag(a,itask,kflag)
      use typre
      implicit none
      class(SolidsProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
      
   end subroutine

    subroutine sld_NULLSUBdtinv(a,dtinv)
      use typre
      implicit none
      class(SolidsProblem) :: a
      real(rp) :: dtinv
      
   end subroutine

   subroutine GetDisplacementArray(a,disp)
      use typre
      implicit none
      class(SolidsProblem), target :: a
      real(rp), pointer :: disp(:,:)
      
      disp => a%disp(:,:,1)
   end subroutine

   subroutine GetVelocityArray(a,veloc)
      use typre
      implicit none
      class(SolidsProblem), target :: a
      real(rp), pointer :: veloc(:,:)
      
      veloc => a%veloc(:,:,1)
   end subroutine

   subroutine SetExternalTractionArray(a,trac)
      use typre
      implicit none
      class(SolidsProblem) :: a
      real(rp), target :: trac(:,:)

      a%trac => trac

   end subroutine

   subroutine GetInternalTractionArray(a,trac)
      use typre
      implicit none
      class(SolidsProblem), target :: a
      real(rp), pointer :: trac(:,:)

      trac => a%btraction_nodal

   end subroutine

   subroutine SetAitkenParam(a,aitken)
      use typre
      implicit none
      class(SolidsProblem) :: a
      real(rp), target :: aitken

      a%sld_dispRelax = aitken

   end subroutine
   
   subroutine sld_GetFoldedElements(a,kfl_folded)
      use typre
      implicit none
      class(SolidsProblem) :: a
      logical::kfl_folded

      kfl_folded = a%kfl_foldedElements
            
   end subroutine

   subroutine sld_GetMatrixOrganization(a,u1,uf,s1,sf,p1,bcstart)
      use typre
      implicit none
      class(SolidsProblem)       :: a
      integer(ip), intent(inout) :: u1,uf,s1,sf,p1,bcstart
      integer(ip)                :: nd,tn

      !Update unknowns
      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      s1 = 0
      sf = 0
      u1 = 1
      uf = nd
      p1 = 0
      bcstart = 0

   end subroutine

   subroutine sld_GetUnkno(a,unkno)
      use typre
      implicit none
      class(SolidsProblem) :: a
      real(rp)             :: unkno(:,:)
      integer(ip)          :: u1,uf,s1,sf,p1
      integer(ip)          :: nd,tn

      !Update unknowns
      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      s1 = 1
      sf = tn
      u1 = sf+1
      uf = sf+nd
      p1 = uf+1
           
      unkno(s1:sf,:) = a%sigma(:,:,1)
      unkno(u1:uf,:) = a%disp(:,:,1)
      unkno(p1,:)    = a%press(:,1)

   end subroutine

end module Mod_Solids




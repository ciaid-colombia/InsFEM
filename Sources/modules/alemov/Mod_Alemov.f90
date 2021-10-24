module Mod_Alemov
   use typre   
   use Mod_PhysicalProblem
   implicit none
   private
   public AlemovProblem,AlemovProblem_Const
   
   type, extends(PhysicalProblem) :: AlemovProblem
	real(rp), allocatable :: Displacement(:,:,:)
	real(rp), allocatable :: Velocity(:,:,:)
	real(rp), allocatable :: restar_displ(:,:,:)
   real(rp), allocatable :: bdisp(:,:,:)
   real(rp), allocatable :: bres(:,:,:)
   real(rp), pointer     :: NSveloc(:,:) => NULL()
   real(rp), pointer     :: sldDisp(:,:) => NULL()

   real(rp) :: sldale_dispRelax = 1.0_rp              !Relaxation param for mesh disp
   real(rp) :: sldale_dispRelax_max = 1.0_rp          !Relaxation param for mesh disp

   integer(ip), allocatable :: kfl_fixno0(:,:)
   
   integer(ip) :: kfl_isFixedMesh = 0
   integer(ip) :: kfl_RemeshingCriteria = 0
   integer(ip) :: DoRemesh = 0

   logical :: kfl_doAitken  = .false.       !Used to know Aitken relaxation should be used
   logical :: kfl_cvgcp     = .true.       !Write header
   logical :: kfl_cvg       = .true.       !Write header

   integer(ip) :: lun_trap               !Used for point tracking
   character(150) :: fil_trap


contains
   procedure :: SetExmod          => ale_SetExmod
   procedure :: SetNdofn          => ale_SetNdofn
   procedure :: SetNdofbc         => ale_SetNdofbc
   procedure :: Elmope            => ale_NULLSUB
   procedure :: Bouope            => ale_NULLSUB
   procedure :: Cvgunk            => ale_cvgunk
   procedure :: Getste            => ale_NULLSUB2
   procedure :: Memall            => ale_memall
   procedure :: SpecificRestart   => ale_restar
   procedure :: Output            => ale_output
   procedure :: Solite            => ale_NULLSUB
   procedure :: BuildLinearSystem => ale_NULLSUB
   procedure :: SpecificBegste    => ale_solite
   procedure :: SpecificBegite    => ale_begite
   procedure :: SpecificBounod    => ale_NULLSUB1
   procedure :: SpecificEndste    => ale_NULLSUB1
   procedure :: SpecificIniunk    => ale_iniunk
   procedure :: SpecificExaerr    => ale_NULLSUB
   procedure :: SpecificEndite    => ale_NULLSUB1
   procedure :: SpecificReaMPI    => ale_reaMPI
   procedure :: SpecificReanut    => ale_reanut
   procedure :: SpecificSolite    => ale_NULLSUB1
   procedure :: SpecificReaphy    => ale_NULLSUB1
   procedure :: SpecificTurnof    => ale_turnof
   procedure :: deferredTurnof    => ale_NULLSUB
   procedure :: SpecificUpdbcs    => ale_updbcs
   procedure :: SpecificTurnon    => ale_turnon
   procedure :: deferredTurnon    => ale_NULLSUB
   procedure :: SpecificOuterr    => ale_NULLSUB
   procedure :: SpecificReaous    => ale_reaous
   procedure :: SpecificReabcs    => ale_reabcs
   procedure :: ALESpecificRefine => ale_NULLSUB4
   procedure :: SpecificReadOnNodes         => ale_NULLSUB
   procedure :: SpecificReadOnBoundaries    => ale_NULLSUB
   procedure :: SpecificReadOnFunctions     => ale_readOnFunctions
   procedure :: SpecificProjectArraysOUM    => ale_ProjectArraysOntoUndeformedMesh
   procedure :: SpecificAdvectArraysOUM     => ale_AdvectArraysOntoUndeformedMesh
   procedure :: SpecificCrankNicolsonEndste => ale_NULLSUB
   procedure :: GetDoRemesh        => ale_GetDoRemesh
   procedure :: SpecificRefine     => ale_refine
   procedure :: ParticularRestart  => ale_NULLSUB1
   procedure :: GetVelocityArray   => ale_GetVelocityArray
   procedure :: SetVelocityArray
   procedure :: SetDisplacementArray
   procedure :: ale_outtpo
   
   end type

interface

   subroutine ale_SetExmod(a)
      use typre
      import AlemovProblem
      implicit none
   
      class(AlemovProblem) :: a
   
   end subroutine

   subroutine ale_SetNdofn(a)
      use typre
      import AlemovProblem
      implicit none
   
      class(AlemovProblem) :: a

   end subroutine

   subroutine ale_SetNdofbc(a)
      use typre
      import AlemovProblem
      implicit none
   
      class(AlemovProblem) :: a

   end subroutine

   subroutine ale_memall(a)
      use typre
      import AlemovProblem
      implicit none
   
      class(AlemovProblem) :: a

   end subroutine
   
   subroutine ale_begite(a)
      use typre
      import AlemovProblem
      implicit none
   
      class(AlemovProblem) :: a

   end subroutine

   subroutine ale_cvgunk (a,itask)
      use typre
      import AlemovProblem
      implicit none
      integer(ip), intent(in) :: itask
      class(AlemovProblem) :: a
   end subroutine

   subroutine ale_output(a,itask)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: itask    
   end subroutine

   subroutine ale_iniunk(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
 
   subroutine ale_reanut(a,itask)
      use typre
      import AlemovProblem
      implicit none 
      integer(ip) :: itask
      class(AlemovProblem) :: a
   end subroutine
   
   subroutine ale_solite(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
   
   subroutine ale_readonfunctions(a,ifunc)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: ifunc
   end subroutine

   subroutine ale_turnof(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine

   subroutine ale_turnon(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine

   subroutine ale_reaous(a,itask)
      use typre
      import AlemovProblem
      implicit none
          
      integer(ip) :: itask
      class(AlemovProblem) :: a
   end subroutine

   subroutine ale_elmope(a,currentbvess)  
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer :: currentbvess
   end subroutine
   
   subroutine ale_ReadOnNodes(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
   
   subroutine ale_BuildLinearSystem(a,currentbvess)  
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer :: currentbvess 
   end subroutine 

   subroutine ale_updbcs(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
   
   subroutine ale_reaMPI(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
   
    subroutine ale_reabcs(a,itask,kflag)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
   end subroutine
   
   subroutine ale_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
      use typre
      use Mod_Mesh
      use Mod_MeshInterpolator
      import 
      implicit none
      class(AlemovProblem) :: a
      type(Interpolator) :: Interp
      integer(ip) :: itask
   end subroutine
   
   subroutine ale_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
      use typre
      use Mod_Mesh
      use Mod_Advector
      import 
      implicit none
      class(AlemovProblem) :: a
      type(Advector) :: Advect
      integer(ip) :: itask
   end subroutine
   
   subroutine ale_refine(a,itask)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      character(6) :: itask
   end subroutine

   subroutine ale_restar(a,itask)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
      integer(ip), intent(in) :: itask
   end subroutine

   subroutine ale_outtpo(a)
      use typre
      import AlemovProblem
      implicit none
      class(AlemovProblem) :: a
   end subroutine
      
end interface

  interface AlemovProblem_Const
      procedure constructor
  end interface AlemovProblem_Const

contains

    function constructor()
        class(AlemovProblem), pointer :: constructor

        allocate(constructor)

    end function constructor

   subroutine ale_NULLSUB(a)
      implicit none
      class(AlemovProblem) :: a
   end subroutine

   subroutine ale_NULLSUB1(a,itask)
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: itask
   end subroutine

   subroutine ale_NULLSUB2(a,dtinv)
      implicit none
      class(AlemovProblem) :: a
      real(rp) :: dtinv
   end subroutine

   subroutine ale_NULLSUB3(a,itask,kflag)
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
   end subroutine
   
   subroutine ale_NULLSUB4(a,itask)
      implicit none
      class(AlemovProblem) :: a
      character(6) :: itask
   end subroutine
   
   subroutine ale_GetDoRemesh(a,DoRemesh)
      implicit none
      class(AlemovProblem) :: a
      integer(ip) :: DoRemesh
      
      DoRemesh = a%DoRemesh
      
   end subroutine

   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(AlemovProblem) :: a
      real(rp), target :: veloc(:,:)
      
      a%NSveloc => veloc
   end subroutine
   
   subroutine ale_GetVelocityArray(a,veloc)
      implicit none
      class(AlemovProblem), target :: a
      real(rp), pointer :: veloc(:,:,:)
      
      veloc => a%Velocity
   end subroutine

   subroutine SetDisplacementArray(a,displ)
      use typre
      implicit none
      class(AlemovProblem) :: a
      real(rp), target :: displ(:,:)
      
      a%sldDisp=>displ 
   end subroutine

end module  

module Mod_LevelSet
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_CutMesh
   use Mod_PhysicalProblem
   use Mod_HeightGauge
   implicit none
   private
   public LevelSetProblem,LevelSetProblem_Const,zelev

   real(rp), parameter :: &
       zelev = epsilon(1.0_rp)

   type, extends(PhysicalProblem) :: LevelSetProblem

      !Numerical treatment
      real(rp)              :: staco(3)              ! Stability constants
      !Physical Problem
      real(rp), allocatable :: level(:,:)

      integer(ip)           :: nReinit=1
      integer(ip)           :: NstepsReinitLevel = 50

      integer(ip)           :: changeTimeIntegrator=0


      ! Lists
      integer(ip), allocatable :: &
         ElementListByLayers(:)                   ! List of level set elements

      !Problem data and flags
      integer(ip)           :: kfl_advec,&       ! Existence of (u.grad)T
                               kfl_ExactLevel,&  !Level set is defined by an analytical field
                               kfl_ReinitLevel,& !Type of reinitialization algorithm
                               kfl_CutCriteria,&
                               kfl_EnrichElem, &
                               kfl_MassCorrection = 0

      !Velocity pointer, obtained from another object (NSI, for instance)
      real(rp), pointer     :: veloc(:,:) => NULL()

      !For Smooth Gradient computation
      integer(ip) :: kfl_SmoothGradient
      real(rp), allocatable :: SmoothGradient(:,:)

      !For ADaptiveRefinement
      integer(ip) :: OutLayers, GeneralRefinementLevels, InterfaceRefinementLevels

      !For Fixed Mesh ALE
      integer(ip) :: kfl_ForceEulerianAdvection = 0

      !Do not redistance
      integer(ip) :: kfl_InitialRedistance = 1

      !For mass loss correction
      integer(ip) :: ipass = 0
      real(rp) :: ReferenceVolume = 0.0_rp

      !For Measuring Height
      integer(ip) :: nHeightGauges = 0
      type(HeightGauge), allocatable :: HeightGauges(:)
      integer(ip) :: lun_outHeight

contains

      procedure :: SetExmod          => lev_SetExmod
      procedure :: SetNdofn          => lev_SetNdofn
      procedure :: SetNdofbc         => lev_SetNdofbc
      procedure :: SpecificReaphy    => lev_reaphy
      procedure :: SpecificReanut    => lev_reanut
      procedure :: SpecificReaous    => lev_reaous
      procedure :: SpecificReaMPI    => lev_reampi
      procedure :: SpecificBounod    => lev_bounod
      procedure :: SpecificReabcs => lev_Reabcs
      procedure :: SpecificReadOnNodes => lev_ReadOnNodes

      procedure :: SpecificOuterr    => lev_outerr
      procedure :: Memall            => lev_memall
      procedure :: SpecificIniunk    => lev_iniunk

      procedure :: Getste            => lev_getste
      procedure :: SpecificBegste    => lev_begste
      procedure :: SpecificBegite    => lev_begite
      procedure :: Solite            => lev_Solite
      procedure :: SpecificSolite    => lev_NULLSUBitask
      procedure :: SpecificEndite    => lev_endite
      procedure :: EnditeElmope      => lev_NULLSUB
      procedure :: SpecificEndste    => lev_endste
      procedure :: SpecificCrankNicolsonEndste => lev_CrankNicolsonEndste
      Procedure :: SpecificUpdbcs    => lev_NULLSUB
      Procedure :: SpecificExaerr    => lev_NULLSUB
      procedure :: SpecificRefine    => lev_refine
      procedure :: GetRefinementCriteria => lev_GetRefinementCriteria

      procedure :: Cvgunk            => lev_cvgunk
      procedure :: Output            => lev_output
      procedure :: Elmope            => lev_elmope
      procedure :: Bouope            => lev_NULLSUB

      procedure :: SpecificTurnon    => lev_NULLSUB
      procedure :: SpecificTurnof    => lev_turnof

      procedure :: Restart           => lev_restar
      procedure :: SpecificRestart   => lev_NULLSUBitaskintentin

      procedure :: GetInitialConditions => lev_GetInitialConditions
      procedure :: SetInitialConditions => lev_SetInitialConditions

      procedure :: SetVelocityArray
      procedure :: SetEnrichedFlag
      procedure :: GetLevelSetArray

      procedure :: SetUnknoA1 => lev_SetUnkno

      procedure :: GetSmoothGradient => lev_getSmoothGradient
      procedure :: SpecificProjectArraysOUM => lev_ProjectArraysOntoUndeformedMesh
      procedure :: SpecificAdvectArraysOUM => lev_AdvectArraysOntoUndeformedMesh

   end type

   interface

      subroutine lev_SetExmod(a)
         import LevelSetProblem
         implicit none

         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_SetNdofn(a)
         import LevelSetProblem
         implicit none

         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_SetNdofbc(a)
         import LevelSetProblem
         implicit none

         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_reaphy(a,itask)
         use typre
         import LevelSetProblem
         implicit none

         integer(ip) :: itask
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_reanut(a,itask)
         use typre
         import LevelSetProblem
         implicit none

         integer(ip) :: itask
         class(LevelSetProblem) :: a

      end subroutine


      subroutine lev_reaous(a,itask)
         use typre
         import LevelSetProblem
         implicit none

         integer(ip) :: itask
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_reampi(a)
         use typre
         import LevelSetProblem
         implicit none

         class(LevelSetProblem) :: a

      end subroutine


      subroutine lev_bounod(a,itask)
         use typre
         import LevelSetProblem
         implicit none

         integer(ip) :: itask
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_outerr(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_memall(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_getste(a,dtinv)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         real(rp) :: dtinv

      end subroutine

      subroutine lev_begste(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_begite(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_endite(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine lev_cvgunk(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip), intent(in) :: itask

      end subroutine

      subroutine lev_output(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: itask


      end subroutine

      subroutine lev_endste(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine lev_iniunk(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_updbcs(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine


      subroutine lev_solite(a)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: itask

      end subroutine

      subroutine lev_elmope(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_turnof(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_CrankNicolsonEndste(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_restar(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip), intent(in) :: itask

         !Restart calculations
      end subroutine


      subroutine lev_Reabcs(a,itask,kflag)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag

      end subroutine

      subroutine lev_ReadOnNodes(a)
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_refine(a,itask)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine lev_GetRefinementCriteria(a,markel)
         use typre
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine lev_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
         use typre
         use Mod_MeshInterpolator
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         type(Interpolator) :: Interp
         integer(ip) :: itask
      end subroutine
      
      subroutine lev_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
         use typre
         use Mod_Advector
         import LevelSetProblem
         implicit none
         class(LevelSetProblem) :: a
         type(Advector) :: Advect
         integer(ip) :: itask
      end subroutine

   end interface

  interface LevelSetProblem_Const
      procedure constructor
  end interface LevelSetProblem_Const

contains

function constructor()
    class(LevelSetProblem), pointer :: constructor

    allocate(constructor)

end function constructor

   subroutine lev_SetUnkno(a,unkno,idofn)
      implicit none
      class(LevelSetProblem), target :: a
      real(rp) :: unkno(:)

      integer(ip) :: idofn,ipoin,npoin

      call a%Mesh%GetNpoin(npoin)

      !First we copy
      do ipoin = 1,npoin
         a%level(ipoin,1) = unkno(ipoin)
      enddo

      if (a%kfl_timei == 0) then
         a%level(:,2) = a%level(:,1)
         a%level(:,3) = a%level(:,1)
      endif
      a%unkno(1,:) = a%level(:,1)

   end subroutine

   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      real(rp), target :: veloc(:,:)

      a%veloc => veloc

   end subroutine

   subroutine SetEnrichedFlag(a,kfl_EnrichElem)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      integer(ip), target :: kfl_EnrichElem

      a%kfl_EnrichElem=kfl_EnrichElem

   end subroutine

   subroutine GetLevelSetArray(a,level)
      use typre
      implicit none
      class(LevelSetProblem), target :: a
      real(rp), pointer :: level(:)

      level => a%level(:,1)
   end subroutine

   subroutine lev_NULLSUB(a)
      implicit none
      class(LevelSetProblem) :: a

   end subroutine

   subroutine lev_NULLSUBitask(a,itask)
      implicit none
      class(LevelSetProblem) :: a
      integer(ip) :: itask
   end subroutine

   subroutine lev_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      integer(ip), intent(in) :: itask

   end subroutine

    subroutine lev_GetInitialConditions(a,unkno)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      real(rp) :: unkno(:,:,:)

      integer(ip) :: npoin,ndime,ispos,ipoin,itwost

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)

      do ipoin = 1,npoin
         unkno(1,ipoin,1) = a%level(ipoin,3)
      enddo

      !Higher order integration schemes
      if (a%ncomp > 3) then
         do itwost = 1,a%ncomp-3
            do ipoin = 1,npoin
               unkno(1,ipoin,itwost+1) = a%level(ipoin,3+itwost)
            enddo
         enddo
      endif

   end subroutine

   subroutine lev_SetInitialConditions(a,unkno)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      real(rp) :: unkno(:,:,:)

      integer(ip) :: npoin,ndime,ispos,ipoin,itwost

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)

      do ipoin = 1,npoin
         a%level(ipoin,3) = unkno(1,ipoin,1)
      enddo

      a%level(:,1) = a%level(:,3)
      a%level(:,2) = a%level(:,2)

      !Higher order integration schemes
      if (a%ncomp > 3) then
         do itwost = 1,a%ncomp -3
            do ipoin = 1,npoin
               a%level(ipoin,3+itwost) = unkno(1,ipoin,itwost+1)
            enddo
         enddo
      endif

   end subroutine

   subroutine lev_getSmoothGradient(a,SmoothGradient)
      use typre
      implicit none
      class(LevelSetProblem) :: a
      real(rp), pointer :: SmoothGradient(:,:)

      SmoothGradient => SmoothGradient

   end subroutine

end module

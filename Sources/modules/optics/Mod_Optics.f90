module Mod_Optics
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Timer
   use Mod_Mesh
   use Mod_PhysicalProblem
   implicit none
   private
   public OpticsProblem,OpticsProblem_Const,LightBeam

   type LightBeam 
      real(rp) :: origin(3)
      real(rp) :: direct(3)
      real(rp) :: length
   end type
   
   type, extends(PhysicalProblem) :: OpticsProblem
   
      integer(ip) :: kfl_dissi         !Dissipation scheme (Smago, etc)
      real(rp)    :: turbu(1)          !Turbulence parameters
      
      real(rp)    :: a_Obukhov
      
      real(rp), allocatable :: cn2(:)  !Cn2 coefficient 
      real(rp), allocatable :: cn2u53(:) !cn2u53 coefficient
      real(rp), allocatable :: ct2(:)  !Ct2 coefficient
      
      integer(ip)                   :: nbeams
      type(LightBeam), allocatable :: beams(:)
      
      real(rp) :: lambda
      
            
      integer(ip) :: kfl_ExternalVelocity, kfl_ExternalPressure,kfl_ExternalTemperature      
      integer(ip) :: kfl_velocity,kfl_pressure, kfl_temperature
            
      real(rp), pointer :: veloc(:,:) => NULL(),press(:) => NULL()
      real(rp), pointer :: tempe(:) => NULL()
      real(rp), pointer :: NSDissipation(:) => NULL(), TempeDissipation(:) => NULL()
      
      !Averages
      real(rp), allocatable :: avg_press(:)
      real(rp), allocatable :: avg_tempe(:)
      real(rp), allocatable :: avg_vnorm(:)
      real(rp), allocatable :: avg_vdiss(:)
      real(rp), allocatable :: avg_tdiss(:)
      
      !The first one corresponds to the average of Cn2 in time
      !The second one corresponds to Cn2 computed from averaged dissipations etc
      real(rp), allocatable :: avg_cn2(:,:)
      real(rp), allocatable :: avg_cn2u53(:,:)
      
            
      real(rp) :: &
         densi, &       !Density
         visco, &       !Fluid Viscosity
         sphea, &       !Specific Heat
         tcond, &       !Thermal Conductivity
         prtur          !Turbulent Prandtl number
      integer(ip) :: &
         units          !Unit system: 0: Meters_Seconds_Pascals_Celsius
         
      integer(ip) :: &
         lun_outres
      
      !For spatially averaging CN2, along 1 dimension
      integer(ip) :: kfl_Avg1DCn2
      integer(ip) :: Avg1DIdime
         
!       class(NavierStokesProblem), pointer :: NavierStokes
!       class(TemperatureProblem), pointer  :: Temperature
   
contains

      
      procedure :: SetNdofn          => opt_NULLSUB
      procedure :: SetNdofbc         => opt_NULLSUB
      procedure :: SpecificBounod    => opt_NULLSUBitask
      procedure :: SpecificReabcs => opt_NULLSUBitaskkflag
      procedure :: SpecificReadOnNodes => opt_NULLSUB
      procedure :: SpecificReadOnBoundaries => opt_NULLSUB
      procedure :: SpecificUpdbcs => opt_NULLSUB
      
      !Some of the usual for the physical problem are NULL for OPTICS
      procedure :: Reabcs             => opt_NULLSUB
      procedure :: Iniunk             => opt_NULLSUB
      procedure :: LinearSystemMemall => opt_NULLSUB
      procedure :: Begste             => opt_NULLSUB
      procedure :: Doiter             => opt_NULLSUB
      
      procedure :: SpecificOuterr     => opt_NULLSUB
      procedure :: SpecificIniunk    => opt_NULLSUB
      
      procedure :: SpecificBegste    => opt_NULLSUB
      procedure :: SpecificBegite    => opt_NULLSUB
      procedure :: SpecificSolite    => opt_NULLSUBitask
      procedure :: SpecificEndite    => opt_NULLSUBitask
      procedure :: EnditeElmope      => opt_NULLSUB
      procedure :: SpecificEndste    => opt_NULLSUBitask
      procedure :: SpecificCrankNicolsonEndste => opt_NULLSUB
      Procedure :: SpecificExaerr    => opt_NULLSUB      
      
      procedure :: Cvgunk            => opt_NULLSUBitaskintentin
      procedure :: Elmope            => opt_NULLSUB
      procedure :: Bouope            => opt_NULLSUB
      
      procedure :: SpecificTurnon    => opt_NULLSUB
      procedure :: SpecificTurnof    => opt_NULLSUB
      
      procedure :: Restart           => opt_NULLSUBitaskintentin
      procedure :: SpecificRestart   => opt_NULLSUBitaskintentin
      
      
      !NonNull procedures
      procedure :: Turnon            => opt_Turnon
      procedure :: OpenFiles         => opt_OpenFiles
      procedure :: CloseFiles        => opt_CloseFiles
      procedure :: SetExmod          => opt_SetExmod
      procedure :: SpecificReaphy    => opt_reaphy
      procedure :: SpecificReanut    => opt_reanut
      procedure :: SpecificReaous    => opt_reaous
      procedure :: SpecificReaMPI    => opt_reampi
      procedure :: Memall            => opt_Memall
      procedure :: Getste            => opt_Getste
      procedure :: Endste            => opt_Endste
      procedure :: Turnof            => opt_Turnof
      procedure :: Output            => opt_output
      procedure :: GetPhysicalParameters => opt_GetPhysicalParameters
   
   
      procedure :: SetDensity
      procedure :: SetViscosity
      procedure :: SetThermalConductivity
      procedure :: SetSpecificHeat
      procedure :: SetVelocityArray
      procedure :: SetPressureArray
      procedure :: SetTemperatureArray
      procedure :: SetNSDissipationArray
      procedure :: SetTempeDissipationArray
      
      procedure :: GetDissipationModel
      
      
      
   
   end type
   
   interface
      subroutine opt_Turnon(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_OpenFiles(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_CloseFiles(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      
      subroutine opt_SetExmod(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_reaphy(a,itask)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine opt_reanut(a,itask)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine opt_reaous(a,itask)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine opt_reampi(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_Memall(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_Getste(a,dtinv)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         real(rp) :: dtinv
      
      end subroutine
      
      subroutine opt_Endste(a,kfl_gotim)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         integer(ip), intent(out) :: kfl_gotim
      
      end subroutine
   
      subroutine opt_turnof(a)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
      
      end subroutine
      
      subroutine opt_output(a,itask)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         integer(ip) :: itask
      
      end subroutine
      
      subroutine opt_GetPhysicalParameters(a,acden,acvis,acsph,actco)
         use typre
         import OpticsProblem
         implicit none
         class(OpticsProblem) :: a
         real(rp) :: acden, acsph,actco,acvis
      end subroutine

   end interface
   
  interface OpticsProblem_Const
      procedure constructor
  end interface OpticsProblem_Const

  contains

    function constructor()
        class(OpticsProblem), pointer :: constructor

        allocate(constructor)

    end function constructor

   
   subroutine opt_NULLSUB(a)
      use typre
      implicit none
      class(OpticsProblem) :: a
   
   end subroutine
   
   subroutine opt_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(OpticsProblem) :: a
      integer(ip) :: itask
      
   end subroutine
   
   subroutine opt_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(OpticsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine
   
    subroutine opt_NULLSUBitaskkflag(a,itask,kflag)
      use typre
      implicit none
      class(OpticsProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
      
   end subroutine

   subroutine SetVelocityArray(a,veloc)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp), target :: veloc(:,:)
      
      a%veloc => veloc
      a%kfl_ExternalVelocity = 1
   end subroutine
   
   subroutine SetPressureArray(a,press)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp), target :: press(:)
      
      a%press => press
      a%kfl_ExternalPressure = 1
   end subroutine
   
   subroutine SetTemperatureArray(a,tempe)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp), target :: tempe(:)
      
      a%tempe => tempe
      a%kfl_ExternalTemperature = 1
   end subroutine
   
   subroutine SetDensity(a,densi)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp) :: densi
      
      a%densi = densi
   
   end subroutine
   
   subroutine SetViscosity(a,visco)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp) :: visco
      
      a%visco = visco
   
   end subroutine
   
   subroutine SetSpecificHeat(a,sphea)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp) :: sphea
      
      a%sphea = sphea
   
   end subroutine
   
   subroutine SetThermalConductivity(a,tcond)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp) :: tcond
      
      a%tcond = tcond
   
   end subroutine 
   
   subroutine SetNSDissipationArray(a,NSDissipation)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp), target :: NSDissipation(:)
      
      a%NSDissipation => NSDissipation
   end subroutine
   
   subroutine SetTempeDissipationArray(a,TempeDissipation)
      use typre
      implicit none
      class(OpticsProblem) :: a
      real(rp), target :: TempeDissipation(:)
      
      a%TempeDissipation => TempeDissipation
   end subroutine
   
   subroutine GetDissipationModel(a,kfl_dissi)
      use typre
      implicit none
      class(OpticsProblem) :: a
      integer(ip) :: kfl_dissi
      
      kfl_dissi = a%kfl_dissi
   end subroutine
   
   
end module   

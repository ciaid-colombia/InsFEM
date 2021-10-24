module Mod_lev_BaseElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_LevelSet
   use Mod_Element
   use Mod_CutMesh   
   use Mod_TemperatureElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_php_Elmdir
   implicit none
   
   class(LevelSetProblem), pointer :: a => NULL()

   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(tem_elmbuv), NOPASS, pointer :: tem_elmbuv => NULL()
      procedure(tem_elmrhu), NOPASS, pointer :: tem_elmrhu => NULL()
      procedure(), NOPASS, pointer :: ComputeTau => NULL()
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
   end type
   type(PPointer) :: ProcPointer
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: Elext => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: Testf => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(), NOPASS, pointer :: PreDirichlet => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()
      procedure(), NOPASS, pointer :: PreAssembly => NULL()
    end type
   type(PHook) :: ProcHook

   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem,igaus,itime
   
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol,dvolt0,dvolt1
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv
   
   integer(ip)           :: ielty0 = 0 !Previous element type
   
   real(rp)              :: elext
   real(rp), allocatable :: ellev(:,:)
   real(rp), allocatable :: gplev(:)
   
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: AGradV(:),gpvel(:),gpadv(:)
   real(rp)              :: gpvno
   
   real(rp), allocatable :: testf(:)
   
   real(rp) :: chale(2),tilev, gprhs
   real(rp) :: acden,acvis,acrcp  

   
   !Coordinates
   real(rp) :: gpcod(3)
   
   real(rp) :: ReferenceDtinv   
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains
   
#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussNonLinear
      implicit none
      
      call e%elmder
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
   end subroutine
   
end module

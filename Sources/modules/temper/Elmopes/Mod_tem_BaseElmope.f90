module Mod_tem_BaseElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_Temperature
   use Mod_Element
   use Mod_TemperatureElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_php_Elmdir
   implicit none
   
   class(TemperatureProblem), pointer :: a => NULL()

   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(tem_elmbuv), NOPASS, pointer :: tem_elmbuv => NULL()
      procedure(tem_elmrhu), NOPASS, pointer :: tem_elmrhu => NULL()
      procedure(), NOPASS, pointer :: ComputeTau => NULL()
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
      procedure(), NOPASS, pointer :: VelocityAndAgradV => NULL()
   end type
   type(PPointer) :: ProcPointer
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: PreLoop => NULL()
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: ElmatsToZero => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: Elext => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: ComputeTestf => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: PreDirichlet => NULL()
      procedure(), NOPASS, pointer :: Assembly => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()
      procedure(), NOPASS, pointer :: PostLoop => NULL()
    end type
   type(PHook) :: ProcHook

   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem,igaus,itime,iboun
   integer(ip) :: igaub,nboun,inodb,inode
   
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol,dvolt0,dvolt1
   real(rp)              :: dsurf
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv
   
   integer(ip)           :: ielty0 = 0 !Previous element type
   
   real(rp)              :: elext
   real(rp), allocatable :: eltem(:,:)
   real(rp), allocatable :: gptem(:)

   
   real(rp), allocatable :: elvel(:,:),grvel(:,:)
   real(rp), allocatable :: AGradV(:),gpvel(:),gpadv(:)
   real(rp)              :: gpvno
   
   real(rp), allocatable :: testf(:)
   
   real(rp) :: chale(2), timom, gprhs
   real(rp) :: acden,acsph,actco,acrea,acsou
   real(rp) :: acvis !contains: actco/acsph (plays the role of viscosity)
   real(rp) :: acrcp !contains: acrea/acsph
   
   !Coordinates
   real(rp) :: gpcod(3)
   
   

   !Tractions for Boundary Operations
   real(rp) :: tract
   
   !Reference Dtinv, for things that have their own time integration scheme (taus, dynamic subscales)
   !It is only necessary for CrankNicolson, where ReferenceDtinv = 2*dtinv. (we want to use backward euler with dt/2 in taus, subscales...)
   real(rp) :: ReferenceDtinv
  
   !Mass Matrix
   real(rp), allocatable :: elmat_mass(:,:,:,:),elrhs_mass(:,:)
   
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains
   
#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine SetPointersAndHooksToNULLSUB
      implicit none
      !External Procedures
      procedure() :: NULLSUB
      
      !Pointers
      ProcPointer%PostGaussElmats   => NULLSUB
      ProcPointer%tem_elmbuv => tem_elmbuv
      ProcPointer%tem_elmrhu => tem_elmrhu
      ProcPointer%ComputeTau        => NULLSUB
      ProcPointer%ExternalForces    => NULLSUB
      ProcPointer%VelocityAndAgradV => NULLSUB
      
      !Hooks
      ProcHook%PreLoop         => NULLSUB
      ProcHook%Initializations => NULLSUB
      ProcHook%ElmatsToZero    => NULLSUB
      ProcHook%Gathers         => NULLSUB
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%PreGauss        => NULLSUB
      ProcHook%Interpolates    => NULLSUB
      ProcHook%Elext           => NULLSUB
      ProcHook%InGauss         => NULLSUB
      ProcHook%InGaussElmats   => NULLSUB
      ProcHook%ComputeTestf    => NULLSUB
      ProcHook%PreDirichlet    => NULLSUB
      ProcHook%Assembly        => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB
      ProcHook%PostLoop        => NULLSUB
   end subroutine
   
end module  
   
 


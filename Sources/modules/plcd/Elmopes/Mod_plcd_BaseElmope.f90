module Mod_plcd_BaseElmope
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_PLCD
   use Mod_Element
   use Mod_plcd_BMatrix
   use Mod_plcd_Material
   use Mod_PointerSetter
   implicit none

   class(PLCDProblem), pointer :: a
   
   type :: PHook
      !Procedure Pointers  
      procedure(), NOPASS, pointer :: Initializations
      procedure(), NOPASS, pointer :: PreGauss
      procedure(), NOPASS, pointer :: InGauss
      procedure(), NOPASS, pointer :: InGaussElmats
      procedure(), NOPASS, pointer :: PostGaussElmats
      procedure(), NOPASS, pointer :: PreDirichlet
      procedure(), NOPASS, pointer :: Finalizations
   
   end type
   type(PHook) :: ProcHook
   
   class(FiniteElement) , pointer     :: e => NULL()
   integer(ip) :: ielem,nelem,igaus

   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol
   
   class(BMatrix), pointer :: BMat => NULL()
   
   class(ElementMaterialData), pointer :: ElementMatData => NULL()
   real(rp), pointer :: C(:,:)
   real(rp), pointer :: NodalForces(:,:) => NULL()
   
   real(rp), pointer :: stress(:) => NULL()
   
   real(rp), allocatable :: GaussElmat(:,:,:,:)
   real(rp), allocatable :: eldisp(:,:,:), elNodalForces(:,:)
   
   real(rp), allocatable, target :: gradDisp(:,:)
   real(rp), pointer :: gradDispHistory(:,:)
   
   real(rp), allocatable :: GaussElInternalForces(:,:), GaussElExternalForces(:,:) ,elInternalForces(:,:), elExternalForces(:,:)

  
   
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"   
   
    subroutine SetPointersAndHooksToNULLSUB
      implicit none
      !External Procedures
      procedure() :: NULLSUB
   

      
      !Hooks
      ProcHook%Initializations    => NULLSUB
      ProcHook%PreGauss           => NULLSUB
      ProcHook%InGauss            => NULLSUB
      ProcHook%InGaussElmats      => NULLSUB
      ProcHook%PostGaussElmats    => NULLSUB
      ProcHook%PreDirichlet       => NULLSUB
      ProcHook%Finalizations      => NULLSUB
   
   end subroutine
   
   
end module

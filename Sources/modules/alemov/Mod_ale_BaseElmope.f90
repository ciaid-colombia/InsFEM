module Mod_ale_BaseElmope
   use typre
   use Mod_Element
   use Mod_Alemov
   implicit none
   
   class(AlemovProblem), pointer :: a => NULL()

   class(FiniteElement), pointer :: e => NULL()
   real(rp) :: dvol
   real(rp),allocatable :: elmat(:,:,:,:), elrhs(:,:)
   integer(ip) :: inode,ielem,nelem,ipoin,npoin
   
   integer(ip) :: currentbvess
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(), NOPASS, pointer :: PreDirichlet => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
    end type
   type(PHook) :: ProcHook
   
!#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"   
   
end module
module Mod_rom_Procedure
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem), pointer :: a => NULL()
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: EigenSystemTurnOf => NULL()
      procedure(), NOPASS, pointer :: Projection => NULL()
   end type
   type(PPointer) :: ProcPointer


   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: TurnOn => NULL()
      procedure(), NOPASS, pointer :: Memall => NULL()
      procedure(), NOPASS, pointer :: Snapsh => NULL()
      procedure(), NOPASS, pointer :: Getste => NULL()
      procedure(), NOPASS, pointer :: Setste => NULL()
      procedure(), NOPASS, pointer :: Begste => NULL()
      procedure(), NOPASS, pointer :: DoIter => NULL()
      procedure(), NOPASS, pointer :: Endste => NULL()
      procedure(), NOPASS, pointer :: Deallo => NULL()
      procedure(), NOPASS, pointer :: TurnOf => NULL()
      procedure(), NOPASS, pointer :: InputD => NULL()
      procedure(), NOPASS, pointer :: Input => NULL()
      procedure(), NOPASS, pointer :: Output => NULL()
      procedure(), NOPASS, pointer :: OutputBkp => NULL()
   end type
   type(PHook) :: ProcHook

   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   

contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine SetPointersAndHooksToNULLSUB
      implicit none
      !External Procedures
      procedure() :: NULLSUB

      !Pointers
      ProcPointer%EigenSystemTurnOf => NULLSUB
      ProcPointer%Projection        => NULLSUB

      !Hooks
      ProcHook%TurnOn => NULLSUB
      ProcHook%Memall => NULLSUB
      ProcHook%Snapsh => NULLSUB
      ProcHook%Getste => NULLSUB
      ProcHook%Setste => NULLSUB
      ProcHook%Begste => NULLSUB
      ProcHook%DoIter => NULLSUB
      ProcHook%Endste => NULLSUB
      ProcHook%Deallo => NULLSUB
      ProcHook%TurnOf => NULLSUB
      ProcHook%InputD => NULLSUB
      ProcHook%Input  => NULLSUB
      ProcHook%Output => NULLSUB
      ProcHook%OutputBkp => NULLSUB
   end subroutine

end module

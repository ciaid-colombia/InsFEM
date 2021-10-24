subroutine rom_SetPointers
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TurnOn
   use Mod_rom_Memall
   use Mod_rom_Snapshots
   use Mod_rom_BuildSystem
   use Mod_rom_TimeStep
   use Mod_rom_Dealloc
   use Mod_rom_TurnOf
   use Mod_rom_IO
   implicit none

   !External Procedures
   procedure() :: NULLSUB

   call ResetProcedureComposition

   !-----------------------------------------------------------
   !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
   call SetPointersAndHooksToNULLSUB

   !------------------------------------------------------------------
   !Initialize the ProcedureFlags
   call SetPointersTurnOn(0)
   call SetPointersMemall(0)
   call SetPointersSnapshots(0)
   call SetPointersBuildSystem(0)
   call SetPointersTimeStep(0)
   call SetPointersDealloc(0)
   call SetPointersTurnOf(0)
   call SetPointersIO(0)

   !-------------------------------------------------------------------
   !Now we set the required pointers 

   call SetPointersTurnOn(1)
   call SetPointersMemall(1)
   call SetPointersSnapshots(1)
   call SetPointersBuildSystem(1)
   call SetPointersTimeStep(1)
   call SetPointersDealloc(1)
   call SetPointersTurnOf(1)
   call SetPointersIO(1)

   !--------------------------------------------------------------------------
   !We deallocate the procedure flags, so that they can be set the next time
   call SetPointersTurnOn(100)
   call SetPointersMemall(100)
   call SetPointersSnapshots(100)
   call SetPointersBuildSystem(100)
   call SetPointersTimeStep(100)
   call SetPointersDealloc(100)
   call SetPointersTurnOf(100)
   call SetPointersIO(100)
end subroutine

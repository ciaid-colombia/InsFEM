module Mod_rom_Dealloc
   use typre
   use Mod_Mesh
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_PhysicalProblem
   implicit none
   private
   public SetPointersDealloc

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersDealloc(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   

      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1

      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            select case (a%kfl_itask)

            case('FOM')
               ProcHook%Deallo => DeallocBuild

            case('BFS')
               ProcHook%Deallo => DeallocBuild

            case('SNP')
               ProcHook%Deallo => DeallocSnaps

            case('INT')
               ProcHook%Deallo => DeallocRun

            case('ROM')
               ProcHook%Deallo => DeallocRun
            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine DeallocBuild
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)

      !The snapshots are local
      call a%Memor%dealloc(a%ndofn,npoin,a%Snapshots,'Snapshots','rom_Memall')
      call a%Memor%dealloc(a%ndofn,npoin,a%SnapMean,'SnapMean','rom_Memall')

   end subroutine

   subroutine DeallocSnaps
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(a%ndofn,npoin,a%Snapshots,'Snapshots','rom_Memall')

   end subroutine

   subroutine DeallocRun
      integer(ip) :: npoin,oldnpoin

      call a%Mesh%GetNpoin(npoin)

      !Allocate space for the Full-order solution
      call a%Memor%dealloc(a%nconv,a%bEnergy,'bEnergy','rom_Memall')
      call a%Memor%dealloc(a%nconv,a%sigma,'sigma','rom_Memall')
      call a%Memor%dealloc(a%ndofn,npoin,a%FOMSolution,'FOMSolution','rom_Memall')
      call a%Memor%dealloc(a%ndofn,npoin,a%ndofr,a%Basis,'Basis','rom_Memall')
      call a%Memor%dealloc(a%ndofn,npoin,a%SnapMean,'SnapMean','rom_Memall')

      if (a%kfl_updateBasis) then
         call a%Problem%OldMesh%GetNpoin(oldnpoin)
         call a%Memor%dealloc(a%ndofn,oldnpoin,a%ndofr,a%OldMeshBasis,'OldMeshBasis','rom_Memall')
         call a%Memor%dealloc(a%ndofn,oldnpoin,a%OldMeshSnapMean,'OldMeshSnapMean','rom_Memall')
      end if

   end subroutine

end module

subroutine rom_Dealloc(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_Dealloc
   implicit none
   class(PodRomProblem),target :: PodRom

   a => PodRom

   call ProcHook%Deallo
end subroutine

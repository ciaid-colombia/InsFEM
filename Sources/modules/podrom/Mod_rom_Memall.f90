module Mod_rom_Memall
   use typre
   use Mod_TimeIntegrator
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_PhysicalProblem
   implicit none
   private
   public SetPointersMemall

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersMemall(itask)
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
               ProcHook%Memall => MemallBuild

            case('BFS')
               ProcHook%Memall => MemallBuild

            case('SNP')
               ProcHook%Memall => MemallSnaps

            case('INT')
               ProcHook%Memall => MemallRun

            case('ROM')
               ProcHook%Memall => MemallRun

            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine MemallBuild
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)

      !The snapshots are local
      call a%Memor%alloc(a%ndofn,npoin,a%Snapshots,'Snapshots','rom_Memall')
      call a%Memor%alloc(a%ndofn,npoin,a%SnapMean,'SnapMean','rom_Memall')

      !Set some counters to zero
      a%isnap = 0
      a%wsnap = 0

   end subroutine

   subroutine MemallSnaps
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(a%ndofn,npoin,a%Snapshots,'Snapshots','rom_Memall')
      a%isnap = 0
      a%wsnap = 0

   end subroutine

   subroutine MemallRun
      integer(ip) :: npoin,oldnpoin

      call a%Mesh%GetNpoin(npoin)
      a%isnap = 0

      !Allocate space for the Full-order solution
      call a%Memor%alloc(a%nconv,a%sigma,'sigma','rom_Memall')
      call a%Memor%alloc(a%ndofn,npoin,a%FOMSolution,'FOMSolution','rom_Memall')
      call a%Memor%alloc(a%ndofn,npoin,a%ndofr,a%Basis,'Basis','rom_Memall')
      call a%Memor%alloc(a%ndofn,npoin,a%SnapMean,'SnapMean','rom_Memall')

      if (a%kfl_updateBasis) then
         call a%Problem%OldMesh%GetNpoin(oldnpoin)
         call a%Memor%alloc(a%ndofn,oldnpoin,a%ndofr,a%OldMeshBasis,'OldMeshBasis','rom_Memall')
         call a%Memor%alloc(a%ndofn,oldnpoin,a%OldMeshSnapMean,'OldMeshSnapMean','rom_Memall')
      end if

   end subroutine

end module

subroutine rom_Memall(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Memall
   use Mod_rom_Procedure
   implicit none
   class(PodRomProblem),target :: PodRom

   a => PodRom

   call ProcHook%Memall
end subroutine

subroutine rom_eigenSystemMemall(a)
   use typre
   use Mod_SVDLibrary
   use Mod_EPSLibrary
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a

#ifdef SLEPC
   select case(a%kfl_eigentype)
   case('EPS')
      a%EigenLibrary => EPSLibraryConst()
      call a%EigenLibrary%Initialize('EPS')
   case('SVD')
      a%EigenLibrary => SVDLibraryConst()
      call a%EigenLibrary%Initialize('SVD')
   end select
   call a%EigenLibrary%CreateSystem(a%EigenSystem,a%Memor)
#endif

end subroutine


module Mod_rom_Snapshots
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_PhysicalProblem
   implicit none
   private
   public SetPointersSnapshots

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersSnapshots(itask)
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
            
            select case(a%kfl_eigentype)
            case('EPS')

            case('SVD')
                select case (a%kfl_itask)
                case('FOM')
                    ProcHook%Snapsh => Snapshots_FOM
                end select
            end select
            
            select case (a%kfl_itask)
            case('SNP')
                ProcHook%Snapsh => Snapshots
            end select
        endif

    case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine Snapshots_FOM
      integer(ip) :: npoinLocal

      call a%Mesh%GetNpoinLocal(npoinLocal)
      !Snapshots of FOM solution
      if (a%kfl_snapshotSpecific) then
         call a%Problem%SpecificGetUnkno(a%Snapshots)
      else
         call a%Problem%GetUnkno(a%Snapshots)
         if (a%ProblemType == 'BOUSS') call a%Problem%SpecificGetUnkno(a%Snapshots)
      endif
      call a%EigenSystem%SetSnapshots(a%Snapshots(:,1:npoinLocal))
   end subroutine

   subroutine Snapshots
      integer(ip) :: npoinLocal

      call a%Mesh%GetNpoinLocal(npoinLocal)
      if (a%kfl_snapshotSpecific) then
         call a%Problem%SpecificGetUnkno(a%Snapshots)
      else
         call a%Problem%GetUnkno(a%Snapshots)
         if (a%ProblemType == 'BOUSS') call a%Problem%SpecificGetUnkno(a%Snapshots)
      endif
   end subroutine

end module

subroutine rom_SnapshotsToSystem(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Snapshots
   use Mod_rom_Procedure
   implicit none
   class(PodRomProblem),target :: PodRom
   integer(ip) :: istep

   a => PodRom

   call a%Problem%GetIstep(istep)

   a%isnap = istep/a%SnapshotsInterval
   if (istep > 0 .and. a%isnap <= a%nsnap .and. mod(istep,a%SnapshotsInterval) == 0) then
      call a%Timer%Assembly%Tic
      call ProcHook%Snapsh
      call a%Timer%Assembly%Toc
      call a%Timer%Output%Tic
      call ProcHook%OutputBkp
      call a%Timer%Output%Toc
   endif
      
end subroutine

module Mod_rom_TurnOf
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   implicit none
   private
   public SetPointersTurnOf

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersTurnOf(itask)
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

            ProcHook%TurnOf => TurnOf

            select case (a%kfl_itask)

            case('FOM')
               call PrependProcedure(ProcHook%TurnOf,TurnofBuild)
               ProcPointer%EigenSystemTurnof => EigenSystemTurnofBuild

            case('BFS')
               call PrependProcedure(ProcHook%TurnOf,TurnofBuild)
               ProcPointer%EigenSystemTurnof => EigenSystemTurnofBuild

            case('SNP')
               call PrependProcedure(ProcHook%TurnOf,TurnofSnaps)

            case('INT')
               ProcPointer%EigenSystemTurnof => EigenSystemTurnofRun

            case('ROM')
               ProcPointer%EigenSystemTurnof => EigenSystemTurnofRun

            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine Turnof

      call a%Output
      call ProcPointer%EigenSystemTurnof
      call a%DeallocAll
   end subroutine

   subroutine TurnofBuild

      a%nsnap = a%isnap

      select case(a%kfl_eigentype)
      case('EPS')

      case('SVD')
      call a%Timer%Assembly%Tic
         if (a%kfl_SubstractMean) then
            call a%EigenSystem%AssemblySnapshots
            call a%EigenSystem%SnapshotsMean(a%nsnap)
         else
            call a%EigenSystem%SnapshotsMeanToZero
         end if

         call a%EigenSystem%AssemblySnapshots
      end select
      call a%Timer%Assembly%Toc

      call a%Timer%Solve%Tic
      if (a%kfl_massMatrix .eqv. .true.) call a%EigenSystem%AssemblySnapshotsMass
      call a%SVD
      if (a%kfl_massMatrix .eqv. .true.) call a%EigenSystem%MassMultBasis
      call a%Timer%Solve%Toc
      
   end subroutine

   subroutine TurnofSnaps

      a%nsnap = a%isnap

   end subroutine

   subroutine EigenSystemTurnofRun

      call a%EigenSystem%DeallocateRun
      call a%EigenLibrary%Finalize
      call a%EigenLibrary%DeallocateSystem(a%EigenSystem,a%Memor)
   end subroutine

   subroutine EigenSystemTurnofBuild

      if (a%kfl_massMatrix .eqv. .true.) call a%EigenSystem%LMassDeallocate
      call a%EigenSystem%Deallocate
      call a%EigenLibrary%Finalize
      call a%EigenLibrary%DeallocateSystem(a%EigenSystem,a%Memor)
   end subroutine

end module

subroutine rom_TurnOf(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TurnOf
   implicit none
   class(PodRomProblem),target :: PodRom

   a => PodRom

   call a%Problem%Timer%Total%Tic
   call a%Problem%Timer%Turnof%Tic
   call ProcHook%Turnof
   call a%OutputTimes('turnof')
   call a%Closefi
   call a%Problem%Timer%Total%Toc
   call a%Problem%Timer%Turnof%Toc
end subroutine

subroutine rom_OutputTimes(a,task)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem)  :: a
   character(6) :: task
   real(rp)     :: cpu_tim(4)

   cpu_tim = 0.0_rp
   !Brief computing times
   select case (task)
   case ('turnon')
   call a%GetTimes(cpu_tim)
      if (a%MPIrank == a%MPIroot) write(a%lun_outro,700) "Input time:              ", (cpu_tim(1))
   case ('turnof')
   call a%GetTimes(cpu_tim)
      if (a%MPIrank == a%MPIroot) write(a%lun_outro,700) "Snapshots assembly time: ", (cpu_tim(2))
      if (a%MPIrank == a%MPIroot) write(a%lun_outro,700) "SVD time:                ", (cpu_tim(3))
      if (a%MPIrank == a%MPIroot) write(a%lun_outro,700) "Output time:             ", (cpu_tim(4))
   end select
   700 format(a25,1(1x,e10.3))

end subroutine

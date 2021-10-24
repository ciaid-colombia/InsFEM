module Mod_rom_TurnOn
   use typre
   use Mod_PETScSystem
   use Mod_ParallelSystemInterface
   use Mod_PodRom
   use Mod_rom_Procedure
   implicit none
   private
   public SetPointersTurnOn

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersTurnOn(itask)
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
               ProcHook%TurnOn => TurnOnBuild

            case('SNP')
               ProcHook%TurnOn => TurnOnSnaps

            case('BFS')
               ProcHook%TurnOn => TurnOnBuild

            case('INT')
               ProcHook%TurnOn => TurnOnRun

            case('ROM')
               ProcHook%TurnOn => TurnOnRun
            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine TurnOnBuild

      call a%InputData
      call a%Memall
      call a%EigenSystem%SetLinearSystem(a%Problem%LinearSystem)
      call a%EigenSystem%SetOrderingPointers
      call a%EigenSystem%InitBuild(a%isnap,a%SnapshotsInterval,a%nsnap)
      if (a%kfl_massMatrix .eqv. .true.) call a%EigenSystem%SetLumpedMass(a%Mesh%vmass)
      if (a%kfl_eigentype .ne. "EPS") call a%Input
   end subroutine

   subroutine TurnOnSnaps

      call a%Memall
   end subroutine

   subroutine TurnOnRun

      call a%InputData
      call a%Memall

      a%Problem%cotol = a%cotol
      a%Problem%maxit = a%nNonLinearIterations
      call a%EigenSystem%SetLinearSystem(a%Problem%LinearSystem)
      call a%EigenSystem%SetMatrixPointers
      call a%EigenSystem%SetRHSPointers
      call a%EigenSystem%SetOrderingPointers
      call a%EigenSystem%InitSystem
      call a%InputBasis

   end subroutine

end module

subroutine rom_TurnOn(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TurnOn
   implicit none
   class(PodRomProblem),target :: PodRom
   integer(ip)    :: gnpoin,npoin,npoinLocal
   character(150) :: optfile
   character(3)   :: exmod

   a => PodRom

   call a%Problem%Timer%Total%Tic
   call a%Problem%Timer%Turnon%Tic
   
   !Open Files
   call a%Openfi

   !Read Data from Files
   if (a%MPIrank == a%MPIroot) then
      call a%Readat
   endif
   call a%ReadatMPI

   !Warning and errors
   call a%outerr

   !Init eigensystem
   call a%EigenSystemMemall

   call a%Mesh%GetGnpoin(gnpoin)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)

   exmod = a%exmod
   optfile = adjustl(trim(a%InputFolder))//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.rom.sol'

   call a%SetEigenSystem
   call a%EigenSystem%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call a%EigenSystem%SetFlush(a%kfl_flush)
   call a%EigenSystem%Init(gnpoin,npoin,npoinLocal,a%ndofn,optfile,a%exmod,a%lun_solro,a%Memor)
   call a%EigenSystem%SetLinearSystem(a%Problem%LinearSystem)
   call a%EigenSystem%SetOrderingPointers
   
   call a%SetPointers
   call ProcHook%TurnOn
   call a%OutputTimes('turnon')

   call a%Problem%Timer%Total%Toc
   call a%Problem%Timer%Turnon%Toc

end subroutine

subroutine rom_InputBasis(a)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a

   call a%Input

   select case(a%Problem%kfl_ProjectionType)
   case(0)
      call a%ComputeLumpedMassMatrix
   case(1)
      call a%ComputeMassMatrix
   end select

end subroutine

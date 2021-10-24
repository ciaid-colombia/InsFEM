module Mod_rom_TimeStep
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   implicit none
   private
   public SetPointersTimeStep

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersTimeStep(itask)
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
               ProcHook%setste => SetsteBuild
               ProcHook%Endste => EndsteBuild

            case('SNP')
               ProcHook%setste => SetsteBuild
               ProcHook%Endste => EndsteBuild

            case('ROM')
               ProcHook%Getste => Getste
               ProcHook%Setste => SetsteRun
               ProcHook%Begste => Begste
               ProcHook%DoIter => doiter
               ProcHook%Endste => EndsteRun

            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine SetsteBuild(dtinv,ctime)
      implicit none
      real(rp) :: dtinv, ctime

      a%ctime = ctime
   end subroutine

   subroutine SetsteRun(dtinv,ctime,timef,nsmax)
      implicit none
      real(rp) :: dtinv, ctime,timef
      integer(ip) :: nsmax
      logical::kfl_coupconv

      a%ctime = ctime
      call a%Problem%Setste(dtinv,ctime,timef,nsmax)
   end subroutine

   subroutine Begste
      implicit none
      a%Problem%cpiter = 0
      call a%Problem%Begste
   end subroutine
   
   subroutine begite
      implicit none
      !Initializations

      call a%Problem%Begite
   end subroutine

   subroutine Getste(dtinv)
      implicit none
      real(rp) :: dtinv

      call a%Problem%getste(dtinv)
   end subroutine

   subroutine doiter
      implicit none

      !Iterate all the non-linear iterations
      call Begite
      do while(a%Problem%kfl_goite==1)
         call TimeStep
         call Endite(0)
         call Endite(1)
      enddo   
      call Endite(2)
      if(a%Problem%kfl_docoupconv) call Endite(4)
   end subroutine

   subroutine TimeStep
      implicit none

      a%Problem%itera = a%Problem%itera + 1

      if (a%MPIrank == a%MPIroot) then
         if(a%Problem%itera==1) write(a%Problem%lun_solve,100) a%Problem%istep
         write(a%Problem%lun_solve,101) a%Problem%itera
      endif

      call a%EigenSystem%InitOperators
      call a%BuildSystem
      call a%Problem%Timer%SolveLinearSystem%Tic
      call a%EigenSystem%SolveProjSystem
      call a%SolutionToPhysicalProblem
      call a%Problem%Timer%SolveLinearSystem%Toc
      call a%EigenSystem%DeallocOperators

      !Formats. 
      100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
      101 format('------------------------------------------------------------', &
         /,'   INNER ITERATION NUMBER: ',i5)
   end subroutine

   subroutine Endite(itask)
      implicit none
      integer(ip) :: itask
      
      select case(itask)
      case(0)
         call a%Problem%Endite(0)
      case(1)
         call a%Problem%Endite(1)
      case(2)
         call a%Problem%Endite(2) 
      case(4)
         call a%Problem%Endite(4) 
         a%Problem%cpiter = a%Problem%cpiter + 1
      end select
   end subroutine

   subroutine EndsteBuild(kfl_gotim)
      implicit none
      integer(ip) :: kfl_gotim
      
      call a%SnapshotsToSystem
   end subroutine

   subroutine EndsteRun(kfl_gotim)
      implicit none
      integer(ip) :: kfl_gotim
      integer(ip) :: istep

      call a%Problem%Endste(kfl_gotim)
      if (a%kfl_outBasis) call a%Postpr
      call a%Problem%GetIstep(istep)
      a%isnap = istep/a%SnapshotsInterval
   end subroutine
   
end module

subroutine rom_Getste(PodRom,dtinv)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TimeStep
   implicit none
   class(PodRomProblem),target :: PodRom
   real(rp) :: dtinv

   a => PodRom

   call ProcHook%Getste
end subroutine

subroutine rom_Setste(PodRom,dtinv,ctime)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TimeStep
   implicit none
   class(PodRomProblem),target :: PodRom
   real(rp) :: dtinv, ctime

   a => PodRom

   call ProcHook%Setste(dtinv,ctime)
end subroutine

subroutine rom_Begste(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TimeStep
   implicit none
   class(PodRomProblem),target :: PodRom

   a => PodRom

   call ProcHook%Begste
end subroutine

subroutine rom_doiter(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TimeStep
   implicit none
   class(PodRomProblem),target :: PodRom

   a => PodRom

   call a%Problem%Timer%Total%Tic
   call a%Problem%Timer%Doiter%Tic
   call ProcHook%DoIter
   call a%Problem%Timer%Total%Toc
   call a%Problem%Timer%Doiter%Toc
end subroutine

subroutine rom_Endste(PodRom,kfl_gotim)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_TimeStep
   implicit none
   class(PodRomProblem),target :: PodRom
   integer(ip) :: kfl_gotim

   a => PodRom

   call ProcHook%Endste(kfl_gotim)
end subroutine

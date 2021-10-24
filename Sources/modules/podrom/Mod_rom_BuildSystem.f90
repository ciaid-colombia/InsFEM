module Mod_rom_BuildSystem 
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   implicit none
   private
   public SetPointersBuildSystem

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersBuildSystem(itask)
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
            
            if ((a%kfl_SubstractMean) .and.(a%kfl_eigentype .ne. "EPS") ) ProcPointer%Projection => MeanToRHS
            
            select case (a%kfl_Precondition)
            case(0)
               call ConcatenateProcedures(ProcPointer%Projection,Galerkin)
            case(1)
               call ConcatenateProcedures(ProcPointer%Projection,PetrovGalerkin)
            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine

   subroutine MeanToRHS
      call a%EigenSystem%MeanToRHS
   end subroutine

   subroutine PetrovGalerkin
      call a%EigenSystem%PetrovGalerkinProj
   end subroutine

   subroutine Galerkin
      call a%EigenSystem%GalerkinProj
   end subroutine

end module

subroutine rom_BuildSystem(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_BuildSystem
   implicit none
   class(PodRomProblem),target :: PodRom
   integer(ip) :: kfl_HangingNodes,kfl_perio
   
   a => PodRom

   call a%Problem%Timer%BuildLinearSystem%Tic
   
   call a%EigenSystem%ToZero
   
   call a%Problem%BuildLinearSystem
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBCToZero(a%ndofn,a%Problem%LinearSystem,a%Memor)

   !Hanging nodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiagToZero(a%ndofn,a%Problem%LinearSystem,a%Memor)

   call a%EigenSystem%SetMatrixPointers
   call a%EigenSystem%SetRHSPointers
   call a%EigenSystem%SetupOperators

   call ProcPointer%Projection
   call a%EigenSystem%MatProjection
   call a%EigenSystem%VecProjection

   call a%Problem%Timer%BuildLinearSystem%Toc

end subroutine

subroutine php_solite(a)
   !-----------------------------------------------------------------------
   !    This routine solves an a%iteration for the Physical Problem Equations
   !-----------------------------------------------------------------------
   use MPI
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a

   !MPI
   integer :: ierr

   integer(ip) :: kfl_perio, kfl_HangingNodes

   !Update inner a%iteration counter and write headings in the solver file.
   a%itera = a%itera + 1
   if (a%MPIrank == a%MPIroot) then
      if(a%itera==1) write(a%lun_solve,100) a%istep
      write(a%lun_solve,101) a%itera
   endif

   call a%SpecificSolite(1)

   if (a%kfl_SkipSystemSolve .eqv. .false.) then
   !   call MPI_BARRIER(a%MPIcomm,ierr)
      call a%Timer%BuildLinearSystem%Tic
      !Construct the system matrix and right-hand-side.
      call a%BuildLinearSystem

   !   call MPI_BARRIER(a%MPIcomm,ierr)
      call a%Timer%BuildLinearSystem%Toc

      !Solve the algebraic system.
      call a%Timer%SolveLinearSystem%Tic
      call a%LinearSystem%Solve(a%unkno)

      !Ghostcommunicate
      call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%unkno)

      !HangingNodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(a%ndofn,a%unkno)

      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(a%ndofn,a%unkno)

   !   call MPI_BARRIER(a%MPIcomm,ierr)
      call a%Timer%SolveLinearSystem%Toc
   endif

   call a%SpecificSolite(2)

   !Formats.
   100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
   101 format('------------------------------------------------------------', &
         /,'   INNER ITERATION NUMBER: ',i5)

end subroutine php_solite

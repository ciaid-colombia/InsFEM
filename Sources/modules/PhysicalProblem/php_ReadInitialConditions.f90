subroutine php_InitializeInitialConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   a%kfl_incnd = 0       ! Initial conditions  
end subroutine

subroutine php_ReadInitialConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   a%kfl_incnd = int(a%Listener%param(1))
   a%picnd = a%Listener%param(2:11)
end subroutine

subroutine php_ScatterInitialConditions(a)
   use Mod_PhysicalProblem
   use MPI
   implicit none
   class(PhysicalProblem) :: a
   integer :: ierr

   CALL MPI_BCAST(a%kfl_incnd, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%picnd, size(a%picnd), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
end subroutine

subroutine php_InitializeConstantBoundaryConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   a%kfl_conbc = 1                                    !Constant boundary conditions
end subroutine

subroutine php_ReadConstantBoundaryConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   !Constant/variable boundary conditions
   if (a%Listener%exists('CONST')) then
      a%kfl_conbc = 1
   else
      a%kfl_conbc = 0
   endif
end subroutine

subroutine php_ScatterConstantBoundaryConditions(a)
   use typre
   use Mod_PhysicalProblem
   use MPI
   implicit none
   class(PhysicalProblem) :: a
   integer :: ierr
   
   !Constant/variable boundary conditions
   CALL MPI_BCAST(a%kfl_conbc, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
end subroutine

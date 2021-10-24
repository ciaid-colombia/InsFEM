subroutine php_Turnon(a)
   !-----------------------------------------------------------------------
   !    This routine performs the following tasks:
   !    - Gets file names and open them.
   !    - Read physical problem data 
   !    - Write some info
   !    - Check if there are errors
   !    - Allocate memory
   !-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   call a%Timer%Total%Tic
   call a%Timer%Turnon%Tic
   
   call a%SetExmod
   call a%SetNdofn
   call a%SetNdofbc
   call a%OpenFiles
   if (a%MPIrank == a%MPIroot) then
      !Read Physical Problem
      call a%Reaphy
      
      !Read Numerical Treatment
      call a%Reanut
      
      !Read Output Strategy
      call a%Reaous
      
   endif
   call a%ReaMPI
   
   !Memory allocation
   call a%Memall

   !Warnings and errors
   call a%Outerr(1)
   
   !Boundary conditions
   call a%Reabcs
   
   !Allocate Memory for solving the linear system
   call a%LinearSystemMemall
   
   !Some Initializations
   call a%Initializations

   call a%SpecificTurnon
   
   !Unknown initialization
   call a%Iniunk
   
   call a%CloseDataFile
   
   !call a%FilePostpr%inipos(a%MPIrank,a%MPIroot)

   !Warnings and errors
   call a%Outerr(2)
   
   call a%Timer%Turnon%Toc
   call a%Timer%Total%Toc
   
end subroutine php_Turnon

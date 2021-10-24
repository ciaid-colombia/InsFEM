subroutine php_LinearSystemMemall(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: npoin
   character(150) :: auxstring
   character(150) :: auxstring2
   
   auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(a%exmod))//'.sol '
   auxstring2 = trim(a%exmod)
   call a%ParallelLibrary%CreateSystem(a%LinearSystem,a%Memor) 
   call a%LinearSystem%SetFlush(a%kfl_flush)
   call a%LinearSystem%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call a%Mesh%InitializeSystem(a%ndofn,0,a%LinearSystem,a%Memor,auxstring,auxstring2,a%lun_solve,a%extramatrix)
   
   call a%Mesh%GetNpoin(npoin)
   call a%Memor%alloc(a%ndofn,npoin,a%unkno,'unkno',trim(a%exmod)//'_LinearSystemMemall')

end subroutine

subroutine php_LinearSystemTurnof(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: npoin
   
   call a%Mesh%GetNpoin(npoin)
   call a%Memor%dealloc(size(a%unkno,1),size(a%unkno,2),a%unkno,'unkno',trim(a%exmod)//'_LinearSystemTurnof')
   call a%LinearSystem%Deallocate
   call a%ParallelLibrary%DeallocateSystem(a%LinearSystem, a%Memor) 

end subroutine


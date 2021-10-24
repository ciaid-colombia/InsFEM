subroutine case_OldDomain(a,itask)
   !-----------------------------------------------------------------------
   !    This routine loads an old mesh
   !-----------------------------------------------------------------------
   use typre
   use Mod_GeneralCase
   use Mod_CaseVariables
   use Mod_DomainVariables
   implicit none
   class(GeneralCase), target :: a
   type(masterVariables), pointer :: m => NULL()
   type(domainVariables), pointer :: d => NULL()

   integer(ip) :: itask

   real(rp)    :: cputim2,cputim1

   m => a%caseVars%masterVars
   d => a%caseVars%domainVars
   
   !Initialize OldMesh
   if (itask == 1) then
   
      call m%cpu_start(2)%Tic
      !OldMesh Input data
      call d%OldMesh%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
      call d%OldMesh%SetReadType(m%ReadTypeString)
      call d%OldMesh%SetInputFile(m%oldnamda)
      call d%OldMesh%SetInputFolder(m%OldDataFolder)
      call d%OldMesh%SetOutputFolder(m%ResultsFolder)
      call d%OldMesh%SetOutputFiles(m%lun_memor,m%lun_outpu)
      call d%OldMesh%SetParallelLibrary(m%ParallelLibrary)
      call d%OldMesh%SetFlush(m%kfl_flush)
      
      !Read data
      call d%OldMesh%Readat
      
      !Initial operations (Graph and Scatter)
      call d%OldMesh%Initialize

      call m%cpu_start(2)%Toc
   
   !Deallocate Read Globals
   elseif (itask == 2) then
      call m%cpu_start(2)%Tic 

      call d%OldMesh%DeallocateReadGlobals
      
      call m%cpu_start(2)%Toc 
   
   !Deallocate Locals
   elseif (itask == 3) then
      call d%OldMesh%Turnof
   
   endif

end subroutine 

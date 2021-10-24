subroutine case_Domain(a,itask)
   use typre
   use Mod_GeneralCase
   use Mod_CaseVariables
   use Mod_int2str
   use Mod_Mesh_AR_Interface
   implicit none
   class(GeneralCase), target :: a
   integer(ip) :: itask
   
   type (domainVariables), pointer :: d => NULL()
   type (masterVariables), pointer :: m => NULL()
   type(adaptiveVariables), pointer :: ad => NULL()
   
   character(150) :: fil_outpu_dom

   m => a%caseVars%masterVars
   d => a%caseVars%domainVars
   
   ad => a%caseVars%adaptiveVars
   
   !Initialize Mesh
   if (itask == 1) then
   
      call m%cpu_start(2)%Tic
      !Mesh Input data
      call d%Mesh%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
      call d%Mesh%SetReadType(m%ReadTypeString)
      call d%Mesh%SetInputFolder(m%BaseDataFolder)
      call d%Mesh%SetInputFile(m%namda)
      call d%Mesh%SetOutputFolder(m%ResultsFolder)
      call d%Mesh%SetOutputFiles(m%lun_memor,m%lun_outpu)
      call d%Mesh%SetParallelLibrary(m%ParallelLibrary)
      call d%Mesh%SetFlush(m%kfl_flush)
      
      !Read data
      call d%Mesh%Readat
   
      !If there is InitialUniformRefinement, we deactivate periodic BC
      !until the refinement has been done
      ! (PBC do not work with adaptive)
      call a%InitialUniformRefinement(0)
   
      !Initial operations (Graph and Scatter)
      call d%Mesh%Initialize
      call d%Mesh%InitialOrdering

      call m%cpu_start(2)%Toc
      
      call m%cpu_start(3)%Tic      
            
      !call d%Mesh%ExportigPointGlobNumber
      
      !Open d%Mesh file 
      call m%FilePostpr%Postpr(d%Mesh)
      
      if (m%kfl_outfo == 0) then
         call m%FilePostpr%bpostpr(d%Mesh%exnor,'exnor',1_ip,0.0_rp,d%Mesh)
         call m%FilePostpr%postpr(d%Mesh%vmass,'vmass',1_ip,0.0_rp,d%Mesh)
      end if

      call m%cpu_start(3)%Toc
   
   !Deallocate Read Globals
   elseif (itask == 2) then
      call d%Mesh%DeallocateReadGlobals
   
   !Deallocate Locals
   elseif (itask == 3) then
      call d%Mesh%Turnof

   endif

end subroutine 

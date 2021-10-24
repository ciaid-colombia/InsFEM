module mod_CaseTurnon
   use Mod_GeneralCase
   use Mod_CaseVariables
   use MPI
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DriverInterface
   use Mod_DCHashCharSize
   implicit none
   
   type(caseVariables),   pointer :: c => NULL()
   type(masterVariables), pointer :: m => NULL()
   
contains
   subroutine LoopTurnonAndSetChannels(myDc)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDc,myDriver)
      
      !We turn on the Drivers
      call myDriver%Turnon(c)
      
      !We set the channels of the driver
      call myDriver%SetOutChannels('ADD')
      call myDriver%SetInChannels
   end subroutine   
      
end module   

subroutine case_Turnon(a)
   use mod_CaseTurnon
   implicit none
   class(GeneralCase), target :: a
   character(DCHashCharSize) :: DoFirstKeys(1)
   integer(ip) :: ierr

   c => a%caseVars
   m => a%caseVars%masterVars

   call MPI_BARRIER(m%MPIcomm,ierr)
   
   !Setup Mesh
   call a%Domain(1)

   !Setup Adaptive
   call a%Adaptive(1)

   call a%InitialUniformRefinement(1)

   call m%cpu_start(1)%Tic()

   call a%SetupIO(0)
   
   !Turnon each driver
   DoFirstKeys(1) = 'ALEPR'
   call a%DriverList%LoopList(LoopTurnonAndSetChannels,DoFirstKeys)

   !Couplings are done externally through the Master
   !Don't need to do anything here
   call m%cpu_start(1)%Toc()

   call m%cpu_start(2)%Tic
   
   !Domain read deallocations
   call a%Domain(2)
   
   !Initial Uniform Refinement
   call a%InitialUniformRefinement(2)

   call a%caseVars%masterVars%FilePostpr%SetDuplicateInterfaceElements(0)

   call m%cpu_start(2)%Toc

end subroutine 

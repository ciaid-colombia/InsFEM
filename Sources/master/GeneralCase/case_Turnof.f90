module Mod_caseTurnof
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   implicit none

   type(caseVariables), pointer  :: c => NULL()
   type(masterVariables), pointer :: m => NULL()    
   
contains
   subroutine LoopTurnof(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()

      call ExtractDriver(myDc,myDriver)
      call myDriver%FinalizeChannelLists  
      call myDriver%Turnof(c)
      call m%MasterMemo%DeallocObj(0,'Driver','case_Turnof',1_ip) 
   end subroutine
   
   subroutine LoopDealloc(myDC)
      class(DistributedContainer), pointer :: myDC 
      class(DriverInterface), pointer :: myDriver => NULL()
         
      call ExtractDriver(myDc,myDriver)
      deallocate(myDriver)
   end subroutine
end module

subroutine case_Turnof(a)
   !-----------------------------------------------------------------------
   !    Stops the run.
   !-----------------------------------------------------------------------
   use Mod_caseTurnof
   implicit none
   class(GeneralCase), target :: a

   c => a%caseVars
   m => a%caseVars%masterVars

   !Modules measure their own times
   call m%cpu_end(1)%Tic
   !We Finalize all the stuff related to the interpolators
   call a%FinalizeInterpolators
   call a%CaseInterpolatorList%Finalize

   !Loops through all the drivers and turns them off
   !HOWEVER, WE DO NOT DEALLOCATE THEM YET, SINCE WE WILL STILL NEED 
   !TO ASK FOR THEIR MEMORY AND TIME
   !**** Deallocation done some lines below
   call a%DriverList%LoopList(LoopTurnof)
  
   !Deallocate Mesh (Locals)
   call a%Domain(3)

   !Deallocate adaptive if necessary
   call a%Adaptive(100)
   !Close postprocess files
   call m%FilePostpr%ClosePostFile
   call m%FilePostpr%CloseGrafFile
   call m%FilePostpr%Finalize
   deallocate(m%FilePostpr)
   
   !Deallocate read and write to disk
   call a%SetupIO(1)
   deallocate(m%Readerpr,m%Writerpr)

   call m%cpu_end(1)%Toc

   !The parallel library and postpro are deallocated afterwards, but count theirs memory now (four Outmem)
   call m%MasterMemo%DeallocObj(0,'FilePostpr','case_Turnof',1_ip)
   call m%MasterMemo%DeallocObj(0,'Readerpr','case_Turnof',1_ip)
   call m%MasterMemo%DeallocObj(0,'Writerpr','case_Turnof',1_ip)
   call m%MasterMemo%DeallocObj(0,'ParallelLibrary','case_Turnof',1_ip) 

   !Writes memory used.
   call a%OutMem
   
   !Write CPU time heading and master's CPU time
   call a%Cputab
   
   !Kill adaptive finally
   call a%Adaptive(101)
   
   !HERE WE DEALLOCATE THE DRIVERS
   call a%DriverList%LoopList(LoopDealloc)
   
   !Free the DriverList
   call a%DriverList%Finalize

   !Finalize the parallel Library 
   call m%ParallelLibrary%Finalize
   deallocate(m%ParallelLibrary)
   

end subroutine 

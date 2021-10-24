subroutine Reapro 
   !This routine reads global problem definition.
   use typre
   use def_master
   use Mod_DistributedContainer
   use Mod_GeneralCase
   use Mod_DC_GeneralCase
   implicit none
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()

   !Timers
   call cpu_start(5)%Tic

   !Arguments and Environment
   call GetArgumentsAndEnvironment

   !Timers
   call cpu_start(5)%Toc
   
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call myCase%Reapro
      
      call CaseList%GetNext(myDC)
   enddo
  
end subroutine 

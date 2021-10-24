subroutine Turnon
   use def_master
   use Mod_DC_GeneralCase
   use Mod_DistributedContainer
   implicit none 
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()

   !Loop the case list and turnon each case
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call myCase%Turnon
      
      call CaseList%GetNext(myDC)
   enddo
   if(MPIrank==MPIroot) write(*,*) 'Case: Turn on ready'

end subroutine Turnon



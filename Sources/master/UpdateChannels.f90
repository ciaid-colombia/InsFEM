subroutine UpdateChannels
   !This routine reads global problem definition.
   use typre
   use def_master
   use Mod_DistributedContainer
   use Mod_GeneralCase
   use Mod_DC_GeneralCase
   implicit none
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()

   !Loop through cases and update the channels
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call myCase%UpdateChannels
      
      call CaseList%GetNext(myDC)
   enddo
   if(MPIrank==MPIroot) write(*,*) '----: Channels updated'
   if(MPIrank==MPIroot) write(*,*) '----: Starting solution'
   if(MPIrank==MPIroot) write(*,*) '--------- For more information see .log files ...'
end subroutine

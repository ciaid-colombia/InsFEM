subroutine case_UpdateChannels(a)
   use Mod_GeneralCase
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   implicit none
   class(GeneralCase), target :: a
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
   
   call a%DriverList%GetFirst(myDC)
   do while(associated(myDC))
      call ExtractDriver(myDC,myDriver)
      
      call myDriver%UpdateChannels
   
      call a%DriverList%GetNext(myDC)
   enddo
end subroutine
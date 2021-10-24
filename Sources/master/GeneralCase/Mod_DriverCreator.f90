module Mod_DriverCreator
   use Mod_DriverInterface
   implicit none
   
   type :: DriverCreatorType
      character(5) :: key
      character(5) :: endword
      character(5) :: nick
      class(DriverInterface), pointer :: Driver => NULL()
   end type

end module

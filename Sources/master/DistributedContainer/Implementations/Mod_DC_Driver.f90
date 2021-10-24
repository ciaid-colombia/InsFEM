module Mod_DC_Driver
   use typre
   use Mod_DistributedContainer
   use Mod_DriverInterface
   implicit none   
   private
   
   public DC_Driver, ExtractDriver, DC_Driver_Const
   
   !The extended types
   type, extends(DistributedContainer) :: DC_Driver
      class(DriverInterface), pointer :: a => NULL()

contains
      procedure :: GetValue
   end type
   
   interface DC_Driver_Const
        procedure constructor
    end interface DC_Driver_Const
    
contains 
   subroutine GetValue(this,value)
      class(DC_Driver) :: this
      class(DriverInterface), pointer :: value

      value => this%a
   end subroutine GetValue

   function constructor(value)
     class(DC_Driver), pointer :: constructor
     class(DriverInterface), target :: value
     
     allocate(constructor)
     constructor%a => value

   end function constructor
   
   subroutine ExtractDriver(myDC,myDriver)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver
      
      select type (myDC)
      type is (DC_Driver)
      
      call myDC%GetValue(myDriver)
      end select
   end subroutine
      



end module

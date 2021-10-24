module Mod_DC_InChannel
   use typre
   use Mod_DistributedContainer
   use Mod_InChannel
   implicit none   
   private
   
   public DC_InChannel, ExtractInChannel, DC_InChannel_const
   
   !The extended types
   type, extends(DistributedContainer) :: DC_InChannel
      type(InChannel), pointer :: a => NULL()

contains
      procedure :: GetValue
   end type
   
   interface DC_InChannel_Const
        procedure constructor
    end interface DC_InChannel_Const
    
contains 
   subroutine GetValue(this,value)
      class(DC_InChannel) :: this
      type(InChannel), pointer :: value

      value => this%a
   end subroutine GetValue

   function constructor(value)
     class(DC_InChannel), pointer :: constructor
     type(InChannel), target :: value
     
     allocate(constructor)
     constructor%a => value

   end function constructor
   
   subroutine ExtractInChannel(myDC,myCase)
      class(DistributedContainer), pointer :: myDC
      type(InChannel), pointer :: myCase
      
      select type (myDC)
      type is (DC_InChannel)
      
      call myDC%GetValue(myCase)
      end select
   end subroutine



end module

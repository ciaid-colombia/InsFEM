module Mod_DC_GeneralCase
   use typre
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   use Mod_GeneralCase
   use Mod_DCHashCharSize
   implicit none   
   private
   
   public DC_GeneralCase, ExtractGeneralCase, DC_GeneralCase_const,GeneralCaseToList
   
   !The extended types
   type, extends(DistributedContainer) :: DC_GeneralCase
      type(GeneralCase), pointer :: a => NULL()

contains
      procedure :: GetValue
      procedure :: SetValue
   end type
   
   interface DC_GeneralCase_Const
        procedure constructor
    end interface DC_GeneralCase_Const
    
contains 
   subroutine GetValue(this,value)
      class(DC_GeneralCase) :: this
      type(GeneralCase), pointer :: value

      value => this%a
   end subroutine GetValue
   
   subroutine SetValue(this,value)
      class(DC_GeneralCase) :: this
      type(GeneralCase), target :: value

      this%a => value 
   end subroutine SetValue

   function constructor(value)
     class(DC_GeneralCase), pointer :: constructor
     type(GeneralCase), target :: value
     
     allocate(constructor)
     constructor%a => value

   end function constructor
   
   subroutine ExtractGeneralCase(myDC,myCase)
      class(DistributedContainer), pointer :: myDC
      type(GeneralCase), pointer :: myCase
      
      select type (myDC)
      type is (DC_GeneralCase)
      
      call myDC%GetValue(myCase)
      end select
   end subroutine
   
   subroutine InsertGeneralCase(myDC,myCase)
      class(DistributedContainer), pointer :: myDC
      type(GeneralCase) :: myCase
      
      select type (myDC)
      type is (DC_GeneralCase)
      
      call myDC%SetValue(myCase)
      end select
   end subroutine
   
   subroutine GeneralCaseToList(myCase,akey,task,mylist)
      type(GeneralCase) :: myCase
      character(*) :: akey, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      character(DCHashCharSize) :: key
      
      key = akey
      
      if (task == 'ADD') then
         myDC => DC_GeneralCase_Const(myCase)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_GeneralCase) 
            call myDC%SetValue(myCase)
         end select
      endif
   end subroutine
      
      



end module

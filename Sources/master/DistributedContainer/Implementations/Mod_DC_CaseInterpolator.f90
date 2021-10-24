module Mod_DC_CaseInterpolator
   use typre
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   use Mod_CaseInterpolator
   use Mod_DCHashCharSize
   implicit none   
   private
   
   public DC_CaseInterpolator, ExtractCaseInterpolator, DC_CaseInterpolator_const,CaseInterpolatorToList
   
   !The extended types
   type, extends(DistributedContainer) :: DC_CaseInterpolator
      type(CaseInterpolator), pointer :: a => NULL()

contains
      procedure :: GetValue
      procedure :: SetValue
   end type
   
   interface DC_CaseInterpolator_Const
        procedure constructor
    end interface DC_CaseInterpolator_Const
    
contains 
   subroutine GetValue(this,value)
      class(DC_CaseInterpolator) :: this
      type(CaseInterpolator), pointer :: value

      value => this%a
   end subroutine GetValue
   
   subroutine SetValue(this,value)
      class(DC_CaseInterpolator) :: this
      type(CaseInterpolator), target :: value

      this%a => value 
   end subroutine SetValue

   function constructor(value)
     class(DC_CaseInterpolator), pointer :: constructor
     type(CaseInterpolator), target :: value
     
     allocate(constructor)
     constructor%a => value

   end function constructor
   
   subroutine ExtractCaseInterpolator(myDC,myCMesh)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator), pointer :: myCMesh
      
      select type (myDC)
      type is (DC_CaseInterpolator)
      
      call myDC%GetValue(myCMesh)
      end select
   end subroutine
   
   subroutine InsertCaseInterpolator(myDC,myCMesh)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator) :: myCMesh
      
      select type (myDC)
      type is (DC_CaseInterpolator)
      
      call myDC%SetValue(myCMesh)
      end select
   end subroutine
   
   subroutine CaseInterpolatorToList(myCMesh,key,task,mylist)
      type(CaseInterpolator) :: myCMesh
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_CaseInterpolator_Const(myCMesh)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_CaseInterpolator) 
            call myDC%SetValue(myCMesh)
         end select
      endif
   end subroutine
   
end module

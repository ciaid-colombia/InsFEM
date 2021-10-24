module Mod_DC_CutMesh
   use typre
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   use Mod_CutMesh
   use Mod_DCHashCharSize
   implicit none   
   private
   
   public DC_CutMesh, ExtractCutMesh, DC_CutMesh_const,CutMeshToList
   
   !The extended types
   type, extends(DistributedContainer) :: DC_CutMesh
      type(CutMesh), pointer :: a => NULL()

contains
      procedure :: GetValue
      procedure :: SetValue
   end type
   
   interface DC_CutMesh_Const
        procedure constructor
    end interface DC_CutMesh_Const
    
contains 
   subroutine GetValue(this,value)
      class(DC_CutMesh) :: this
      type(CutMesh), pointer :: value

      value => this%a
   end subroutine GetValue
   
   subroutine SetValue(this,value)
      class(DC_CutMesh) :: this
      type(CutMesh), target :: value

      this%a => value 
   end subroutine SetValue

   function constructor(value)
     class(DC_CutMesh), pointer :: constructor
     type(CutMesh), target :: value
     
     allocate(constructor)
     constructor%a => value

   end function constructor
   
   subroutine ExtractCutMesh(myDC,myCMesh)
      class(DistributedContainer), pointer :: myDC
      type(CutMesh), pointer :: myCMesh
      
      select type (myDC)
      type is (DC_CutMesh)
      
      call myDC%GetValue(myCMesh)
      end select
   end subroutine
   
   subroutine InsertCutMesh(myDC,myCMesh)
      class(DistributedContainer), pointer :: myDC
      type(CutMesh) :: myCMesh
      
      select type (myDC)
      type is (DC_CutMesh)
      
      call myDC%SetValue(myCMesh)
      end select
   end subroutine
   
   subroutine CutMeshToList(myCMesh,key,task,mylist)
      type(CutMesh) :: myCMesh
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_CutMesh_Const(myCMesh)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_CutMesh) 
            call myDC%SetValue(myCMesh)
         end select
      endif
   end subroutine
      
      



end module

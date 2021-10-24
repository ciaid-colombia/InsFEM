module Mod_DC_rp
   use typre
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   implicit none   
   private
   
   public DC_rp, ExtractRp, DC_rp_const,RpToList
   
   !The extended types
   type, extends(DistributedContainer) :: DC_rp
      real(rp), pointer :: a3(:,:,:)=>NULL()
      real(rp), pointer :: a2(:,:)=>NULL()
      real(rp), pointer :: a1(:)=>NULL()
      real(rp), pointer :: a0=>NULL()
      
      type(r3p), pointer :: a_r3p(:)=>NULL()
      type(r2p), pointer :: a_r2p(:)=>NULL()

contains
      procedure :: GetValue0
      procedure :: GetValue1
      procedure :: GetValue2
      procedure :: GetValue3
      procedure :: GetValue_r3p
      procedure :: GetValue_r2p
      
      procedure :: SetValue0
      procedure :: SetValue1
      procedure :: SetValue2
      procedure :: SetValue3
      procedure :: SetValue_r3p
      procedure :: SetValue_r2p
      
      
      generic :: GetValue => GetValue0, GetValue1, GetValue2, GetValue3, GetValue_r3p, GetValue_r2p
      
      generic :: SetValue => SetValue0, SetValue1, SetValue2, SetValue3, SetValue_r3p, SetValue_r2p
   end type
   
   interface DC_rp_const
        procedure constructor1, constructor2, constructor3, constructor0, constructor_r3p, constructor_r2p
        
   end interface 
   
   interface ExtractRp
      procedure ExtractRp0, ExtractRp1, ExtractRp2, ExtractRp3, ExtractRp_r3p, ExtractRp_r2p
   end interface
   
   interface RPToList
      procedure RPToList0, RPToList1, RPToList2, RPToList3, RPToList_r3p, RPToList_r2p
   end interface
    
contains 

   !Gets
   subroutine GetValue0(this,value)
      class(DC_rp) :: this
      real(rp), pointer :: value
      
      value => this%a0
   end subroutine
   
   subroutine GetValue1(this,value)
      class(DC_rp) :: this
      real(rp), pointer :: value(:)
      
      value => this%a1
   end subroutine
   
   subroutine GetValue2(this,value)
      class(DC_rp) :: this
      real(rp), pointer :: value(:,:)
      
      value => this%a2
   end subroutine
    
   subroutine GetValue3(this,value)
      class(DC_rp) :: this
      real(rp), pointer :: value(:,:,:)
      
      value => this%a3
   end subroutine 
   
   subroutine GetValue_r3p(this,value)
      class(DC_rp) :: this
      type(r3p), pointer :: value(:)
      
      value => this%a_r3p
   end subroutine 

   subroutine GetValue_r2p(this,value)
      class(DC_rp) :: this
      type(r2p), pointer :: value(:)

      value => this%a_r2p
   end subroutine 
   
    !Sets
   subroutine SetValue0(this,value)
      class(DC_rp) :: this
      real(rp), target :: value

      this%a0 => value
   end subroutine
   
   subroutine SetValue1(this,value)
      class(DC_rp) :: this
      real(rp), target :: value(:)

      this%a1 => value
   end subroutine
   
   subroutine SetValue2(this,value)
      class(DC_rp) :: this
      real(rp), target :: value(:,:)

      this%a2 => value
   end subroutine
    
   subroutine SetValue3(this,value)
      class(DC_rp) :: this
      real(rp), target :: value(:,:,:)

      this%a3 => value
   end subroutine 
   
   subroutine SetValue_r3p(this,value)
      class(DC_rp) :: this
      type(r3p), target :: value(:)

      this%a_r3p => value
   end subroutine 

   subroutine SetValue_r2p(this,value)
      class(DC_rp) :: this
      type(r2p), target :: value(:)

      this%a_r2p => value
   end subroutine 

   !Constructors
   function constructor0(value)
     class(DC_rp), pointer :: constructor0
     real(rp), target :: value
     
     allocate(constructor0)
     constructor0%a0 => value

   end function

   function constructor1(value)
     class(DC_rp), pointer :: constructor1
     real(rp), target :: value(:)
     
     allocate(constructor1)
     constructor1%a1 => value

   end function

   function constructor2(value)
     class(DC_rp), pointer :: constructor2
     real(rp), target :: value(:,:)
     
     allocate(constructor2)
     constructor2%a2 => value

   end function

   function constructor3(value)
     class(DC_rp), pointer :: constructor3
     real(rp), target :: value(:,:,:)
     
     allocate(constructor3)
     constructor3%a3 => value
   end function

   function constructor_r3p(value)
     class(DC_rp), pointer :: constructor_r3p
     type(r3p), target :: value(:)
     
     allocate(constructor_r3p)
     constructor_r3p%a_r3p => value
   end function
   
   function constructor_r2p(value)
     class(DC_rp), pointer :: constructor_r2p
     type(r2p), target :: value(:)
     
     allocate(constructor_r2p)
     constructor_r2p%a_r2p => value
   end function
   
   !Extractors
   subroutine ExtractRp0(myDC,myRp)
      class(DistributedContainer), pointer :: myDC
      real(rp), pointer :: myRp
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine

   subroutine ExtractRp1(myDC,myRp)
      class(DistributedContainer), pointer :: myDC
      real(rp), pointer :: myRp(:)
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine
   
   subroutine ExtractRp2(myDC,myRp)
      class(DistributedContainer), pointer :: myDC 
      real(rp), pointer :: myRp(:,:)
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine
   
   subroutine ExtractRp3(myDC,myRp)
      class(DistributedContainer), pointer :: myDC
      real(rp), pointer :: myRp(:,:,:)
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine
   
   subroutine ExtractRp_r3p(myDC,myRp)
      class(DistributedContainer), pointer :: myDC
      type(r3p), pointer :: myRp(:)
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine
   
   subroutine ExtractRp_r2p(myDC,myRp)
      class(DistributedContainer), pointer :: myDC
      type(r2p), pointer :: myRp(:)
      
      select type (myDC)
      type is (DC_rp)
      
      call myDC%GetValue(myRp)
      end select
   end subroutine
   
   !RpToList
   subroutine RPToList0(value,key,task,mylist)
      real(rp) :: value
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine RPToList1(value,key,task,mylist)
      real(rp) :: value(:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine RPToList2(value,key,task,mylist)
      real(rp) :: value(:,:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine RPToList3(value,key,task,mylist)
      real(rp) :: value(:,:,:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine RPToList_r2p(value,key,task,mylist)
      type(r2p) :: value(:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine RPToList_r3p(value,key,task,mylist)
      type(r3p) :: value(:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_RP_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_RP) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
  


end module

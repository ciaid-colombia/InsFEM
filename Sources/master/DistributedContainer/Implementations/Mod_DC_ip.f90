module Mod_DC_ip
   use typre
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   implicit none   
   private
   
   public DC_ip, Extractip, DC_ip_Const,ipToList
   
   !The extended types
   type, extends(DistributedContainer) :: DC_ip
      integer(ip), pointer :: a3(:,:,:) => NULL()
      integer(ip), pointer :: a2(:,:) => NULL()
      integer(ip), pointer :: a1(:) => NULL()
      integer(ip), pointer :: a0 => NULL()
      
      type(i2p), pointer :: a_i2p(:) => NULL()

contains
      procedure :: GetValue0
      procedure :: GetValue1
      procedure :: GetValue2
      procedure :: GetValue3
      procedure :: GetValue_i2p
      
      procedure :: SetValue0
      procedure :: SetValue1
      procedure :: SetValue2
      procedure :: SetValue3
      procedure :: SetValue_i2p
      
      
      generic :: GetValue => GetValue0, GetValue1, GetValue2, GetValue3, GetValue_i2p
      
      generic :: SetValue => SetValue0, SetValue1, SetValue2, SetValue3, SetValue_i2p
   end type
   
   interface DC_ip_const
        procedure constructor1, constructor2, constructor3, constructor0, constructor_i2p
        
   end interface 
   
   interface Extractip
      procedure Extractip0, Extractip1, Extractip2, Extractip3,  Extractip_i2p
   end interface
   
   interface ipToList
      procedure ipToList0, ipToList1, ipToList2, ipToList3,  ipToList_i2p
   end interface
    
contains 

   !Gets
   subroutine GetValue0(this,value)
      class(DC_ip) :: this
      integer(ip), pointer :: value
      
      value => this%a0
   end subroutine
   
   subroutine GetValue1(this,value)
      class(DC_ip) :: this
      integer(ip), pointer :: value(:)
      
      value => this%a1
   end subroutine
   
   subroutine GetValue2(this,value)
      class(DC_ip) :: this
      integer(ip), pointer :: value(:,:)
      
      value => this%a2
   end subroutine
    
   subroutine GetValue3(this,value)
      class(DC_ip) :: this
      integer(ip), pointer :: value(:,:,:)
      
      value => this%a3
   end subroutine 
   
   subroutine GetValue_i2p(this,value)
      class(DC_ip) :: this
      type(i2p), pointer :: value(:)
      
      value => this%a_i2p
   end subroutine 
   
    !Sets
   subroutine SetValue0(this,value)
      class(DC_ip) :: this
      integer(ip), target :: value

      this%a0 => value
   end subroutine
   
   subroutine SetValue1(this,value)
      class(DC_ip) :: this
      integer(ip), target :: value(:)

      this%a1 => value
   end subroutine
   
   subroutine SetValue2(this,value)
      class(DC_ip) :: this
      integer(ip), target :: value(:,:)

      this%a2 => value
   end subroutine
    
   subroutine SetValue3(this,value)
      class(DC_ip) :: this
      integer(ip), target :: value(:,:,:)

      this%a3 => value
   end subroutine 
   
   subroutine SetValue_i2p(this,value)
      class(DC_ip) :: this
      type(i2p), target :: value(:)

      this%a_i2p => value
   end subroutine 

   !Constructors
   function constructor0(value)
     class(DC_ip), pointer :: constructor0
     integer(ip), target :: value
     
     allocate(constructor0)
     constructor0%a0 => value

   end function

   function constructor1(value)
     class(DC_ip), pointer :: constructor1
     integer(ip), target :: value(:)
     
     allocate(constructor1)
     constructor1%a1 => value

   end function

   function constructor2(value)
     class(DC_ip), pointer :: constructor2
     integer(ip), target :: value(:,:)
     
     allocate(constructor2)
     constructor2%a2 => value

   end function

   function constructor3(value)
     class(DC_ip), pointer :: constructor3
     integer(ip), target :: value(:,:,:)
     
     allocate(constructor3)
     constructor3%a3 => value
   end function

   function constructor_i2p(value)
     class(DC_ip), pointer :: constructor_i2p
     type(i2p), target :: value(:)
     
     allocate(constructor_i2p)
     constructor_i2p%a_i2p => value
   end function
   
   !Extractors
   subroutine Extractip0(myDC,myip)
      class(DistributedContainer), pointer :: myDC
      integer(ip), pointer :: myip
      
      select type (myDC)
      type is (DC_ip)
      
      call myDC%GetValue(myip)
      end select
   end subroutine

   subroutine Extractip1(myDC,myip)
      class(DistributedContainer), pointer :: myDC
      integer(ip), pointer :: myip(:)
      
      select type (myDC)
      type is (DC_ip)
      
      call myDC%GetValue(myip)
      end select
   end subroutine
   
   subroutine Extractip2(myDC,myip)
      class(DistributedContainer), pointer :: myDC 
      integer(ip), pointer :: myip(:,:)
      
      select type (myDC)
      type is (DC_ip)
      
      call myDC%GetValue(myip)
      end select
   end subroutine
   
   subroutine Extractip3(myDC,myip)
      class(DistributedContainer), pointer :: myDC
      integer(ip), pointer :: myip(:,:,:)
      
      select type (myDC)
      type is (DC_ip)
      
      call myDC%GetValue(myip)
      end select
   end subroutine
   
   subroutine Extractip_i2p(myDC,myip)
      class(DistributedContainer), pointer :: myDC
      type(i2p), pointer :: myip(:)
      
      select type (myDC)
      type is (DC_ip)
      
      call myDC%GetValue(myip)
      end select
   end subroutine
   
   !ipToList
   subroutine ipToList0(value,key,task,mylist)
      integer(ip) :: value
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_ip_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_ip) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine ipToList1(value,key,task,mylist)
      integer(ip) :: value(:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_ip_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_ip) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine ipToList2(value,key,task,mylist)
      integer(ip) :: value(:,:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_ip_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_ip) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine ipToList3(value,key,task,mylist)
      integer(ip) :: value(:,:,:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_ip_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_ip) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   subroutine ipToList_i2p(value,key,task,mylist)
      type(i2p) :: value(:)
      character(*) :: key, task
      type(DistributedContainerList) :: myList
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      if (task == 'ADD') then
         myDC => DC_ip_Const(value)
         call myDC%SetKey(key)
         call mylist%Add(myDC)
      elseif (task == 'UPDATE') then
         call mylist%GetFromKey(key,myDC)
         if (.not. associated(mydC)) call runend('DCToList, update, DC not found')
         select type (myDC)
         type is (DC_ip) 
            call myDC%SetValue(value)
         end select
      endif
   end subroutine
   
   


end module

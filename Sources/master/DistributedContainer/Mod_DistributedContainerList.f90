module Mod_DistributedContainerList
    use typre
    use Mod_DCHashCharSize
    use Mod_DistributedContainer
    use Mod_DChash
    implicit none

    private

    public :: DistributedContainerList
    
    type :: DCListElement
      class(DistributedContainer), pointer :: DC => NULL()
      class(DCListElement), pointer :: next => null()
    end type
    
    interface DCListElement_const
      procedure constructor_DCListElement
    end interface


    type :: DistributedContainerList

        private

        class(DCListElement), pointer :: first_eDC => null()
        class(DCListElement), pointer :: last_eDC => null()
        class(DCListElement), pointer :: curr_eDC=> null()  ! list iterator
        integer ::  numDC = 0
        type(DChash) :: hash
        
        logical :: isInit = .false.

    contains

        procedure , non_overridable :: Initialize
        procedure , non_overridable :: Add
        procedure , non_overridable :: GoToFirst        ! Rewind DistributedContainerList iterator
        procedure , non_overridable :: GetNext
        procedure , non_overridable :: GetFirst
        procedure , non_overridable :: GetNumberOfDC
        procedure , non_overridable :: GetFromKey
        procedure , non_overridable :: Finalize
        procedure , non_overridable :: LoopList

    end type DistributedContainerList

  
  

contains

   function constructor_DCListElement(DC)
      class(DCListElement), pointer :: constructor_DCListElement
      class(DistributedContainer), target :: DC
      
      allocate(constructor_DCListElement)
      constructor_DCListElement%DC => DC
   end function

   subroutine Initialize(a)
      class(DistributedContainerList) :: a
      
      a%first_eDC => null()
      a%last_eDC  => null()
      a%curr_eDC  => null()
      
      !Initialize the hash table
      call a%hash%Init
      
      a%isInit = .true.
      
   end subroutine

   subroutine Add(a,DC)
     class(DistributedContainerList) :: a
     class(DistributedContainer), target :: DC
     
     class(DCListElement), pointer :: eDC => NULL()
     
     if (a%isInit .eqv. .false.) call runend('DCList not initialized')
     
     eDC => DCListElement_Const(DC)
     
     a%numDC = a%numDC + 1
     if (.not. associated(a%first_eDC)) then
         a%first_eDC => eDC
     else
         a%last_eDC%next => eDC
     end if
     a%last_eDC => eDC
     
     if (DC%IsSetKey) then
       call a%hash%put(DC%key,DC) 
       if (DC%IsSetNick) then
           call a%hash%put(DC%nick,DC) 
       endif
   endif

   end subroutine Add

   subroutine GetNumberOfDC(a, numberDC)
      class(DistributedContainerList) :: a
      integer :: numberDC
      if (a%isInit .eqv. .false.) call runend('DCList not initialized')
      
      NumberDC = a%numDC
   end subroutine

   subroutine GoToFirst(a)
      class(DistributedContainerList) :: a
   
      if (a%isInit .eqv. .false.) call runend('DCList not initialized')
      a%curr_eDC => a%first_eDC
   end subroutine GoToFirst

   subroutine GetNext(a,Next)
      class(DistributedContainerList) :: a
      class(DistributedContainer), pointer :: Next
      
      if (a%isInit .eqv. .false.) call runend('DCList not initialized')
      if (associated(a%curr_eDC)) then
         Next => a%curr_eDC%DC   
      else
         Next => null()
      endif
      if (associated(a%curr_eDC)) then
         a%curr_eDC => a%curr_eDC%next
      endif

   end subroutine
   
   subroutine GetFirst(a,Next)
      class(DistributedContainerList) :: a
      class(DistributedContainer), pointer :: Next
      
      call a%GoToFirst
      call a%GetNext(Next)
   end subroutine   
   
   
   
   subroutine GetFromKey(a,key,DC)
      use typre
      class(DistributedContainerList) :: a
      class(DistributedContainer), pointer :: DC
      character(*) :: key
      
      character(DCHashCharSize) :: akey
      
      if (a%isInit .eqv. .false.) call runend('DCList not initialized')      
      akey = key

      
      call a%hash%get(akey,DC) 
   end subroutine
   
   
   
   subroutine Finalize(a, eonlylist)
      class(DistributedContainerList) :: a
      logical, optional :: eonlylist
      
      logical :: onlylist
      
      !if onlylist is true, it kills only the list, not the containers themself
      
      class(DCListElement), pointer :: aux_eDC => NULL()
      
      if (.not. present(eonlylist)) then
          onlylist = .false.
      else
         onlylist = eonlylist
      endif

      if (a%isInit .eqv. .false.) call runend('DCList not initialized')
      if (associated(a%first_eDC)) then
         a%curr_eDC => a%first_eDC
         do while (associated(a%curr_eDC))

            !Set auxpointer
            aux_eDC => a%curr_eDC
            
            !Deallocate the containter if not onlylist
            if (onlylist .eqv. .false.) then
               if (associated(a%curr_eDC%DC)) deallocate(a%curr_eDC%DC)
            endif

            !Move on if possible
            if (associated(a%curr_eDC%next)) then
               a%curr_eDC => a%curr_eDC%next  
            else
               a%curr_eDC => null()
            endif         
            
            !deallocate
            deallocate(aux_eDC)
            
            
         enddo
      endif
      
      !Finalize the hash table
      call a%hash%free
      
      a%isInit = .false.
   end subroutine
   
   !Additional useful subroutines
   subroutine AddDCLists(list1,list2)
       class(DistributedContainerList) :: list1
       class(DistributedContainerList) :: list2
       
       class(DistributedContainer), pointer :: theirsDC => NULL()
       
       call list2%GoToFirst
       call list2%GetNext(theirsDC)
       do while (associated(theirsDC))
         call list1%Add(theirsDC)
         call list2%GetNext(theirsDC)
       enddo
   end subroutine
   
   subroutine LoopList(a,ExternalSubroutine,DoFirstKeys)
      class(DistributedContainerList) :: a
      external :: ExternalSubroutine
      character(DCHashCharSize), optional :: DoFirstKeys(:)
      
      interface
         subroutine ExternalSubroutine(myDC)
            use Mod_DistributedContainer
            implicit none
            class(DistributedContainer), pointer :: myDC
         end subroutine
      end interface
      
      class(DistributedContainer), pointer :: myDC
      integer(ip) :: i
      logical :: dokey
      character(DCHashCharSize) :: key
      
      !If a list of do first keys is passed, then do these keys first
      if (present(DoFirstKeys)) then
         do i = 1,size(DoFirstKeys)
            call a%GetFromKey(DoFirstKeys(i),myDC)
            !Do what needs to be done
            if (associated(myDC)) call ExternalSubroutine(myDC)
         enddo
      endif
            
      dokey = .true.
      call a%GoToFirst
      call a%GetNext(myDC)
      do while (associated(myDC))
         
         !If the list of dokeys list is present, ensure that we do not do them twice
         if (present(DoFirstKeys)) then
            dokey = .true.
            call myDC%GetKey(key)
            do i = 1,size(DoFirstKeys)
               if (DoFirstKeys(i) == key) dokey = .false.
            enddo
         endif
         
         !Do what needs to be done      
         if (dokey .eqv. .true.) call ExternalSubroutine(myDC)
         
         call a%GetNext(myDC)
      enddo
   end subroutine
    
      

   

end module Mod_DistributedContainerList

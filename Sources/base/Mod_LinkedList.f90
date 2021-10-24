module Mod_LinkedList
   use Mod_Memor
   use typre
   private
   public LinkedListHeader, LinkedListStorage
   
   type :: LinkedListHeader
      integer(ip) :: first = -1
      integer(ip) :: last = -1
      integer(ip) :: nelem = 0

contains
      procedure :: Initialize => HeaderInitialize
      procedure :: GetNelem
      
   end type
   
   type :: LinkedListStorage
      integer(ip) :: maxnelem
      
      integer(ip) :: CanRemove
      integer(ip) :: CanRepeat
      
      !next is used for 2 things: when the element is in the list, it 
      !stores the next element in the list
      !when the element is free, it stores the next free space
      integer(ip), allocatable :: next(:)
      integer(ip), allocatable :: prev(:)
      integer(ip), allocatable :: value(:)
      
      integer(ip) :: nextfree   !If we free an element we add it to the list of free elements
            
      !For moving through the list
      integer(ip) :: nextItem
      
      type(MemoryMan), pointer :: Memor => NULL()
   
      !Depending on the options these point to a subroutine or another
      procedure(GetNextInterface), pointer :: GetNext => NULL()
      procedure(GetValueInterface), private, pointer :: GetValue => NULL()
      procedure(RemoveInterface), pointer :: Remove => NULL()
      procedure(RemoveAllInterface), pointer :: RemoveAll => NULL()
   
contains

      procedure :: Initialize
      procedure :: Add
      procedure :: GoToFirstOfList
      procedure :: ListToArray
      procedure :: Dealloc

      
   
   end type
   
   abstract interface
      subroutine GetNextInterface(a,value)
         use typre
         import LinkedListStorage
         implicit none
         class(LinkedListStorage) :: a
         integer(ip) :: value
         
      end subroutine
      
      subroutine GetValueInterface(a,ipos,value)
         use typre
         import LinkedListStorage
         implicit none
         class(LinkedListStorage) :: a
         integer(ip) :: ipos,value
         
      end subroutine
      
      subroutine RemoveInterface(a,Header,ielem,ierr)
         use typre
         import LinkedListStorage
         import LinkedListHeader
         implicit none
         class(LinkedListStorage) :: a
         type(LinkedListHeader) :: Header
         integer(ip) :: ielem,ierr
      end subroutine
      
      subroutine RemoveAllInterface(a,Header)
         use typre
         import LinkedListStorage
         import LinkedListHeader
         implicit none
         class(LinkedListStorage) :: a
         type(LinkedListHeader) :: Header
      end subroutine
      
   end interface
  
contains

   subroutine Initialize(a,maxnelem,CanRemove,CanRepeat,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(LinkedListStorage) :: a
      integer(ip) :: maxnelem
      character(6) :: CanRemove, CanRepeat
      type(MemoryMan), target :: Memor
      !CanRemove values: 'Removeable' or 'NonRemoveable'
      !CanRepeat values: 'Repeatable' or 'NonRepeatable'
      
      integer(ip) :: ielem
      
      a%maxnelem = maxnelem
      
      a%Memor => Memor

      call a%Memor%alloc(maxnelem,a%next,'next','LinkedList')
      
      !The first free element is 1
      a%nextfree = 1
      
      !CanRepeat
      if (CanRepeat == 'Repeat') then
         a%CanRepeat = 1 !Elements can belong to more than one list
      else
         a%CanRepeat = 0
      endif
      if (a%CanRepeat == 1) then
         call a%Memor%alloc(maxnelem,a%value,'value','LinkedList')
         a%GetNext => GetNextCanRepeat
         a%GetValue => GetValueCanRepeat
         
         !Default is the next position availabe is ipos + 1
         forall (ielem = 1:maxnelem)
            a%next(ielem) = ielem+1
         end forall
      else  
         a%GetNext => GetNextNoRepeat
         a%GetValue => GetValueNoRepeat
         
         !Elements are in no list
         a%next = 0
      endif
      
      !CanRemove
      if (CanRemove == 'Remova') then
         a%CanRemove = 1 !We can Remove elements
      else
         a%CanRemove = 0
      endif
      if (a%CanRemove == 1) then
         call a%Memor%alloc(maxnelem,a%prev,'prev','LinkedList')
         if (a%CanRepeat == 1) then
            a%Remove => RemoveCanRemove
         else
            a%Remove => RemoveFromIpos
         endif
         a%RemoveAll => RemoveAll
      else
         a%Remove => RemoveNoRemove
         a%RemoveAll => RemoveAllNoRemove
      endif
         
      
   end subroutine
   
   subroutine Dealloc(a)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      
      !Memory deallocations
      
      call a%Memor%dealloc(a%maxnelem,a%next,'next','LinkedList')
      if (a%CanRepeat == 1) call a%Memor%dealloc(a%maxnelem,a%value,'value','LinkedList')
      if (a%CanRemove == 1) call a%Memor%dealloc(a%maxnelem,a%prev,'prev','LinkedList')
   end subroutine
   
   subroutine Add(a,Header,ielem)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      integer(ip) :: ielem
      type(LinkedListHeader) :: Header
      
      integer(ip) :: ipos
      
      integer(ip) :: new_maxnelem
      integer(ip), allocatable :: aux_next(:), aux_prev(:), aux_value(:)
      
      !a%nextfree contains the next free space where I am going to add the element
      if (a%CanRepeat == 1) then
         ipos = a%nextfree
      else
         ipos = ielem
      endif
      
      !Check if the total number of elements is too large, reallocate if necessary
      if (ipos > a%maxnelem) then
         if (a%CanRepeat == 1) then
            new_maxnelem = 2_ip*a%maxnelem
         else
            new_maxnelem = ipos
         endif
         
         call a%Memor%realloc(new_maxnelem,a%next,'next','AddLinkedList')
         if (a%CanRemove == 1) call a%Memor%realloc(new_maxnelem,a%prev,'prev','AddLinkedList')
         if (a%CanRepeat == 1) then
            call a%Memor%realloc(new_maxnelem,a%value,'value','AddLinkedList')
         
            !Default is the next position availabe is ipos + 1
            forall (ielem = a%maxnelem+1:new_maxnelem)
               a%next(ielem) = ielem+1
            end forall
         else
            !New elements belong to no list
            a%next(a%maxnelem+1:new_maxnelem) = 0
         endif
         
         a%maxnelem = new_maxnelem
      endif
         
      
      !Add
      if (a%CanRepeat == 1) then
         a%value(ipos) = ielem
      else
         if (a%next(ipos) /= 0) call runend('LinkedList: Element Already assigned to a list, and list is of non-repeateable type')
      endif
      if (a%CanRemove == 1) a%prev(ipos) = Header%last
      
      if (Header%last /= -1) then
         a%next(Header%last) = ipos
      else
         Header%first = ipos
      endif
      Header%last = ipos
      
      !a%next has the next free position (when the element was in no list)
      if (a%CanRepeat == 1) a%nextfree = a%next(ipos)
      a%next(ipos) = -1
      
      !The list has an additional element
      Header%nelem = Header%nelem + 1
      
   end subroutine
   
   subroutine RemoveFromIpos(a,Header,ipos,ierr)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      type(LinkedListHeader) :: Header
      integer(ip) :: ipos,ierr
      
      integer(ip) :: ifree,inext,iprev
      
      ierr = 1
      
      iprev = a%prev(ipos)
      inext = a%next(ipos)
      
      !We skip the element, set the nexts in the previous element
      if (iprev /= -1) then
         a%next(iprev) = inext
      else
         Header%first = inext
      endif
   
      !We skip the element, set the prevs in the next element
      if (inext /= -1) then
         a%prev(inext) = iprev
      else
         Header%last = iprev
      endif
      
      !Now we add the element to the list of free elements
      if (a%CanRepeat == 1) then
         ifree = a%nextfree
         a%next(ipos) = ifree
         a%nextfree = ipos
      else
         a%next(ipos) = 0
      endif   
      
      !The list has one fewer element
      Header%nelem = Header%nelem - 1
      
      ierr = 0
   end subroutine
   
   subroutine RemoveCanRemove(a,Header,ielem,ierr)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      type(LinkedListHeader) :: Header
      integer(ip) :: ielem,ierr
      
      integer(ip) :: ifree,ipos,inext,iprev,value
      
      ierr = 1
      
      !First we look for the element
      ipos = Header%first
      do while (ipos /= -1) 
         call a%GetValue(ipos,value)
         if (value == ielem) then
            !Now straightly remove
            call RemoveFromIpos(a,Header,ipos,ierr)
            
            !We are done
            ipos = -1
         
         !If it is not the next one, continue through the list
         else
            ipos = a%next(ipos)
         endif
         
      enddo
      
      if (ierr /= 0) then
         call runend('LinkedList: element to remove was not found')
      endif
      
   end subroutine
   
   subroutine RemoveNoRemove(a,Header,ielem,ierr)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      type(LinkedListHeader) :: Header
      integer(ip) :: ielem,ierr
      
      call runend('LinkedList: Tried to remove an element but the list is of non-removeable type')
   end subroutine   

   subroutine GoToFirstOfList(a,Header)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      class(LinkedListHeader) :: Header
      
      a%nextItem = Header%first
      
   end subroutine
   
   subroutine RemoveAllNoRemove(a,Header)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      type(LinkedListHeader) :: Header
   
      if (a%CanRemove == 0) call runend('LinkedList: Tried to remove all but the list is of non-removable type')
   end subroutine   
   
   subroutine RemoveAll(a,Header)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      type(LinkedListHeader) :: Header   
      integer(ip) :: ipos, ierr, aux_next
      
      !First we look for the element
      ipos = Header%first
      do while (ipos /= -1) 
         !Now straightly remove
         aux_next = a%next(ipos)
         call RemoveFromIpos(a,Header,ipos,ierr)
      
         ipos = aux_next
      enddo
   end subroutine
   
   subroutine GetNextCanRepeat(a,value)
      use typre
      class(LinkedListStorage) :: a
      integer(ip) :: value
      
      integer(ip) :: inext
      
      if (a%nextItem /= -1) then
         value = a%value(a%nextItem)
         a%nextItem = a%next(a%nextItem)
      else
         value = a%nextItem
      endif
   end subroutine
   
   subroutine GetNextNoRepeat(a,value)
      use typre
      class(LinkedListStorage) :: a
      integer(ip) :: value
      
      value = a%nextItem
      if (a%nextItem /= -1) then
         a%nextItem=a%next(a%nextItem)
      endif
   end subroutine
   
   subroutine GetValueCanRepeat(a,ipos,value)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      integer(ip) :: ipos,value
      
      value = a%value(ipos)
      
   end subroutine
   
   subroutine GetValueNoRepeat(a,ipos,value)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      integer(ip) :: ipos,value
      
      value = ipos
      
   end subroutine
  
   subroutine ListToArray(a,Header,array)
      use typre
      implicit none
      class(LinkedListStorage) :: a
      class(LinkedListHeader) :: Header
      integer(ip) :: array(*)
      
      integer(ip) :: nz,ivalue
      
      call a%GoToFirstOfList(Header)
      call a%GetNext(ivalue)
      
      nz = 0
      do while (ivalue /= -1)
         nz = nz+1
         array(nz) = ivalue
         
         call a%GetNext(ivalue)
      enddo
   end subroutine
  
   subroutine GetNelem(a,nelem)
      use typre
      class(LinkedListHeader) :: a
      integer(ip) :: nelem
      
      nelem = a%nelem
   end subroutine

     
   
   subroutine HeaderInitialize(a)
      use typre
      class(LinkedListHeader) :: a
      
      a%nelem = 0
      a%first = -1
      a%last = -1
   end subroutine


end module

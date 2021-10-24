subroutine LinkedListTest
   use typre
   use Mod_LinkedList
   use Mod_Memor
   implicit none
   
   type(LinkedListHeader), allocatable  :: MyListHeader(:)
   type(LinkedListStorage) :: MyListStorage
   type(MemoryMan)  :: Memor
   
   integer(ip) :: ielem, ilist,ierr
 
   
   
   write(*,*) 'Testing Removable, NonRepeateable'
   allocate(MyListHeader(10))
   call MyListStorage%Initialize(50,'Removable','NonRepeatable',Memor)
   call testlistAdd
   call testlistRemove   
   call MyListStorage%Dealloc
   deallocate(MyListHeader)
   
   write(*,*) 'Testing Removable, Repeateable'
   allocate(MyListHeader(10))
   call MyListStorage%Initialize(50,'Removable','Repeatable',Memor)
   call testlistAdd 
   call testlistRemove  
   call MyListStorage%Dealloc
   deallocate(MyListHeader)
   
   write(*,*) 'Testing NonRemovable, Repeateable'
   allocate(MyListHeader(10))
   call MyListStorage%Initialize(50,'NonRemovable','Repeatable',Memor)
   call testlistAdd  
   call testlistRemove 
   call MyListStorage%Dealloc
   deallocate(MyListHeader)
   
   write(*,*) 'Testing NonRemovable, NonRepeateable'
   allocate(MyListHeader(10))
   call MyListStorage%Initialize(50,'NonRemovable','NonRepeatable',Memor)
   call testlistAdd  
   call testlistRemove 
   call MyListStorage%Dealloc
   deallocate(MyListHeader)
   
   
   
   
contains 
   
   subroutine testlistAdd
   
      ilist = 0
      do ielem = 1,100
         ilist = ilist + 1
         if (ilist > 10) ilist = 1
         
         call MyListStorage%Add(MyListHeader(ilist),ielem)
      enddo
      
      do ilist = 1,10
         call MyListStorage%GoToFirstOfList(MyListHeader(ilist))
         
         call MyListStorage%GetNext(ielem)
         do while (ielem /= -1)
            write(*,*) ilist, ielem
         
         
            call MyListStorage%GetNext(ielem)
         enddo
      enddo
   end subroutine
   
   subroutine testlistRemove   
      
      do ielem = 20,33
         ilist = mod(ielem,10)
         if (ilist == 0) ilist = 10
         call MyListStorage%Remove(MyListHeader(ilist),ielem,ierr)
      enddo
      
      call MyListStorage%Remove(MyListHeader(1),11,ierr)
      call MyListStorage%Add(MyListHeader(10),11)
      call MyListStorage%Add(MyListHeader(3),11)
      call MyListStorage%RemoveAll(MyListHeader(7))

      do ilist = 1,10
         call MyListStorage%GoToFirstOfList(MyListHeader(ilist))
         
         call MyListStorage%GetNext(ielem)
         do while (ielem /= -1)
            write(*,*) ilist, ielem
         
         
            call MyListStorage%GetNext(ielem)
         enddo
      enddo
   
   end subroutine   
   
end subroutine

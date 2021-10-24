module Mod_HeapSortExternal
   use typre
   implicit none
   
   type AbstractHeapSort
      
contains
      procedure :: IsLowerThan
   
   end type
   
contains

   subroutine IsLowerThan(this,a,b,TestResult)
      class(AbstractHeapSort) :: this
      integer(ip) :: a,b
      logical TestResult
      
      call runend('Islowerthan not ready')
   
   
   end subroutine

   subroutine HeapSortExternalObject(n,a,AHeapSort)
      use typre
      implicit none
      integer(ip) :: n
      integer(ip), intent(in out) :: a(0:)
      class(AbstractHeapSort) :: AHeapSort

      integer :: start, bottom
      integer :: temp
   
      do start = (n - 2) / 2, 0, -1
      call siftdownExternalObject(a, start, n,AHeapSort);
      end do
   
      do bottom = n - 1, 1, -1
         temp = a(0)
         a(0) = a(bottom)
         a(bottom) = temp;
         call siftdownExternalObject(a, 0, bottom,AHeapSort)
      end do
   
   end subroutine heapsortExternalObject
   
   subroutine siftdownExternalObject(a, start, bottom, AHeapSort)
      use typre
      implicit none
   
      integer(ip), intent(in out) :: a(0:)
      class(AbstractHeapSort) :: AHeapSort
      integer(ip), intent(in) :: start, bottom
      integer(ip) :: child, root
      integer(ip) :: temp
      
      logical :: TestResult

      
      root = start
      do while(root*2 + 1 < bottom)
         child = root * 2 + 1
      
         if (child + 1 < bottom) then
            call AHeapSort%IsLowerThan(a(child),a(child+1),TestResult)
            if (TestResult) child = child + 1
         end if
      
         call AHeapSort%IsLowerThan(a(root),a(child),TestResult)
         if (TestResult) then
            temp = a(child)
            a(child) = a (root)
            a(root) = temp
            root = child
         else
            return
         end if  
      end do      
   
   end subroutine siftdownExternalObject

   subroutine HeapSortExternal(n,a,IsLowerThan)
      use typre
      implicit none
      integer(ip) :: n
      integer(ip), intent(in out) :: a(0:)
      logical, external :: IsLowerThan

      integer :: start, bottom
      integer :: temp
   
      do start = (n - 2) / 2, 0, -1
      call siftdownExternal(a, start, n,IsLowerThan);
      end do
   
      do bottom = n - 1, 1, -1
         temp = a(0)
         a(0) = a(bottom)
         a(bottom) = temp;
         call siftdownExternal(a, 0, bottom,IsLowerThan)
      end do
   
   end subroutine heapsortExternal
   
   subroutine siftdownExternal(a, start, bottom, IsLowerThan)
      use typre
      implicit none
   
      integer(ip), intent(in out) :: a(0:)
      logical, external :: IsLowerThan
      integer(ip), intent(in) :: start, bottom
      integer(ip) :: child, root
      integer(ip) :: temp
      
      logical :: TestResult

      
      root = start
      do while(root*2 + 1 < bottom)
         child = root * 2 + 1
      
         if (child + 1 < bottom) then
            TestResult = IsLowerThan(a(child),a(child+1))
            if (TestResult) child = child + 1
         end if
      
         TestResult = IsLowerThan(a(root),a(child))
         if (TestResult) then
            temp = a(child)
            a(child) = a (root)
            a(root) = temp
            root = child
         else
            return
         end if  
      end do      
   
   end subroutine siftdownExternal
   
   subroutine HeapSort(n,a)
      use typre
      implicit none
      integer(ip) :: n
      integer(ip), intent(in out) :: a(0:)
      logical, external :: IsLowerThan

      integer :: start, bottom
      integer :: temp
   
      do start = (n - 2) / 2, 0, -1
      call siftdown(a, start, n);
      end do
   
      do bottom = n - 1, 1, -1
         temp = a(0)
         a(0) = a(bottom)
         a(bottom) = temp;
         call siftdown(a, 0, bottom)
      end do
   
   end subroutine heapsort
   
   subroutine siftdown(a, start, bottom)
      use typre
      implicit none
   
      integer(ip), intent(in out) :: a(0:)
      logical, external :: IsLowerThan
      integer(ip), intent(in) :: start, bottom
      integer(ip) :: child, root
      integer(ip) :: temp
      
      logical :: TestResult

      
      root = start
      do while(root*2 + 1 < bottom)
         child = root * 2 + 1
      
         if (child + 1 < bottom) then
            if (a(child)<a(child+1)) child = child + 1
         end if
         if (a(root)<a(child)) then
            temp = a(child)
            a(child) = a (root)
            a(root) = temp
            root = child
         else
            return
         end if  
      end do      
   
   end subroutine siftdown
   
   function binarySearch_I (a, value)
      use typre
      implicit none
      integer(ip)                  :: binarySearch_I
      integer(ip), intent(in), target :: a(:)
      integer(ip), intent(in)         :: value
      integer(ip), pointer            :: p(:) => NULL()
      integer(ip)                  :: mid, offset
   
      p => a
      binarySearch_I = 0
      offset = 0
      do while (size(p) > 0)
         mid = size(p)/2 + 1
         if (p(mid) > value) then
               p => p(:mid-1)
         else if (p(mid) < value) then
               offset = offset + mid
               p => p(mid+1:)
         else
               binarySearch_I = offset + mid    ! SUCCESS!!
               return
         end if
      end do
   end function binarySearch_I
   
   subroutine HeapSortAuxArray(n,a,auxarray)
      use typre
      implicit none
      integer(ip) :: n
      integer(ip), intent(in out) :: a(0:)
      integer(ip) :: auxarray(0:)
      logical, external :: IsLowerThan

      integer :: start, bottom
      integer :: temp,ipoin
      
      do ipoin = 1,n
         auxarray(ipoin-1) = ipoin-1
      enddo
   
      do start = (n - 2) / 2, 0, -1
      call siftdownAuxArray(a, start, n,auxarray);
      end do
   
      do bottom = n - 1, 1, -1
         temp = auxarray(0)
         auxarray(0) = auxarray(bottom)
         auxarray(bottom) = temp;
         call siftdownAuxArray(a, 0, bottom,auxarray)
      end do
      
      do ipoin = 1,n
         auxarray(ipoin-1) = auxarray(ipoin-1)+1
      enddo
   
   end subroutine 
   
   subroutine siftdownauxArray(a, start, bottom,auxarray)
      use typre
      implicit none
   
      integer(ip), intent(in out) :: a(0:)
      integer(ip) :: auxarray(0:)   
      integer(ip), intent(in) :: start, bottom
      integer(ip) :: child, root
      integer(ip) :: temp
      
      logical :: TestResult

      
      root = start
      do while(root*2 + 1 < bottom)
         child = root * 2 + 1
      
         if (child + 1 < bottom) then
            if (a(auxarray(child))<a(auxarray(child+1))) child = child + 1
         end if
         if (a(auxarray(root))<a(auxarray(child))) then
            temp = auxarray(child)
            auxarray(child) = auxarray(root)
            auxarray(root) = temp
            root = child
         else
            return
         end if  
      end do      
   
   end subroutine 
   
end module

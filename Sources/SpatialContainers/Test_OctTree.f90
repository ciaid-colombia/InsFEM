module Mod_MyRangeGiver
   use Mod_Octree

   !This is a RangeGiver implementation
   type, extends(RangeGiver) :: MyRangeGiver
      
      integer(ip) :: npoinx, npoiny,npoinz
      
contains      
      procedure :: GiveMeRange => MyGiveMeRange
      procedure :: Initialize
   end type

contains

   subroutine MyGiveMeRange(a,ielem,range)
      use typre
      implicit none
      class(MyRangeGiver) :: a
      integer(ip) :: ielem
      real(rp) :: range(*)
      
      integer(ip) :: auxipoin
      integer(ip) :: resultpoint(3)

      
      
      auxipoin = ielem
      resultpoint(3) = (auxipoin-1)/(a%npoinx*a%npoiny)+1
      auxipoin = auxipoin - (resultpoint(3)-1)*(a%npoinx*a%npoiny)
      resultpoint(2) = (auxipoin-1)/a%npoinx+1
      auxipoin = auxipoin - (resultpoint(2)-1)*a%npoinx
      resultpoint(1) = auxipoin
      
      range(1) = resultpoint(1)
      range(2) = resultpoint(1)+1
      range(3) = resultpoint(2)
      range(4) = resultpoint(2)+1
      range(5) = resultpoint(3)
      range(6) = resultpoint(3)+1
   end subroutine
   
   subroutine Initialize(a,npoinx,npoiny,npoinz)
      use typre
      implicit none
      class(MyRangeGiver) :: a
      integer(ip) :: npoinx, npoiny, npoinz
      
      a%npoinx = npoinx
      a%npoiny = npoiny
      a%npoinz = npoinz
   end subroutine   
   
end module




subroutine Test_OctTree
   use typre
   use Mod_Tree
   use Mod_Octree
   use Mod_MyRangeGiver
   implicit none
   
   integer(ip) :: npoin 
   integer(ip) :: npoinx = 5
   integer(ip) :: npoiny = 5
   integer(ip) :: npoinz = 5
   real(rp), pointer :: coord(:,:,:,:) => NULL()

   type(Octree) :: MyTree
   type(MemoryMan) :: Memor
   
   integer(ip) :: ipoinx,ipoiny,ipoinz

   class(LinkedListHeader), pointer :: Head => NULL()
   class(LinkedListStorage), pointer :: Store => NULL()
   
   integer(ip) :: mypoint(3), resultpoint(3)
   integer(ip) :: maxpoin
   
   integer(ip) :: auxipoin
   integer(ip) :: check,foundpoin,idime,ipoin
   
   real(rp) :: time1,time2
   
   type(MyRangeGiver) :: MacGyver
   
   real(rp) :: range(2,3)
   
   !-----------------------------------------------------
   !Test Octree for points
   maxpoin = 1
   
   npoin = npoinx*npoiny*npoinz
   allocate(coord(3,npoinx,npoiny,npoinz))
   
   
   do ipoinx = 1,npoinx
      do ipoiny = 1,npoiny
         do ipoinz = 1,npoinz
            coord(1,ipoinx,ipoiny,ipoinz) = ipoinx
            coord(2,ipoinx,ipoiny,ipoinz) = ipoiny
            coord(3,ipoinx,ipoiny,ipoinz) = ipoinz
         enddo
      enddo
   enddo
   
   call cpu_time(time1)
   
   call MyTree%InitializePointOctree(3_ip,npoin,coord,maxpoin,Memor)
   call cpu_time(time2)
   !write(*,*) 'Time for building Octree: ', time2- time1
   
   call cpu_time(time1)
   do ipoinx = 1,npoinx
      do ipoiny = 1,npoiny
         do ipoinz = 1,npoinz
         
            mypoint(1) = ipoinx
            mypoint(2) = ipoiny
            mypoint(3) = ipoinz
         
            call MyTree%GetListForCoord(coord(:,mypoint(1),mypoint(2),mypoint(3)),Head,Store)
            call Store%GoToFirstOfList(Head)
            
            foundpoin = -1
            call Store%GetNext(ipoin)
            do while (ipoin /= -1)
               !ipoin to ipoinx, ipoiny, ipoinz
               auxipoin = ipoin
               resultpoint(3) = (auxipoin-1)/(npoinx*npoiny)+1
               auxipoin = auxipoin - (resultpoint(3)-1)*(npoinx*npoiny)
               resultpoint(2) = (auxipoin-1)/npoinx+1
               auxipoin = auxipoin - (resultpoint(2)-1)*npoinx
               resultpoint(1) = auxipoin
               
               check = 0
               do idime = 1,3
                  if (coord(idime,resultpoint(1),resultpoint(2),resultpoint(3)) == coord(idime,mypoint(1),mypoint(2),mypoint(3)) ) then
                     check = check + 1
                  endif
               enddo
               
               if (check == 3) then
                  foundpoin = ipoin
                  write(*,*) 'Found the point. x: ',resultpoint(1), ' y: ',resultpoint(2), ' z: ',resultpoint(3)
                  ipoin = -1
               else
                  call Store%GetNext(ipoin)
               endif
            enddo
            
            if (foundpoin == -1) then
               write(*,*) 'Did not find the point. x: ',resultpoint(1), ' y: ',resultpoint(2), ' z: ',resultpoint(3)
               call runend('One of the points was not found in Octree')
            endif
            
         enddo
      enddo
   enddo
   
   call cpu_time(time2)
   !write(*,*) 'Time for looking for all the points: ', time2- time1
   
   
   !----------------------------------------------------------------------
   !Now I test for getting a range
   range(1,1) = 2
   range(2,1) = 5
   range(1,2) = 3
   range(2,2) = 7
   range(1,3) = 2
   range(2,3) = 3
   
   call MyTree%GetListForRange(range,Head,Store)
   call Store%GoToFirstOfList(Head)
            
   foundpoin = -1
   call Store%GetNext(ipoin)
   do while (ipoin /= -1)
      !ipoin to ipoinx, ipoiny, ipoinz
      auxipoin = ipoin
      resultpoint(3) = (auxipoin-1)/(npoinx*npoiny)+1
      auxipoin = auxipoin - (resultpoint(3)-1)*(npoinx*npoiny)
      resultpoint(2) = (auxipoin-1)/npoinx+1
      auxipoin = auxipoin - (resultpoint(2)-1)*npoinx
      resultpoint(1) = auxipoin
      
      write(*,*) (resultpoint(idime),idime=1,3)
      
      call Store%GetNext(ipoin)
   enddo

   
   
   
   !Deallocate Everything
   call MyTree%Dealloc(Memor)
   
   write(*,*) 'after Mytreedealloc'
   
   !------------------------------------------------------------------------
   !Test Octree for elements
   
   maxpoin = 10
   
   call cpu_time(time1)
   !Initialize MacGyver
   call MacGyver%Initialize(npoinx,npoiny,npoinz)
   
   call MyTree%InitializeElementOctree(3,npoin,MacGyver,maxpoin,Memor)
   call cpu_time(time2)
   !write(*,*) 'Time for initializing element Octree: ', time2-time1
      
   call cpu_time(time1)
   do ipoinx = 1,npoinx
      do ipoiny = 1,npoiny
         do ipoinz = 1,npoinz
         
            mypoint(1) = ipoinx
            mypoint(2) = ipoiny
            mypoint(3) = ipoinz
            write(*,*) 'Element list for point: ',(mypoint(idime),idime=1,3)
         
            call MyTree%GetListForCoord(coord(:,mypoint(1),mypoint(2),mypoint(3)),Head,Store)
            call Store%GoToFirstOfList(Head)
            
            foundpoin = -1
            call Store%GetNext(ipoin)
            do while (ipoin /= -1)
               !ipoin to ipoinx, ipoiny, ipoinz
               auxipoin = ipoin
               resultpoint(3) = (auxipoin-1)/(npoinx*npoiny)+1
               auxipoin = auxipoin - (resultpoint(3)-1)*(npoinx*npoiny)
               resultpoint(2) = (auxipoin-1)/npoinx+1
               auxipoin = auxipoin - (resultpoint(2)-1)*npoinx
               resultpoint(1) = auxipoin
               
               write(*,*) (resultpoint(idime),idime=1,3)
               
               call Store%GetNext(ipoin)
            enddo
         enddo
      enddo
   enddo  
   call cpu_time(time2)
   !write(*,*) 'Time for getting the list of all elements: ', time2-time1
   
   call MyTree%Dealloc(Memor)
      
   stop




end subroutine
module Mod_BuildGraph
use MPI
use Mod_ToCSR
contains

   subroutine BuildElpo(pelpo,lelpo,npoin,nelem,pnods,lnods,nelty,nnode,Memor)
      use typre 
      use Mod_memor
      implicit none

      type(MemoryMan) :: Memor
      integer(ip) :: npoin,nelem,pnods(*),lnods(*),nelty,nnode(*)
      integer(ip) :: kfl_ValuesOutOfRange
      integer(ip), allocatable :: pelpo(:), lelpo(:)
      
      integer(ip) :: ielem, inode,ipoin, ispos,pnode
      integer(ip), allocatable :: peaux(:)
      integer :: ierr

      !Create pelpo
      call Memor%alloc(npoin+1,pelpo,'pelpo','BuildElpo')
      pelpo=0
      pelpo(1)=1
      do ielem=1,nelem
         call getlnods(ielem,pnode,ispos,pnods,lnods,nelty,nnode)
         do inode=1,pnode
            ipoin=lnods(ispos+inode)
            pelpo(ipoin+1)=pelpo(ipoin+1)+1
         end do
      end do
      
      !Compress pelpo
      do ipoin=1,npoin
         pelpo(ipoin+1)=pelpo(ipoin+1)+pelpo(ipoin)
      end do
      
      !Create lelpo
      call Memor%alloc(npoin,peaux,'peaux','BuildElpo')
      peaux=-1
      call Memor%alloc(pelpo(npoin+1)-1,lelpo,'lelpo','BuildElpo')

      do ielem=1,nelem
         call getlnods(ielem,pnode,ispos,pnods,lnods,nelty,nnode)
         do inode=1,pnode
            ipoin=lnods(ispos+inode)
            peaux(ipoin)=peaux(ipoin)+1
            lelpo(pelpo(ipoin)+peaux(ipoin))=ielem
         end do
      end do
      call Memor%dealloc(npoin,peaux,'peaux','BuildElpo')

   end subroutine 


   subroutine BuildGraph(ia,ja,pelpo,lelpo,npoin,nelem,pnods,lnods,nelty,nnode,Memor)
      use typre
      use Mod_memor
      implicit none
      
      integer(ip), allocatable :: ia(:),ja(:)
      integer(ip) :: lnods(*),pnods(*),pelpo(*), lelpo(*),nelty,nnode(*),npoin,nelem
      type(MemoryMan) :: Memor
      
      integer(ip), allocatable :: iwa(:),iwa2(:)
      integer(ip), pointer     :: ilnods(:) => NULL()
      type(i1p), allocatable   :: i1wa(:)
      
      integer(ip) :: icoun,elemi,iaux1,iaux2,ielem,ipoin,jpoin,jnode,pnode,ispos,totalcoun
      

      call Memor%alloc(npoin,iwa,'iwa','BuildGraph')
      call Memor%alloc(npoin,iwa2,'iwa2','BuildGraph')
      call Memor%alloc(npoin,i1wa,'i1wa','BuildGraph')

      
      !first to count and fill
      iwa = 0
      totalcoun = 0
      do ipoin = 1,npoin
         icoun = 0
         iwa(ipoin) = ipoin
         do elemi = pelpo(ipoin),pelpo(ipoin+1)-1
            ielem = lelpo(elemi)
            call getlnods(ielem,pnode,ispos,pnods,lnods,nelty,nnode)
            do jnode = 1,pnode
               jpoin = lnods(ispos+jnode)
               if (iwa(jpoin) /= ipoin) then
                  icoun = icoun+1
                  totalcoun = totalcoun+1
                  iwa(jpoin) = ipoin
                  iwa2(icoun) = jpoin
               endif
            enddo
         enddo
         allocate(i1wa(ipoin)%l(icoun))
         i1wa(ipoin)%l = iwa2(1:icoun)
      enddo
      call Memor%allocObj(0_ip,'i1wa%l','BuildGraph',totalcoun*ip)
      
      call Memor%alloc(npoin+1,ia,'ia','BuildGraph')
      call Memor%alloc(totalcoun,ja,'ja','BuildGraph')
      call i1p2CSR(npoin,i1wa,ia,ja,Memor,'i1wa%l')
      
      call Memor%dealloc(npoin,iwa,'iwa','BuildGraph')
      call Memor%dealloc(npoin,iwa2,'iwa2','BuildGraph')
      call Memor%dealloc(npoin,i1wa,'i1wa','BuildGraph')
      
   end subroutine
   
end module

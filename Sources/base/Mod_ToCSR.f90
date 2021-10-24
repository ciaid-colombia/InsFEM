module Mod_ToCSR
   implicit none
   
   
contains
   
   subroutine i1p2CSR(npoin,list,ia,ja,Memor,vanam)
      use typre
      use Mod_Memor
      implicit none
      
      integer(ip) :: npoin,ia(*),ja(*)
      type(i1p)   :: list(*)
      type(MemoryMan) :: Memor
      character(*) :: vanam
      integer(ip) :: iaux,iaux2,ipoin
      
      iaux2 = 1
      do ipoin = 1,npoin
         ia(ipoin) = iaux2
         if (associated(list(ipoin)%l)) then
            iaux = size(list(ipoin)%l)
            ja(iaux2:iaux2+iaux-1) = list(ipoin)%l
            deallocate(list(ipoin)%l)
            iaux2 = iaux2+iaux
         endif   
      enddo
      ia(npoin+1) = iaux2
      call Memor%deallocObj(0,trim(vanam),'i1p2CSR',(iaux2-1)*ip)
      
   end subroutine
   
   subroutine r1p2CSR(npoin,list,ia,ja,Memor,vanam)
      use typre
      use Mod_Memor
      implicit none
      
      integer(ip) :: npoin,ia(*)
      real(rp) :: ja(*)
      type(r1p)   :: list(*)
      type(MemoryMan) :: Memor
      character(*) :: vanam
      integer(ip) :: iaux,iaux2,ipoin
      
      iaux2 = 1
      do ipoin = 1,npoin
         ia(ipoin) = iaux2
          if (associated(list(ipoin)%a)) then
            iaux = size(list(ipoin)%a)
            ja(iaux2:iaux2+iaux-1) = list(ipoin)%a
            deallocate(list(ipoin)%a)
            iaux2 = iaux2+iaux
         endif
      enddo
      ia(npoin+1) = iaux2
      call Memor%deallocObj(0,trim(vanam),'r1p2CSR',(iaux2-1)*rp)
      
   end subroutine
         
      
   





end module
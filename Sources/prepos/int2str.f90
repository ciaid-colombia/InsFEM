module Mod_int2str
use typre
contains

   function int2str(i) result(res)
      !-------------------------------------
      !
      !  Convert an integer(ip) to a string
      !
      !-------------------------------------
       character(:),allocatable :: res
       integer,intent(in) :: i
       character(range(i)+2) :: tmp
       write(tmp,'(i0)') i
       res = trim(tmp)
   end function
   
   function real2str(integ)
      !-------------------------------------
      !
      !  Convert an integer(ip) to a string
      !
      !-------------------------------------
      implicit none
      real(rp)   :: integ
      character(20) :: real2str
      
      write(real2str,"(F10.2)") integ
      real2str=adjustl(real2str)
      
   end function real2str

   function str2int(str)
      !-------------------------------------
      !
      !  Convert an integer(ip) to a string
      !
      !-------------------------------------
      implicit none
      character(len=*)   :: str
      integer(ip)     :: str2int 
      
      read(str,*) str2int
      
   end function str2int 
end module

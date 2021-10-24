module Mod_Chdir
   use, intrinsic :: ISO_C_Binding 
   implicit none

   public :: getcwd

contains

   subroutine getcwd(str,stat)
      character(*),intent(out) :: str
      integer,intent(out) :: stat
      interface   
         subroutine c_getcwd(str,stat) bind(C, name="getcwd")
            import
            character(kind=c_char) :: str(*)
            integer(kind=c_int)    :: stat
         end subroutine
      end interface

      call c_getcwd(str,stat)
   end subroutine

end module Mod_Chdir

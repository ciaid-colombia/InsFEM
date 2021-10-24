module Mod_OneToAllPLCDReadOnBoundaries
   use typre
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllReadOnBoundaries
   use Mod_Listen
   implicit none
   private
   public OneToAllPLCDReadOnBoundaries
   
   
   type, extends(OneToAllReadOnBoundaries) :: OneToAllPLCDReadOnBoundaries
      
   
contains
      procedure :: SpecificRootAddToBuffer => RootAddToBuffer
   
   end type
   
contains
   
  
   
   
   
   subroutine RootAddToBuffer(a)
      implicit none
      class(OneToAllPLCDReadOnBoundaries) :: a
      
      integer(ip) :: cproc    !Number of processors to which this boundary belongs
      
      integer(ip) :: inode,pnode, ipoin, irank
      
      if (a%Listener%words(1) == 'STAGE') then
         !The change of stage is denoted as -1
         a%Listener%param(2) = a%Listener%param(1)
         a%Listener%param(1) = -1
         
         do irank = 0,a%MPIsize-1
            !We add Listener to the corresponding buffer
            call a%Buffer%AddToBuffer(irank,2,a%Listener%param)
         enddo
      elseif (a%Listener%words(1) == '     ') then
         !We call the parent procedure
         call a%OneToAllReadOnBoundaries%SpecificRootAddToBuffer
      endif
   end subroutine
   
end module

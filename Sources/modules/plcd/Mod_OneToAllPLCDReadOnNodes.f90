module Mod_OneToAllPLCDReadOnNodes
   use typre
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllReadOnNodes
   use Mod_Listen 
   implicit none
   private
   public OneToAllPLCDReadOnNodes
   
   type, extends(OneToAllReadOnNodes) :: OneToAllPLCDReadOnNodes
      
      
   
contains
      
      procedure :: SpecificRootAddToBuffer => RootAddToBufferOnNodes
   
   end type
   
contains
   
   subroutine RootAddToBufferOnNodes(a)
      implicit none
      class(OneToAllPLCDReadOnNodes) :: a  
      
      integer(ip) :: irank
      
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
         call a%OneToAllReadOnNodes%SpecificRootAddToBuffer
      endif

   end subroutine
   
   
end module
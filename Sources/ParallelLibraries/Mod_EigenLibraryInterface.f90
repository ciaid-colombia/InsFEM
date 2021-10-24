module Mod_EigenLibraryInterface
   use typre
   use Mod_Memor
   use Mod_EigenSystemInterface
   use Mod_ParallelSystemInterface
   implicit none
   private
   public EigenLibraryInterface
   
   type :: EigenLibraryInterface
   private
   
contains

   procedure :: Initialize
   procedure :: Finalize
   procedure :: CreateSystem
   procedure :: DeallocateSystem
   procedure :: CreateSystemProj
   
   end type
   
contains

   subroutine Initialize(a,namda)
      implicit none
      class(EigenLibraryInterface) :: a
      character(*) :: namda
      
      call runend('EigenLibrary, Initialize procedure not implemented')
   end subroutine
   
   subroutine Finalize(a)
      implicit none
      class(EigenLibraryInterface) :: a
      
      call runend('EigenLibrary, Finalize procedure not implemented')
   end subroutine
      
   subroutine CreateSystem(a,EigenSystem,Memor)   !FACTORY 
      implicit none
      class(EigenLibraryInterface) :: a
      class(EigenSystemInterface), pointer  :: EigenSystem
      type(MemoryMan) :: Memor
      
      call runend('EigenLibrary, CreateSystem procedure not implemented')
   end subroutine

   subroutine DeallocateSystem(a,EigenSystem,Memor)   !FACTORY 
      implicit none
      class(EigenLibraryInterface) :: a
      class(EigenSystemInterface), pointer  :: EigenSystem
      type(MemoryMan) :: Memor
      
      call runend('EigenLibrary, DeallocateSystem procedure not implemented')
   end subroutine
   
   subroutine CreateSystemProj(a,ParallelSystem,Memor)  !FACTORY
      implicit none
      class(EigenLibraryInterface) :: a
      class(ParallelSystemInterface), pointer :: ParallelSystem
      type(MemoryMan) :: Memor
      
      call runend('EigenLibrary, CreateSystemProj procedure not implemented')
   end subroutine
end module

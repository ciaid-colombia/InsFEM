module Mod_EPSLibrary
   use typre
   use Mod_Memor
   use Mod_SLEPcSystem
   use Mod_EPSSystem
   use Mod_MyCommunicator
   use Mod_SLEPcLibrary
   use Mod_EigenLibraryInterface
   use Mod_EigenSystemInterface
   implicit none
#ifdef SLEPC
   private
   public EPSLibrary,EPSLibraryConst

   type, extends(SLEPcLibrary) :: EPSLibrary


   contains   

      procedure :: CreateSystem => CreateSystemEPS

   end type

   interface EPSLibraryConst
      procedure constructor
   end interface EPSLibraryConst

contains

   function constructor()
      class(EPSLibrary), pointer :: constructor

      allocate(constructor)

   end function constructor

   subroutine CreateSystemEPS(a, EigenSystem, Memor)  !FACTORY
      implicit none
      class(EPSLibrary) :: a
      class(EigenSystemInterface), pointer :: EigenSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      class(EPSSystem), pointer :: iEigenSystem => NULL()

      allocate(iEigenSystem,stat=istat)   
      call Memor%allocObj(istat,'EigenSystem','CreateSystemEPS',1)

      EigenSystem => iEigenSystem

   end subroutine
#endif

end module Mod_EPSLibrary

module Mod_SVDLibrary
   use typre
   use Mod_Memor
   use Mod_MyCommunicator
   use Mod_SVDSystem
   use Mod_SLEPcSystem
   use Mod_EigenLibraryInterface
   use Mod_EigenSystemInterface
   use Mod_SLEPcLibrary
   implicit none
#ifdef SLEPC
   private
   public SVDLibrary,SVDLibraryConst

   type, extends(SLEPcLibrary) :: SVDLibrary

   contains   

      procedure :: CreateSystem         => CreateSystemSVD

   end type

   interface SVDLibraryConst
      procedure constructor
   end interface SVDLibraryConst

contains

   function constructor()
      class(SVDLibrary), pointer :: constructor
      allocate(constructor)
   end function constructor

   subroutine CreateSystemSVD(a,EigenSystem,Memor)  !FACTORY
      implicit none
      class(SVDLibrary) :: a
      class(EigenSystemInterface), pointer :: EigenSystem
      type(MemoryMan) :: Memor
      integer(ip) :: istat

      class(SVDSystem), pointer :: iEigenSystem => NULL()

      allocate(iEigenSystem,stat=istat)   
      call Memor%allocObj(istat,'EigenSystem','CreateSystemSVD',1)

      EigenSystem => iEigenSystem

   end subroutine
#endif

end module Mod_SVDLibrary

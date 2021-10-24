module Mod_Mesh_AR_Interface
   use typre
   use Mod_AdaptiveInterface
   use Mod_Mesh, only : FemMesh
   implicit none
   private
   public Mesh_IniData
   
   type, extends(AdaptiveRefinerInitialData) :: Mesh_IniData
   class(FemMesh), pointer :: Mesh => NULL()
   
   
contains
      procedure :: GetLnode 
      procedure :: Local2Global
      procedure :: GetProcessor
      
      procedure :: Initialize
   end type

contains
   subroutine Initialize(a,Mesh)
      use typre
      implicit none
      class(Mesh_IniData) :: a
      class(FemMesh), target :: Mesh
      
      a%Mesh => Mesh
      call Mesh%GetNdime(a%ndime)
      call Mesh%GetNpoin(a%npoin)
      call Mesh%GetGnpoin(a%gnpoin)
      call Mesh%GetNpoinLocal(a%npoinLocal)
      call Mesh%GetNelem(a%nelem)
   end subroutine

   subroutine GetLnode(a,ielem,pnode,lnode)
      use typre
      implicit none
      class(Mesh_IniData) :: a
      integer(ip) :: ielem,pnode
      integer(ip), pointer :: lnode(:)
      
      call a%Mesh%GetLnode(ielem,pnode,lnode)
   end subroutine
   
   subroutine Local2Global(a,npoin,LocalNumbering,GlobalNumbering)
      use typre
      implicit none
      class(Mesh_IniData) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalNumbering(npoin),GlobalNumbering(npoin)
      
      call a%Mesh%Local2Global(npoin,LocalNumbering,GlobalNumbering)
   end subroutine
   
   subroutine GetProcessor(a,npoin,LocalNumbering,ProcessorList)
      use typre
      implicit none
      class(Mesh_IniData) :: a
      integer(ip) :: npoin
      integer(ip) :: LocalNumbering(npoin),ProcessorList(npoin)
      
      call a%Mesh%LocalOrdering%GetProc(npoin,LocalNumbering,ProcessorList)
   end subroutine
end module

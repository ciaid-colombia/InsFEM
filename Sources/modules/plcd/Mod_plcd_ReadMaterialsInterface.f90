module Mod_plcd_ReadMaterials

interface
 

   subroutine ReadMaterials_Root(Listener,NumberOfMaterials,Materials,AuxRead_MaterialTypeList,MPIcomm,MPIsize,MPIroot,MPIrank,ndime,LargeStrainsflag,Memor)
      use typre
      use Mod_Listen
      use Mod_Mesh
      use Mod_plcd
      use Mod_Memor

      implicit none
      type(ListenFile) :: Listener
      integer(ip) :: NumberOfMaterials
      type(MatArray), allocatable :: Materials(:)
      character(5), allocatable :: AuxRead_MaterialTypeList(:)
      integer(ip) :: MPIsize,MPIroot,MPIrank,ndime,MPIcomm,LargeStrainsflag
      type(MemoryMan) :: Memor
      
   
   end subroutine
   
   subroutine ReadMaterials_Scatter(NumberOfMaterials,Materials,AuxRead_MaterialTypeList,MPIcomm,MPIsize,MPIroot,MPIrank,ndime,LargeStrainsflag,Memor)
      use MPI
      use typre
      use Mod_Listen
      use Mod_Mesh
      use Mod_plcd
      use Mod_Memor
      implicit none
      type(ListenFile) :: Listener
      integer(ip) :: NumberOfMaterials
      type(MatArray), allocatable :: Materials(:)
      character(5), allocatable :: AuxRead_MaterialTypeList(:)
      integer(ip) :: MPIsize,MPIroot,MPIrank,ndime,MPIcomm,LargeStrainsflag
      type(MemoryMan) :: Memor
    
   end subroutine
end interface

end module
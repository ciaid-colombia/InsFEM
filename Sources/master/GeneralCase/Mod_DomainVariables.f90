module Mod_DomainVariables
  use typre
  use Mod_Memor
  use Mod_Mesh
  use Mod_MeshInterpolator
  use Mod_Timer
  use Mod_Advector
  implicit none

  type :: DomainVariables   
     type (FemMesh) :: Mesh
     type (FemMesh) :: OldMesh
     integer(8)     :: DomCurrentMemo,DomMaxMemo !Memo counters, modules
     
     !For FixedMeshALE
     type(Interpolator) :: FMALEInterpolator
     type(Advector)     :: FMALEAdvector
     
     type(Timer) ::  cpu_MeshProjections
  end type
end module 

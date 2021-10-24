module Mod_plcd_BMatrixFactory
   use typre
   use Mod_Memor
   use Mod_PLCD
   use Mod_plcd_BMatrix
   use Mod_plcd_BMatrix_InfinitesimalStrain_3D
   use Mod_plcd_BMatrix_InfinitesimalStrain_2D
   use Mod_plcd_BMatrix_TotalLargeStrains_3D
   use Mod_plcd_BMatrix_TotalLargeStrains_2D
   implicit none

   
contains
   
   subroutine CreateBMatrix(a,Memor,BMat) 
      class(PLCDProblem) :: a
      type(MemoryMan) :: Memor
      class(BMatrix), pointer :: BMat
      
      integer(ip) :: ndime
      
      
      !Here we need to put all the ifs
      !2d, 3d,finite strains, etc, shells
      call a%Mesh%GetNdime(ndime)
      if (ndime == 3) then
         if (a%kfl_LargeStrains == 2) then
             allocate(BMatrix_LSIS3D::BMat)
         else
            allocate(BMatrix_IS3D::BMat)
         endif
      elseif (ndime == 2) then
         if (a%kfl_LargeStrains == 2) then 
            allocate(BMatrix_LSIS2D::BMat)
         else
            allocate(BMatrix_IS2D::BMat)
         endif
      else
         call runend('Bmat for this case not ready')
      endif
      call Memor%allocObj(0,'BMat','CreateBMatrix',1)
   end subroutine
  
   subroutine DestroyBMatrix(a,Memor,BMat)
      class(PLCDProblem) :: a
      type(MemoryMan) :: Memor      
      class(BMatrix), pointer :: BMat
      
      deallocate(BMat)
      call Memor%deallocObj(0,'BMat','CreateBMatrix',1)
   end subroutine
      

end module
   
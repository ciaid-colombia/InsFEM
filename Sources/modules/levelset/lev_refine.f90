subroutine lev_refine(a,itask)
   use typre
   use Mod_LevelSet
   use Mod_phpRefineArrays
   implicit none
   class(LevelSetProblem), target :: a
   character(6) :: itask
   
   real(rp), pointer, save :: auxsmoothgradient(:,:) => NULL()

   interface
     subroutine lev_ComputeCuts(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_CutElementsAndListByLayers(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine
   end interface


   integer(ip) :: nelem,npoin,ndime

   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)


   call php_RefineArrays(a,itask,a%level,'level')

   call a%Memor%realloc(nelem,a%ElementListByLayers,'ElementListByLayers','lev_refine')

   !Dealloc the cut element
   call a%CutMesh%deallocCutElement

   call a%CutMesh%deallocCutMesh
   call a%CutMesh%allocCutMesh(a%Memor,a%Mesh)

   if (a%kfl_ForceEulerianAdvection == 0) then
      !First we deallocate the previous one
      call a%CutMesh%deallocCutElement

      !now we recompute the cuts
      call lev_CutElementsAndListByLayers(a)
   endif

   if (a%kfl_SmoothGradient == 1) then
      call php_RefineArrays_b(a,itask,a%SmoothGradient,'SmoothGradient')
      !call a%Memor%realloc(ndime,npoin,a%SmoothGradient,'SmoothGradient','lev_memall')
      auxsmoothgradient => a%SmoothGradient
   endif

end subroutine




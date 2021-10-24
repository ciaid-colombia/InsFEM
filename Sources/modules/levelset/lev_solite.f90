subroutine lev_solite(a)
   use typre
   use Mod_LevelSet
   use def_parame
   implicit none
   class(LevelSetProblem) :: a
   
   interface
      subroutine php_solite(a)
         use Mod_PhysicalProblem
         class(PhysicalProblem) :: a
      end subroutine
   end interface
   
   
   logical :: isALE
   
   call a%Mesh%GetALE(isALE)
   !If we are doing ALE, do not advect the level set function
   if (isALE .and. a%kfl_ForceEulerianAdvection == 0) then
      a%itera = a%maxit+1 
      return
   else
      call php_solite(a)
   endif
   
   
end subroutine
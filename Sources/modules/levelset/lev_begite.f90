subroutine lev_begite(a)
   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_TimeIntegrator   
   use Mod_CutMesh
   implicit none
   class(LevelSetProblem) :: a
   
   interface
     subroutine lev_ComputeCuts(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine  
   end interface
   
   !Not for the first outer iteration since it has been computed in begste
   if (a%OutIiter > 1) call lev_ComputeCuts(a)
   
   

end subroutine
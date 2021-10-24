subroutine lev_reabcs(a,itask,kflag)
   use Mod_Memor
   use typre
  
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag

   integer(ip) :: npoin
   
   !Initializations
   if (itask == 0) then
      
    
   !Header
   elseif(itask == 1) then
      
      
   !Finalization
   elseif(itask == 100) then
      
   
   endif
   
end subroutine


subroutine lev_ReadOnNodes(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem) :: a  

end subroutine   

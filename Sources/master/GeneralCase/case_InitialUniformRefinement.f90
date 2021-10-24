subroutine case_InitialUniformRefinement(a,itask)
   !Deals with the adaptive process at the beginning of each time step
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_Mesh_AR_Interface   
   use Mod_Refiners
   implicit none
   class(GeneralCase), target :: a
   type(adaptiveVariables), pointer :: ad => NULL()
   type(domainVariables), pointer :: d => NULL()

   integer(ip) :: kfl_AdaptiveRefinement, iref, itask
   character(5) :: RefinementLeader
   
   ad => a%caseVars%adaptiveVars

   if (itask == 0) then
      !if initial uniformrefinement and periodic boundary conditions
      !we deactivate periodic boundary conditions till the last refinement is done
      !(PBC do not work with adaptivity so far)
      if (ad%NumberOfInitialUniformRefinementSteps > 0) then
         d => a%caseVars%domainVars
         
         call d%Mesh%GetPerio(ad%kfl_perio)
         if (ad%kfl_perio == 1) then
            call d%Mesh%SetPerio(0)
         endif
      endif
   
   elseif (itask == 1) then

   !Initial Uniform Refinement
      kfl_AdaptiveRefinement = ad%kfl_AdaptiveRefinement
      if (ad%NumberOfInitialUniformRefinementSteps > 0) then
         if (ad%kfl_AdaptiveRefinement == 0) then
            ad%kfl_AdaptiveRefinement = 1
            call a%Adaptive(1)
         endif
      endif
      ad%kfl_AdaptiveRefinement = kfl_AdaptiveRefinement

   elseif (itask == 2) then
      
      if (ad%NumberOfInitialUniformRefinementSteps > 0) then
         RefinementLeader = ad%RefinementLeader
         kfl_AdaptiveRefinement = ad%kfl_AdaptiveRefinement
         ad%RefinementLeader = 'ALLRE'
         ad%kfl_AdaptiveRefinement = 1
         !All steps minus one (because periodic boundary conditions need to be reactivated)
         do iref = 1,ad%NumberOfInitialUniformRefinementSteps-1
            call a%Adaptive(2)
         enddo  
         
         !In case of periodic boundary conditions we reactivate them here
         d => a%caseVars%domainVars
         call d%Mesh%SetPerio(ad%kfl_perio)
         
         !The last one (with pBC reactivated)
         call a%Adaptive(2)
      
         ad%RefinementLeader = RefinementLeader
         ad%kfl_AdaptiveRefinement = kfl_AdaptiveRefinement

         if (ad%kfl_AdaptiveRefinement == 0) then
            call ad%Refiner%Dealloc
         endif
      endif

   endif
end subroutine

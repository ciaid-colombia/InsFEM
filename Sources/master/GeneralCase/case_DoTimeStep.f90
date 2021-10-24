subroutine case_DoTimeStep(a)
   use typre
   use Mod_masterVariables
   use Mod_GeneralCase
   implicit none
   class(GeneralCase), target :: a
   
   type (masterVariables), pointer :: m => NULL()
   m => a%caseVars%masterVars
   
   !Iterate the case
   m%kfl_gocou = 1
   m%iiter = 1

   !This is only local coupling, not case coupling
   !-----local to the container------
   local_coupling: do while (m%kfl_gocou==1)

       !Mesh Projections for FMALE
       call a%MeshProjections
       call a%Doiter
       call a%LocalCouplingConvergence
       
   end do local_coupling

end subroutine 

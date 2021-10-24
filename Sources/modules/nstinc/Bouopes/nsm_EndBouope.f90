module Mod_nsm_EndBouope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope
   use Mod_nsm_Forces
   use Mod_nsm_Output
   !use Mod_nsm_ViscosityBouope
   use Mod_nsm_PhysicalProperties
   use Mod_nsm_BaseBouopeRoutines
   use Mod_NavierStokes
   implicit none

contains
   subroutine SetPointers
      implicit none

      call ResetProcedureComposition
      call SetPointersAndHooksToNULLSUB
      call SetPointersForces%Initialize
      call SetPointersOutput%Initialize
      call SetPointersPhysicalProperties%Initialize 
      !call SetPointersViscosity%Initialize
      call SetPointersForces%Set
      call SetPointersPhysicalProperties%SetTask('Bouope')
      call SetPointersPhysicalProperties%Set
      !call SetPointersViscosity%Set
      call SetPointersOutput%Set
      call SetPointersForces%Finalize
      call SetPointersPhysicalProperties%Finalize 
      !call SetPointersViscosity%Finalize
      call SetPointersOutput%Finalize

      call ConcatenateProcedures(ProcHook_AllocateArrays,AllocateBaseBouopeArrays)
      call ConcatenateProcedures(ProcHook_DeallocateArrays,DeallocateBaseBouopeArrays)
      call ConcatenateProcedures(ProcHook_Interpolates,BoundaryInterpolates)
      call ConcatenateProcedures(ProcHook_Gathers,BoundaryGathers)
   end subroutine

end module

subroutine nsm_EndBouope(NSProblem,task)
!-----------------------------------------------------------------------
!****f* Nstinc/nsm_endbouope
! NAME 
!    nsm_endbouope
! DESCRIPTION
!    Navier-Stokes boundary elemental operations:
!***
!-----------------------------------------------------------------------
   use Mod_nsm_EndBouope

   implicit none
   class(NavierStokesProblem), target :: NSProblem
   character(6) :: task
   integer(ip)  :: aux_logic
   
   !SetPointers for NavierStokesProblem
   a=>NSProblem
   itask = task
   
   if (itask .eq. 'Endite') then
      aux_logic = 0 
      if (a%kfl_computeTractions == 1) aux_logic =1
      if (a%kfl_bousg == 1) aux_logic =1
      if (aux_logic == 0) return
   elseif (itask .eq. 'Endste') then
      aux_logic = 0
      if (a%kfl_outfm == 1) aux_logic =1
      if (a%kfl_postBtract==1) aux_logic =1
      if (a%kfl_computeTractions == 1) aux_logic =1
      if (a%kfl_bousg == 1) aux_logic =1
      if (aux_logic == 0) return
   endif

   !Pointers
   call SetPointers

   call nsm_bouLoop
   
end subroutine

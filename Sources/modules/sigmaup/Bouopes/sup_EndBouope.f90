module Mod_sup_EndBouope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope
   use Mod_sup_Forces
   use Mod_nsm_Output
   use Mod_nsm_BaseBouopeRoutines
   use Mod_nsm_Viscosity
   use Mod_NavierStokes
   implicit none

contains
   subroutine SetPointers
      implicit none

      call ResetProcedureComposition
      call SetPointersAndHooksToNULLSUB
      call SetPointersForces%Initialize
      call SetPointersOutput%Initialize
      call SetPointersForces%Set
      call SetPointersOutput%Set
      call SetPointersForces%Finalize
      call SetPointersOutput%Finalize

      call ConcatenateProcedures(ProcHook_PhysicalProp,PhysicalProp)
      call ConcatenateProcedures(ProcHook_AllocateArrays,AllocateBaseBouopeArrays)
      call ConcatenateProcedures(ProcHook_DeallocateArrays,DeallocateBaseBouopeArrays)
      call ConcatenateProcedures(ProcHook_AllocateArrays,AllocateBaseBouopeArraysSUP)
      call ConcatenateProcedures(ProcHook_DeallocateArrays,DeallocateBaseBouopeArraysSUP)
      call ConcatenateProcedures(ProcHook_InGauss,InGauss)
      call ConcatenateProcedures(ProcHook_Interpolates,BoundaryInterpolates)
      call ConcatenateProcedures(ProcHook_Interpolates,BoundaryInterpolatesSUP)
      call ConcatenateProcedures(ProcHook_Gathers,BoundaryGathers)
      call ConcatenateProcedures(ProcHook_Gathers,Gathers)

   end subroutine

  subroutine PhysicalProp
      implicit none

      call a%GetPhysicalParameters(imat,acden,acvis)

      if(a%MatProp(imat)%lawvi<0)then
          beta=a%MatProp(imat)%LawViParam(2)
          acvis=acvis*beta         
      end if
  end subroutine

  subroutine Gathers
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      !Gathers in boundary elements 
      if (sup%LogFormulation==0) call e%gatherb(tn,bosig,sup%sigma(:,:,1))
      if (sup%LogFormulation==1) call e%gatherb(tn,bosig,sup%sigmaold(:,:))

  end subroutine

  subroutine InGauss

            !Non-Newtonian viscosity
            if(sup%MatProp(imat)%lawvi > 0)then
               call nsi_vislaw(e,grbve,sup%MatProp(imat)%lawvi,sup%MatProp(imat)%LawViParam,acvis)
            end if  
            
            !Coupling with temperature, WLF model
            if (sup%kfl_cotem ==1) then
               call e%gather(1,eltem,a%tempe)
                if (sup%kfl_cotem_WLF==1) call sup_templaw_WLF(e, eltem, sup%ReferenceTemp, sup%c1_WLF, sup%c2_WLF, sup%MatProp(imat)%LawViParam, acvis, lambda)
                if (sup%kfl_cotem_Arrhenius==1) call sup_templaw_Arrhenius(e, eltem, sup%ReferenceTemp, sup%alpha_Arrhenius, sup%MatProp(imat)%LawViParam, acvis, lambda)
            endif  
            
            !Smagorinsky
            if(sup%kfl_cotur == -1) then
               call nsm_smago(e,grbve,acden,sup%turbu(1),vista)
               acvis = acvis + vista
            endif

  end subroutine

end module

subroutine sup_EndBouope(TFProblem,task)
!-----------------------------------------------------------------------
!    sup_forces
! DESCRIPTION
!    This routine computes forces and moments at a given time step.
!
!-----------------------------------------------------------------------
   use Mod_sup_EndBouope

   implicit none
   class(ThreeFieldNSProblem), target :: TFProblem   
   character(6) :: task
   integer(ip)  :: aux_logic
  
   a   => TFProblem
   sup => TFProblem

   if (task .eq. 'Endite') then
      aux_logic = 0 
      if (a%kfl_computeTractions == 1) aux_logic =1
      if (aux_logic == 0) return
   elseif (task .eq. 'Endste') then
      aux_logic = 0
      if (a%kfl_outfm == 1) aux_logic =1
      if (a%kfl_postBtract==1) aux_logic =1
      if (a%kfl_computeTractions == 1) aux_logic =1
      if (aux_logic == 0) return
   endif
   
   !Pointers
   call SetPointers

   call nsm_bouLoop

end subroutine

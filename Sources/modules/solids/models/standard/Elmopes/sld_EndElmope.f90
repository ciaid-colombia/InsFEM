module Mod_sld_EndElmope
   use typre
   use Mod_sld_BaseElmope
   use Mod_sld_BaseElmopeRoutines
   use Mod_sld_ExternalForces
   use Mod_sld_NonLinearDerivatives
   use Mod_sld_calculateVmass
   use Mod_sld_calculateSigma
   use Mod_sld_calculatePress
   implicit none

contains       
    
   subroutine SetPointers
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      !Reset Procedure Composition
      call ResetProcedureComposition
      !Set All Pointers To NULLSUB (in Mod_sld_BaseElmope)
      call SetPointersAndHooksToNULLSUB

      !Non-linear elements
      call SetPointersNonLinearDerivatives(0)
      call SetPointersCalculateVmass(0)
      call SetPointersCalculateSigma(0)
      call SetPointersCalculatePress(0)
      
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateSolidBase)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateSolidBase)

      call ConcatenateProcedures(ProcHook%Gathers,Gathers)

      call ConcatenateProcedures(ProcHook%InGauss,InGauss)

      call ConcatenateProcedures(ProcHook%PostGauss,PostGauss)

      if(a%sld_type== 'NONLI' ) then
          call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUP)
          call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

          call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateNonLinearSolidArrays)
          call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateNonLinearSolidArrays)

          call SetPointersCalculateVmass(1)
          call SetPointersCalculateSigma(1)
          call SetPointersCalculatePress(1)

      endif

      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
          call SetPointersNonLinearDerivatives(1)
          !If more than one element type, then pointers should be reset when ielty changes!
          call a%Mesh%GetNelty(nelty)
          if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

          call ConcatenateProcedures(ProcHook%PreGauss,InGaussVolumesNonLinear)
      end if

      call SetPointersCalculateVmass(100)
      call SetPointersCalculateSigma(100)
      call SetPointersCalculatePress(100)
      call SetPointersNonLinearDerivatives(100)

   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine

   subroutine Gathers
       implicit none

       call displacementGather

   end subroutine

   subroutine InGauss
       implicit none
       integer(ip) :: ndime

       call a%Mesh%GetNdime(ndime)

       !reset
       strain = 0.0_rp
       stress = 0.0_rp
       gradDisp = 0.0_rp
       elext=0.0_rp

       dvol = 0.0_rp
       dvol = e%weigp(e%igaus)*e%detjm
       dvolt0 = dvol + dvolt0

       call calculateGradientsAndDeter

       call calculateStress(ndime,ielem)

       a%strain_g(ielem)%a(:,igaus) = strain
       a%stress_g(ielem)%a(:,igaus) = stress

   end subroutine

   subroutine PostGauss
      implicit none
      integer(ip) :: ndime,sz

      if(a%kfl_saveStrainsALE) then

          call a%Mesh%GetNdime(ndime)
          sz  = (ndime*(ndime+1))/2

          !this is allocated by the ale module
          a%pStrain(ielem)%a(:) = 0.0_rp
          a%pStrain(ielem)%a(:) = getPrincipalStrains(ndime,sz,strain)

      end if

  end subroutine

end module

subroutine sld_EndElmope(SldProblem,task)
   use Mod_Solids
   use Mod_sld_EndElmope
   implicit none
   class(SolidsProblem), target :: SldProblem
   character(6) :: task
   logical(lg) :: aux_logic

   itask = task

   a=>SldProblem

   if (itask .eq. 'Endite') then

   elseif (itask .eq. 'Endste') then

      call SetPointers

      call sld_elemLoop

   end if
   
end subroutine

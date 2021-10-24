module Mod_sldsup_elmope_lin
  use Mod_CauchyElement
  use Mod_SUPSolids_lin
  use Mod_sld_BaseElmopeRoutines
  implicit none   
   
  contains

  subroutine SetPointersElmopeSUP_lin
     use typre
     implicit none
     procedure() :: NULLSUB
     integer(ip) :: kfl_nonlinear,nelty
           
     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateAssemblyArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateAssemblyArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUP)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUPLinear)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUPLinear)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUP_SGS)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP_SGS)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateSUPAssemblyMatrices)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateSUPAssemblyMatrices)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays)
     call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays_SGS)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetAssemblyMatrices)

     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpDisplacements)
     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpSigma)
     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)

     call SetPointers_Linear

     if(a%kfl_linea == 0 ) then      !Standard RHS linearization
         call SetRHSPointers
     else if(a%kfl_linea == 2 ) then !Newton-Raphson linearization
         !call SetNRPointers
         write(*,*) "sldsup: NR linearization not ready, setting RHS"
         call SetRHSPointers
     endif

     ProcPointer%getTauParameters => getTauParameters_linear

     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP)
     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP_SGS)
     call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleRhsSUP)

     call ConcatenateProcedures(ProcHook%InGauss,InGaussResetExternal)

     !--------------------------------------------------------------------------
     if(sup_lin%kfl_timei==1) then
         !Dynamic allocations
         call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateDynamicArray)
         call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateDynamicArray)

         !Dynamic resets
         call ConcatenateProcedures(ProcHook%ResetArrays,ResetDynamicMatrix)

         !Dynamic elemental build
         call ConcatenateProcedures(ProcHook%DynamicMass,buildMassMat)
     
         !Dynamic build to LHS
         call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleMassMatrix)

         call ConcatenateProcedures(ProcHook%Initializations,AllocVelandAccel)
         call ConcatenateProcedures(ProcHook%Initializations,InitTimeIntegration)
         call ConcatenateProcedures(ProcHook%Gathers        ,GatherDynamic)
         call ConcatenateProcedures(ProcHook%Finalizations  ,DeallocVelandAccel)

         !Dynamic matrix routines
         call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicForce)
         call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicMass)

     endif

     !Non-linear elements
     call sup_lin%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
     if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%InGaussElmats,ProcPointer%PostGauss)
         ProcPointer%PostGauss => NULLSUB      
     end if

     call ConcatenateProcedures(ProcHook%ToLinearSystem,AssembleLinearSystem)
     !-------------------------------------------------------
     !If more than one element type, then pointers should be reset when ielty changes!
     call a%Mesh%GetNelty(nelty)
     if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

    !Mass matrix
    call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%Dynamic)
    !K_mat and K_geo, and internal forces
    call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussElmats)
    call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussVec)

  end subroutine SetPointersElmopeSUP_lin
     
  subroutine SetPointers_Linear
     use typre
     implicit none

     !Set the Specific Pointers for this subroutine
     ProcPointer%PostGauss  => PostGaussElmats_linear
     ProcPointer%sld_elmrhu => sld_elmrhu

     call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp)

     call ConcatenateProcedures(ProcHook%Gathers,sigmaGather)
     call ConcatenateProcedures(ProcHook%Gathers,displacementGather)
     call ConcatenateProcedures(ProcHook%Gathers,pressGather)

  end subroutine SetPointers_Linear

  !----------------------------------------------------------
  !PostGauss Matrices
  subroutine PostGaussElmats_linear
     implicit none
     integer(ip) :: nd,tn

     call sup_lin%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call ProcHook%PhysicalProp

     call sup_elmVS(e,nd,tn,dvol,P,     elmsv)
     call sup_elmVP(e,nd,   dvol,       elmpv)
     call sup_elmEU(e,nd,tn,dvol,P,     elmut)
     call sup_elmES(e,nd,tn,dvol,D,     elmst)
     call sup_elmQU(e,nd,   dvol,       elmuq)
     call sup_elmQP(e,nd,   dvol,Dv_scalar,elmpq)

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_ES(e,nd,tn,dvol,elm_usgs_ES)
     call sup_elm_usgs_EP(e,nd,tn,dvol,elm_usgs_EP)
     call sup_elm_usgs_QS(e,nd,tn,dvol,elm_usgs_QS)
     call sup_elm_usgs_QP(e,nd,tn,dvol,elm_usgs_QP)

     !-------------------Tau_s matrix calculation----------
     call sup_elm_ssgs_VU(e,nd,tn,dvol,C_dev,elm_ssgs_VU)
     call sup_elm_ssgs_VS(e,nd,tn,dvol,P    ,elm_ssgs_VS)
     call sup_elm_ssgs_EU(e,nd,tn,dvol,P    ,elm_ssgs_EU)
     call sup_elm_ssgs_ES(e,nd,tn,dvol,D    ,elm_ssgs_ES)

     !-------------------Tau_p matrix calculation----------
     call sup_elm_psgs_VU(e,nd,   dvol,K,G      ,elm_psgs_VU)
     call sup_elm_psgs_VP(e,nd,   dvol          ,elm_psgs_VP)
     call sup_elm_psgs_QU(e,nd,   dvol          ,elm_psgs_QU)
     call sup_elm_psgs_QP(e,nd,   dvol,Dv_scalar,elm_psgs_QP)
  end subroutine

  subroutine PhysicalProp
      implicit none
      integer(ip) :: nd,tn

      call sup_lin%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call sup_lin%SUPGetPhysicalParameters(nd,tn,densi,K,G,C_dev,D,D_dev,Dv_scalar,P)
  end subroutine

  !Gauss matrices
  subroutine InGaussElmats_linear
      implicit none
      
  end subroutine

  subroutine OnIeltyChange
     implicit none
      
     if (e%ielty /= ielty0) then
        call SetPointersElmopeSUP_lin
        ielty0 = e%ielty
     endif
  end subroutine

end module   

!Solids_ELMOPE subroutine   
subroutine sldsup_Elmope_lin(SldProblem)
   use Mod_sldsup_elmope_lin

   implicit none
   class(SUPSolidsProblem_lin), target :: SldProblem
   
   !Setting sld problem pointer in sld_BaseElmope
   up      => SldProblem
   sup_lin => SldProblem
   sup     => SldProblem
   up      => SldProblem
   a       => SldProblem

   call SetPointersGeneral
   call SetPointersPointForces
   call SetPointersElmopeSUP_lin
   
   call sld_elemLoop

end subroutine sldsup_Elmope_lin

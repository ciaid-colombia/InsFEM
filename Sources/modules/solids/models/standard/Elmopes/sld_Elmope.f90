module Mod_sld_elmope
   use Mod_CauchyElement
   use Mod_sld_BaseElmopeRoutines
   use Mod_Solids
   implicit none   
   
   contains

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine

   subroutine SetPointers
      use typre
      implicit none
      procedure() :: NULLSUB
      integer(ip) :: kfl_nonlinear,nelty
            
      !Allocates
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateAssemblyArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateAssemblyArrays)

      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateSolidBase)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateSolidBase)

      !Resets
      call ConcatenateProcedures(ProcHook%ResetArrays,ResetSolidBaseArrays)
      call ConcatenateProcedures(ProcHook%ResetArrays,ResetAssemblyMatrices)

      !Gathers
      call ConcatenateProcedures(ProcHook%Gathers,Gathers)

      !Interpolates
      call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpDisplacements)

      !Physical properties
      call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp)

      !Final assembly
      call ConcatenateProcedures(ProcHook%ToLinearSystem,AssembleLinearSystem)

      if(a%sld_type== 'NONLI' ) then
          call SetPointers_nonLinear
      else
          call SetPointers_Linear
      endif

      !Linearization scheme
      if(a%kfl_linea == 0 ) then      !Standard RHS linearization
          call SetRHSPointers
      else if(a%kfl_linea == 2 ) then !Newton-Raphson linearization
          call SetNRPointers
      endif

      !Elementat assemblies 
      call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrix)
      call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleRhs)

      !--------------------------------------------------------------------------
      if(a%kfl_timei==1) then
          !Dynamic allocations
          call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateDynamicArray)
          call ConcatenateProcedures(ProcHook%Initializations,AllocVelandAccel)

          !Dynamic initialization
          call ConcatenateProcedures(ProcHook%Initializations,InitTimeIntegration)

          !Dynamic resets
          call ConcatenateProcedures(ProcHook%ResetArrays,ResetDynamicMatrix)

          !Dynamic gather
          call ConcatenateProcedures(ProcHook%Gathers,GatherDynamic)

          !Dynamic elemental assembly
          call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleMassMatrix)

          !Dynamic matrix routines
          call ConcatenateProcedures(ProcHook%DynamicMass, buildMassMat)
          call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicForce)
          call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicMass)

          !Dynamic deallocs
          call ConcatenateProcedures(ProcHook%Finalizations,DeallocVelandAccel)
          call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateDynamicArray)
      endif

      !Non linear elements
      !call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      !if (kfl_nonlinear == 1) then
      !    call ConcatenateProcedures(ProcHook%InGaussElmats,ProcPointer%PostGauss)
      !    ProcPointer%PostGauss => NULLSUB      
      !end if

      call ConcatenateProcedures(ProcHook%InGauss,InGaussResets)
      call ConcatenateProcedures(ProcHook%InGauss,InGaussResetExternal)
      call ConcatenateProcedures(ProcHook%InGauss,GradientsAndDeterminant)

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

      !Pointers to pointers have to be in the end
      !Mass matrix
      call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%Dynamic)
      !Stiffness matrix and RHS
      call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussElmats)
      call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussVec)

   end subroutine SetPointers
      
   subroutine SetPointers_nonLinear
      use typre
      implicit none

      !Set the Specific Pointers for this subroutine
      ProcPointer%sld_elmk           => sld_elmk
      ProcPointer%sld_elmgeo         => sld_elmgeo
      ProcPointer%sld_elmrhu         => sld_elmrhu
      !ProcPointer%PostGauss          => PostGaussElmats
      ProcHook%EndGauss              => PostGaussElmats

      !------------------------------------------------------------------
      call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateNonLinearSolidArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateNonLinearSolidArrays)

      call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleGeometricMatrix)
      call ConcatenateProcedures(ProcHook%ResetArrays,ResetNonLinearSolidArrays)

   end subroutine SetPointers_nonLinear

   subroutine SetPointers_Linear
      use typre
      implicit none

      !Set the Specific Pointers for this subroutine
      ProcPointer%sld_elmk        => sld_elmk
      ProcPointer%sld_elmrhu      => sld_elmrhu
      !ProcPointer%PostGauss       => PostGaussElmats_linear
      ProcHook%EndGauss              => PostGaussElmats_linear

      !Modal analysis
      if(a%kfl_timei==1) then
          if(a%kfl_eigensolve .and. a%kfl_eigenter) then
              call ConcatenateProcedures(ProcHook%preAssemblyLhs,AssemblyEigenSystem)
              call ConcatenateProcedures(ProcHook%Finalizations,CloseEigenProblem)
          endif
      endif

   end subroutine SetPointers_Linear

   subroutine SetNRPointers

     !Allocates
     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateNRSolidArrays)

     !Resets
     call ConcatenateProcedures(ProcHook%ResetArrays,ResetNRSolidArrays)

     !Residual
     call ConcatenateProcedures(ProcHook%calculateResidual,staticResidual)
     if(a%kfl_timei==1) then
         call ConcatenateProcedures(ProcHook%calculateResidual,dynamicResidual)
     endif

     !Elemental assemblies
     call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmats)
     call ConcatenateProcedures(ProcHook%preAssemblyLhs,residual_and_jacobian)
     call ConcatenateProcedures(ProcHook%preAssemblyRhs,completeForcingTerm)

     !Deallocs
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateNRSolidArrays)

   end subroutine SetNRPointers

   !----------------------------------------------------------
   subroutine PhysicalProp
      implicit none
      integer(ip) :: ndime,sz

      call a%Mesh%GetNdime(ndime)
      sz  = (ndime*(ndime+1))/2
      
      call a%GetPhysicalParameters(det,sz,densi,c_elas,ielem)
  end subroutine

   subroutine PostGaussElmats_linear
      implicit none
      integer(ip) :: ndime,sz

      call a%Mesh%GetNdime(ndime)
      sz  = (ndime*(ndime+1))/2

      call ProcPointer%sld_elmk(e,ndime,sz,c_elas,dvol,kmat)

   end subroutine

   subroutine PostGaussElmats
      implicit none
      integer(ip) :: ndime,sz

      call a%Mesh%GetNdime(ndime)
      sz  = (ndime*(ndime+1))/2

      call calculateGradientsAndDeter

      !recalculate elastic tensor for K_mat
      call ProcHook%PhysicalProp
      call ProcPointer%sld_elmk(e,ndime,sz,c_elas,dvol,kmat)

      !calculate stress for K_geo
      call calculateStress(ndime)
      call ProcPointer%sld_elmgeo(e,ndime,sz,stress,dvol,geomat)

  end subroutine

   subroutine InGaussElmats
      implicit none
      integer(ip) :: ndime,sz

      call a%Mesh%GetNdime(ndime)
      sz  = (ndime*(ndime+1))/2

      !Forces due to current state of stress
      call calculateStress(ndime)
      call getInternalForce(e,ndime,sz,stress,dvol,intForce)

  end subroutine

   subroutine residual_and_jacobian
      implicit none

      !calculate dynamic or static residual
      call ProcHook%calculateResidual

      if(a%kfl_timei==1) then
          !We complete the mass matrix Jacobian
          massmat = tsch_deriv*a%dtinv*a%dtinv*massmat
      endif

   end subroutine

   subroutine AssemblyEigenSystem
      implicit none
      !Dirichlet Boundary Conditions, set residual for dirichlet nodes
      call sld_elmdir(a,e,massmat,elrhs)
      call sld_elmdir(a,e,kmat,elrhs)

      !Generalized eigen value problem [K-lam*M]x=0
      call a%LinearSystem%AssemblyEigenMat(e,-massmat,kmat)
   end subroutine

   subroutine completeForcingTerm
      implicit none
      real(rp) :: auxvec(e%ndime,e%pnode)
      !Very useful formulation for ROM, instead of calculating residual 
      !we calculate the variable
      !by solving A*U_i+1 = r + A*U_i, ------> U_i+1 - U_i = d_u
      !Finally, the problem to solve is : U_i+1 = A^-1*(r + A*U_i)
      !Usual Newton Raphson assembly is: d_u = A^-1 r

      !We add the previous value so we do not calculate residual 
      !but the total variable

      call matvecmult(e%ndime,e%pnode,elmat,eldisp(:,:,1),auxvec)
      elrhu = elrhu + auxvec

   end subroutine

   subroutine InGaussResets
       implicit none

       !reset
       strain = 0.0_rp
       stress = 0.0_rp
       gradDisp = 0.0_rp

   end subroutine

   subroutine GradientsAndDeterminant
       implicit none

       call calculateGradientsAndDeter

   end subroutine

   subroutine Gathers
       implicit none

       call displacementGather

   end subroutine

end module   

subroutine sld_elmope(SldProblem)
   use Mod_sld_elmope

   implicit none
   class(SolidsProblem), target :: SldProblem
   
   a=>SldProblem

   call SetPointersGeneral
   call SetPointersPointForces
   call SetPointers
   
   call sld_elemLoop

end subroutine sld_elmope

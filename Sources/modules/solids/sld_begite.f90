module Mod_sld_begiteElmope
   use Mod_Solids
   use Mod_sld_BaseElmopeRoutines
   use Mod_sld_ExternalForces
   use Mod_sld_NonLinearDerivatives
   implicit none   
   
   contains

   subroutine SetPointersBegite
      use typre
      implicit none
      integer(ip) :: kfl_nonlinear,nelty
            
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_sld_BaseElmope)
      call ResetProcedureComposition
      call SetPointersAndHooksToNULLSUB
   
      call SetPointersExternalForces(0)
      call SetPointersNonLinearDerivatives(0)
   
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateAssemblyArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateAssemblyArrays)
   
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)
   
      call ConcatenateProcedures(ProcHook%ResetArrays,ResetAssemblyMatrices)
      call ConcatenateProcedures(ProcHook%ResetArrays,ResetForceVector)
   
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
          call SetPointersNonLinearDerivatives(1)
          call ConcatenateProcedures(ProcHook%PreGauss,InGaussVolumesNonLinear)
      end if

      densi = a%densi

      ProcPointer%sld_elmrhu => sld_elmrhu

      call SetPointersExternalForces(1)
   
      call ConcatenateProcedures(ProcHook%InGauss,calculateVolume)
      call ConcatenateProcedures(ProcHook%InGauss,InGaussResetExternal)
      call ConcatenateProcedures(ProcHook%InGaussVec,calculateForces)
   
      !-------------------------------------------------------
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChangeBegite
   
     !K_mat and K_geo, and internal forces
     call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussVec)
   
     call SetPointersNonLinearDerivatives(100)
     call SetPointersExternalForces(100)
   
   end subroutine SetPointersBegite

   subroutine SetPointersBegiteIrreducible
      use typre
      implicit none
      integer(ip) :: kfl_nonlinear,nelty
            
      call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleForceRhs)
   
   end subroutine SetPointersBegiteIrreducible
   
   subroutine OnIeltyChangeBegite
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointersBegite
         ielty0 = e%ielty
      endif
   end subroutine
   
   subroutine ResetForceVector
         implicit none
   
          a%extForce(ielem)%a(:,:)=0.0_rp
   
   end subroutine
   
   subroutine calculateForces
         implicit none
   
         !Compute vector of external forces
         call ProcPointer%ExternalForces   
   
         !Compute contributions to RHS :
         call ProcPointer%sld_elmrhu(e,dvol,elext,elrhu)
   
   end subroutine
   
   subroutine AssembleForceRhs
         implicit none
         integer(ip) :: u1,uf,s1,sf,p1,bc

         call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

         a%extForce(ielem)%a(u1:uf,:) = elrhu 

   end subroutine

end module Mod_sld_begiteElmope

subroutine sld_begite(sldProblem)
!-----------------------------------------------------------------------
! NAME 
!    sld_begite
! DESCRIPTION
!    This routine starts an internal iteration for the elastic solid prob. 
!-----------------------------------------------------------------------
    use typre
    use Mod_Solids
    use Mod_sld_begiteElmope
    implicit none
    class(SolidsProblem), target :: sldProblem

   a=>sldProblem

   call SetPointersBegite
   call SetPointersBegiteIrreducible

   call sld_elemLoop

end subroutine sld_begite

module Mod_plcd_LargeStrainsOperations
   use typre  
   use Mod_plcd_BaseElmope
   use Mod_PointerSetter
   implicit none
   private
   public SetPointersLargeStrains, SetPointersLargeStrainsEndite 
   
   type, extends(PointerSetter) :: SPLargeStrains
contains
      procedure :: SpecificSet => SpecificSetPointersLargeStrains
   end type 
   type(SPLargeStrains) :: SetPointersLargeStrains
   
   type, extends(PointerSetter) :: SPLargeStrainsEndite
contains
      procedure :: SpecificSet => SpecificSetPointersLargeStrainsEndite
   end type 
   type(SPLargeStrainsEndite) :: SetPointersLargeStrainsEndite
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SpecificSetPointersLargeStrains(d)
      implicit none
      class(SPLargeStrains) :: d
            !-------------------------------------------------------
            !LargeStrains
            if (a%kfl_LargeStrains /= 0) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyGeometricStiffness)
            endif
            if (a%kfl_LargeStrains == 2) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocDisplacementGradient)
               call ConcatenateProcedures(ProcHook%InGauss,SetDisplacementGradienttoBMatrix)
               call ConcatenateProcedures(ProcHook%Finalizations,DeAllocDisplacementGradient)
            endif
   end subroutine
   
   subroutine ComputeAndAssemblyGeometricStiffness
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental geometric matrix 
      !    B_t*K_geo*B  , see Belytschko for notation
      !-----------------------------------------------------------------------
      implicit none

      integer(ip)                :: idime
      real(rp)                   :: Hmat(e%pnode,e%pnode) 
      
      call ElementMatData%GetStressTensorPointer(e%igaus,stress)
      
      Hmat = 0.0_rp
      
      call BMat%Betat_Times_Vector_Times_Beta(stress,Hmat)
      Hmat = Hmat*dvol
      
      forall(idime=1:e%ndime)
         elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + Hmat(1:e%pnode,1:e%pnode)
      end forall

   end subroutine ComputeAndAssemblyGeometricStiffness
   
   subroutine AllocDisplacementGradient
      call a%Memor%alloc(e%ndime,e%ndime,gradDisp,'gradDisp','plcd_LargeStrainsOperations')
   
   end subroutine AllocDisplacementGradient
      
   subroutine DeAllocDisplacementGradient
      call a%Memor%dealloc(e%ndime,e%ndime,gradDisp,'gradDisp','plcd_LargeStrainsOperations')
   
   end subroutine DeAllocDisplacementGradient
   
   subroutine SetDisplacementGradienttoBMatrix
      !ComputeDisplacementGradients
      gradDisp =0.0_rp
      call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)
      call BMat%SetDisplacementGradient(gradDisp)
      
   end subroutine SetDisplacementGradienttoBMatrix
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
   !ENDITE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   subroutine SpecificSetPointersLargeStrainsEndite(d)
      implicit none
      class(SPLargeStrainsEndite) :: d
            !-------------------------------------------------------
            if (a%kfl_LargeStrains == 2) then
               call ConcatenateProcedures(ProcHook%InGauss,SetDisplacementGradienttoBMatrixEndite)
            endif
   end subroutine   
   
   subroutine SetDisplacementGradienttoBMatrixEndite
      call BMat%SetDisplacementGradient(gradDisp)
   
   end subroutine SetDisplacementGradienttoBMatrixEndite

end module

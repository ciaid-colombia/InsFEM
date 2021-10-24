module Mod_nsm_NonLinearDerivatives
   use typre
   use Mod_PointerSetter
   use mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersNonLinearDerivatives

   type, extends(PointerSetter) :: SPNonLinearDerivatives
contains
      procedure :: SpecificSet => SpecificSetNonLinearDerivatives
   end type
   type(SPNonLinearDerivatives) :: SetPointersNonLinearDerivatives
   
contains
   
   subroutine SpecificSetNonLinearDerivatives(d)
      implicit none
      class(SPNonLinearDerivatives) :: d
      integer(ip) :: kfl_nonlinear
         
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook_InGauss,InGaussNonLinearDerivatives)
   end subroutine
   
   !Computation Subroutine
   subroutine InGaussNonLinearDerivatives
      implicit none
      
      call e%elmder
      call e%elmhes
   end subroutine

end module

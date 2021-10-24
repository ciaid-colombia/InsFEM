  module Mod_nsm_TurbulenceModel
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_Viscosity
   implicit none
   private
   public SetPointersTurbulenceModel

   type, extends(PointerSetter) :: SPTurbulenceModel
contains
      procedure :: SpecificSet => SpecificSetTurbulenceModel
   end type
   type(SPTurbulenceModel) :: SetPointersTurbulenceModel
   
   real(rp) :: vista
   
contains
   
   subroutine SpecificSetTurbulenceModel(d)
      implicit none
      class(SPTurbulenceModel) :: d
         
      !Turbulence model
      if (a%kfl_cotur < 0 ) then
         !Interpolate gradients
         call SetPointersInterpolateGradients%Set
      endif
      !Smagorinsky
      if (a%kfl_cotur == -1) then
         !Modify the viscosity
         call ConcatenateProcedures(ProcHook_PhysicalProp,Smagorinsky)
      !WALE model
      elseif (a%kfl_cotur == -2) then
         !Modify the viscosity
         call ConcatenateProcedures(ProcHook_PhysicalProp,WaleModel)
      endif
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine Smagorinsky
      implicit none
      call nsm_smago(e,grvel,acden,a%turbu(1),vista)
      acvis = acvis + vista
   end subroutine  
   
   subroutine WaleModel
      implicit none
      call nsm_wale(e,grvel,acden,a%turbu(1),vista)
      acvis = acvis + vista
   end subroutine
   
end module


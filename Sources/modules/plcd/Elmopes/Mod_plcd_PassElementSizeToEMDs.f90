module Mod_plcd_PassElementSizeToEMDs
   use typre
   use Mod_plcd_BaseElmope
   implicit none
   private
   public SetPointersPassElementSizeToEMDs
   
   type, extends(PointerSetter) :: SPPassElementSizeToEMDs
contains
      procedure :: SpecificSet => SpecificSetPassElementSizeToEMDs
   end type
   
   type(SPPassElementSizeToEMDs) :: SetPointersPassElementSizeToEMDs
   
   integer(ip), allocatable :: kfl_IsSet
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetPassElementSizeToEMDs(d)
      class(SPPassElementSizeToEMDs) :: d
      !Element size required by EMDS
      if (a%IsElementSizeRequiredByEMDs .eqv. .true.) then
         call ConcatenateProcedures(ProcHook%PreGauss,PassElementSizeToEMDs)
      endif
   end subroutine   
   
   
   !For actual computations
   !------------------------------------------------------
   !Hanging nodes
   subroutine PassElementSizeToEMDs
      implicit none
      
      real(rp) :: hleng
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      !Element length at center of gravity
      call e%elmlen
      hleng = e%hleng(1)
      !If smoothing, damage occurs along multiple (2) elements (less localization)
      !As a consequence hleng needs to be larger
      if (a%UseSmoothedDisplacementGradient) hleng = 2.0_rp*hleng
      call a%ElementMaterialsData(ielem)%p%SetElementSize(hleng)
   end subroutine
end module

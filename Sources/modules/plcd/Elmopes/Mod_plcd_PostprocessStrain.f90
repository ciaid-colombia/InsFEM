module Mod_plcd_PostprocessStrain
   use typre
   use Mod_plcd_BaseElmope
   implicit none
   private
   public SetPointersPostprocessStrain
   
   type, extends(PointerSetter) :: SPPostprocessStrain
contains
      procedure :: SpecificSet => SpecificSetPostprocessStrain
   end type
   type(SPPostprocessStrain) :: SetPointersPostprocessStrain
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetPostprocessStrain(d)
      implicit none
      class(SPPostprocessStrain) :: d
   
      !Post Process Strain
      if(a%PostprocessStrain) then
         call ConcatenateProcedures(ProcHook%InGaussElmats,SaveStrainToPostprocess)
      endif

   end subroutine   
   
   !Save Strain to PostProcess
   subroutine SaveStrainToPostprocess
      implicit none

      class(PLCDMaterial), pointer :: Mater => NULL()
      
      call a%ElementMaterialsData(ielem)%p%GetMaterialPointer(Mater)
      call Mater%CT%GetStrain(gradDisp,a%Strain(ielem)%a(:,igaus))
   end subroutine
    
end module 

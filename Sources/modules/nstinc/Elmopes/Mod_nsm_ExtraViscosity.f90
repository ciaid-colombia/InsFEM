module Mod_nsm_ExtraViscosity
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersExtraViscosity
   
   !For elmopes
   type, extends(PointerSetter) :: SPExtraViscosity
contains
      procedure :: SpecificSet => SpecificSetExtraViscosity
   end type
   type(SPExtraViscosity) :: SetPointersExtraViscosity
   
contains
   subroutine SpecificSetExtraViscosity(d)
      implicit none
      class(SPExtraViscosity) :: d
      
      if (a%kfl_ExtraInitialViscosity == 1) then
         call PrependProcedure(ProcPointer_PostGaussElmats,ExtraViscosityMatrices)
      endif
   end subroutine
   
   subroutine ExtraViscosityMatrices
      implicit none
      real(rp) :: ExtraVis
      
      if (a%istep < a%EIV_nsteps) then
         Extravis = (1.0_rp-real(a%istep)/real(a%EIV_nsteps))*a%EIV_Visco
         
         !BLOCK U,P
         ! Viscosity terms : we only consider mu*(grad v, grad u)
         call elmvis(e,dvolt0,ExtraVis,wrmat1)

         ! If you want the complete term div(e·grad u) uncomment next line 
         if ( a%fvins > zensi ) then
            call nsm_elmvis_div(e,dvolt0,Extravis,elmuv)
         endif
      endif
   
   end subroutine
   
end module

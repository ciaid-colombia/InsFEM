module Mod_nsm_ViscosityBouope
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   use Mod_nsm_Viscosity
   implicit none
   private
   public SetPointersViscosity

   type, extends(PointerSetter) :: SPViscosity
contains
      procedure :: SpecificSet => SpecificSetViscosity
   end type
   type(SPViscosity) :: SetPointersViscosity

contains

   subroutine SpecificSetViscosity(d)
      implicit none
      class(SPViscosity) :: d
   
      !Smagorinsky
      if (a%kfl_cotur< 0.or.a%MatProp(imat)%lawvi/=0.or.a%kfl_ExchangeLocationWallLaw==1) then
         call ConcatenateProcedures(ProcHook_Gathers,GatherSmagorinsky)
      endif
      if (a%kfl_cotur<0) call ConcatenateProcedures(ProcHook_InGauss,InGaussSmagorinsky)
      !Non-newtonian viscosity
      if(a%MatProp(imat)%lawvi /= 0) call ConcatenateProcedures(ProcHook_InGauss,InGaussNonNewtonianViscosity)
   end subroutine

   !*************!
   ! Smagorinsky !
   !*************!
   subroutine GatherSmagorinsky
      implicit none
      call e%gather(e%ndime,elvel,a%veloc(:,:,1))
   end subroutine

   subroutine InGaussSmagorinsky 
      implicit none
      call e%elmlen
      if (a%kfl_cotur == -1) then   
         call nsm_smago(e,grbve,acden,a%turbu(1),vista)
      elseif (a%kfl_cotur == -2) then
         call nsm_wale(e,grbve,acden,a%turbu(1),vista)
      endif
      acvis = acvis + vista
   end subroutine

   !*************************!
   ! Non-Newtonian Viscosity !
   !*************************!
   subroutine InGaussNonNewtonianViscosity
      implicit none
      call nsi_vislaw(e,grbve,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis)
   end subroutine

end module

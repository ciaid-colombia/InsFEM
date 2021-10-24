module Mod_nsm_PhysicalProperties
   use typre
   use Mod_PointerSetter
   use Mod_nsm_Viscosity
   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope
   use Mod_nsm_InterpolateGradients
   implicit none
   private
   public SetPointersPhysicalProperties
   
   type, extends(PointerSetter) :: SPPhysicalProperties
contains
      procedure :: SpecificSet => SpecificSetPhysicalProperties
      procedure :: SetTask => SpecificSetPhysicalPropertiesTask
   end type
   type(SPPhysicalProperties) :: SetPointersPhysicalProperties

   character(6) :: task
   
contains
   
   subroutine SpecificSetPhysicalProperties(d)
      implicit none
      class(SPPhysicalProperties) :: d
         
      !Viscosity law
      if (a%MatProp(imat)%lawvi /= 0) then
         if (task .eq. 'Elmope') then
            call SetPointersInterpolateGradients%Set
            call ConcatenateProcedures(ProcHook_PhysicalProp,ViscosityLawElmope)
         elseif (task .eq. 'EndElm') then   
            call ConcatenateProcedures(ProcHook_PhysicalProp,ViscosityLawEndElmope)
         elseif (task .eq. 'Bouope') then   
            call ConcatenateProcedures(ProcHook_PhysicalProp,ViscosityLawBouope)
         endif
      endif 
   end subroutine
 
   subroutine SpecificSetPhysicalPropertiesTask(d,intask)
      implicit none
      class(SPPhysicalProperties) :: d
      character(6) :: intask
      task = intask
   end subroutine   
   
   !---------------------------------------------------------------------------
   subroutine ViscosityLawEndElmope 
      implicit none 
      ! is recomended use the same viscosity compute in elmope
       
      real(rp) :: auxelvel(e%ndime,e%pnode)
      real(rp) :: grvelpre(e%ndime,e%ndime)
      
      auxelvel = a%veloc(:,e%lnods(1:e%pnode),2)
      call e%gradient(e%ndime,auxelvel,grvelpre) 
      call nsi_vislaw(e,grvelpre,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis) 
   end subroutine 
   
   subroutine ViscosityLawElmope 
      implicit none   
      call nsi_vislaw(e,grvel,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis) 
      if ((a%kfl_repro >= 1).or.(a%npp_stepi(5)==1)) then
         a%viscarray(ielem)%a(e%igaus) = acvis
      endif
   end subroutine
   
   subroutine ViscosityLawBouope 
      implicit none   
      call nsi_vislaw(e,grbve,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis) 
   end subroutine

end module

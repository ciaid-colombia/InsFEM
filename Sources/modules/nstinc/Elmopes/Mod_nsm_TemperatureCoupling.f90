module Mod_nsm_TemperatureCoupling
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_ExternalForces
   use Mod_nsi_BoussinesqForces
   implicit none
   private
   public :: SetPointersTemperatureCoupling

   type, extends(PointerSetter) :: SPTemperatureCoupling
contains
      procedure :: SpecificSet => SpecificSetTemperatureCoupling
   end type
   type(SPTemperatureCoupling) :: SetPointersTemperatureCoupling

   real(rp), allocatable :: eltem(:)
   real(rp) :: gptem(1)

contains

   subroutine SpecificSetTemperatureCoupling(d)
      implicit none
      class(SPTemperatureCoupling) :: d
         
      if (a%kfl_cotem == 1) then
         call SetPointersExternalForces%Set
         !Boussinesq Coupling
         call ConcatenateProcedures(ProcHook_Initializations,AllocTem)
         call ConcatenateProcedures(ProcHook_Gathers,GatherTem)
         call ConcatenateProcedures(ProcHook_Interpolates,InterpolateTem)
         if (a%kfl_ExternalTemperatureSGS == 1) then
            call ConcatenateProcedures(ProcHook_Interpolates,InterpolateTemSGS)
         endif
         call ConcatenateProcedures(ProcPointer_ExternalForces,BoussinesqForces)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocTem)
      endif
   end subroutine
   
   !-------------------------------------------------------------------
   subroutine AllocTem
      implicit none
      call a%Memor%alloc(e%mnode,eltem,'eltem','nsm_elmope')
   end subroutine
   
   subroutine DeallocTem
      implicit none
      call a%Memor%dealloc(e%mnode,eltem,'eltem','nsm_elmope')
   end subroutine
   
   subroutine GatherTem
      implicit none
      call e%gather(1,eltem,a%tempe)
   end subroutine
   
   subroutine InterpolateTem
      implicit none
      call e%interpg(1,eltem,gptem(1))
   end subroutine   
   
   subroutine InterpolateTemSGS
      implicit none
      gptem = gptem + a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
   
   subroutine BoussinesqForces
      implicit none
      call nsi_BoussinesqForces(e,acden,a%bougr,a%boube,a%boutr,gptem(1),elext)
   end subroutine  
   
end module

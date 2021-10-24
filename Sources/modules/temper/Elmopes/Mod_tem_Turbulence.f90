module Mod_tem_Turbulence
   use typre
   use mod_tem_BaseElmope
   use Mod_tem_VelocityGradient
   use Mod_nsm_Viscosity
   implicit none
   private
   public SetPointersTurbulence

   integer(ip), allocatable :: kfl_IsSet
   real(rp), allocatable    :: testf_static(:)
   
contains
   
   !Set Pointers
   subroutine SetPointersTurbulence(itask)
      implicit none
      integer(ip) :: itask
      
      integer(ip) :: kfl_nonlinear
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            !------------------------------------------------------------------
            !LES Turbulence models
            if (a%kfl_cotur < 0) then
               call SetPointersVelocityGradient(1)
               
               !Smagorinsky
               if (a%kfl_cotur == -1) then
                  call ConcatenateProcedures(ProcHook%PhysicalProp,Smagorinsky)
            
               !Wale
               elseif (a%kfl_cotur == -2) then
                  call ConcatenateProcedures(ProcHook%PhysicalProp,Wale)
               endif
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
   !-------------------------------------------------------------------
   !Turbulence modelling
   subroutine Smagorinsky
      implicit none
      real(rp) :: vista
      
      call nsm_smago(e,grvel,acden,a%turbu,vista)
      !Prandtl turbulent number
      vista = vista/a%prtur
      acvis = acvis + vista
   end subroutine
   
   subroutine Wale
      implicit none
      real(rp) :: vista
      
      call nsm_wale(e,grvel,acden,a%turbu,vista)
      !Prandtl turbulent number
      vista = vista/a%prtur
      acvis = acvis + vista
   end subroutine
end module

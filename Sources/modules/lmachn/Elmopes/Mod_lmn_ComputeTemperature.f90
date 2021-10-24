module Mod_lmn_ComputeTemperature
   use typre
   use Mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersTemperature
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersTemperature(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            ProcPointer%ComputeTemperature => ComputeTemperature
            
            !NonLinear Subscales
            if (a%kfl_nolsg == 1) then
               call ConcatenateProcedures(ProcPointer%ComputeTemperature,NonLinearSGSTemperature)
            endif
            if (a%kfl_eqnst == 1) then
               call ConcatenateProcedures(ProcPointer%ComputeTemperature,ComputeThermalExpansion)
            end if
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Actual Computations
   !Temperature
   subroutine ComputeTemperature
      implicit none      
      gtemp = gptem(1) 
   end subroutine
   
   subroutine NonLinearSGSTemperature
      implicit none 
      gtemp = gtemp + a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
   
   subroutine ComputeThermalExpansion
      implicit none      
      actex = 1.0_rp /gtemp/acsph 
   end subroutine

end module

module Mod_supm_TemporalDerivatives
   use Mod_supm_BaseElmope
   use Mod_ThreeFieldElement
   use Mod_LogarithmicElement
   implicit none
   private
   public SetPointersTemporalDerivatives
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersTemporalDerivatives(itask)
      integer(ip) :: itask
      !procedure() :: NULL()
      
      select case (itask)   
         case(0)
            allocate(kfl_IsSet)
            call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
            kfl_IsSet = -1
      
         case(1)
            if (kfl_IsSet == -1) then
               kfl_IsSet = 1
                if (a%LogFormulation==0) then
                  ProcPointer%TemporalDerivatives => sup_temporalpart
                else if (a%LogFormulation==1) then
                  ProcPointer%TemporalDerivatives => sup_temporalpart_LCR
                end if
            endif
            
         case(100)
            deallocate(kfl_IsSet)
            call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
       end select   
   end subroutine  
       
   
   subroutine sup_temporalpart
      implicit none
      call sup_TimeIntegrationToElext(e,auxtens,integrator,auxVE,a%dtinv,gpsig,elextS)  
   end subroutine
   
   subroutine sup_temporalpart_LCR
      implicit none
      call sup_LCR2_TimeIntegrationToElext(e,auxtens,nsteps,integrator,a%dtinv,lambda,lambda0,gpsig,elextS)
   end subroutine
   
end module

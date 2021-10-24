 module Mod_tem_ComputeTaus
   use typre
   use Mod_tem_BaseElmope
   implicit none
   private
   public SetPointersComputeTaus, timom_static
   
   !Tau Smoothing
   real(rp), allocatable :: elSmoothedTau(:,:)
   
   !TransientSubgridScales
   real(rp) :: timom_static
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersComputeTaus(itask)
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
         
            !Taus are computed as usual
            !if (a%kfl_Tausm == 0 .or. (a%istep == 1)) then
               ProcPointer%ComputeTau => ComputeStaticTau
            
            !Taus are interpolated from an elementary array
            !elseif (a%kfl_Tausm == 1) then
            !   call ConcatenateProcedures(ProcHook%Initializations,AllocSmoothedTau)
            !   call ConcatenateProcedures(ProcHook%Finalizations,DeallocSmoothedTau)
            !   call ConcatenateProcedures(ProcHook%Gathers,GatherTau)
            !   call ConcatenateProcedures(ProcHook%ComputeTaus,InterpolateTau)
            !endif
            
            !Dynamic subgrid scales
            if (a%kfl_tacsg == 1) then
               ProcPointer%ComputeTau => TransientTau
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
  !--------------------------------------------------------------------
   !For computing Tau
   subroutine ComputeStaticTau
      call ComputeTauCDR(e,acden,acvis,acrcp,gpvno,a%staco,chale,timom)
   end subroutine   
   
   subroutine TransientTau
      implicit none
      call ComputeTauCDR(e,acden,acvis,acrcp,gpvno,a%staco,chale,timom)
      !timom is the transient one, the static one goes to timom_static
      timom_static = timom
      call ComputeTransientTau(ReferenceDtinv,acden,timom_static,timom)
   end subroutine
end module

 

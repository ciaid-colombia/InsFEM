 module Mod_nsc_ComputeTaus
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersComputeTaus, timom_static
   
   !External Procedures
   procedure() :: NULLSUB

   !TransientSubgridScales
   real(rp) :: timom_static(3)
   
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
            ProcHook_nsc_ComputeTaus => ComputeTaus
            
            if (a%kfl_tacsg == 1) then
               call ConcatenateProcedures(ProcHook_nsc_ComputeTaus,TransientTaus)
            endif

            !Sources Exists
            if (a%kfl_sourc == 1) then
               call ConcatenateProcedures(ProcHook_nsc_ComputeTaus,SourcesTaus)
            endif

         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !Compute Tau values
   subroutine ComputeTaus
     implicit none
      
     call nsc_ComputeTau(e,gpspd,acvis*invgpd,actco*invgpd/accph,gpmno*invgpd,a%staco,chale,timom)
    
   end subroutine
   
   !Computes the transient stabilization parameter
   subroutine TransientTaus
     !timom is the transient one, the static one goes to timom_static
     timom_static = timom
     call ComputeTransientTau(ReferenceDtinv,1.0_rp,timom_static(1),timom(1))
     call ComputeTransientTau(ReferenceDtinv,1.0_rp,timom_static(2),timom(2))
     call ComputeTransientTau(ReferenceDtinv,1.0_rp,timom_static(3),timom(3))
   end subroutine

   subroutine SourcesTaus

     call nsc_ComputeSourcesTau(e,gpspd,acvis*invgpd,actco*invgpd/accph,gpmno*invgpd,a%grnor,a%srce,a%staco,chale,timom)

   end subroutine
end module

 

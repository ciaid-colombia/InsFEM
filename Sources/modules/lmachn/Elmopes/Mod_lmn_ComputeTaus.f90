module Mod_lmn_ComputeTaus
   use typre
   use Mod_lmn_BaseElmope
!   use Mod_lmn_ComputeTauSmoothing
   implicit none
   private
   public SetPointersComputeTaus, timom_static, tiene_static
   
   !Tau Smoothing
   real(rp), allocatable :: elSmoothedTau(:,:)
   
   !TransientSubgridScales
   real(rp) :: timom_static, tiene_static
   
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
!            if (a%kfl_Tausm == 0 .or. (a%istep == 1)) then
               call ConcatenateProcedures(ProcHook%ComputeTaus,ComputeTaus_mom)
               call ConcatenateProcedures(ProcHook%ComputeTaus,ComputeTaus_ene)
            
            !Taus are interpolated from an elementary array
!            elseif (a%kfl_Tausm == 1) then
!               call ConcatenateProcedures(ProcHook%Initializations,AllocSmoothedTau)
!               call ConcatenateProcedures(ProcHook%Finalizations,DeallocSmoothedTau)
!               call ConcatenateProcedures(ProcHook%Gathers,GatherTau)
!               call ConcatenateProcedures(ProcHook%ComputeTaus,InterpolateTau)
!            endif
            
            if (a%kfl_tacsg == 1) then
               call ConcatenateProcedures(ProcHook%ComputeTaus,TransientTaus_mom)
               call ConcatenateProcedures(ProcHook%ComputeTaus,TransientTaus_ene)
               if (a%kfl_nolsg == 1) then
                  select case (a%kfl_nolsgScheme)
                     case (1)
                     case (2)
                     call ConcatenateProcedures(ProcHook%ComputeTaus_seg,ComputeTaus_ene)
                     call ConcatenateProcedures(ProcHook%ComputeTaus_seg,TransientTaus_ene)
                     case (3)
                     call ConcatenateProcedures(ProcHook%ComputeTaus_seg,ComputeTaus_ene)
                     call ConcatenateProcedures(ProcHook%ComputeTaus_seg,TransientTaus_ene)
                     case default
                  end select
               end if
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !Computation of Taus
   subroutine ComputeTaus_mom
      implicit none
      !Tau momento
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)
      !Tau continuity
      call ComputeTauCon(e,acden,a%staco(1),chale,timom,ticon)
   end subroutine
   
   subroutine ComputeTaus_ene
      implicit none
      !Tau energy
      call ComputeTau(e,acden,actco,gpvno,a%staco,chale,tiene)
!     call ComputeTauCDR(e,acden,actco,acrea,gpvno,a%staco,chale,tiene)
   end subroutine
   
!   !To interpolate the tau values from smoothed array Tausmo
!   subroutine AllocSmoothedTau
!      implicit none
!      
!      call a%Memor%alloc(2,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
!   end subroutine
!   
!   subroutine DeallocSmoothedTau
!      implicit none
!      
!      call a%Memor%dealloc(2,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
!   end subroutine
!   
!   subroutine GatherTau
!      implicit none
!      
!      call e%gather(2,elSmoothedTau,a%Tausmo)
!   end subroutine
!   
!   subroutine InterpolateTau
!      implicit none
!      real(rp) :: gptau(2)
!      
!      call e%interpg(2,elSmoothedTau,gptau)
!      timom = gptau(1)
!      tidiv = gptau(2)
!   end subroutine
   
   !Computes the transient stabilization parameter
   subroutine TransientTaus_mom
      !timom is the transient one, the static one goes to timom_static
      timom_static = timom
      call ComputeTransientTau(ReferenceDtinv,acden,timom_static,timom)
   end subroutine 

   subroutine TransientTaus_ene
      !timom is the transient one, the static one goes to timom_static
      tiene_static = tiene
      call ComputeTransientTau(ReferenceDtinv,acden,tiene_static,tiene)
   end subroutine 
end module

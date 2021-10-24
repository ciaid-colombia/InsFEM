 module Mod_nsm_ComputeTaus
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_ComputeTauSmoothing
   implicit none
   private
   public SetPointersComputeTaus, timom_static, tidiv_static
   
   type, extends(PointerSetter) :: SPComputeTaus
contains
      procedure :: SpecificSet => SpecificSetComputeTaus
   end type
   type(SPComputeTaus) :: SetPointersComputeTaus
   
   !Tau Smoothing
   real(rp), allocatable :: elSmoothedTau(:,:)
   !TransientSubgridScales
   real(rp) :: timom_static, tidiv_static
   !For ALE
   logical  :: isALE
   real(rp) :: gpvelnor
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetComputeTaus(d)
      implicit none
      class(SPComputeTaus) :: d
            
      call a%Mesh%GetALE(isALE)
      !Taus are computed as usual
      if (a%kfl_Tausm == 0 .or. (a%istep == 1)) then
         ProcHook_ComputeTaus => ComputeTaus
      !Taus are interpolated from an elementary array
      elseif (a%kfl_Tausm >= 1) then
         call ConcatenateProcedures(ProcHook_Initializations,AllocSmoothedTau)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocSmoothedTau)
         call ConcatenateProcedures(ProcHook_Gathers,GatherTau)
         call ConcatenateProcedures(ProcHook_ComputeTaus,InterpolateTau)
      endif
      if (a%kfl_tacsg == 1) then
         call ConcatenateProcedures(ProcHook_ComputeTaus,TransientTaus)
      endif
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !Compute Tau values
   subroutine ComputeTaus
      implicit none
      
      !If Fixed -MeshALE or ALE, compute tau with the real velocity (u, instead of u-vmesh)
      if (isALE) then
         call vecnor(gpvel,e%ndime,gpvelnor,2)
         call ComputeTau(e,acden,acvis,gpvelnor,a%staco,chale,timom)
      else
         call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)
      endif
      if (a%kfl_advco == 1) then 
         call vecnor(a%advco,e%ndime,gpvelnor,2)
         call ComputeTau(e,acden,acvis,gpvelnor,a%staco,chale,timom)
      endif
      call ComputeTauDiv(e,a%staco(4),chale,timom,tidiv)
   end subroutine
   
   !To interpolate the tau values from smoothed array Tausmo
   subroutine AllocSmoothedTau
      implicit none
      
      call a%Memor%alloc(2,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
   end subroutine
   
   subroutine DeallocSmoothedTau
      implicit none
      
      call a%Memor%dealloc(2,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
   end subroutine
   
   subroutine GatherTau
      implicit none
      
      call e%gather(2,elSmoothedTau,a%Tausmo)
   end subroutine
   
   subroutine InterpolateTau
      implicit none
      real(rp) :: gptau(2)
      
      call e%interpg(2,elSmoothedTau,gptau)
      timom = gptau(1)
      tidiv = gptau(2)
   end subroutine
   
   !Computes the transient stabilization parameter
   subroutine TransientTaus
      !timom is the transient one, the static one goes to timom_static
      timom_static = timom
      call ComputeTransientTau(ReferenceDtinv,acden,timom_static,timom)
   end subroutine
   
end module

 

module Mod_nsm_DynSubsElmope
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersDynSubsElmope
   
   type, extends(PointerSetter) :: SPDynSubsElmope
contains
      procedure :: SpecificSet => SpecificSetDynSubsElmope
   end type
   type(SPDynSubsElmope) :: SetPointersDynSubsElmope
   
   real(rp), allocatable :: testf_static(:)

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetDynSubsElmope(d)
      implicit none
      class(SPDynSubsElmope) :: d
         
      !Tracking of the subscales
      !Dynamic subscales
      if (a%kfl_tacsg == 1) then
         call ConcatenateProcedures(ProcHook_Initializations,AllocDSS)
         call ConcatenateProcedures(ProcHook_Finalizations, DeallocDSS)
         
         if (a%kfl_repro == 1) then
            !This adds the subscale at the previous time step to the RHS
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_OSS)
         
         elseif (a%kfl_repro == 0) then
            !If the subscales are not orthogonal, we also need to add the term
            !(v,\partial_t \tilde(u))
            !See the documentation folder
            call ConcatenateProcedures(ProcHook_ComputeTestf,DynamicSgsTestf_ASGS)
            
            !This adds the terms involving the subscale at the previous time step to the RHS
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_ASGS)
         elseif (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
         
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS)
         else
            call runend('nsm_dynsubsElmope not ready for kfl_repro not 0, 1, 2 or 3')
            
         endif
      endif
   end subroutine   

   !-------------------------------------------------------------------
   !Computations subroutines
   !Dynamic subscales
   subroutine AllocDSS
      implicit none
      
      call a%Memor%alloc(e%mnode,testf_static,'testf_static','nsm_elmope')
   end subroutine
   
   subroutine DeAllocDSS
      implicit none
      
      call a%Memor%dealloc(e%mnode,testf_static,'testf_static','nsm_elmope')
   end subroutine
   
   subroutine InGaussElmatsDSS_OSS
      use Mod_nsm_ComputeTestf
      use Mod_nsm_ComputeTaus
      use typre
      implicit none
      integer(ip) :: inode
      
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      !Compute contributions to RHS : Block U
      call nsm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,testf,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhu)
            
      !Compute contributions to RHS : Block P
      call nsm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhp)
   end subroutine
   
   subroutine InGaussElmatsDSS_ASGS
      use typre
      use Mod_nsm_ComputeTestf
      use Mod_nsm_ComputeTaus
      implicit none
      real(rp) :: aux_testf(e%pnode)
      integer(ip) :: inode
      !tau_t*(-L*(v_h)+v_h*tau_k^-1*tau_t,ro tessgs_n/dtime
      
      aux_testf(1:e%pnode) = testf_static(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*timom/timom_static
      
      !We add the contribution of the subscales at the previous time step 
     
      !Compute contributions to RHS : Block U
      call nsm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,aux_testf,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)
      
      !Compute contributions to RHS : Block P
      call nsm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs(ielem)%a(:,2,e%igaus),elrhp)
   end subroutine
   
    subroutine InGaussElmatsDSS_SOSS
      use Mod_nsm_ComputeTestf
      use Mod_nsm_ComputeTaus
      use typre
      implicit none
      integer(ip) :: inode
      
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      
      !Compute contributions to RHS : Block U
      !first subscale of velocity (convective component)
      call nsm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,testf,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhu)

            
      !Compute contributions to RHS : Block P
      !second subscale of velocity (pressure gradient component)
      call nsm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs2(ielem)%a(1:e%ndime,2,e%igaus),elrhp)
   end subroutine
   
   subroutine DynamicSgsTestf_ASGS
      use Mod_nsm_ComputeTestf
      use Mod_nsm_ComputeTaus
      !Here we modify the test function so that it takes into account
      !tau_t*(-L*(v_h))-v_h[1-tau_t*tau_k^-1],-R)

      testf_static(1:e%pnode) = testf(1:e%pnode)
      testf(1:e%pnode) = testf(1:e%pnode) - e%shape(1:e%pnode,e%igaus)*(1-timom/timom_static)
   end subroutine  
   
end module

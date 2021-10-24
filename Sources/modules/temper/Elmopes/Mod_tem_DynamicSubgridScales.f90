module Mod_tem_DynamicSubgridScales
   use typre
   use mod_tem_BaseElmope
   use Mod_tem_ComputeTaus
   implicit none
   private
   public SetPointersDynamicSubgridScales

   integer(ip), allocatable :: kfl_IsSet
   
   real(rp), allocatable :: testf_static(:)
   
contains
   
   !Set Pointers
   subroutine SetPointersDynamicSubgridScales(itask)
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
         
             !---------------------------------------------------------------
            !Dynamic Subgrid Scales
            !See the documentation folder
            if (a%kfl_tacsg == 1) then
            
               call ConcatenateProcedures(ProcHook%Initializations,AllocDSS)
               call ConcatenateProcedures(ProcHook%Finalizations, DeallocDSS)
               
               !Tau is changed to the transient one
               !this is done in the computeTaus module
               call SetPointersComputeTaus(1)
               
               if (a%kfl_repro == 1) then
                  !This adds the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsDSS_OSS)
               
               elseif (a%kfl_repro == 0) then
                  !If the subscales are not orthogonal, we also need to add the term
                  !(v,\partial_t \tilde(u))
                  !See the documentation folder
                  call ConcatenateProcedures(ProcHook%ComputeTestf,DynamicSgsTestf_ASGS)
                  
                  !This adds the terms involving the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsDSS_ASGS)
                  
               endif
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
    !-------------------------------------------------------------------
   !TRANSIENT SUBGRID SCALES
   subroutine AllocDSS
      implicit none
      
      call a%Memor%alloc(e%mnode,testf_static,'testf_static','tem_elmope')
   end subroutine
   
   subroutine DeAllocDSS
      implicit none
      
      call a%Memor%dealloc(e%mnode,testf_static,'testf_static','tem_elmope')
   end subroutine
   
   subroutine InGaussElmatsDss_OSS
      use typre
      implicit none
      real(rp) :: aux
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      
      !We add the contribution of the subscales at the previous time step 
      !Compute contributions to RHS 
      call tem_elmrhs_dss(e,dvol,acden,ReferenceDtinv,testf,a%tesgs(ielem)%a(2,e%igaus),elrhs)
      
   end subroutine
   
   subroutine InGaussElmatsDss_ASGS
      use typre
      implicit none
      real(rp) :: aux_testf(e%pnode)

      !tau_t*(-L*(v_h)+v_h*tau_k^-1*tau_t,ro tessgs_n/dtime
      
      aux_testf(1:e%pnode) = testf_static(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*timom/timom_static
      !We add the contribution of the subscales at the previous time step 
      !Compute contributions to RHS 
      call tem_elmrhs_dss(e,dvol,acden,ReferenceDtinv,aux_testf,a%tesgs(ielem)%a(2,e%igaus),elrhs)
      
   end subroutine
      
   subroutine DynamicSgsTestf_ASGS
      !Here we modify the test function so that it takes into account
      !tau_t*(-L*(v_h))-v_h[1-tau_t*tau_k^-1],-R)

      testf_static(1:e%pnode) = testf(1:e%pnode)
      testf(1:e%pnode) = testf(1:e%pnode) - e%shape(1:e%pnode,e%igaus)*(1-timom/timom_static)
   end subroutine

end module

module Mod_nsc_DynSubsElmope
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersDynSubsElmope,nsc_elmrhd_dss,nsc_elmrhm_dss,nsc_elmrhe_dss
   
   integer(ip), allocatable :: kfl_IsSet
   
   !Transient subscales RHS 
   real(rp)    :: SGSrhd,SGSrhe
   real(rp), allocatable :: SGSrhm(:)
   !Test functions
   !Adjoint Matrix coefficients
   real(rp), allocatable :: LTdd_static(:), LTmm_static(:,:,:),   LTee_static(:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersDynSubsElmope(itask)
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
         
            !Tracking of the subscales
            !See the documentation folder
            !Dynamic subscales
            if (a%kfl_tacsg == 1) then
            
               call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocDSS)
               call ConcatenateProcedures(ProcHook_nsc_Finalizations, DeallocDSS)
               
               if (a%kfl_repro == 1) then
                  !This adds the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,InGaussElmatsDSS_OSS)
               
               elseif (a%kfl_repro == 0) then
                  !If the subscales are not orthogonal, we also need to add the term
                  !(v,\partial_t \tilde(u))
                  !See the documentation folder
                  call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,DynamicSgsTestf_ASGS)
                  
                  !This adds the terms involving the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,InGaussElmatsDSS_ASGS)
               else
                  call runend('dynsubs not ready for kfl_repro not 0 or 1')
               endif
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   

   !-------------------------------------------------------------------
   !Computations subroutines
   !Dynamic subscales
   subroutine AllocDSS
      implicit none
      
      call a%Memor%alloc(e%ndime,SGSrhm,'SGSrhm','nsc_elmope')
      call a%Memor%alloc(e%mnode,LTdd_static,'LTdd_static','nsc_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,LTmm_static,'LTmm_static','nsc_elmope')
      call a%Memor%alloc(e%mnode,LTee_static,'LTee_static','nsc_elmope')

   end subroutine
   
   subroutine DeAllocDSS
      implicit none
      
      call a%Memor%dealloc(e%ndime,SGSrhm,'SGSrhm','nsc_elmope')
      call a%Memor%dealloc(e%mnode,LTdd_static,'LTdd_static','nsc_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,LTmm_static,'LTmm_static','nsc_elmope')
      call a%Memor%dealloc(e%mnode,LTee_static,'LTee_static','nsc_elmope')
   end subroutine
   
   subroutine InGaussElmatsDss_OSS
      use typre
      implicit none
      real(rp)    :: auxarray(e%ndime)

      auxarray=0.0_rp
      
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),Ao Usgs_n/dtime)

      !Compute contributions to RHS : Continuity
       call nsc_elmrhd_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTdd,LTdm,LTde,elrhd)
            
      !Compute contributions to RHS : Momentum
       call nsc_elmrhm_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTmd,LTmm,LTme,elrhm)

      !Compute contributions to RHS : Energy
       call nsc_elmrhe_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTed,LTem,LTee,elrhe)

   end subroutine
   
   subroutine InGaussElmatsDss_ASGS
      use typre
      use Mod_nsc_ComputeTaus
      implicit none
      integer(ip) :: idime
      real(rp)    :: auxarray(e%ndime)

      auxarray=0.0_rp

      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h)+v_h*tau_k^-1*tau_t,Ao Usgs_n/dtime)
      
      LTdd_static(1:e%pnode) = LTdd_static(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*timom(1)/timom_static(1)
      do idime=1,e%ndime
         LTmm_static(idime,idime,:) = LTmm_static(idime,idime,:) &
                           + e%shape(1:e%pnode,e%igaus)*timom(2)/timom_static(2)
      end do
      LTee_static(1:e%pnode) = LTee_static(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*timom(3)/timom_static(3)
     
      !Compute contributions to RHS : Continuity
       call nsc_elmrhd_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTdd_static,LTdm,LTde,elrhd)
            
      !Compute contributions to RHS : Momentum
       call nsc_elmrhm_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTmd,LTmm_static,LTme,elrhm)

      !Compute contributions to RHS : Energy
       call nsc_elmrhe_dss(e,dvol,ReferenceDtinv,a%cosgs(ielem)%a(2,e%igaus),a%mosgs(ielem)%a(1:e%ndime,2,e%igaus),a%ensgs(ielem)%a(2,e%igaus),LTed,LTem,LTee_static,elrhe)

   end subroutine
      
   subroutine DynamicSgsTestf_ASGS
      use typre
      use Mod_nsc_ComputeTaus
      implicit none
      integer(ip)                :: idime
      !Here we modify the test function so that it takes into account
      !tau_t*(-L*(v_h))-v_h[1-tau_t*tau_k^-1],-R)

      LTdd_static = LTdd !LTee(p)
      LTmm_static = LTmm !LTmm(i,d,p)
      LTee_static = LTee !LTee(p)
      
      LTdd(1:e%pnode) = LTdd(1:e%pnode) - e%shape(1:e%pnode,e%igaus)*(1-timom(1)/timom_static(1))
      do idime=1,e%ndime
         LTmm(idime,idime,:) = LTmm(idime,idime,:) &
                       - e%shape(1:e%pnode,e%igaus)*(1-timom(2)/timom_static(2))
      end do
      LTee(1:e%pnode) = LTee(1:e%pnode) - e%shape(1:e%pnode,e%igaus)*(1-timom(3)/timom_static(3))

   end subroutine  
   
   
   subroutine nsc_elmrhd_dss(e,dvolu,dtinv,SGSrhd,SGSrhm,SGSrhe,LTdd,LTdm,LTde,elrhd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for mass equation for DSS
      !    + Taut_d(q·L*dd, Atdd·preSGS_n/dt + Atde·temSGS_n/dt) 
      !    + Taut_m(q·L*dm, Atmd·preSGS_n/dt + Atmm·velSGS_n/dt + Atme·temSGS_n/dt)
      !    + Taut_e(q·L*de, Ated·preSGS_n/dt + Atem·velSGS_n/dt + Atee·temSGS_n/dt)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: SGSrhd,SGSrhe
      real(rp),    intent(in)    :: SGSrhm(e%ndime)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elrhd(e%mnode)
      real(rp)                   :: aux(e%mnode)

      aux =  matmul(SGSrhm,LTdm) 

      elrhd(1:e%pnode) = elrhd(1:e%pnode) + &
       (LTdd(1:e%pnode)*SGSrhd + aux(1:e%pnode) + LTde(1:e%pnode)*SGSrhe)*dvolu*dtinv
   end subroutine
   
   subroutine nsc_elmrhm_dss(e,dvolu,dtinv,SGSrhd,SGSrhm,SGSrhe,LTmd,LTmm,LTme,elrhm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for momentum equation for DSS
      !    + Taut_d(n·L*md, Atdd·preSGS_n/dt + Atde·temSGS_n/dt)
      !    + Taut_m(n·L*mm, Atmd·preSGS_n/dt + Atmm·velSGS_n/dt + Atme·temSGS_n/dt) 
      !    + Taut_e(n·L*me, Ated·preSGS_n/dt + Atem·velSGS_n/dt + Atee·temSGS_n/dt) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: SGSrhd,SGSrhe
      real(rp),    intent(in)    :: SGSrhm(e%ndime)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elrhm(e%ndime,e%mnode)

      integer(ip)                :: idime
      real(rp)                   :: aux(e%mnode)

      do idime=1,e%ndime
         aux(:) =  matmul(transpose(LTmm(idime,:,:)),SGSrhm) 
         elrhm(idime,1:e%pnode) = elrhm(idime,1:e%pnode) + &
          (LTmd(idime,1:e%pnode)*SGSrhd + LTme(idime,1:e%pnode)*SGSrhe + &
           aux(1:e%pnode))*dvolu*dtinv 
      end do

   end subroutine 

   subroutine nsc_elmrhe_dss(e,dvolu,dtinv,SGSrhd,SGSrhm,SGSrhe,LTed,LTem,LTee,elrhe)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for energy for DSS
      !    + Taut_d(g·L*ed, Atdd·preSGS_n/dt + Atde·temSGS_n/dt)  
      !    + Taut_m(g·L*em, Atmd·preSGS_n/dt + Atmm·velSGS_n/dt + Atme·temSGS_n/dt)
      !    + Taut_e(g·L*ee, Ated·preSGS_n/dt + Atem·velSGS_n/dt + Atee·temSGS_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: SGSrhd,SGSrhe
      real(rp),    intent(in)    :: SGSrhm(e%ndime)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu,dtinv
      real(rp),    intent(inout) :: elrhe(e%mnode)
      real(rp)                   :: aux(e%mnode)

      aux =  matmul(SGSrhm,LTem) 

      elrhe(1:e%pnode) = elrhe(1:e%pnode) + &
       (LTee(1:e%pnode)*SGSrhe + LTed(1:e%pnode)*SGSrhd + aux(1:e%pnode))*dvolu*dtinv

   end subroutine 
   
end module

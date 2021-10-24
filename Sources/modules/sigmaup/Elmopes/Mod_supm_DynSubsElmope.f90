module Mod_supm_DynSubsElmope
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersDynSubsElmopeTF
   
   integer(ip), allocatable   :: kfl_IsSet
   real(rp), allocatable      :: testf_static(:)
   real(rp)                   :: auxVEpol
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersDynSubsElmopeTF(itask)
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
            
               call ConcatenateProcedures(ProcHook_Initializations,AllocDSS)
               call ConcatenateProcedures(ProcHook_Finalizations, DeallocDSS)
                
               if (a%kfl_repro == 1) then
                  !This adds the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_OSS) 
                  if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_OSS_VE) 

               
               elseif (a%kfl_repro == 0) then
                  !If the subscales are not orthogonal, we also need to add the term
                  !(v,\partial_t \tilde(u))
                  !See the documentation folder
                  call ConcatenateProcedures(ProcHook_ComputeTestf,DynamicSgsTestf_ASGS)
                  
                  !This adds the terms involving the subscale at the previous time step to the RHS

                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_ASGS) 
                   if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_ASGS_VE) 
                  
               elseif (a%kfl_repro == 2 .or. a%kfl_repro==3) then
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS_mixed)
                  if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS_mixed_VE)  
               
               elseif (a%kfl_repro==4) then
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS)
                  if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS_VE)
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
      call a%Memor%alloc(e%mnode,testf_static,'testf_static','nsm_elmope')
   end subroutine
   
   
   subroutine DeAllocDSS
      implicit none     
      call a%Memor%dealloc(e%mnode,testf_static,'testf_static','nsm_elmope')
   end subroutine
   

   subroutine InGaussElmatsDSS_OSS
      use Mod_supm_ComputeTestf
      use Mod_supm_ComputeTaus
      use typre
      implicit none
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      
      !Equations implemented in nstinc module and inherited here
      !Compute contributions to RHS : Block U
      call ProcPointer%supm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,testf,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhu)
   
      !Compute contributions to RHS : Block P
      call ProcPointer%supm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhp)
      
      !Compute contributions to RHS : Block S
      call supm_elmrhc_dss(e,ReferenceDtinv,dvol,acden,beta*(1-a%LogFormulation),timom,auxtens,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhc)

   end subroutine
   
   
   subroutine InGaussElmatsDSS_OSS_VE
      implicit none
      real(rp)    :: auxVEnew
       auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)
      
      !Compute contributions to RHS : Block U
      call supm_elmrhuVES_dss(e,auxVEnew,tisig,dvol,a%sisgs(ielem)%a(1:auxtens,2,e%igaus),auxtens,ReferenceDtinv,elrhu) 
      !Viscoelastic case, remainder terms
      call supm_elmrhcVES_dss(e,auxVEnew,beta,tisig,acvis,dvol,auxtens,AgradV_VE,grvel,a%sisgs(ielem)%a(1:auxtens,2,e%igaus),ReferenceDtinv,a%LogFormulation,elrhc)
      
   end subroutine    
   
   
   subroutine InGaussElmatsDSS_ASGS
      use typre
      use Mod_supm_ComputeTestf
      use Mod_supm_ComputeTaus
      implicit none
      real(rp) :: aux_testf(e%pnode), auxVEnew
      integer(ip) :: inode
      !tau_t*(-L*(v_h)+v_h*tau_k^-1*tau_t,ro tessgs_n/dtime
      
      aux_testf(1:e%pnode) = testf_static(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*timom/timom_static
      !We add the contribution of the subscales at the previous time step 
     
      !Compute contributions to RHS : Block U
      call ProcPointer%supm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,aux_testf,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)
      
      !Compute contributions to RHS : Block P
      call ProcPointer%supm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs(ielem)%a(:,2,e%igaus),elrhp)
      

      call supm_elmrhc_dss(e,ReferenceDtinv,dvol,acden,beta*(1-a%LogFormulation),timom,auxtens,a%vesgs(ielem)%a(:,2,e%igaus),elrhc)

   end subroutine
   
   subroutine InGaussElmatsDSS_ASGS_VE
      implicit none
      real(rp) :: aux_testf(e%pnode), auxVEnew
      
      auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)
      
      !Compute contributions to RHS : Block U
      call supm_elmrhuVES_dss(e,auxVEnew,tisig,dvol,a%sisgs(ielem)%a(:,2,e%igaus),auxtens,ReferenceDtinv,elrhu) 
      !Viscoelastic case, remainder terms
      call supm_elmrhcVES_dss(e,auxVEnew,beta,tisig,acvis,dvol,auxtens,AgradV_VE,grvel,a%sisgs(ielem)%a(:,2,e%igaus),ReferenceDtinv,a%LogFormulation,elrhc)
      
   end subroutine   
   
   
   
   subroutine InGaussElmatsDSS_SOSS_mixed
      use Mod_supm_ComputeTestf
      use Mod_supm_ComputeTaus
      use typre
      implicit none
      integer(ip) :: inode
 
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      !Compute contributions to RHS : Block U
      !first subscale of velocity (convective component)
      call ProcPointer%supm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,testf,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)

      !Compute contributions to RHS : Block P
      !second subscale of velocity (pressure gradient component)
      call ProcPointer%supm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs2(ielem)%a(:,2,e%igaus),elrhp)

      !Compute contributions to RHS : Block S
      !third subscale of velocity (sigma divergence component)
      call supm_elmrhc_dss(e,ReferenceDtinv,dvol,acden,beta*(1-a%LogFormulation),timom,auxtens,a%vesgs3(ielem)%a(:,2,e%igaus),elrhc)
       
   end subroutine
   
   subroutine InGaussElmatsDSS_SOSS_mixed_VE
      implicit none
      real(rp)    :: auxVEnew
      auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)

      !Compute contributions to RHS : Block U
      call supm_elmrhuVES_dss(e,auxVEnew,tisig,dvol,a%sisgs(ielem)%a(1:auxtens,2,e%igaus),auxtens,ReferenceDtinv,elrhu) 

      !Viscoelastic case, remainder terms
      call supm_elmrhcVES_dss(e,auxVEnew,beta,tisig,acvis,dvol,auxtens,AgradV_VE,grvel,a%sisgs(ielem)%a(1:auxtens,2,e%igaus),ReferenceDtinv,a%LogFormulation,elrhc) 
   end subroutine   
   
   
   
   subroutine InGaussElmatsDSS_SOSS
      use Mod_supm_ComputeTestf
      use Mod_supm_ComputeTaus
      use typre
      implicit none
      integer(ip) :: inode

      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      !Compute contributions to RHS : Block U
      !first subscale of velocity (convective component)
      call ProcPointer%supm_elmrhu_dss(e,dvol,ReferenceDtinv,acden,testf,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)

      !Compute contributions to RHS : Block P
      !second subscale of velocity (pressure gradient component)
      call ProcPointer%supm_elmrhp_dss(e,dvol,ReferenceDtinv,acden,timom,a%vesgs2(ielem)%a(:,2,e%igaus),elrhp)
      
      !Compute contributions to RHS : Block S
      !third subscale of velocity (sigma divergence component)
      call supm_elmrhc_dss(e,ReferenceDtinv,dvol,acden,beta*(1-a%LogFormulation),timom,auxtens,a%vesgs3(ielem)%a(:,2,e%igaus),elrhc)

   end subroutine
   
   
   
   subroutine InGaussElmatsDSS_SOSS_VE
      implicit none
      real(rp)    :: auxVEnew
      auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)

      !Compute contributions to RHS : Block U
      !first subscale of stress (symmetric gradient of velocity)
      call supm_elmrhuVES_dss(e,auxVEnew,tisig,dvol,a%sisgs(ielem)%a(:,2,e%igaus),auxtens,ReferenceDtinv,elrhu) 
      !second subscale of stress (convective stress term)
      call supm_elmrhcVES_dss_split2(e,auxVEnew,tisig,dvol,auxtens,AgradV_VE,a%sisgs2(ielem)%a(:,2,e%igaus),ReferenceDtinv,elrhc)
      !third subscale of stress (deformation terms)
      call supm_elmrhcVES_dss_split3(e,auxVEnew,tisig,dvol,auxtens,grvel,a%sisgs3(ielem)%a(:,2,e%igaus),ReferenceDtinv,elrhc)
      
   end subroutine   
   
   
   subroutine DynamicSgsTestf_ASGS
      use Mod_supm_ComputeTestf
      use Mod_supm_ComputeTaus
      !Here we modify the test function so that it takes into account
      !tau_t*(-L*(v_h))-v_h[1-tau_t*tau_k^-1],-R)

      testf_static(1:e%pnode) = testf(1:e%pnode)
      testf(1:e%pnode) = testf(1:e%pnode) - e%shape(1:e%pnode,e%igaus)*(1-timom/timom_static)
   end subroutine  
   
end module

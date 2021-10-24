module Mod_nsf_DynSubsElmope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_DynSubsElmope
   implicit none

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   
   !NSF_Elmope1rst
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersDynSubsElmope1rst(itask)
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
               if (a%kfl_repro >= 1) then !OSS and Split OSS
                  !This adds the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_OSS_nsf_Elmope1rst)
               elseif (a%kfl_repro == 0) then
                  call runend('Dynsubs not ready for ASGS in fractional step')
               endif
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !Computation Subroutines
    subroutine InGaussElmatsDSS_OSS_nsf_Elmope1rst
      use Mod_nsm_ComputeTestf
      use typre
      implicit none
      integer(ip) :: inode
      
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      !Compute contributions to RHS : Block U
      call nsm_elmrhu_dss(e,dvol,a%dtinv,acden,testf,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhs)
            
      
   end subroutine
   
   !NSF_Elmope2nd
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersDynSubsElmope2nd(itask)
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
               if (a%kfl_repro == 1) then ! OSS 
                  !This adds the subscale at the previous time step to the RHS
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_OSS_nsf_Elmope2nd)
	       elseif (a%kfl_repro == 2 .or. a%kfl_repro == 3)	then ! Split OSS
	          call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsDSS_SOSS_nsf_Elmope2nd)
               elseif (a%kfl_repro == 0) then
                  call runend('Dynsubs not ready for ASGS in fractional step')
               endif
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !Computation Subroutines
    subroutine InGaussElmatsDSS_OSS_nsf_Elmope2nd
      use Mod_nsm_ComputeTestf
      use typre
      implicit none
      integer(ip) :: inode
      
      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h),ro tesgs_n/dtime)
      !Compute contributions to RHS : Block U
            
      !Compute contributions to RHS : Block P
      call nsm_elmrhp_dss(e,dvol,a%dtinv,acden,timom,a%vesgs(ielem)%a(1:e%ndime,2,e%igaus),elrhs)
      
   end subroutine
   
	subroutine InGaussElmatsDSS_SOSS_nsf_Elmope2nd !split OSS
      use Mod_nsm_ComputeTestf
      use typre
      implicit none
      integer(ip) :: inode
      
      !We add the contribution of the subscales at the previous time step 
            
      !Compute contributions to RHS : Block P (2nd problem)
      call nsm_elmrhp_dss(e,dvol,a%dtinv,acden,timom,a%vesgs2(ielem)%a(1:e%ndime,2,e%igaus),elrhs)
      
   end subroutine
   
end module

module Mod_nsc_ComputeSubscales
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_SubgridSpaceResidual
   use Mod_nsc_ComputeTaus
   implicit none
   private
   public SetPointersComputeSubscales,GPSGS, SetPointersGetSubscales

   !SubgridScales
   real(rp) :: GpSGS(5)
   
   integer(ip), allocatable :: kfl_IsSet, kfl_IsSetGetSubscales
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubscales(itask)
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
            
            !We need to compute the residual and project it 
            !to the subgrid scale space
            call SetPointersComputeSubgridSpaceResidual(1) 
            
            if (a%kfl_tacsg == 0) then
               
               !We need to compute the Tau values
               call SetPointersComputeTaus(1)
               
               call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,ComputeSubgridScaleQSS)
               
            elseif ( a%kfl_tacsg == 1) then
            
               !We need to modify the tau so that it includes 
               !the time step, but we keep the old timom in timom_static
               call SetPointersComputeTaus(1)
               
               !We compute the transient subscales
               call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,ComputeSubgridScaleDSS)
               
            endif
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersGetSubscales(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSetGetSubscales)
         call a%Memor%allocObj(0,'kfl_IsSetGetSubscales','InitProcedurePointer',1)
         kfl_IsSetGetSubscales = -1
      
      case(1)
      
         if (kfl_IsSetGetSubscales == -1) then
            kfl_IsSetGetSubscales = 1
            
            
            if (a%kfl_trasg == 0) then
               call SetPointersComputeSubscales(1)
               
            else 

               call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,GetSubgridScale)

            endif
         endif  
      case(100)
         deallocate(kfl_IsSetGetSubscales)
         call a%Memor%deallocObj(0,'kfl_IsSetGetSubscales','InitProcedurePointer',1)
      end select  
         
   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute 

   !Tracking of Subscales
   !Dynamic subscales
   subroutine ComputeSubgridScaleDSS
      implicit none

      !We add the contribution of the subscales at the previous time step 
      !tau_t*(-L*(v_h)+v_h*tau_k^-1*tau_t,Ao Usgs_n/dtime)
      
      !Uses a Backward Euler scheme
      a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = timom(2)*(gpmSGSpaceResidual(1:e%ndime) + a%mosgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      a%cosgs(ielem)%a(1,e%igaus) = timom(1)*(gpdSGSpaceResidual(1) + a%cosgs(ielem)%a(2,e%igaus)*ReferenceDtinv)
      a%ensgs(ielem)%a(1,e%igaus) = timom(3)*(gpeSGSpaceResidual(1) + a%ensgs(ielem)%a(2,e%igaus)*ReferenceDtinv)
      
      GpSGS(1:e%ndime) = a%mosgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(4) = a%cosgs(ielem)%a(1,e%igaus)
      GpSGS(5) = a%ensgs(ielem)%a(1,e%igaus)
   end subroutine
   
    subroutine GetSubgridScale
      implicit none
      
      GpSGS(1:e%ndime) = a%mosgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(4) = a%cosgs(ielem)%a(1,e%igaus)
      GpSGS(5) = a%ensgs(ielem)%a(1,e%igaus)
   end subroutine
      
   !Static subscales   
   subroutine ComputeSubgridScaleQSS
      implicit none

      GpSGS(1:e%ndime)   = timom(2)*gpmSGSpaceResidual(1:e%ndime)
      GpSGS(4)           = timom(1)*gpdSGSpaceResidual(1)
      GpSGS(5)           = timom(3)*gpeSGSpaceResidual(1)
      if (a%kfl_trasg /= 0) then
         a%cosgs(ielem)%a(1,e%igaus)           = GpSGS(4)          
         a%mosgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime) 
         a%ensgs(ielem)%a(1,e%igaus)           = GpSGS(5)  
      endif
   end subroutine
end module 

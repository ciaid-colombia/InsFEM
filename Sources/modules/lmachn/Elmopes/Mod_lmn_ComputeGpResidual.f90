  module Mod_lmn_ComputeGpResidual
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_InterpolateGradients
   implicit none
   private
   public SetPointersComputeGpResidual, gpres
   
   real(rp), allocatable :: gpres(:)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SetPointersComputeGpResidual(itask)
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
         
            !We need the gradients in the gauss point
            call SetPointersInterpolateGradients(1)
         
            call ConcatenateProcedures(ProcHook%Initializations,AllocGpRes)
            call ConcatenateProcedures(ProcHook%Finalizations,DeallocGpRes)
               
            if (a%kfl_nolsg == 0) then
               if (a%kfl_repro == 0) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_adv)
               elseif (a%kfl_repro == 1) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_OSS)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_adv)
                  call ConcatenateProcedures(ProcHook%ComputeResidual,ComputeGpRes_OSS)
                  call ConcatenateProcedures(ProcHook%ComputeResidual,ComputeGpRes_adv)
               end if
            elseif (a%kfl_nolsg == 1) then
               if (a%kfl_repro == 0) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_adv_nln)
               elseif (a%kfl_repro == 1) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_OSS_nln)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_adv_nln)
                  call ConcatenateProcedures(ProcHook%ComputeResidual,ComputeGpRes_OSS)
                  call ConcatenateProcedures(ProcHook%ComputeResidual,ComputeGpRes_adv)
               else
                  call runend('lmn_EndsteElmope: Dissipation option not implemented?')
               endif
            endif

            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpNonLinearRes)
 
            !Store residual for postprocess
            if (a%npp_stepi(11) /= 0) then
               call ConcatenateProcedures(ProcHook%ComputeResidual,StoreGpRes)
            endif
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
  
 
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------
   !Residual computation
   subroutine AllocGpRes
      implicit none
      
      call a%Memor%alloc(a%ndofn,gpres,'gpres','lmn_EndElmope')
   end subroutine
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(a%ndofn,gpres,'gpres','lmn_EndElmope')
   end subroutine
   
   subroutine ComputeGpRes_adv
      implicit none

      call lmn_elmrfe_mom_adv(e,acden,gpadv,grvel,gpres(1:e%ndime))
      call lmn_elmrfe_ene_adv(e,acden,gpadv,grtem,gpres(e%ndime+2))

   end subroutine

   subroutine ComputeGpRes_adv_nln
      implicit none

      call lmn_elmrfe_mom_adv(e,acden,gpvel(:,1),grvel,gpres(1:e%ndime))
      call lmn_elmrfe_ene_adv(e,acden,gpvel(:,1),grtem,gpres(e%ndime+2))

   end subroutine

   subroutine ComputeGpRes
      implicit none

      gpres = 0.0_rp 
      call lmn_elmrfe_pre(e,Integrator,a%dtinv,actex,gpadv,gpden,grtem,divvel,gpres(e%ndime+1))
      call lmn_elmrfe_mom(e,LHSdtinv,acden,gpvel(:,1),grpre,elext_mom,gpres(1:e%ndime))
      call lmn_elmrfe_ene(e,Integrator,LHSdtinv,a%dtinv,acden,actex,gtemp,gptem(1),acpth(1:nsteps),elext_ene,gpres(e%ndime+2))

   end subroutine

   subroutine ComputeGpRes_OSS
      implicit none

      gpres = 0.0_rp 
      call lmn_elmrfe_pre(e,Integrator,a%dtinv,actex,gpvel(:,1),gpden,grtem,divvel,gpres(e%ndime+1))
      call lmn_elmrfe_mom_trm(e,grpre,elext_mom,gpres(1:e%ndime))
      call lmn_elmrfe_ene_trm(e,Integrator,a%dtinv,actex,gtemp,acpth(1:nsteps),elext_ene,gpres(e%ndime+2))

   end subroutine
   
   subroutine ComputeGpRes_OSS_nln
      implicit none

      gpres = 0.0_rp 
      call lmn_elmrfe_pre(e,Integrator,a%dtinv,actex,gpadv,gpden,grtem,divvel,gpres(e%ndime+1))
      call lmn_elmrfe_mom_trm(e,grpre,elext_mom,gpres(1:e%ndime))
      call lmn_elmrfe_ene_trm(e,Integrator,a%dtinv,actex,gtemp,acpth(1:nsteps),elext_ene,gpres(e%ndime+2))

   end subroutine
   
   subroutine ComputeGpNonLinearRes
      implicit none
   
      call lmn_elmrfe_mom_nonlinear(e,acvis,elvel(:,:,1),gpres(1:e%ndime))
      call lmn_elmrfe_ene_nonlinear(e,actco,eltem(:,1),gpres(e%ndime+2))

   end subroutine
   
   subroutine StoreGpRes
      implicit none
      
      a%residualU(ielem)%a(:,e%igaus) = gpres(1:e%ndime)
      a%residualP(ielem)%a(e%igaus)   = gpres(e%ndime+1)
      a%residualT(ielem)%a(e%igaus)   = gpres(e%ndime+2)
   end subroutine
   
end module

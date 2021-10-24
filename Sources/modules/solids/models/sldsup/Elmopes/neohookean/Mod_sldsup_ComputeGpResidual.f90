  module Mod_sldsup_ComputeGpResidual
   use typre
   use Mod_sldsup_calculateAU
   use Mod_CauchyElement
   use Mod_sld_BaseElmope
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
            call concatenateprocedures(prochook%InGauss,calculateGradients)
            call ConcatenateProcedures(ProcHook%Initializations,AllocGpRes)
            call ConcatenateProcedures(ProcHook%Finalizations,DeallocGpRes)
               
            if (sup%kfl_nolsg == 0) then
               if (sup%kfl_repro == 0) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_ASGS)
               elseif (sup%kfl_repro == 1) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_OSS)
                  !call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpRes_ASGS)
               end if
            elseif (sup%kfl_nolsg == 1) then
                  call runend('sldsup_EndsteElmope: nonlinear sgs not implemented')
            endif

            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            !if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpNonLinearRes)
 
            !Store residual for postprocess
            if (sup%kfl_printResiduals) then
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
      
      call a%Memor%alloc(a%ndofn,gpres,'gpres','sldsup_EndElmope')
   end subroutine
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(a%ndofn,gpres,'gpres','sldsup_EndElmope')
   end subroutine

   subroutine ComputeGpRes_ASGS
      implicit none
      integer(ip)          :: nd,tn,u1,uf,s1,sf,p1,bc

      call sup%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      gpres = 0.0_rp 

      if(a%kfl_timei==1) call sup_residualU_dyn(e,nd,a%accel,gpres(u1:uf))
      call sup_residualU(e,nd,grdpre,divstr,elext,gpres(u1:uf))
      call sup_residualS(e,nd,tn,G,det,B,gpsigma,gpres(s1:sf))
      call sup_residualP(e,nd,G,lam,det,B,gppress,gpres(p1))

   end subroutine
   
   subroutine ComputeGpRes_OSS
      implicit none
      integer(ip)  :: nd,tn,u1,uf,s1,sf,p1,bc
      real(rp)     :: zero(e%ndime)

      call sup%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      gpres = 0.0_rp 
      zero  = 0.0_rp 

      call sup_residualU(e,nd,grdpre,divstr,zero,gpres(u1:uf))
      call sup_residualS_OSS(e,nd,tn,B,gpres(s1:sf))
      call sup_residualP_OSS(e,nd,G,lam,det,B,gpres(p1))

  end subroutine
   
   subroutine StoreGpRes
      implicit none
      integer(ip)          :: nd,tn,u1,uf,s1,sf,p1,bc

      call sup%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)
      
      sup%residualU(ielem)%a(:,e%igaus) = gpres(u1:uf)
      sup%residualS(ielem)%a(:,e%igaus) = gpres(s1:sf)
      sup%residualP(ielem)%a(e%igaus)   = gpres(p1)
   end subroutine
   
end module

module Mod_sldup_ComputeSubscales
   use typre
   use Mod_sld_BaseElmope
   use Mod_sldup_SubgridSpaceResidual
   implicit none
   private
   public SetPointersComputeSubscales

   !SubgridScales
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

            if (up%kfl_tacsg == 0) then

               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleQSS)

            elseif ( up%kfl_tacsg == 1) then

               call ConcatenateProcedures(ProcHook%Initializations,InitSGSTimeIntegrator)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleDSS) 

            endif
         endif

      case(100)

         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine   

   !------------------------------------------------------------------
   !Tracking of Subscales
   !Dynamic subscales
   subroutine InitSGSTimeIntegrator

       call php_SetTimeIntegrator(a,IntegratorSGS,LHSdtinv2SGS,nsteps)
       if (a%kfl_tsche_2nd_current == 'NEWMA') then
           call php_SetNewmarkCoefficients(a,IntegratorSGS,a%beta,a%gamma)
       endif
       call php_GetTimeSchemeDerivative(a,IntegratorSGS,tsch_deriv)

   end subroutine

   !Computes the transient stabilization parameter
   subroutine ComputeSubgridScaleDSS
      implicit none
      integer(ip) :: nd,tn,u1,uf,s1,sf,p1,bc,aux(e%ndime)
      real(rp)    :: tau_u,tau_s,tau_p
      real(rp)    :: gprhs(e%ndime)

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      call IntegratorSGS%GetRHS(e%ndime,up%u_sgs(ielem)%a(:,:,e%igaus),gprhs)
      aux = gprhs*a%densi*a%dtinv2

      !-------Tau_u------------
      up%u_sgs(ielem)%a(:,1,e%igaus) = tau_u*(gpSubscaleSpaceResidual(u1:uf) + aux)
       
      !-------Tau_p------------
      up%p_sgs(ielem)%a(e%igaus) = tau_p*gpSubscaleSpaceResidual(p1)

   end subroutine

   !Static subscales
   subroutine ComputeSubgridScaleQSS
      implicit none
      integer(ip) :: nd,tn,u1,uf,s1,sf,p1,bc,aux(e%ndime)
      real(rp)    :: tau_u,tau_s,tau_p

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      !-------Tau_u------------
      up%u_sgs(ielem)%a(:,1,e%igaus) = tau_u*gpSubscaleSpaceResidual(u1:uf)

      !-------Tau_p------------
      up%p_sgs(ielem)%a(e%igaus)     = tau_p*gpSubscaleSpaceResidual(p1)

   end subroutine

end module 

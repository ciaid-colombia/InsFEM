 module Mod_tem_Dissipation
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_GaussPointSGS
   use Mod_tem_ComputeGpSubscaleSpaceResidual
   use Mod_tem_GaussPointSGS
   use Mod_tem_TempeGradient
   implicit none
   private
   public SetPointersComputeDissipation
   
   integer(ip), allocatable :: kfl_IsSet
   logical  :: logical_aux
   
   !For Dissipation
   real(rp), allocatable :: eldis(:)
   real(rp)              :: gpdisTest,gpdis
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersComputeDissipation(itask)
      implicit none
      integer(ip) :: itask

      logical_aux = .false.
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            if (a%kfl_dispa > 1 ) logical_aux = .true.
            if (a%npp_stepi(4) /= 0) then   
               if (mod(a%istep,a%npp_stepi(4))==0) then
                  logical_aux = .true.
               endif
            endif

            !Dissipation
            if (logical_aux)  then
               call ConcatenateProcedures(ProcHook%PreLoop,PreLoopDis)
               call ConcatenateProcedures(ProcHook%Initializations,InitializationsDis)
               call ConcatenateProcedures(ProcHook%PreGauss,ElmatsToZeroDis)
               call SetPointersComputeTempeGradient(1) 
               
               call ConcatenateProcedures(ProcHook%Interpolates,GpDisToZero)
               
               !Smagorinsky or wale
               if (a%kfl_cotur < 0) then
               
               !Numerical Dissipation
               else
                  if (a%kfl_repro == 0) then
                     !-L*(u_h)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpDisTestASGS)
                  else 
                     call SetPointersComputeGpSubscaleSpaceResidual(1)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpDisTestOSS)
                  endif
                  
                  !Sets all the pointers so that the subgrid scale at the gauss point is computed
                  call SetPointersComputeGpSGS(1)
                     
                  !Once the subscales are computed, I can add the corresponding numerical dissipation
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpDisNumerical)
                  
                  !Transient Subscales
                  if (a%kfl_tacsg == 1) then
                     !If subscales are not orthogonal, I need to add the term
                     !(u_h,\partial \tilde u)
                     if (a%kfl_repro == 0) then
                        !In order to do this, I add (u_h,\partial \tilde u) to the computation of the dissipation
                        call ConcatenateProcedures(ProcHook%Interpolates,AdditionalDissipation_DSS_ASGS)
                     endif
                  endif
                  
               endif   
               
               !Adds the molecular dissipation (including turbulent if Smagorinsky has been used)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGpDis)
               
               call ConcatenateProcedures(ProcHook%InGaussElmats,GpdisAssemblyToElDis)
               
               call ConcatenateProcedures(ProcHook%Assembly,AssemblyDissipation)
               call ConcatenateProcedures(ProcHook%Finalizations,FinalizationsDis)
               call ConcatenateProcedures(ProcHook%PostLoop,ProjectDissipation)
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !FOR DISSIPATION
   subroutine PreLoopDis
      implicit none
      
      a%dissipation = 0.0_rp
   end subroutine
   
   subroutine InitializationsDis
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%mnode,eldis,'eldis','tem_EndsteElmope')
      gpdis = 0.0_rp
   end subroutine

   subroutine FinalizationsDis
      implicit none
      
      !Matrices dealloc
      call a%Memor%dealloc(e%mnode,eldis,'eldis','tem_EndsteElmope')
   end subroutine

   subroutine ElmatsToZeroDis
      implicit none
      
      eldis = 0.0_rp
   end subroutine
      
   subroutine ComputeGpDisTestASGS
      implicit none
      
      integer(ip) :: idime,jdime
      
      GpDisTest = -dot_product(testf(1:e%pnode),eltem(1:e%pnode,1))/timom
   end subroutine
   
   subroutine ComputeGpDisTestOSS
      implicit none
      
      integer(ip) :: idime,jdime
      
      !Here we are cheating and we reuse the projection of the residual as gptest functions
      
      GpDisTest = -gpSubscaleSpaceResidual
   end subroutine
   
   subroutine GpDisToZero
      implicit none
      
      gpdis = 0.0_rp
   end subroutine
   
   subroutine ComputeGpDisNumerical
      implicit none
      
      !Numerical Dissipation
      gpdis = gpdis + gpdisTest*GpTempeSgs(1)
   end subroutine   
   
   subroutine ComputeGpDis
      implicit none
      
      !Molecular and turbulent (acvis smago) dissipation
      gpdis = gpdis + acvis*dot_product(grtem(1:e%ndime),grtem(1:e%ndime))
   end subroutine
   
   subroutine GpdisAssemblyToElDis
      implicit none
      
      eldis(1:e%pnode) = eldis(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*gpdis*dvol
   end subroutine
   
   subroutine AssemblyDissipation
      implicit none
      
      !a%dissipation(e%lnods(1:e%pnode)) = a%dissipation(e%lnods(1:e%pnode)) + eldis(1:e%pnode)
      call a%Mesh%AssemblyToArray(e,1_ip,eldis,a%dissipation)
   end subroutine

   subroutine ProjectDissipation
      implicit none
      call a%Project(1,a%Dissipation) 
   end subroutine
   
   !Dissipation due to transient subscales if ASGS, (u_h,\partial \tilde u)
   subroutine AdditionalDissipation_DSS_ASGS
      implicit none
      
      !First order subscales
      gpdis = gpdis + gptem(1)*(a%tesgs(ielem)%a(1,e%igaus)-a%tesgs(ielem)%a(2,e%igaus))*ReferenceDtinv
   end subroutine
  
   
end module

 

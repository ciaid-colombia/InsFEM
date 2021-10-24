module Mod_tem_EndsteElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_Temperature
   use Mod_Element
   use Mod_TemperatureElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_tem_BaseElmope
   use Mod_tem_ComputeTaus
   use Mod_tem_ComputeTestf
   use Mod_tem_NonLinearDerivatives
   use Mod_tem_VelocityGradient
   use Mod_tem_Turbulence
   use Mod_tem_Advection
   use Mod_tem_TempeGradient
   use Mod_tem_ComputeGpResidual
   use Mod_tem_ComputeGpSubscaleSpaceResidual
   use Mod_tem_InterpolateResidualProjection
   use Mod_tem_GaussPointSGS
   use Mod_tem_Dissipation
   use Mod_tem_DustTransport
   use Mod_tem_ComputeGradientProjection
   implicit none
   
contains

   !SetPointers
   subroutine SetPointers
      use typre
      implicit none
      
      integer(ip) :: kfl_nonlinear,nelty
      
      !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
      
      call SetPointersAndHooksToNULLSUB
      
      call SetPointersComputeTaus(0)
      call SetPointersComputeTestf(0)
      call SetPointersNonLinearDerivatives(0)
      call SetPointersVelocityGradient(0)
      call SetPointersTurbulence(0)
      call SetPointersAdvectionVelocity(0)
      call SetPointersComputeTempeGradient(0)
      call SetPointersComputeGpResidual(0)
      call SetPointersComputeGpSubscaleSpaceResidual(0)
      call SetPointersComputeGpSGS(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeDissipation(0)
      call SetPointersDustTransport(0)
      call SetPointersComputeGradientProjection(0) 

      
      
      call SetPointersComputeTaus(1)
      call SetPointersComputeTestf(1)
      call SetPointersNonLinearDerivatives(1)
      call SetPointersAdvectionVelocity(1)
      
      !-------------------------------------------------------------------
      !Tracking of the subgrid scales
      !or transient subgrid scales
      if (a%kfl_tacsg == 1 .or. a%kfl_trasg == 1) then
         call SetPointersComputeGpSGS(1)
      endif
      
      !-----------------------------------------------------------
      !For adaptivity error indicators and subscales on the element boundaries
      if (a%kfl_adapsgs == 1) call SetPointersComputeGradientProjection(1)

      !-------------------------------------------------------------------
      call SetPointersComputeDissipation(1)
      
      call SetPointersTurbulence(1)
      
      call SetPointersDustTransport(1)
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      
      call SetPointersComputeTaus(100)
      call SetPointersComputeTestf(100)
      call SetPointersNonLinearDerivatives(100)
      call SetPointersVelocityGradient(100)
      call SetPointersTurbulence(100)
      call SetPointersAdvectionVelocity(100)
      call SetPointersComputeTempeGradient(100)
      call SetPointersComputeGpResidual(100)
      call SetPointersComputeGpSubscaleSpaceResidual(100)
      call SetPointersComputeDissipation(100)
      call SetPointersComputeGpSGS(100)
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersDustTransport(100)
      call SetPointersComputeGradientProjection(100) 
   end subroutine
 
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module
 
 
subroutine tem_EndsteElmope(TempeProblem) 
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_Temperature
   use def_parame
   use Mod_Tem_EndsteElmope
   implicit none
   class(TemperatureProblem), target :: TempeProblem
   a => TempeProblem
   
   
   !Return if there is nothing to be done
   if (a%kfl_dispa == 0 .and. a%kfl_trasg == 0 .and. a%kfl_tacsg == 0) return
   
   if (a%kfl_sourc==3) then
      call runend('tem_EndsteElmope not ready for user defined source terms')
   endif
   
   if (a%kfl_advec == 9) then
      call runend('tem_EndsteElmope not ready for Burgers equation')
   endif
   
   !If only dissipation
   if (a%kfl_trasg == 0 .and. a%kfl_tacsg == 0) then
      if (a%kfl_dispa > 1) then
         !Compute the dissipation always, as this is needed for something else (optics module)
      elseif (a%istep>=a%npp_inits .and.a%npp_stepi(4)>0 ) then
         if(mod(a%istep,a%npp_stepi(4))/=0) then
            return
         endif
      endif
   endif
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call ProcHook%PreLoop
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   ReferenceDtinv = a%dtinv
   if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = a%dtinv*2
   
   call a%Mesh%GetNelem(nelem)
   
   !Element Allocation
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'tem_EndsteElmope')
      
   !Other Arrays Alloc
   call a%Memor%alloc(e%mnode,a%ncomp-1,eltem,'eltem','tem_elmope')
   call a%Memor%alloc(a%ncomp-1,gptem,'gptem','tem_elmope')

   !Advection
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','tem_elmope')
   gpvno = 0.0_rp
   
   !Physical Parameters
   call a%GetPhysicalParameters(1,acden,acsph,actco,acrea,acsou)
   acvis = actco/acsph
   acrcp = acrea/acsph
   
   !Hook
   call ProcHook%Initializations
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      if (a%NumberofMaterials > 1) then
         !Physical Parameters
         call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
         acvis = actco/acsph
         acrcp = acrea/acsph
      endif
   
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
   
      !Elmats to Zero
      call ProcHook%ElmatsToZero
     
      !Gathers
      call e%gather(1,eltem(:,1),a%tempe(:,1))
      call e%gather(1,eltem(:,2),a%tempe(:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(1,eltem(:,itime),a%tempe(:,itime+1)) 
      enddo
      
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call e%interpg(1,eltem(:,1),gptem(1))
         call e%interpg(1,eltem(:,2),gptem(2))
         do itime = 3,nsteps ! Time bdf2 and others
            call e%interpg(1,eltem(:,itime),gptem(itime))
         enddo
         !Hook
         call ProcHook%Interpolates
         
         !Pointer 
         call ProcPointer%VelocityAndAgradV
         
         !PhysicalProperties
         acvis = actco/acsph
         acrcp = acrea/acsph
         !Hook
         call ProcHook%PhysicalProp
         
         !Compute the stability parameters
         call ProcPointer%ComputeTau
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         !Hook
         call ProcHook%ComputeTestf
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Compute vector of external forces (source term)
         call tem_ComputeExternalForces(e,acsou,acsph,elext)
         !Time integration
         call tem_TimeIntegrationToElext(e,Integrator,acden,a%dtinv,gptem,elext)
         !Hook
         call ProcHook%Elext
         
         !InGaussElmats
         !Hook
         call ProcHook%InGaussElmats
   
      enddo gauss_points
   
      !Assembly Endste
      !Hook
      call ProcHook%Assembly
   
   enddo elements
   
   !Hook
   call ProcHook%Finalizations
   
   !Other Arrays Alloc
   call a%Memor%dealloc(e%mnode,a%ncomp-1,eltem,'eltem','tem_elmope')
   call a%Memor%dealloc(a%ncomp-1,gptem,'gptem','tem_elmope')
   !Advection
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
   
   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook%PostLoop
   
   !Element Deallocation
   call a%Mesh%ElementDeAlloc(e,a%Memor,a%EndLoopQuadrature,'tem_EndsteElmope')
   
   
end subroutine

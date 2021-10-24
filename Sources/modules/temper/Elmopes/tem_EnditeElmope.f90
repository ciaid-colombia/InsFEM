module Mod_tem_EnditeElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_NavierStokes
   use Mod_Element
   use Mod_ConvectiveElement
   use Mod_TemperatureElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_PhysicalProblem
   use Mod_Temperature
   use Mod_tem_BaseElmope
   use Mod_tem_NonLinearDerivatives
   use Mod_tem_Advection
   use Mod_tem_Turbulence
   use Mod_tem_ALE
   use Mod_tem_TempeGradient
   use Mod_tem_ComputeGpResidual
   use Mod_tem_ComputeResidualProjection
   use Mod_tem_TempeGradientProjection
   use Mod_tem_DustTransport
   use Mod_tem_VelocityGradient
   use Mod_tem_SigmaupCoupling
   implicit none
   
   
 
contains

   !SetPointers
   subroutine SetPointers
      use typre
      implicit none
      
      integer(ip) :: nelty
      

      
      !External Procedures
      procedure() :: NULLSUB

      
      call ResetProcedureComposition
      
      call SetPointersAndHooksToNULLSUB
      
      call SetPointersNonLinearDerivatives(0)
      call SetPointersAdvectionVelocity(0)
      call SetPointersTurbulence(0)
      call SetPointersALE(0)
      call SetPointersComputeTempeGradient(0)
      call SetPointersComputeGpResidual(0)
      call SetPointersComputeResidualProjection(0)
      call SetPointersComputeTempeGradientProjection(0)
      call SetPointersDustTransport(0)
      call SetPointersVelocityGradient(0)
      call SetPointersSigmaupCoupling(0,'Endite')
      
      
      
      call SetPointersNonLinearDerivatives(1)
      call SetPointersAdvectionVelocity(1)
      if (a%kfl_repro == 1 .or. a%kfl_shock == 1 .or. a%kfl_shock == 2) then
         call SetPointersComputeResidualProjection(1)
      endif   
      !-------------------------------------------------------------------
      !Shock Capturing
      if (a%kfl_shock == 3 .or. a%kfl_shock == 4) then
         call SetPointersComputeTempeGradientProjection(1)
      endif
      call SetPointersTurbulence(1)
      call SetPointersALE(1)
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      
      call SetPointersDustTransport(1)
      
      call SetPointersSigmaupCoupling(1,'Endite')      
      
      call SetPointersNonLinearDerivatives(100)
      call SetPointersAdvectionVelocity(100)
      call SetPointersTurbulence(100)
      call SetPointersALE(100)
      call SetPointersComputeTempeGradient(100)
      call SetPointersComputeGpResidual(100)
      call SetPointersComputeResidualProjection(100)
      call SetPointersComputeTempeGradientProjection(100)
      call SetPointersDustTransport(100)
      call SetPointersVelocityGradient(100)
      call SetPointersSigmaupCoupling(100,'Endite')
      
   end subroutine
 
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
  
   
 
end module
 
 
 
 
 
subroutine tem_enditeElmope(TempeProblem) 
   use Mod_Tem_EnditeElmope
   use typre
   use def_parame
   use Mod_PhysicalProblem
   use Mod_Temperature


   implicit none
   class(TemperatureProblem), target :: TempeProblem
   a => TempeProblem
   
   !Return if there is nothing to be done
   if (a%kfl_repro == 0 .and. a%kfl_shock == 0 .and. .true.) return

   if (a%kfl_sourc==3) then
      call runend('tem_EnditeElmope not ready for user defined source terms')
   endif
   
   if (a%kfl_advec == 9) then
      call runend('tem_EnditeElmope not ready for Burgers equation')
   endif
  
   !Element Allocation
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'tem_EnditeElmope')
  
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call ProcHook%PreLoop
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   call a%Mesh%GetNelem(nelem)
   
  
   
   !Other Arrays Alloc
   call a%Memor%alloc(e%mnode,a%ncomp-1,eltem,'eltem','tem_elmope')
   call a%Memor%alloc(a%ncomp-1,gptem,'gptem','tem_elmope')
   call a%Memor%alloc(e%mnode,testf,'testf','tem_elmope')
   !Advection
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','tem_elmope')
   gpvno = 0.0_rp
   
   !Physical Parameters
   call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
   acvis = actco/acsph
   acrcp = acrea/acsph
   
   !Hook
   call ProcHook%Initializations
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      if (a%NumberofMaterials > 1) then
         !Physical Parameters
         call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
         acvis = actco/acsph
         acrcp = acrea/acsph
      endif
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
   
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
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Compute vector of external forces (source term)
         call ProcPointer%ExternalForces
         !Time integration
         call tem_TimeIntegrationToElext(e,Integrator,acden,a%dtinv,gptem,elext)
         !Hook
         call ProcHook%Elext
         
         !InGaussElmats
         !Hook
         call ProcHook%InGaussElmats
   
      enddo gauss_points
   
      !Assembly Endite
      !Hook
      call ProcHook%Assembly
   
   enddo elements
   
   !Hook
   call ProcHook%Finalizations
   
   !Other Arrays Alloc
   call a%Memor%dealloc(e%mnode,a%ncomp-1,eltem,'eltem','tem_elmope')
   call a%Memor%dealloc(a%ncomp-1,gptem,'gptem','tem_elmope')
   call a%Memor%dealloc(e%mnode,testf,'testf','tem_elmope')
   !Advection
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
   
   
   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook%PostLoop
   
   !Element Deallocation
   call a%Mesh%ElementDeAlloc(e,a%Memor,a%EndLoopQuadrature,'nsm_EnditeElmope')
   
   
end subroutine

module Mod_tem_elmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_Temperature
   use Mod_Element
   use Mod_TemperatureElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_php_Elmdir
   use Mod_tem_BaseElmope
   use Mod_tem_Advection
   use Mod_tem_VelocityGradient
   use Mod_tem_ExternalForces
   use Mod_tem_BaseElmope
   use Mod_tem_NonLinearDerivatives
   use Mod_tem_ComputeTestf
   use Mod_tem_InterpolateResidualProjection
   use Mod_tem_DynamicSubgridScales
   use Mod_tem_ComputeTaus
   use Mod_tem_Turbulence
   use Mod_tem_ALE
   use Mod_tem_HangingNodes
   use Mod_tem_ShockCapturing
   use Mod_tem_TempeGradient
   use Mod_tem_DustTransport
   use Mod_tem_SigmaupCoupling
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
      
      call SetPointersAndHooksTONULLSUB
      
      !Defaults for this subroutine
      ProcPointer%PostGaussElmats => PostGaussElmats
      ProcPointer%tem_elmbuv => tem_elmbuv
      ProcPointer%tem_elmrhu => tem_elmrhu

      
      
      !InitializePointers
      call SetPointersAdvectionVelocity(0)
      call SetPointersVelocityGradient(0)
      call SetPointersExternalForces(0)
      call SetPointersNonLinearDerivatives(0)
      call SetPointersComputeTestf(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersDynamicSubgridScales(0)
      call SetPointersComputeTaus(0)
      call SetPointersTurbulence(0)
      call SetPointersALE(0)
      call SetPointersHangingNodes(0)
      call SetPointersShockCapturing(0)
      call SetPointersComputeTempeGradient(0)
      call SetPointersDustTransport(0)
      call SetPointersSigmaupCoupling(0,'Elmope')
      
      !Set the pointers
      call SetPointersAdvectionVelocity(1)
      call SetPointersExternalForces(1)
      call SetPointersComputeTaus(1)
     
      !-------------------------------------------------------------------
      !Non-Linear Elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersNonLinearDerivatives(1)
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
         call ConcatenateProcedures(ProcHook%InGaussElmats,ProcPointer%PostGaussElmats)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsNonLinear)
         ProcPointer%PostGaussElmats => NULLSUB
      endif
            
      !-----------------------------------------------------------
      !ResidualProjection
      if (a%kfl_repro == 1) then
         call SetPointersInterpolateResidualProjection(1)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsRep)
         ProcPointer%tem_elmbuv => tem_elmbuv_rep
         ProcPointer%tem_elmrhu => tem_elmrhu_rep
      endif
      
      call SetPointersComputeTestf(1)
      
      if (a%kfl_tacsg /= 0) then
         call SetPointersDynamicSubgridScales(1)
      endif   
                 
      !-------------------------------------------------------------------
      !Shock Capturing
      call SetPointersShockCapturing(1)
      
      call SetPointersTurbulence(1)

      !------------------------------------------------------
      !ALE
      call SetPointersALE(1)
      
      !-------------------------------------------------------
      !Hanging Nodes
      call SetPointersHangingNodes(1)
      
      !DustTransport
      call SetPointersDustTransport(1)
      call SetPointersSigmaupCoupling(1,'Elmope')
      
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      
      !Finalize SetPointers
      call SetPointersAdvectionVelocity(100)
      call SetPointersVelocityGradient(100)
      call SetPointersExternalForces(100)
      call SetPointersNonLinearDerivatives(100)
      call SetPointersComputeTestf(100)
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersDynamicSubgridScales(100)
      call SetPointersComputeTaus(100)
      call SetPointersTurbulence(100)
      call SetPointersALE(100)
      call SetPointersHangingNodes(100)
      call SetPointersShockCapturing(100)
      call SetPointersComputeTempeGradient(100)
      call SetPointersDustTransport(100)
      call SetPointersSigmaupCoupling(100,'Elmope')
      
   end subroutine
  

   !--------------------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none
      
      !BLOCK U,P
      ! Viscosity terms : we only consider mu*(grad v, grad u)
      call elmvis(e,dvolt0,acvis,elmat)
   end subroutine
      
   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussNonLinear
      implicit none
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
   end subroutine
   
   subroutine InGaussElmatsNonLinear
      implicit none

      call tem_elmbuv_lap(e,dvol,actco,testf,elmat)
   end subroutine
   
   !Residual Projection
   subroutine InGaussElmatsRep
      implicit none
      
      !Compute contributions to RHS 
      call tem_elmrhs_oss(e,dvol,testf,gprep(1),elrhs)
   end subroutine
   
   !Element Change
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
  
end module


subroutine tem_elmope(TempeProblem)
   use typre
   use Mod_Temperature
   use Mod_Tem_elmope
   use Mod_TemperatureElement
   implicit none
   class(TemperatureProblem), target :: TempeProblem
   integer(ip) :: checkpoints(20)
   
   a=>TempeProblem

   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','tem_elmope')
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   ReferenceDtinv = a%dtinv
   if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = a%dtinv*2
      
   !Matrices Alloc
   call a%Memor%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','tem_elmope')
   call a%Memor%alloc(1_ip,e%mnode,elrhs,'elrhs','tem_elmope')
   
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

   !do itest = 1,100
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
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      
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
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      
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
         !Hook
         call ProcHook%ComputeTestf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Compute vector of external forces (source term)
         call ProcPointer%ExternalForces
         !Time integration
         call tem_TimeIntegrationToElext(e,Integrator,acden,a%dtinv,gptem,elext)
         !Hook
         call ProcHook%Elext
         
         !InGaussElmats
         !Compute contributions to RHS : Block U
         call ProcPointer%tem_elmrhu(e,dvol,testf,elext,elrhs)
               
         !Compute contributions to elemental matrix : Block U,V
         call ProcPointer%tem_elmbuv(e,dvol,acden,acrcp,LHSdtinv,AGradV,testf,elmat(1,:,1,:))

         !Hook
         call ProcHook%InGaussElmats
    
       enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer%PostGaussElmats
      
      call ProcHook%PreDirichlet
      
      !Dirichlet Boundary Conditions
      call php_elmdir(a,e,1_ip,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
     
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   enddo elements
   
   !Finalizations
   !Hook
   call ProcHook%Finalizations
   
   !Matrices Alloc
   call a%Memor%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','tem_elmope')
   call a%Memor%dealloc(1_ip,e%mnode,elrhs,'elrhs','tem_elmope')
   
   !Other Arrays Alloc
   call a%Memor%dealloc(e%mnode,a%ncomp-1,eltem,'eltem','tem_elmope')
   call a%Memor%dealloc(a%ncomp-1,gptem,'gptem','tem_elmope')
   !Advection
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','tem_elmope')

end subroutine

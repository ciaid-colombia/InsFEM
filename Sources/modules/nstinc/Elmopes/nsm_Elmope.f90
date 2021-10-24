module Mod_nsm_Elmope
   use Mod_Element
   use Mod_ConvectiveElement
   use Mod_NavierStokesElement
   use Mod_NavierStokes
   use Mod_php_SetTimeIntegrator
   use Mod_nsm_elmdir
   use Mod_NsiExacso    
   
   !For setting pointers
   use Mod_nsm_TemperatureCoupling
   use Mod_nsm_ComputeAdvectionVelocity
   use Mod_nsm_NonLinearDerivatives
   use Mod_nsm_PhysicalProperties
   use Mod_nsm_ExternalForces
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_TurbulenceModel
   use Mod_nsm_ComputeGpResidual
   use Mod_nsm_ComputeResidualProjection
   use Mod_nsm_InterpolateResidualProjection
   use Mod_nsm_SubgridSpaceResidual
   use Mod_nsm_ComputeTauSmoothing
   use Mod_nsm_ComputeTaus
   use Mod_nsm_ComputeTestf
   use Mod_nsm_ComputeSubscales
   use Mod_nsm_ComputeNonLinearSubscales
   use Mod_nsm_DynSubsElmope
   use Mod_nsm_NewtonRaphson
   use Mod_nsm_SwitchOff
   use Mod_nsm_LevelSetCoupling
   use Mod_nsm_FreeSurfaceMatrices   
   use Mod_nsm_EnrichElement
   use Mod_nsm_HangingNodes
   use Mod_nsm_FreeSurfaceStabilization
   use Mod_nsm_OrthogonalSubscales
   use Mod_nsm_BoundarySubscales
   use Mod_nsm_ElasticBoundaryDir
   use Mod_nsm_ExtraViscosity
   use Mod_nsm_PorousMedia
   use Mod_nsm_PressurePenalty
   use Mod_nsm_FORAxesRotation
   use Mod_nsm_PressTempSubscale

   implicit none   
   
contains

   !----------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      use typre
      implicit none

      !External Procedures
      procedure() :: NULLSUB

      integer(ip) :: kfl_nonlinear,nelty
            
      call ResetProcedureComposition
      
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Set the Specific Pointers for this subroutine
      ProcPointer_PreAssembly     => NULLSUB
      ProcPointer_PostGaussElmats => PostGaussElmats
      ProcPointer_nsm_elmbuv      => nsm_elmbuv
      ProcPointer_nsm_elmrhu      => nsm_elmrhu
      ProcPointer_nsm_elmrhp      => nsm_elmrhp
      ProcPointer_nsm_elmbuq      => nsm_elmbuq
      ProcPointer_nsm_elmbpv      => nsm_elmbpv
      
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersAdvectionVelocity%Initialize 
      call SetPointersExternalForces%Initialize 
      call SetPointersInterpolateGradients%Initialize 
      call SetPointersComputeTaus%Initialize 
      call SetPointersComputeTestf%Initialize 
      call SetPointersComputeGpResidual%Initialize 
      call SetPointersInterpolateResidualProjection%Initialize
      call SetPointersComputeSubgridSpaceResidual%Initialize 
      call SetPointersNonLinearDerivatives%Initialize 
      call SetPointersPhysicalProperties%Initialize 
      call SetPointersTemperatureCoupling%Initialize 
      call SetPointersTurbulenceModel%Initialize 
      call SetPointersDynSubsElmope%Initialize
      call SetPointersNewtonRaphson%Initialize
      call SetPointersLevelSetCoupling%Initialize
      call SetPointersFreeSurfaceMatrices%Initialize 
      call SetPointersEnrichElement%Initialize
      call SetPointersHangingNodes%Initialize
      call SetPointersSwitchOff%Initialize
      call SetPointersFSS%Initialize
      call SetPointersOrthogonalSubscales%Initialize
      !call SetPointersBoundarySubscales%Initialize
      call SetPointersElasticBoundaryDir%Initialize
      call SetPointersPorousMedia%Initialize
      call SetPointersFORAxesRotation%Initialize
      call SetPointersPressTempSubscale%Initialize
      call SetPointersExtraViscosity%Initialize
      call SetPointersPressurePenalty%Initialize

      !-------------------------------------------------------------------
      !Now we set the required pointers 
      
      call SetPointersAdvectionVelocity%Set
      
      !PenaltyTerm for the pressure
      call SetPointersPressurePenalty%Set
            
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersNonLinearDerivatives%Set
         if (a%kfl_repro .le. 1) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsNonLinear)
      end if
      
      if(kfl_nonlinear==1 .or. a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook_InGauss,InGaussVolumesNonLinear)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer_PostGaussElmats)
         ProcPointer_PostGaussElmats => NULLSUB      
      end if
      
      call SetPointersComputeTaus%Set
      call SetPointersOrthogonalSubscales%Set
      !call SetPointersBoundarySubscales%Set
      call SetPointersTurbulenceModel%Set
      call SetPointersExternalForces%Set
      call SetPointersTemperatureCoupling%Set
      call SetPointersComputeTaus%Set
      call SetPointersComputeTestf%Set
      call SetPointersHangingNodes%Set
      call SetPointersPorousMedia%Set
      call SetPointersFORAxesRotation%Set
      
      !Tracking of the subscales                      
      !Dynamic subscales
      if (a%kfl_tacsg == 1) call SetPointersDynSubsElmope%Set
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange

      !-----------------------------------------------------------------------
      if (a%kfl_linea==2) then
         call SetPointersInterpolateGradients%Set 
         call SetPointersNewtonRaphson%Set
      endif

      !-----------------------------------------------------------------------
      !Disconnection of elements
      if (a%kfl_SwitchOff == 1) call SetPointersSwitchOff%Set
      
      !We do it twice to ensure that postGauss will be free for 
      !Surface Tension
      if(kfl_nonlinear==1 .or. a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer_PostGaussElmats)
         ProcPointer_PostGaussElmats => NULLSUB      
      end if
      call SetPointersLevelSetCoupling%Set
      call SetPointersFreeSurfaceMatrices%Set 
      call SetPointersEnrichElement%Set 
      call SetPointersPhysicalProperties%SetTask('Elmope')
      call SetPointersPhysicalProperties%Set
      
      !Free surface stabilization
      call SetPointersFSS%Set
      call SetPointersElasticBoundaryDir%Set
      call SetPointersExtraViscosity%Set
      
      call SetPointersPressTempSubscale%Set
      
      !--------------------------------------------------------------------------
      !We deallocate the procedure flags, so that they can be set the next time
      call SetPointersAdvectionVelocity%Finalize     
      call SetPointersExternalForces%Finalize 
      call SetPointersInterpolateGradients%Finalize 
      call SetPointersComputeTaus%Finalize 
      call SetPointersComputeTestf%Finalize 
      call SetPointersComputeGpResidual%Finalize 
      call SetPointersInterpolateResidualProjection%Finalize
      call SetPointersComputeSubgridSpaceResidual%Finalize 
      call SetPointersNonLinearDerivatives%Finalize 
      call SetPointersPhysicalProperties%Finalize 
      call SetPointersTemperatureCoupling%Finalize 
      call SetPointersTurbulenceModel%Finalize 
      call SetPointersDynSubsElmope%Finalize
      call SetPointersNewtonRaphson%Finalize
      call SetPointersLevelSetCoupling%Finalize
      call SetPointersFreeSurfaceMatrices%Finalize 
      call SetPointersEnrichElement%Finalize 
      call SetPointersHangingNodes%Finalize
      call SetPointersSwitchOff%Finalize
      call SetPointersFSS%Finalize
      call SetPointersOrthogonalSubscales%Finalize
      !call SetPointersBoundarySubscales%Finalize
      call SetPointersElasticBoundaryDir%Finalize
      call SetPointersPorousMedia%Finalize
      call SetPointersFORAxesRotation%Finalize
      call SetPointersExtraViscosity%Finalize
      call SetPointersPressurePenalty%Finalize
      call SetPointersPressTempSubscale%Finalize
      
   end subroutine

   !----------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none
      
      !BLOCK U,P
      ! Viscosity terms : we only consider mu*(grad v, grad u)
      call elmvis(e,dvolt0,acvis,wrmat1)

      ! If you want the complete term div(e·grad u) uncomment next line 
      if ( a%fvins > zensi ) then
         call nsm_elmvis_div(e,dvolt0,acvis,elmuv)
      endif
      
      ! tau2*(div v, div u)
      call nsm_elmdiv(e,dvolt2,elmuv)
      
      !BLOCK P,Q : tau1*(graq q, grad p)
      call nsm_elmbpq(e,dvolt1,elmpq)
   end subroutine
   
   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
   end subroutine

   subroutine InGaussElmatsNonLinear
      implicit none
      call nsm_elmbuv_lap(e,dvol,acvis,testf,wrmat1)
      call nsm_elmbuq_lap(e,dvol,acvis,timom,elmuq)
   end subroutine

   subroutine OnIeltyChange
      implicit none
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module   


!NSM_ELMOPE subroutine   
subroutine nsm_Elmope(NSProblem)
   use Mod_nsm_Elmope
   implicit none
   class(NavierStokesProblem), target :: NSProblem
   
   a=>NSProblem
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsm_EnditeElmope')
   
   call ProcHook_PreAllocate
   
   !Allocate Arrays in BaseElmope
   call AllocateBaseElmopeArrays

   !Matrices Alloc
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_elmope')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','nsm_elmope')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsm_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','nsm_elmope')
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','nsm_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','nsm_elmope')
   call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','nsm_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,elrhu,'elrhu','nsm_elmope')
   call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','nsm_elmope')
   
   
   !Initialize Statistics
   call a%InitStats
   
   !Hook
   call ProcHook_Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  

      !Hook
      call ProcHook_OnIeltyChange
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      elmuv=0.0_rp
      elrhu=0.0_rp
      wrmat1=0.0_rp
      elmpq=0.0_rp
      elmpv=0.0_rp
      elmuq=0.0_rp
      elrhp=0.0_rp
      
      !Gathers
      call VelocityAndPressureGathers
      !Hook
      call ProcHook_Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp

      !Hook
      call ProcHook_PreGauss
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook_InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate         
         call InterpolateGpVelocities
         !Hook
         call ProcHook_Interpolates
         
         !PostHook 
         call ProcHook_PostInterpolate
         
         !Physical Properties
         call a%GetPhysicalParameters(imat,acden,acvis)
         call ProcHook_PhysicalProp
         
         !Advection velocity      
         call ProcPointer_ComputeAdvectionVelocity
         
         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
         
         !Compute a·grad(V)
         call ComputeAGradV(e,gpadv,AGradV)
         
         !Compute the stability parameters 
         call ProcHook_ComputeTaus
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call ProcHook_ComputeTestf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0   ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1   ! w(gp)*detjm*tau1
         dvolt2=dvol*tidiv + dvolt2   ! w(gp)*detjm*tau2
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Compute vector of external forces
         call ProcPointer_ExternalForces   
         
         !Time integration
         eltemp = 0.0_rp
         call nsi_TimeIntegrationToElTemp(e,Integrator,acden,a%dtinv,gpvel,eltemp)
         
         !InGaussElmats
         !Compute contributions to RHS : Block U
         call ProcPointer_nsm_elmrhu(e,dvol,testf,elext,eltemp,elrhu)
         
         !Compute contributions to elemental matrix : Block U,V
         call ProcPointer_nsm_elmbuv(e,dvol,acden,LHSdtinv,AGradV,testf,wrmat1)
         
         !Compute contributions to elemental matrix : Block V,P
         call ProcPointer_nsm_elmbpv(e,dvol,testf,elmpv)
            
         !Compute contributions to RHS : Block P       
         call ProcPointer_nsm_elmrhp(e,timom,dvol,elext,eltemp,elrhp)
         
         !Compute contributions to elemental matrix : Block U,Q
         call ProcPointer_nsm_elmbuq(e,timom,dvol,acden,LHSdtinv,AGradV,elmuq)     
         
         !Hook
         call ProcHook_InGaussElmats
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
      
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer_PostGaussElmats
      
      !Matrix composition
      ! Assembly wrmat1  in elmuv
      forall (idime = 1:e%ndime)
         elmuv(idime,1:e%pnode,idime,1:e%pnode) = elmuv(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      end forall   
   
      ! Assembly elmuv to elmat
      elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) + elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      
      ! Assembly elrhu to elrhs
      elrhs(1:e%ndime,1:e%pnode) = elrhu(1:e%ndime,1:e%pnode) + elrhs(1:e%ndime,1:e%pnode)
   
      ! Assembly elmuq to elmat
      elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) + elmuq(1,1:e%pnode,1:e%ndime,1:e%pnode)

      ! Assembly elmpv & elmpq to elmat
      elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) + elmpv(1:e%ndime,1:e%pnode,1,1:e%pnode)
      elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) + elmpq(1,1:e%pnode,1,1:e%pnode)

      ! Assembly  elrhp to elrhs
      elrhs(e%ndime+1,1:e%pnode) = elrhp(1,1:e%pnode) + elrhs(e%ndime+1,1:e%pnode)
      
      !PreDirichlet
      call ProcHook_PreDirichlet
      
      !Dirichlet Boundary Conditions
      call nsm_elmdir(a,e,elmat,elrhs)
            
      !Assembly
      call ProcHook_Assembly
      
   enddo elements
   !enddo
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook_Finalizations
   
   !Matrix Deallocations
   call a%Memor%dealloc(a%ndofn,size(elmat,2),a%ndofn,size(elmat,4),elmat,'elmat','nsm_elmope')
   call a%Memor%dealloc(a%ndofn,size(elrhs,2),elrhs,'elrhs','nsm_elmope')
   call a%Memor%dealloc(size(wrmat1,1),size(wrmat1,2),wrmat1,'wrmat1','nsm_elmope')
   call a%Memor%dealloc(e%ndime,size(elmuv,2),e%ndime,size(elmuv,4),elmuv,'elmuv','nsm_elmope')
   call a%Memor%dealloc(1,size(elmpq,2),1,size(elmpq,4),elmpq,'elmpq','nsm_elmope')
   call a%Memor%dealloc(e%ndime,size(elmpv,2),1,size(elmpv,4),elmpv,'elmpv','nsm_elmope')
   call a%Memor%dealloc(1,size(elmuq,2),e%ndime,size(elmuq,4),elmuq,'elmuq','nsm_elmope')
   call a%Memor%dealloc(e%ndime,size(elrhu,2),elrhu,'elrhu','nsm_elmope')
   call a%Memor%dealloc(1,size(elrhp,2),elrhp,'elrhp','nsm_elmope')  
   
   !Arrays Deallocations
   call DeallocateBaseElmopeArrays
   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsm_Elmope')
   
end subroutine

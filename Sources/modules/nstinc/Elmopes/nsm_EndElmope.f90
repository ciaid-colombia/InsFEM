module Mod_nsm_EndElmope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_TemperatureCoupling
   use Mod_nsm_LevelSetCoupling   
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
   use Mod_nsm_ComputeDissipation
   use Mod_nsm_ComputeVorticity
   use Mod_nsm_ComputeQfactor
   use Mod_nsm_SwitchOff
   use Mod_nsm_FreeSurfaceStabilization
   use Mod_nsm_ComputeDivergence
   use Mod_nsm_ComputeGradientProjection
   implicit none 
   
contains       
    
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      procedure() :: NULLSUB      
      integer(ip) :: nelty
      
      !Reset Procedure Composition
      call ResetProcedureComposition
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Set the Specific Pointers for this subroutine
      ProcPointer_PreAssembly     => NULLSUB      
      
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
      call SetPointersGetSubscales%Initialize 
      call SetPointersSwitchOff%Initialize 
      call SetPointersComputeDivergence%Initialize
      
      !Also for the specific EndElmope Procedures
      call SetPointersComputeDissipation%Initialize 
      call SetPointersNonLinearSubscales%Initialize
      call SetPointersComputeResidualProjection%Initialize 
      call SetPointersComputeGradientProjection%Initialize 
      call SetPointersComputeTauSmoothing%Initialize 
      call SetPointersComputeSubscales%Initialize 
      call SetPointersComputeVorticity%Initialize
      call SetPointersComputeQFactor%Initialize
      
      !LevelSet Coupling
      call SetPointersLevelSetCoupling%Initialize
      call SetPointersFSS_Endite%Initialize

      !----------------------------------------------------------------------------------------         
      !Now we set the required pointers         
      call SetPointersAdvectionVelocity%Set
      call SetPointersNonLinearDerivatives%Set
      call SetPointersExternalForces%Set
      call SetPointersTurbulenceModel%set
      call SetPointersTemperatureCoupling%Set
      call SetPointersLevelSetCoupling%Set
      call SetPointersFSS_Endite%Set
      call SetPointersPhysicalProperties%SetTask('EndElmope')
      call SetPointersPhysicalProperties%Set
      call SetPointersComputeDivergence%Set
      
      !-----------------------------------------------------------------------
      !Disconnection of elements
      if (a%kfl_SwitchOff == 1) call SetPointersSwitchOff%Set
     
      !-----------------------------------------------------------
      !Residual Projection
      if (itask .eq. 'Endite' .and. a%kfl_repro >= 1) call SetPointersComputeResidualProjection%Set
      
      !-----------------------------------------------------------
      !For adaptivity error indicators and subscales on the element boundaries
      if (itask .eq. 'Endste' .and. a%kfl_adapsgs == 1) call SetPointersComputeGradientProjection%Set
      
      !---------------------------------------------------------------
      !Dynamic, non-linear subscales
      !For Endste, only if tracking or transient, do not do it if non-linear
      if (itask .eq. 'Endste' .and. a%kfl_nolsg == 0 .and. (a%kfl_trasg == 1 .or. a%kfl_tacsg == 1)) then
         call SetPointersComputeSubscales%Set
      endif
      
      !For Endite, if non-linear subscales
      if (itask .eq. 'Endite' .and. a%kfl_nolsg == 1) call SetPointersNonLinearSubscales%Set
     
      !Tau smoothing
      if (itask .eq. 'Endite' .and. a%kfl_Tausm >= 1) call SetPointersComputeTauSmoothing%Set
         
      if (itask .eq. 'Endste') then
         !Vorticity
         call SetPointersComputeVorticity%Set
         !Q-factor
         call SetPointersComputeQFactor%Set
         !Dissipation computation
         call SetPointersComputeDissipation%Set
      endif

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
     
      !Deallocate the procedure flags, so that they can be used in the next run
      !Initialize the ProcedureFlags
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
      call SetPointersGetSubscales%Finalize 
      call SetPointersSwitchOff%Finalize
      
      !Also for the specific EndElmope Procedures
      call SetPointersComputeDissipation%Finalize
      call SetPointersNonLinearSubscales%Finalize
      call SetPointersComputeResidualProjection%Finalize 
      call SetPointersComputeGradientProjection%Finalize 
      call SetPointersComputeTauSmoothing%Finalize 
      call SetPointersComputeSubscales%Finalize
      call SetPointersComputeVorticity%Finalize
      call SetPointersComputeQFactor%Finalize
      call SetPointersComputeDivergence%Finalize

      !LevelSet
      call SetPointersLevelSetCoupling%Finalize 
      call SetPointersFSS_Endite%Finalize

   end subroutine
   
   !--------------------------------------------------------------------
   !Multiple type of elements
   subroutine OnIeltyChange
      implicit none
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module

subroutine nsm_EndElmope(NSProblem,task)
   use Mod_nsm_BaseElmope
   use Mod_nsm_EndElmope
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem), target :: NSProblem
   character(6) :: task
   integer(ip)  :: aux_logic
   logical      :: aux_logic2
   integer(ip)  :: itime
   
   a=>NSProblem
   itask = task
   
   !Things to be done if endite
   if (itask .eq. 'Endite') then
      !Return if there is nothing to be done
      aux_logic = 0 
      
      if (a%kfl_repro >= 1) aux_logic = 1
      if (a%kfl_nolsg == 1) aux_logic = 1
      if (a%kfl_Tausm == 1) aux_logic = 1
      if (a%kfl_StabilizeFreeSurface == 1) aux_logic = 1
      !Divergence postprocess
      if (a%npp_stepi(11) /= 0) then   
         if (mod(a%istep,a%npp_stepi(11))==0) then
            aux_logic = 1
         endif
      endif
      if (aux_logic == 0) return
   
   !Things to be done if endste
   elseif (itask .eq. 'Endste') then
      !Return if there is nothing to be done
      aux_logic = 0
      if (a%kfl_dispa /= 0) aux_logic = aux_logic + 1
      if (a%kfl_trasg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      if (a%kfl_tacsg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      !Logics for deciding if we need to compute vorticity
      aux_logic2 = .false.
      if (a%npp_stepi(10) /= 0) then   
        if (mod(a%istep,a%npp_stepi(10))==0) then
           aux_logic2 = .true.
        endif
      endif

      !Logics for deciding if we need to compute Q-factor
      if (a%npp_stepi(21) /= 0) then   
         if (mod(a%istep,a%npp_stepi(21))==0) then
            aux_logic2 = .true.
         endif
      endif      

      if (aux_logic2 .eqv. .true.) aux_logic = aux_logic + 1
     
      if (aux_logic == 0) return
      
      !If only dissipation for postprocessing and not a postprocessing step return
      if (a%kfl_dispa == 1 .and. aux_logic == 1) then
         if (a%istep < a%npp_inits) return
         if (a%npp_stepi(15)>0 ) then
            if (mod(a%istep,a%npp_stepi(15))/=0) then
               return
            endif
         endif
      endif
   endif
   
   !We force closed rule for smoothing
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'nsm_EnditeElmope')
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Hook
   call ProcHook_PreLoop
   
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_EnditeElmope')

   !AllocateArrays in BaseElmope
   call AllocateBaseElmopeArrays
   
   !Hook
   call ProcHook_Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      
      !Hook
      call ProcHook_OnIeltyChange
      
      !Elmats to Zero
      call ProcHook_ElmatsToZero
      
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
         
         !Physical Properties
         call a%GetPhysicalParameters(imat,acden,acvis)
         !Hook
         call ProcHook_PhysicalProp          
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Compute vector of external forces
         call ProcPointer_ExternalForces  
         
         !Time integration
         eltemp = 0.0_rp
         call nsi_TimeIntegrationToElTemp(e,Integrator,acden,a%dtinv,gpvel,eltemp)
         
         !Default is just one iteration
         kfl_GoIteInGauss = 1
         do while (kfl_GoIteInGauss > 0)
         
            !Advection velocity      
            call ProcPointer_ComputeAdvectionVelocity
         
            !Advection velocity norm
            call vecnor(gpadv,e%ndime,gpvno,2)
            
            !Compute a·grad(V)
            call ComputeAGradV(e,gpadv,AGradV)
            
            !Compute the stability parameters, only if required
            call ProcHook_ComputeTaus
            
            !Compute the stabilization test function, only if required
            call ProcHook_ComputeTestf
            
            !InGaussElmats
            !Hook
            call ProcHook_InGaussElmats
            
            kfl_GoIteInGauss = kfl_GoIteInGauss - 1
         enddo
         
         !InGaussElmats Assembly
         call ProcHook_InGaussElmatsAssembly
         
      enddo gauss_points      
  
      
      !PreDirichlet
      call ProcHook_PreDirichlet      
      
      !Assembly Endite
      !Hook
      call ProcHook_AssemblyEndite
      
   enddo elements
   !enddo
   
   !Hook
   call ProcHook_Finalizations
   
   !Deallocate arrays in BaseElmope
   call DeallocateBaseElmopeArrays

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_EnditeElmope')

   !Element
   call a%Mesh%ElementDeAlloc(e,a%Memor,a%EndLoopQuadrature,'nsm_EnditeElmope') 

   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook_PostLoop
   
end subroutine

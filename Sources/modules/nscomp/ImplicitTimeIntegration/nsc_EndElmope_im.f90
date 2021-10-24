module Mod_nsc_EndElmope_im
   use typre

   use Mod_nsc_BaseElmope
   use Mod_nsc_ComputeVariables
   use Mod_nsc_ComputeConvectionCoefficients
   use Mod_nsc_ComputeConvectionJacobianGradientCoefficients
   use Mod_nsc_ComputeDiffusionCoefficients
   use Mod_nsc_ComputeNonLinearDiffusionCoefficients
   use Mod_nsc_ComputeReactionCoefficients
   use Mod_nsc_ComputeTaus
   use Mod_nsc_ComputeTestf
   use Mod_nsc_NonLinearDerivatives
   use Mod_nsc_ExternalForces

   use Mod_nsc_ComputeGpResidual
   use Mod_nsc_ComputeResidualProjection
   use Mod_nsc_InterpolateResidualProjection
   use Mod_nsc_ComputeGradientProjection
   use Mod_nsc_ArbitraryLagrangianEulerian
   use Mod_nsc_SubgridSpaceResidual
   use Mod_nsc_ComputeSubscales
   use Mod_nsc_ComputeNonLinearSubscales
   implicit none
   
  
   
contains       
    
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      !External Procedures
      procedure() :: NULLSUB      
      
      integer(ip) :: kfl_nonlinear, nelty
      
      !Reset Procedure Composition
      call ResetProcedureComposition

      !Set All Pointers To NULLSUB (in Mod_nsc_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Set the Specific Pointers for this subroutine
      ProcPointer_nsc_PreAssembly     => NULLSUB      
      
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersVariables(0)     
      call SetPointersConvectionCoefficients(0)     
      call SetPointersConvectionJacobianGradientCoefficients(0)     
      call SetPointersDiffusionCoefficients(0)     
      call SetPointersNonLinearDiffusionCoefficients(0)     
      call SetPointersReactionCoefficients(0)
      call SetPointersExternalForces(0)        
     !call SetPointersInterpolateGradients(0)         
      call SetPointersComputeTaus(0)                  
      call SetPointersComputeTestf(0)                 

      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeSubgridSpaceResidual(0)  
      call SetPointersNonLinearDerivatives(0)         
      call SetPointersArbitraryLagrangianEulerian(0)
      call SetPointersGetSubscales(0)   
      
      !Also for the specific EndElmope Procedures
      call SetPointersComputeGpResidual(0)            
      call SetPointersComputeResidualProjection(0)     
      call SetPointersComputeGradientProjection(0)     
      call SetPointersComputeNonLinearSubscales(0)     
      !call SetPointersComputeTauSmoothing(0)           
      call SetPointersComputeSubscales(0)              
      
      !-----------------------------------------------------------------------
      !Now we set the required pointers         
          
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersNonLinearDerivatives(1)
      endif
      call SetPointersVariables(1)
      
      !Convection
      call SetPointersConvectionCoefficients(1)
          
       !Convection Jacobian Gradient Matrix Coefficients
      if (a%kfl_jacgr == 1) then
         call SetPointersConvectionJacobianGradientCoefficients(1)
      endif

      !Diffusion Exists
      if (a%kfl_visco == 1) then
          call SetPointersDiffusionCoefficients(1)
          !Non-linear elements
          if (kfl_nonlinear == 1) then
             call SetPointersNonLinearDiffusionCoefficients(1)
         endif
      endif

      !Reaction (transformed sources) Exists
      if (a%kfl_react == 1) then
         call SetPointersReactionCoefficients(1)
      endif

      !External forces 
      call SetPointersExternalForces(1)

      !-----------------------------------------------------------
      !Shock Capturing
      if (a%kfl_shock /= 0 ) then
        if (a%kfl_shock == 1) then
             call SetPointersComputeGpResidual(1)            
         else if (a%kfl_shock == 2) then
             call SetPointersComputeGradientProjection(1)            
         else if (a%kfl_shock > 2) then
             call runend('Nsc_elmope: Other shock capturing methods not ready')
        endif
      endif
      !-----------------------------------------------------------
      !Residual Projection
      if (itask .eq. 'Endite' .and. a%kfl_repro >= 1) call SetPointersComputeResidualProjection(1)
      
      if ((itask .eq. 'Endste') .and. (a%npp_stepi(13)>0)) call SetPointersComputeGpResidual(1)            
      
      !---------------------------------------------------------------
      !Dynamic, non-linear subscales
      !For Endste, only if tracking or transient, do not do it if non-linear
      if (itask .eq. 'Endste' .and. a%kfl_nolsg == 0 .and. (a%kfl_trasg == 1 .or. a%kfl_tacsg == 1)) then
         call SetPointersComputeSubscales(1)
      endif
      
      !For Endite, if non-linear subscales
      if (itask .eq. 'Endite' .and. a%kfl_nolsg == 1) then
         call SetPointersComputeNonLinearSubscales(1)
      endif
     
      !-----------------------------------------------------------
      !For adaptivity error indicators and subscales on the element boundaries
      if (itask .eq. 'Endste' .and. a%ErrorEstimatorTypeOfSubscales == 1) then
         call SetPointersComputeGradientProjection(1)            
      endif

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_nsc_OnIeltyChange => OnIeltyChange
      
      !Deallocate the procedure flags, so that they can be used in the next run
      !Initialize the ProcedureFlags
      call SetPointersVariables(100)     
      call SetPointersConvectionCoefficients(100)     
      call SetPointersConvectionJacobianGradientCoefficients(100)     
      call SetPointersDiffusionCoefficients(100)     
      call SetPointersNonLinearDiffusionCoefficients(100)     
      call SetPointersReactionCoefficients(100)
      call SetPointersExternalForces(100)        
      !call SetPointersInterpolateGradients(100)         
      call SetPointersComputeTaus(100)                  
      call SetPointersComputeTestf(100)                 
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeSubgridSpaceResidual(100)  
      call SetPointersNonLinearDerivatives(100)         
      !call SetPointersTurbulenceModel(100)              
      call SetPointersGetSubscales(100)   
      call SetPointersArbitraryLagrangianEulerian(100)
      
      !Also for the specific EndElmope Procedures
      call SetPointersComputeGpResidual(100)            
      call SetPointersComputeResidualProjection(100)     
      call SetPointersComputeGradientProjection(100)     
      call SetPointersComputeNonLinearSubscales(100)     
      call SetPointersComputeSubscales(100)
      
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


subroutine nsc_EndElmope_im(NSCompImplicitProblem,task)
   use Mod_nsc_BaseElmope
   use Mod_nsc_EndElmope_im
   use Mod_NSCompressibleImplicit
   implicit none
   class(NSCompressibleImplicitProblem), target :: NSCompImplicitProblem
   character(6) :: task
   integer(ip) :: aux_logic
   
   a=>NSCompImplicitProblem
   itask = task
   

   !Things to be done if endite
   if (itask .eq. 'Endite') then
      !Return if there is nothing to be done
      aux_logic = 0 
      
      if (a%kfl_repro >= 1) aux_logic = 1
      if (a%kfl_nolsg == 1) aux_logic = 1
      if (a%kfl_shock /= 0 ) aux_logic = 1

      if (aux_logic == 0) return
      
   
   !Things to be done if endste
   elseif (itask .eq. 'Endste') then
      !Return if there is nothing to be don
      aux_logic = 0
      if (a%kfl_trasg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      if (a%kfl_tacsg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1

      if (aux_logic == 0) return
      
      !If only dissipation for postprocessing and not a postprocessing step return
!      if (a%kfl_dispa == 1 .and. aux_logic == 1) then
!         if (a%istep < a%npp_inits) return
!         if (a%npp_stepi(15)>0 ) then
!            if(mod(a%istep,a%npp_stepi(15))/=0) then
!               return
!            endif
!         endif
!      endif
   endif
   

   !We force closed rule for smoothing
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'nsc_EnditeElmope')
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Hook
   call ProcHook_nsc_PreLoop
   


   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_EndElmope')

   !AllocateArrays in BaseElmope
   call AllocateBaseElmopeArrays
   
   !Hook
   call ProcHook_nsc_Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      
      !Hook
      call ProcHook_nsc_OnIeltyChange
      
      !Elmats to Zero
      call ProcHook_nsc_ElmatsToZero
      elmat=0.0_rp
      
      !Gathers
      call GatherBase
      !Hook
      call ProcHook_nsc_Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      call CalculateElmchl
      
      !Hook
      call ProcHook_nsc_PreGauss      
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook_nsc_InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call InterpolateBase
         !Hook
         call ProcHook_nsc_Interpolates
         
         !Physical Properties
         !Hook
         call ProcHook_nsc_PhysicalProp          
         
         !Linear transport variables         
         call ComputeLinearVariables

         !Nonlinear transport variables      
         call ProcPointer_nsc_ComputeVariables

         !Sound speed
         call nsc_ComputeSoundSpeed(accph,accvh,gptem,gpspd)
         
         !Compute gradients 
         call ComputeBaseGradients
         
         !Advection matrix      
         call ProcPointer_nsc_ComputeConvectionCoefficients

         !Diffusion matrix      
         call ProcPointer_nsc_ComputeDiffusionCoefficients

         !Other transport matrices      
         call ProcPointer_nsc_ComputeTransportCoefficients
         
         !Compute the stability parameters, only if required
         call ProcHook_nsc_ComputeTaus
         
         !Compute the stabilization test function, only if required
         call ProcHook_nsc_ComputeTestf
            
         !Compute Elext, Temporal Derivatives
         elexd = 0.0_rp
         elexm = 0.0_rp
         elexe = 0.0_rp
         !Compute vector of external forces
         call ProcPointer_nsc_ExternalForces  
         
         !Time integration
         eltemd = 0.0_rp
         eltemm = 0.0_rp
         elteme = 0.0_rp
         call nsc_TimeIntegrationToElTemp(e,1_ip,Integrator,a%dtinv,gpden,eltemd)
         call nsc_TimeIntegrationToElTemp(e,e%ndime,Integrator,a%dtinv,gpmom,eltemm)
         call nsc_TimeIntegrationToElTemp(e,1_ip,Integrator,a%dtinv,gpene,elteme)
            
         !InGaussElmats
         !Hook
         call ProcHook_nsc_InGaussElmats
            
         !InGaussElmats Assembly
         call ProcHook_nsc_InGaussElmatsAssembly
         
      enddo gauss_points      
  
      
      !PreDirichlet
      call ProcHook_nsc_PreDirichlet      
      
      !Assembly Endite
      !Hook
      call ProcHook_nsc_AssemblyEndite
      
      
   enddo elements
   !enddo
   
   !Hook
   call ProcHook_nsc_Finalizations
   
   !Deallocate arrays in BaseElmope
   call DeallocateBaseElmopeArrays

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_EndElmope')

   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook_nsc_PostLoop
   
   !Element
   call a%Mesh%ElementDeAlloc(e,a%Memor,a%EndLoopQuadrature,'nsc_EnditeElmope') 
   
end subroutine












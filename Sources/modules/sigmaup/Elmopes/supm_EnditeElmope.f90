module Mod_supm_EnditeElmope
   use Mod_supm_BaseElmope
   use Mod_supm_TemperatureCoupling
   use Mod_supm_ExternalForces
   use Mod_supm_ComputeAdvectionVelocity
   use Mod_supm_PhysicalProperties
   use Mod_supm_TemporalDerivatives
   use Mod_supm_LevelSetCoupling
   use Mod_supm_InterpolateGradients
   use Mod_supm_DiscontCapturing
   use Mod_supm_ElemVElastic
   use Mod_supm_NonLinearElem
   use Mod_supm_ComputeSplitTermsProjection
   use Mod_supm_ComputeGpSplitTerms
   use Mod_supm_ElemLogarithmic
   use Mod_supm_LogarithmicProblem
   use Mod_supm_ElemNoVElastic
   use Mod_supm_ComputeTaus
   use Mod_supm_ComputeTauSmoothing
   use Mod_supm_ComputeResidualProjection
   use Mod_supm_ComputeGpResidual
   use Mod_supm_InterpolateResidualProjection
   use Mod_supm_SubgridSpaceResidual
   use Mod_supm_ComputeSubscales
   use Mod_supm_InterpolateSplitTermsProjection
   use Mod_supm_SubgridSpaceSplitTerms
   use Mod_supm_ComputeTestf
   use Mod_supm_ComputeVorticity
   implicit none
    
    
contains

   !*********************
   !SetPointers
   !*********************
   subroutine SetPointers
      implicit none
      integer(ip) :: kfl_nonlinear, nelty
      logical     :: aux_logic
      !External Procedures
      call ResetProcedureComposition
      !Set All Pointers To NULL() (in Mod_supm_BaseElmope)
      call SetPointersAndHooksToNULL
      !PreAssembly to use in free-surface and in enriched elements
      ProcPointer%PreAssembly_sup => NULLSUB
      imat=1
      
      !----------------------------------------------------- 
      !Initialize the ProcedureFlags
      !-----------------------------------------------------
      call SetPointersNonLinearElem(0,'Endite')
      call setPointersElemVElastic(0,'Endite')
      call SetPointersDiscontCapturing(0,'Endite')
      
      call SetPointersLevelSetCoupling(0,'Endite')
      call SetPointersPhysicalPropertiesSUP(0,'Endite')
      call SetPointersAdvectionVelocity(0)
      call SetPointersComputeTaus(0)
      call SetPointersComputeTauSmoothing(0)
      call SetPointersExternalForces(0)
      call SetPointersTemperatureCoupling(0) 
      call SetPointersNoVElastic(0,'Endite')
      
      call SetPointersComputeGpResidual(0)         
      call SetPointersComputeResidualProjection(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeSubgridSpaceResidual(0)
      call SetPointersComputeSubscales_sup(0)   
      
      call SetPointersComputeGpSplitTerms(0)
      call SetPointersComputeSplitTermsProjection(0)
      call SetPointersInterpolateSplitTermsProjection(0)
      call SetPointersComputeSubgridSpaceSplitTerms(0)

      call SetPointersInterpolateGradients(0)
      call SetPointersTemporalDerivatives(0)
      
      call SetPointersLogarithmicComponents(0)
      call SetPointersElemLCR(0,'Endite')
      call SetPointersComputeTestf(0)
      call SetPointersComputeVorticity(0)
      
      !Element allocation
	   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope')  
      call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'supm_EnditeElmope') 
      
      call SetPointersNoVElastic(1,'Endite')
      call SetPointersElemVElastic(1,'Endite')

      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) call SetPointersNonLinearElem(1,'Endite')
      !Temporal Derivatives
      if((itask .eq. 'Endite') .and. a%MatProp(imat)%lawvi < 0)then
        call SetPointersTemporalDerivatives(1)        
      end if

      auxdim=e%ndime
      
      if (itask .eq. 'Endite') then
         !Residual Projection
         call SetPointersComputeResidualProjection(1)
         !Split-Oss
         call SetPointersComputeSplitTermsProjection(1)
      end if   
      
      !DynamicSS ------------------------------------
      !For Endste, only if tracking or transient, do not do it if non-linear
      if (itask .eq. 'Endste' .and. a%kfl_nolsg == 0 .and. (a%kfl_trasg == 1 .or. a%kfl_tacsg == 1)) then
         call SetPointersComputeSubscales_sup(1)
      endif
      
      !Tau smoothing --------------------------------
      if (itask .eq. 'Endite') call SetPointersComputeTauSmoothing(1)
      
      
      !Vorticity ------------------------------------
      if (itask .eq. 'Endste') then
 
         !Logics for deciding if we need to compute vorticity
         aux_logic = .false.
         if (a%npp_stepi(10) /= 0) then   
           if (mod(a%istep,a%npp_stepi(10))==0) then
              aux_logic = .true.
           endif
         endif
        
         !If vorticity needs to be computed
         if (aux_logic .eqv. .true.) then   
            call SetPointersComputeVorticity(1)
         endif
      endif
      
      
      !Shock Capturing -------------------------------
      call SetPointersDiscontCapturing(1,'Endite')
           
      !Physical Properties
      call SetPointersPhysicalPropertiesSUP(1,'Endite')

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange

      !Level Set
      call SetPointersLevelSetCoupling(1,'Endite')
      
      !Advection Velocity
      call SetPointersAdvectionVelocity(1)
      
      !Compute External Forces
      call SetPointersExternalForces(1)
      
      !Temperature Coupling
      call SetPointersTemperatureCoupling(1) 
      
      !Logarithm conformation tensor model
      call SetPointersLogarithmicComponents(1)
      call SetPointersElemLCR(1,'Endite')
      
      !--------------------------------------------------------------------------
      !We deallocate the procedure flags, so that they can be set the next time
      !--------------------------------------------------------------------------
      call SetPointersNoVElastic(100,'Elmope')
      !call SetPointersSplitOss(100,'Endite')
      call SetPointersComputeSplitTermsProjection(100)
      call SetPointersComputeSubgridSpaceSplitTerms(100)
      call SetPointersComputeGpSplitTerms(100)
      call SetPointersNonLinearElem(100,'Endite')
      call setPointersElemVElastic(100,'Endite')
      call SetPointersDiscontCapturing(100,'Endite')
      call SetPointersInterpolateGradients(100)
      call SetPointersComputeResidualProjection(100)
      call SetPointersLevelSetCoupling(100,'Endite')
      call SetPointersTemporalDerivatives(100)
      call SetPointersPhysicalPropertiesSUP(100,'Endite')
      call SetPointersAdvectionVelocity(100)
      call SetPointersComputeTaus(100)
      call SetPointersComputeTauSmoothing(100)
      call SetPointersExternalForces(100)
      call SetPointersTemperatureCoupling(100) 
      
      call SetPointersComputeGpResidual(100)            
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeSubgridSpaceResidual(100)
      call SetPointersComputeSubscales_sup(100)   
      call SetPointersInterpolateSplitTermsProjection(100)
      
      call SetPointersLogarithmicComponents(100)
      call SetPointersElemLCR(100,'Endite')
      call SetPointersComputeTestf(100)
      call SetPointersComputeVorticity(100)
      
   end subroutine
   
   !Change of type of element
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
   !***********************************************
   !General Subroutines called by SUPM_ENDITEELMOPE
   !************************************************
   
   subroutine ComputeElmlen
      call e%elmlen
   end subroutine
   
   subroutine ComputeCGDerivatives
      call e%elmdcg
   end subroutine
   
   subroutine ComputeDerivatives
      call e%elmder
   end subroutine  
     
end module



!*****************************
!SUPM_ENDITEELMOPE subroutine  
!*****************************

subroutine supm_EnditeElmope(TFNSProblem,task)
   use Mod_supm_EnditeElmope
   !use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem), target :: TFNSProblem
   character(6) :: task
   integer(ip) :: aux_logic
   a=>TFNSProblem
   
   itask=task
   
   !Reset
   !Things to be done if endite
   if (itask .eq. 'Endite') then
      !Return if there is nothing to be done
      aux_logic = 0 
      
      if (a%kfl_repro >= 1) aux_logic = 1
      if (a%kfl_nolsg == 1) aux_logic = 1
      if (a%kfl_Tausm == 1) aux_logic = 1
      if (a%kfl_shock == 1) aux_logic = 1
      if (aux_logic == 0) return
   !Things to be done if endste
   elseif (itask .eq. 'Endste') then
      !Return if there is nothing to be don
      aux_logic = 0
      if (a%kfl_dispa /= 0) aux_logic = aux_logic + 1 !dissipation
      if (a%kfl_trasg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      if (a%kfl_tacsg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      if (aux_logic == 0) return
   endif
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sub_elem')

   !call Initialize
    ielem = 1

   !Set the Pointers for execution
   call SetPointers
   
   !Allocate Matrices and Arrays
   call AllocateBaseElmopeArrays
   call ProcHook_PreLoop
   
   !Initializations
   call ProcHook_Initializations
   
   call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','supm_EnditeElmope')
   
   call a%Mesh%GetNelem(nelem)

   do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)
      
      call ComputeCGDerivatives
      call ComputeElmlen
      !Hook
      call ProcHook_OnIeltyChange
      
      !Hook
      call ProcHook_PreGauss

      !Elmats to Zero
      call ProcHook_ElmatsToZero
      
      !Gathers
      call VelocityAndSigmaGathers
   
      call e%gather(1_ip   ,elpre(1:e%pnode,1),a%press(1:e%pnode,1))
      
      !Hook
      call ProcHook_Gathers
      
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)

      do igaus = 1,e%pgaus
         e%igaus = igaus
         !Hook
         call ProcHook_InGauss
         
         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call InterpolateGpVelocityAndSigma
         
         !Hook
         call ProcHook_Interpolates
         
         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)
         
         !vista = 0.0_rp
         !Hook
         call ProcHook_PhysicalProp          
         
         !Advection velocity      
         call ProcPointer%ComputeAdvectionVelocity_sup
         
         !Advection velocity norm
         call ProcPointer%ComputeAdvectionVelocityNorm_sup
         
         !Compute a·grad(V)
         call ProcPointer%ComputeAGradV_sup
         
         !Compute estabilization terms if these are dynamics
         call ProcPointer%TauConstitutive
         
         !Compute the stabilization test function, only if required
         call ProcHook_ComputeTestf

         elext = 0.0_rp
         elextC= 0.0_rp
         elextS= 0.0_rp
         elextSEstab=0.0_rp
         elextSEstab2=0.0_rp
         elextSEstab3=0.0_rp
         elextSEstab4=0.0_rp
         elextSMat=0.0_rp
         elextSEstabMat=0.0_rp
         elextEstab=0.0_rp
         elextEstab2=0.0_rp
         elextEstab3=0.0_rp
         elextEstab4=0.0_rp
         elextEstab5=0.0_rp
         
         !Time integration
         call  nsi_TimeIntegrationToEltemp(e,Integrator,acden,a%dtinv,gpvel,elext)
         call  ProcPointer%TemporalDerivatives  
         
         !Compute vector of external forces
         call  ProcPointer%ExternalForces_sup                  
                     
         !InGaussElmats
         call ProcHook_InGaussElmats
      enddo

      !Pre Assembly Modifications
      !Pointer
      call  ProcPointer%PreAssembly_sup       
            
      !Assembly Endite
      !Hook
      call ProcHook_AssemblyEndite
   enddo
   
   call a%Memor%dealloc(1_ip   ,e%ndime,grpre,'grpre','supm_EnditeElmope')
   
   call ProcHook_Finalizations
   
   !Deallocate Matrices and Arrays
   call DeallocateBaseElmopeArrays
   
   call ProcHook_PostLoop
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','sub_elem')
      
end subroutine

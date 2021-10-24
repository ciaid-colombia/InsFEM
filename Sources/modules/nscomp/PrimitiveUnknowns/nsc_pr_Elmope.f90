module Mod_nsc_pr_elmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_NSCompressiblePrimitive
   use Mod_Element
   use Mod_php_SetTimeIntegrator
!   use Mod_NSCompressiblePrimitiveElement
   use Mod_NSCompressibleSubroutines
   use Mod_ConvectiveElement
   use Mod_NSCompressibleImplicitElement
   use Mod_TimeIntegrator
   use Mod_NscExacso    
   use Mod_nsc_pr_elmdir


   !For setting pointers
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_ComputeVariables
   use Mod_nsc_pr_MatricesContribution
   use Mod_nsc_pr_ComputeTransientCoefficients
   use Mod_nsc_pr_ComputeConvectionCoefficients
   use Mod_nsc_pr_ComputeConvectionJacobianGradientCoefficients
   use Mod_nsc_pr_ComputeDiffusionCoefficients
   use Mod_nsc_pr_ComputeNonLinearDiffusionCoefficients
   use Mod_nsc_pr_ComputeDiffusiveTerm
   use Mod_nsc_pr_ComputeShockCapturing
   use Mod_nsc_pr_ComputeReactionCoefficients
   use Mod_nsc_pr_ComputeTaus
   use Mod_nsc_pr_ComputeTestf
   use Mod_nsc_pr_NonLinearDerivatives
   use Mod_nsc_pr_ExternalForces
   use Mod_nsc_pr_ComputeGpResidual
   use Mod_nsc_pr_InterpolateGradientProjection
   !use Mod_nsc_pr_InterpolateGradients
   !use Mod_nsc_pr_TurbulenceModel
   use Mod_nsc_pr_ComputeResidualProjection
   use Mod_nsc_pr_InterpolateResidualProjection
   use Mod_nsc_pr_SubgridSpaceResidual
   !use Mod_nsm_pr_ComputeTauSmoothing
   use Mod_nsc_pr_ComputeSubscales
   !use Mod_nsc_pr_ComputeNonLinearSubscales
   !use Mod_nsc_pr_ComputeDissipation
   use Mod_nsc_pr_DynSubsElmope
   !use Mod_nsc_pr_NewtonRaphson
   use Mod_nsc_pr_HangingNodes
   use Mod_nsc_pr_ArbitraryLagrangianEulerian
   use Mod_nsc_pr_OrthogonalSubscales
   use Mod_nsc_pr_ConfinedFLow
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
            
      !Reset Procedure Composition
      call ResetProcedureComposition
      
      !Set All Pointers To NULLSUB (in Mod_nsc_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Set the Specific Pointers for this subroutine
      ProcPointer_nsc_pr_PreAssembly     => NULLSUB
      ProcPointer_nsc_pr_PostGaussElmats => NULLSUB!PostGaussElmats
      ProcPointer_nsc_pr_elmbdq      => nsc_pr_elmbdq
      ProcPointer_nsc_pr_elmbdn      => nsc_pr_elmbdn
      ProcPointer_nsc_pr_elmbdg      => nsc_pr_elmbdg
      ProcPointer_nsc_pr_elmbmq      => nsc_pr_elmbmq
      ProcPointer_nsc_pr_elmbmn      => nsc_pr_elmbmn
      ProcPointer_nsc_pr_elmbmg      => nsc_pr_elmbmg
      ProcPointer_nsc_pr_elmbeq      => nsc_pr_elmbeq
      ProcPointer_nsc_pr_elmben      => nsc_pr_elmben
      ProcPointer_nsc_pr_elmbeg      => nsc_pr_elmbeg
      ProcPointer_nsc_pr_elmrhd      => nsc_pr_elmrhd
      ProcPointer_nsc_pr_elmrhm      => nsc_pr_elmrhm
      ProcPointer_nsc_pr_elmrhe      => nsc_pr_elmrhe
                             
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersVariables(0)     
      call SetPointersMatricesContribution(0)     
      call SetPointersTransientCoefficients(0)     
      call SetPointersConvectionCoefficients(0)     
      call SetPointersConvectionJacobianGradientCoefficients(0)     
      call SetPointersDiffusionCoefficients(0)     
      call SetPointersNonLinearDiffusionCoefficients(0)     
      call SetPointersDiffusiveTerm(0)     
      call SetPointersShockCapturing(0)
      call SetPointersReactionCoefficients(0)
      call SetPointersExternalForces(0)        
     !call SetPointersInterpolateGradients(0)         
      call SetPointersComputeTaus(0)                  
      call SetPointersComputeTestf(0)                 
      call SetPointersComputeGpResidual(0)            
      call SetPointersInterpolateGradientProjection(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeSubgridSpaceResidual(0)  
      call SetPointersNonLinearDerivatives(0)         
      call SetPointersGetSubscales(0)   
      call SetPointersDynSubsElmope(0)
      !call SetPointersNewtonRaphson(0)
      call SetPointersHangingNodes(0)
      call SetPointersArbitraryLagrangianEulerian(0)
      call SetPointersOrthogonalSubscales(0)
      call SetPointersConfinedFlow(0)
      
      !-------------------------------------------------------------------
      !Now we set the required pointers 

      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
          call SetPointersNonLinearDerivatives(1)
      endif

      call SetPointersVariables(1)
      call SetPointersMatricesContribution(1)     
      
      !Transient term
      call SetPointersTransientCoefficients(1)

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
             call ConcatenateProcedures(ProcHook_nsc_pr_InGauss,InGaussVolumesNonLinear)
         endif
         call SetPointersDiffusiveTerm(1)     
         call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,InGaussElmatsDiffusion)
      endif

      !Reaction Exists
      if (a%kfl_react == 1) then
         call SetPointersReactionCoefficients(1)
         call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,InGaussElmatsReaction)
      endif

      !External forces 
      call SetPointersExternalForces(1)
      
      call SetPointersComputeTaus(1)
      
      call SetPointersOrthogonalSubscales(1)
      
                    
      !call SetPointersTurbulenceModel(1)
      call SetPointersComputeTestf(1)
      call SetPointersHangingNodes(1)
      
      !Dynamic subscales
      if (a%kfl_tacsg == 1) call SetPointersDynSubsElmope(1)
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_nsc_pr_OnIeltyChange => OnIeltyChange

      call SetPointersConfinedFlow(1)
      
      !-----------------------------------------------------------------------
      !We deallocate the procedure flags, so that they can be set the next time
      call SetPointersVariables(100)     
      call SetPointersMatricesContribution(100)     
      call SetPointersTransientCoefficients(100)     
      call SetPointersConvectionCoefficients(100)     
      call SetPointersConvectionJacobianGradientCoefficients(100)     
      call SetPointersDiffusionCoefficients(100)     
      call SetPointersNonLinearDiffusionCoefficients(100)     
      call SetPointersDiffusiveTerm(100)     
      call SetPointersShockCapturing(100)
      call SetPointersReactionCoefficients(100)
      call SetPointersExternalForces(100)        
      !call SetPointersInterpolateGradients(100)         
      call SetPointersComputeTaus(100)                  
      call SetPointersComputeTestf(100)                 
      call SetPointersComputeGpResidual(100)            
      call SetPointersInterpolateGradientProjection(100)
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeSubgridSpaceResidual(100)  
      call SetPointersNonLinearDerivatives(100)         
      !call SetPointersTurbulenceModel(100)              
      call SetPointersGetSubscales(100)   
      call SetPointersDynSubsElmope(100)
      !call SetPointersNewtonRaphson(100)
      call SetPointersHangingNodes(100)
      call SetPointersArbitraryLagrangianEulerian(100)
      call SetPointersOrthogonalSubscales(100)
      call SetPointersConfinedFlow(100)
      
   end subroutine
   
  
   !----------------------------------------------------------

   !InGauss Matrices

   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
   end subroutine

   
   subroutine InGaussElmatsDiffusion
      implicit none

         !Compute contributions to elemental matrix : Block q,rho 
         call nsc_pr_elmbdq_diff(e,dvol,KGdd,elmdq) !Radial Damping Diffusion
         !Compute contributions to elemental matrix : Block q,vel 
         call nsc_elmbmq_diff(e,dvol,Kmm,Kem,LTdm,LTde,elmmq)
         !Compute contributions to elemental matrix : Block q,tem 
         call nsc_elmbeq_diff(e,dvol,Kee,LTde,elmeq)
         !Compute contributions to elemental matrix : Block n,vel
         call nsc_elmbmn_diff(e,dvol,KGmm,Kmm,Kem,LTmm,LTme,elmmn)
         !Compute contributions to elemental matrix : Block n,tem 
         call nsc_elmben_diff(e,dvol,Kee,LTme,elmen)
         !Compute contributions to elemental matrix : Block g,vel 
         call nsc_elmbmg_diff(e,dvol,KGem,Kmm,Kem,LTem,LTee,elmmg)
         !Compute contributions to elemental matrix : Block g,tem 
         call nsc_elmbeg_diff(e,dvol,KGee,Kee,LTee,elmeg)

   end subroutine

   subroutine InGaussElmatsReaction
      implicit none

         !Compute contributions to elemental matrix : Block q,rho
         call nsc_elmbdq_reac(e,dvol,Sdd,Smd,LTdd,LTdm,elmdq)
         !Compute contributions to elemental matrix : Block q,mom
         call nsc_elmbmq_reac(e,dvol,Smm,LTdm,elmmq)
         !Compute contributions to elemental matrix : Block q,ene
         call nsc_elmbeq_reac(e,dvol,See,LTde,elmeq)
         !Compute contributions to elemental matrix : Block n,rho
         call nsc_elmbdn_reac(e,dvol,Sdd,Smd,Sed,LTmd,LTmm,LTme,elmdn)
         !Compute contributions to elemental matrix : Block n,mom
         call nsc_elmbmn_reac(e,dvol,Smm,Sem,LTmm,LTme,elmmn)
         !Compute contributions to elemental matrix : Block n,ene
         call nsc_elmben_reac(e,dvol,See,LTme,elmen)
         !Compute contributions to elemental matrix : Block g,rho
         call nsc_elmbdg_reac(e,dvol,Sdd,Smd,Sed,LTed,LTem,LTee,elmdg)
         !Compute contributions to elemental matrix : Block g,mom 
         call nsc_elmbmg_reac(e,dvol,Smm,Sem,LTem,LTee,elmmg)
         !Compute contributions to elemental matrix : Block g,ene 
         call nsc_elmbeg_reac(e,dvol,See,LTee,elmeg)

   end subroutine
   !----------------------------------------------------------

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module   
   
   
   
   
   
   
   
!NSC_ELMOPE subroutine   
subroutine nsc_pr_elmope(NSCompPrimitiveProblem)
   use Mod_nsc_pr_elmope
   implicit none
   class(NSCompressiblePrimitiveProblem), target :: NSCompPrimitiveProblem
   a=>NSCompPrimitiveProblem
  

   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_elmope')
   
   call ProcHook_nsc_pr_PreAllocate
   
   !Allocate Arrays in BaseElmope
   call AllocateBaseElmopeArrays

   !Matrices Alloc
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_pr_elmope')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%mnode,elmdq,'elmdq','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%ndime,e%mnode,elmmq,'elmmq','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%mnode,elmeq,'elmeq','nsc_pr_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,e%mnode,elmdn,'elmdn','nsc_pr_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmmn,'elmmn','nsc_pr_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,e%mnode,elmen,'elmen','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%mnode,elmdg,'elmdg','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%ndime,e%mnode,elmmg,'elmmg','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,e%mnode,elmeg,'elmeg','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,elrhd,'elrhd','nsc_pr_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,elrhm,'elrhm','nsc_pr_elmope')
   call a%Memor%alloc(e%mnode,elrhe,'elrhe','nsc_pr_elmope')
   
   !Hook
   call ProcHook_nsc_pr_Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook_nsc_pr_OnIeltyChange
      
      !ElmatsToZero
      call ProcHook_nsc_pr_ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      elmdq=0.0_rp
      elmdn=0.0_rp
      elmdg=0.0_rp
      elrhd=0.0_rp
      elmmq=0.0_rp
      elmmn=0.0_rp
      elmmg=0.0_rp
      elrhm=0.0_rp
      elmeq=0.0_rp
      elmen=0.0_rp
      elmeg=0.0_rp
      elrhe=0.0_rp
      
      !Gathers
      call GatherBase
      !Hook
      call ProcHook_nsc_pr_Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elvel,chale)

      !Hook
      call ProcHook_nsc_pr_PreGauss
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook_nsc_pr_InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate         
         call InterpolateBase

         !Hook
         call ProcHook_nsc_pr_Interpolates
         
         !PostHook 
         call ProcHook_nsc_pr_PostInterpolate
         
         !Hook
         call ProcHook_nsc_pr_PhysicalProp    
         
         !Linear transport variables         
         call ComputeLinearVariables

         !Nonlinear transport variables      
         call ProcPointer_nsc_pr_ComputeVariables

         !Sound speed
         call nsc_ComputeSoundSpeed(accph,accvh,gpadt,gpspd)
         
         !Compute gradients 
         call ComputeBaseGradients
         
         !Transient matrix 
         call ProcPointer_nsc_pr_ComputeTransientCoefficients

         !Advection matrix      
         call ProcPointer_nsc_pr_ComputeConvectionCoefficients

         !Diffusion matrix      
         call ProcPointer_nsc_pr_ComputeDiffusionCoefficients

         !Other transport matrices      
         call ProcPointer_nsc_pr_ComputeTransportCoefficients
         
         !Compute the stability parameters 
         call ProcHook_nsc_pr_ComputeTaus
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call ProcHook_nsc_pr_ComputeTestf
         
         !Compute vector of external forces (exact sol)
         elexp = 0.0_rp
         elexv = 0.0_rp
         elext = 0.0_rp
         call ProcPointer_nsc_pr_ExternalForces 
         
         !Time integration
         eltemp = 0.0_rp
         eltemv = 0.0_rp
         eltemt = 0.0_rp
         call nsc_TimeIntegrationToElTemp(e,1_ip,Integrator,a%dtinv,gppre,eltemp)
         call nsc_TimeIntegrationToElTemp(e,e%ndime,Integrator,a%dtinv,gpvel,eltemv)
         call nsc_TimeIntegrationToElTemp(e,1_ip,Integrator,a%dtinv,gptem,eltemt)
         
         !InGaussElmats
         !Euler contribution
         call ComputeEulerContribution
         !Compute contributions to elemental matrix : Block q,pre
         call ProcPointer_nsc_pr_elmbdq(e,dvol,Edd,Emd,Eed,LTdd,LTdm,LTde,elmdq)
         !Compute contributions to elemental matrix : Block q,vel 
         call ProcPointer_nsc_pr_elmbmq(e,dvol,Edm,Emm,Eem,LTdd,LTdm,LTde,elmmq)
         !Compute contributions to elemental matrix : Block q,tem 
         call ProcPointer_nsc_pr_elmbeq(e,dvol,Ede,Eme,Eee,LTdd,LTdm,LTde,elmeq)
         !Compute contributions to elemental matrix : Block n,pre
         call ProcPointer_nsc_pr_elmbdn(e,dvol,Edd,Emd,Eed,LTmd,LTmm,LTme,elmdn)
         !Compute contributions to elemental matrix : Block n,vel
         call ProcPointer_nsc_pr_elmbmn(e,dvol,Edm,Emm,Eem,LTmd,LTmm,LTme,elmmn)
         !Compute contributions to elemental matrix : Block n,tem 
         call ProcPointer_nsc_pr_elmben(e,dvol,Ede,Eme,Eee,LTmd,LTmm,LTme,elmen)
         !Compute contributions to elemental matrix : Block g,pre
         call ProcPointer_nsc_pr_elmbdg(e,dvol,Edd,Emd,Eed,LTed,LTem,LTee,elmdg)
         !Compute contributions to elemental matrix : Block g,vel 
         call ProcPointer_nsc_pr_elmbmg(e,dvol,Edm,Emm,Eem,LTed,LTem,LTee,elmmg)
         !Compute contributions to elemental matrix : Block g,tem 
         call ProcPointer_nsc_pr_elmbeg(e,dvol,Ede,Eme,Eee,LTed,LTem,LTee,elmeg)
         !Compute RHS contribution
         call ComputeRHSContribution(e,eltemp(1),eltemv,eltemt(1),elexp(1),elexv,elext(1),Erhd,Erhm,Erhe)
        !Compute contributions to RHS : Block pre
         call ProcPointer_nsc_pr_elmrhd(e,dvol,Erhd,Erhm,Erhe,LTdd,LTdm,LTde,elrhd)
         !Compute contributions to RHS : Block vel
         call ProcPointer_nsc_pr_elmrhm(e,dvol,Erhd,Erhm,Erhe,LTmd,LTmm,LTme,elrhm)
         !Compute contributions to RHS : Block tem
         call ProcPointer_nsc_pr_elmrhe(e,dvol,Erhd,Erhm,Erhe,LTed,LTem,LTee,elrhe)

         !Hook
         call ProcHook_nsc_pr_InGaussElmats
         
      enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer_nsc_pr_PostGaussElmats

      !Matrix composition
   
      ! Assembly elmdq to elmat
      elmat(1,1:e%pnode,1,1:e%pnode) = elmat(1,1:e%pnode,1,1:e%pnode) + elmdq(1:e%pnode,1:e%pnode)

      ! Assembly elmmq to elmat
      elmat(1,1:e%pnode,2:e%ndime+1,1:e%pnode) = elmat(1,1:e%pnode,2:e%ndime+1,1:e%pnode) + elmmq(1:e%pnode,1:e%ndime,1:e%pnode)

      ! Assembly elmeq to elmat
      elmat(1,1:e%pnode,e%ndime+2,1:e%pnode) = elmat(1,1:e%pnode,e%ndime+2,1:e%pnode) + elmeq(1:e%pnode,1:e%pnode)

      ! Assembly elmdn to elmat
      elmat(2:e%ndime+1,1:e%pnode,1,1:e%pnode) = elmat(2:e%ndime+1,1:e%pnode,1,1:e%pnode) + elmdn(1:e%ndime,1:e%pnode,1:e%pnode)

      ! Assembly elmmn to elmat
      elmat(2:e%ndime+1,1:e%pnode,2:e%ndime+1,1:e%pnode) = elmat(2:e%ndime+1,1:e%pnode,2:e%ndime+1,1:e%pnode) + elmmn(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      
      ! Assembly elmen to elmat
      elmat(2:e%ndime+1,1:e%pnode,e%ndime+2,1:e%pnode) = elmat(2:e%ndime+1,1:e%pnode,e%ndime+2,1:e%pnode) + elmen(1:e%ndime,1:e%pnode,1:e%pnode)

      ! Assembly elmdg to elmat
      elmat(e%ndime+2,1:e%pnode,1,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,1,1:e%pnode) + elmdg(1:e%pnode,1:e%pnode)

      ! Assembly elmmg to elmat
      elmat(e%ndime+2,1:e%pnode,2:e%ndime+1,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,2:e%ndime+1,1:e%pnode) + elmmg(1:e%pnode,1:e%ndime,1:e%pnode)

      ! Assembly elmeg to elmat
      elmat(e%ndime+2,1:e%pnode,e%ndime+2,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,e%ndime+2,1:e%pnode) + elmeg(1:e%pnode,1:e%pnode)

      ! Assembly  elrhd to elrhs
      elrhs(1,1:e%pnode) = elrhd(1:e%pnode) + elrhs(1,1:e%pnode)

      ! Assembly elrhm to elrhs
      elrhs(2:e%ndime+1,1:e%pnode) = elrhm(1:e%ndime,1:e%pnode) + elrhs(2:e%ndime+1,1:e%pnode)
   
      ! Assembly  elrhp to elrhs
      elrhs(e%ndime+2,1:e%pnode) = elrhe(1:e%pnode) + elrhs(e%ndime+2,1:e%pnode)
   
      !PreDirichlet
      call ProcHook_nsc_pr_PreDirichlet
      
      !Dirichlet Boundary Conditions
      call nsc_pr_elmdir(a,e,elmat,elrhs)
            
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)

   enddo elements

   !Hook
   call ProcHook_nsc_pr_Finalizations
   
   call a%Memor%dealloc(a%ndofn,size(elmat,2),a%ndofn,size(elmat,4),elmat,'elmat','nsc_pr_elmope')
   call a%Memor%dealloc(a%ndofn,size(elrhs,2),elrhs,'elrhs','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmdq,1),size(elmdq,2),elmdq,'elmdq','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmmq,1),e%ndime,size(elmmq,3),elmmq,'elmmq','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmeq,1),size(elmeq,2),elmeq,'elmeq','nsc_pr_elmope')
   call a%Memor%dealloc(e%ndime,size(elmdn,2),size(elmdn,3),elmdn,'elmdn','nsc_pr_elmope')
   call a%Memor%dealloc(e%ndime,size(elmmn,2),e%ndime,size(elmmn,4),elmmn,'elmmn','nsc_pr_elmope')
   call a%Memor%dealloc(e%ndime,size(elmen,2),size(elmen,3),elmen,'elmen','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmdg,1),size(elmdg,2),elmdg,'elmdg','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmmg,1),e%ndime,size(elmmg,3),elmmg,'elmmg','nsc_pr_elmope')
   call a%Memor%dealloc(size(elmeg,1),size(elmeg,2),elmeg,'elmeg','nsc_pr_elmope')
   call a%Memor%dealloc(size(elrhd,1),elrhd,'elrhd','nsc_pr_elmope')
   call a%Memor%dealloc(e%ndime,size(elrhm,2),elrhm,'elrhm','nsc_pr_elmope')
   call a%Memor%dealloc(size(elrhe,1),elrhe,'elrhe','nsc_pr_elmope')
   !Arrays Deallocations
   call DeallocateBaseElmopeArrays
   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsc_pr_elmope')
   
!   deb_PostprocessMatrix = 1
end subroutine

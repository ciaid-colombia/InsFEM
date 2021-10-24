module Mod_lmn_EndElmope
   use MPI
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_ComputeAdvectionVelocity
   use Mod_lmn_ComputeTemperature
   use Mod_lmn_ComputeDensity
   use Mod_lmn_NonLinearDerivatives
   use Mod_lmn_ExternalForces
   use Mod_lmn_InterpolateGradients
   use Mod_lmn_ComputeGpResidual
   use Mod_lmn_ComputeResidualProjection
   use Mod_lmn_InterpolateResidualProjection
   use Mod_lmn_SubgridSpaceResidual
   use Mod_lmn_ComputeTaus
   use Mod_lmn_ComputeTestf
   use Mod_lmn_ComputeSubscales
   use Mod_lmn_ComputeNonLinearSubscales
   use Mod_lmn_HangingNodes
   use Mod_lmn_ComputeGradientProjection
   implicit none
  
contains

   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      integer(ip) :: nelty
      
      call ResetProcedureComposition
      !Set All Pointers To NULLSUB (in Mod_lmn_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Initialize the ProcedureFlags
      call SetPointersAdvectionVelocity(0)
      call SetPointersTemperature(0)
      call SetPointersDensity(0)
      call SetPointersExternalForces(0)
      call SetPointersInterpolateGradients(0)
      call SetPointersComputeTaus(0)
      call SetPointersComputeTestf(0)
      call SetPointersComputeGpResidual(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeSubgridSpaceResidual(0)
      call SetPointersNonLinearDerivatives(0)
      call SetPointersGetSubscales(0)
      call SetPointersHangingNodes(0)
 
      !Also for the specific EndElmope Procedures
      call SetPointersComputeNonLinearSubscales(0)     
      call SetPointersComputeResidualProjection(0)     
      call SetPointersComputeGradientProjection(0) 
      call SetPointersComputeSubscales(0)              
               
      !Now we set the required pointers         
      ProcPointer%TimeIntegrationToElext => lmn_TimeIntegrationToElext
      call SetPointersAdvectionVelocity(1)
      call SetPointersTemperature(1)     
      call SetPointersDensity(1)     
      call SetPointersNonLinearDerivatives(1)
      call SetPointersExternalForces(1)
      call SetPointersHangingNodes(1)


      !-----------------------------------------------------------
      !Residual Projection
      if (itask .eq. 'Endite' .and. a%kfl_repro == 1) call SetPointersComputeResidualProjection(1)

      !-----------------------------------------------------------
      !For adaptivity error indicators and subscales on the element boundaries
      if (itask .eq. 'Endste' .and. a%kfl_adapsgs == 1) call SetPointersComputeGradientProjection(1)

      !---------------------------------------------------------------
      !Dynamic, non-linear subscales
      !For Endste, only if tracking or transient, do not do it if non-linear
      if (itask .eq. 'Endste' .and. a%kfl_nolsg == 0 .and. (a%kfl_trasg == 1 .or. a%kfl_tacsg == 1)) then
         call SetPointersComputeSubscales(1)
      endif
      
      !For Endite, if non-linear subscales
      if (itask .eq. 'Endite' .and. a%kfl_nolsg == 1) call SetPointersComputeNonLinearSubscales(1)

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

      !Deallocate the procedure flags, so that they can be used in the next run
      !Initialize the ProcedureFlags
      call SetPointersAdvectionVelocity(100)
      call SetPointersTemperature(100)     
      call SetPointersDensity(100)     
      call SetPointersExternalForces(100)        
      call SetPointersInterpolateGradients(100)         
      call SetPointersComputeTaus(100)                  
      call SetPointersComputeTestf(100)                 
      call SetPointersComputeGpResidual(100)            
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeSubgridSpaceResidual(100)  
      call SetPointersNonLinearDerivatives(100)         
      call SetPointersGetSubscales(100)   
      call SetPointersHangingNodes(100)
      
      !Also for the specific EndElmope Procedures
      call SetPointersComputeNonLinearSubscales(100)     
      call SetPointersComputeResidualProjection(100)     
      call SetPointersComputeGradientProjection(100) 
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


subroutine lmn_EndElmope(LMProblem,task)
   use Mod_lmn_BaseElmope
   use Mod_lmn_EndElmope
   use Mod_LowMach
   implicit none
   class(LowMachProblem), target :: LMProblem
   character(6) :: task
   integer(ip) :: aux_logic
   
   a=>LMProblem
   itask = task

   !Things to be done if endite
   if (itask .eq. 'Endite') then
      !Return if there is nothing to be done
      aux_logic = 0
      
      if (a%kfl_repro == 1) aux_logic = 1
      if (a%kfl_nolsg == 1) aux_logic = 1
      
      if (aux_logic == 0) return
   
   !Things to be done if endste
   elseif (itask .eq. 'Endste') then
      !Return if there is nothing to be done
      aux_logic = 0
!      if (a%kfl_dispa /= 0) aux_logic = aux_logic + 1
      if (a%kfl_trasg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      if (a%kfl_tacsg /= 0 .and. a%kfl_nolsg == 0) aux_logic = aux_logic + 1
      
      if (aux_logic == 0) return
  
   end if
    
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
  
   call ProcHook%PreLoop
   
   !We force closed rule for smoothing
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'lmn_EndElmope')
   !AllocateArrays in BaseElmope
   call AllocateBaseElmopeArrays

   !Hook
   call ProcHook%Initializations
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      !Elmats to Zero
      call ProcHook%ElmatsToZero
     
      call ElementGathers 
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
         
         call ElementInterpolates 
         call ProcPointer%ComputeDensity
         call vecnor(gpden,1,acden,2)

         !Hook
         call ProcHook%Interpolates
         
         !Compute Elext, Temporal Derivatives
         elext_mom=0.0_rp
         elext_ene=0.0_rp

         !Compute vector of external forces
         call ProcPointer%ExternalForces   

         !Time integration
         call ProcPointer%TimeIntegrationToElext

         !Default is just one iteration
         kfl_GoIteInGauss = 1

         do while (kfl_GoIteInGauss > 0)
         
            !Advection velocity      
            call ProcPointer%ComputeAdvectionVelocity
        
            !Temperature
            call ProcPointer%ComputeTemperature
    
            !Advection velocity norm
            call vecnor(gpadv,e%ndime,gpvno,2)
            
            !Compute the stability parameters, only if required
            call ProcHook%ComputeTaus
            
            !InGaussElmats
            !Hook
            call ProcHook%InGaussElmats

            !segregated 
            call ProcPointer%ComputeAdvectionVelocity
            call vecnor(gpadv,e%ndime,gpvno,2)
            call ProcHook%ComputeTaus_seg
            call ProcHook%InGaussElmats_seg

            kfl_GoIteInGauss = kfl_GoIteInGauss - 1
 
         enddo
         
         !InGaussElmats Assembly
         call ProcHook%ComputeResidual
         call ProcHook%InGaussElmatsAssembly

      enddo gauss_points
      
      !Assembly Endite
      !Hook
      call ProcHook%AssemblyEndite
       
   enddo elements
  
   !Hook
   call ProcHook%Finalizations
   
   !Deallocate arrays in BaseElmope
   call DeallocateBaseElmopeArrays
   
   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook%PostLoop

   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,a%EndLoopQuadrature,'lmn_EndElmope')

   
end subroutine

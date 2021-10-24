module Mod_nsf_elmope_1stNew
   use typre
   use Mod_nsFractionalStep
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_NavierStokesElement
   use Mod_NSF_Element
   use Mod_ConvectiveElement
   use Mod_nsm_Elmdir
   use Mod_php_Elmdir
   use Mod_php_SetTimeIntegrator
   
   !For setting pointers
   use Mod_nsm_BaseElmope
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
   use Mod_nsm_ComputeDissipation
   use Mod_nsf_DynSubsElmope
   use Mod_nsm_SwitchOff
   use Mod_nsm_LevelSetCoupling   
   use Mod_nsm_FreeSurfaceMatrices 
   use Mod_nsm_HangingNodes
   use Mod_nsm_FreeSurfaceStabilization
   use Mod_nsm_OrthogonalSubscales
   use Mod_nsm_ElasticBoundaryDir
   use Mod_nsm_ExtraViscosity
   use Mod_nsm_PorousMedia
   implicit none   
   
   
   class(NSFractionalStepProblem), pointer :: b => NULL()
   integer(ip) :: inode,ipoin
 
contains


   !--------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
      
      !------------------------------------------------------------
      !Defaults
      
      !Pointers
      ProcPointer_PostGaussElmats => PostGaussElmats
      ProcPointer_nsm_elmbuv  => nsm_elmbuv
      ProcPointer_nsm_elmrhu  => nsm_elmrhu
      
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersNonLinearDerivatives%Initialize
      call SetPointersAdvectionVelocity%Initialize 
      call SetPointersComputeTaus%Initialize
      call SetPointersComputeTestf%Initialize
      call SetPointersPhysicalProperties%Initialize 
      call SetPointersTurbulenceModel%Initialize 
      call SetPointersTemperatureCoupling%Initialize 
      call SetPointersDynSubsElmope1rst(0)
      call SetPointersInterpolateResidualProjection%Initialize
      call SetPointersExternalForces%Initialize
      call SetPointersInterpolateGradients%Initialize
      call SetPointersSwitchOff%Initialize
      call SetPointersLevelSetCoupling%Initialize
      call SetPointersFreeSurfaceMatrices%Initialize 
      call SetPointersHangingNodes%Initialize
      call SetPointersFSS%Initialize
      call SetPointersElasticBoundaryDir%Initialize
      call SetPointersPorousMedia%Initialize

      !Now we set the pointers
      !------------------------------------------------------------------
      call SetPointersAdvectionVelocity%Set
      call SetPointersComputeTaus%Set
      call SetPointersComputeTestf%Set
      call SetPointersTurbulenceModel%Set
      call SetPointersExternalForces%Set
      call SetPointersTemperatureCoupling%Set
      call SetPointersLevelSetCoupling%Set
      call SetPointersFreeSurfaceMatrices%Set 
      call SetPointersPhysicalProperties%SetTask('Elmope') 
      call SetPointersPhysicalProperties%Set
      call SetPointersFSS%Set
      call SetPointersElasticBoundaryDir%Set
      call SetPointersPorousMedia%Set
      call SetPointersExtraViscosity%Initialize
      
      if(a%kfl_colev==1)then
         !kfl_repro_SkipFE==2 skip only the temporal term
         if (a%kfl_repro_SkipFE == 1) a%kfl_repro_SkipFE = 2
      endif 
      
      if (a%kfl_repro == 1) then
         call SetPointersInterpolateResidualProjection%Set
         ! Hook for including RHS terms involving pressure
         call ConcatenateProcedures(ProcHook_InGaussElmats, InGaussElmrhPressureTerms)
         if (a%kfl_repro_SkipFE == 1) then
            ProcPointer_nsm_elmbuv => nsm_elmbuv_trm
            ProcPointer_nsm_elmrhu => nsm_elmrhu_trm
         elseif (a%kfl_repro_SkipFE == 2) then
            ProcPointer_nsm_elmbuv => nsm_elmbuv_trm
            ProcPointer_nsm_elmrhu => nsm_elmrhu_skipt
         end if
      
      !Include split OSS 
      elseif (a%kfl_repro == 2) then
         call SetPointersInterpolateResidualProjection%Set
         call ConcatenateProcedures(ProcHook_InGaussElmats, InGaussElmrhPressureTermsSplit)
         ! For block U,V and its rhs
         ProcPointer_nsm_elmbuv => nsm_elmbuv_split
         ProcPointer_nsm_elmrhu => nsm_elmrhu_split
      else
         call runend('ASGS not possible in fractional step')
      endif
      
     !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         if ( a%fvins > zensi ) then
            !call runend('Only Laplacian form ready for nonlinear elements')
         endif
         call SetPointersNonLinearDerivatives%Set
         call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsNonLinear)
      endif
      if (kfl_nonlinear == 1 .or. a%kfl_colev==1) then
         call ConcatenateProcedures(ProcHook_InGauss,InGaussVolumesNonLinear)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer_PostGaussElmats)
         ProcPointer_PostGaussElmats => NULLSUB
      endif       
              
      call SetPointersExtraViscosity%Set
      
      !Dynamic subgrid scales
      if (a%kfl_tacsg /= 0) call SetPointersDynSubsElmope1rst(1)
         
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange

      !-----------------------------------------------------------------------
      !Disconnection of elements
      if (a%kfl_SwitchOff == 1) call SetPointersSwitchOff%Set
      
      !HangingNodes
      call SetPointersHangingNodes%Set
         
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
      
      call SetPointersNonLinearDerivatives%Finalize
      call SetPointersAdvectionVelocity%Finalize 
      call SetPointersComputeTaus%Finalize
      call SetPointersComputeTestf%Finalize
      call SetPointersPhysicalProperties%Finalize 
      call SetPointersTurbulenceModel%Finalize 
      call SetPointersExternalForces%Finalize
      call SetPointersTemperatureCoupling%Finalize 
      call SetPointersDynSubsElmope1rst(100)
      call SetPointersInterpolateResidualProjection%Finalize
      call SetPointersInterpolateGradients%Finalize
      call SetPointersSwitchOff%Finalize
      call SetPointersLevelSetCoupling%Finalize 
      call SetPointersFreeSurfaceMatrices%Finalize 
      call SetPointersHangingNodes%Finalize
      call SetPointersFSS%Finalize
      call SetPointersElasticBoundaryDir%Finalize
      call SetPointersPorousMedia%Finalize
      call SetPointersExtraViscosity%Finalize
   end subroutine

   !---------------------------------------------------------------------
   !PostGauss Elmats
   subroutine PostGaussElmats
      implicit none
      
      ! Viscosity terms : we only consider mu*(grad v, grad u)
      call elmvis(e,dvolt0,acvis,wrmat1)

      ! If you want the complete term div(e·grad u) uncomment next line 
      if ( a%fvins > zensi ) then
            call nsm_elmvis_div(e,dvolt0,acvis,elmat)
      endif
      
      ! tau2*(div v, div u)
      call nsm_elmdiv(e,dvolt2,elmat)
   end subroutine
   
   
   !InGaussElmats
 
   
   !---------------------------------------------------------
   !Pressure terms in RHS
   subroutine InGaussElmrhPressureTerms
      implicit none
      call nsm_elmrhuf(e,dvol,testf,grpre,gppre(1),elext,eltemp,elrhs)
   end subroutine
   
   subroutine InGaussElmrhPressureTermsSplit
      implicit none
      call nsm_elmrhuf_split(e,dvol,gppre(1),elrhs)
   end subroutine

   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none

      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
   end subroutine

   subroutine InGaussElmatsNonLinear
      implicit none
      
      call nsm_elmbuv_lap(e,dvol,acvis,testf,wrmat1)
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module

subroutine nsf_elmope_1stNew(NSProblem)
   use Mod_nsf_elmope_1stNew
   use Mod_Int2str
   implicit none
   class(NSFractionalStepProblem), target :: NSProblem
   character(150) :: names
   
   a=>NSProblem%NavierStokesProblem
   b=>NSProblem
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
      
   !Allocations
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsf_elmope_1st')
   
   call AllocateBaseElmopeArrays
   
   !Matrix Allocations
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsf_elmope_1st')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','nsf_elmope_1st')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_1st')
   
   !Statistics
   call a%InitStats
   
   !Hook
   call ProcHook_Initializations

   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook_OnIeltyChange
      
      !Initializations
      elmat=0.0_rp
      elrhs=0.0_rp
      wrmat1 = 0.0_rp
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      !Gathering operations
      call VelocityAndPressureGathers
      !Hook
      call ProcHook_Gathers
      
      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)
     
      !Hook
      call ProcHook_PreGauss      
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp

      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook_InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call InterpolateGpVelocities
         call e%interpg(1,elpre(:,1),gppre)
         call e%gradient(1,elpre(:,1),grpre)
         !Hook
         call ProcHook_Interpolates
         
         !Physical Properties
         call a%GetPhysicalParameters(imat,acden,acvis)
         !Hook
         call ProcHook_PhysicalProp           
      
         !Advection velocity      
         call ProcPointer_ComputeAdvectionVelocity
         
         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
      
         !Compute a·grad(V)
         call ComputeAGradV(e,gpadv,AGradV)
         
         !Compute the stability parameters 
         !Pointer
         call ProcHook_ComputeTaus
         
         !Adjoint Test Function
         call ProcHook_ComputeTestf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
         dvolt2=dvol*tidiv + dvolt2                ! w(gp)*detjm*tau2
                  
         !Compute Elext, Temporal Derivatives, Repro...
         elext=0.0_rp
         !Compute vector of external forces
         call ProcPointer_ExternalForces        

         !Time integration
         eltemp=0.0_rp
         call nsi_TimeIntegrationToEltemp(e,Integrator,acden,a%dtinv,gpvel,eltemp)

         !Compute contributions to elemental matrix : Block U,V
         call ProcPointer_nsm_elmbuv(e,dvol,acden,LHSdtinv,aGradV,testf,wrmat1)
         
         !Compute contribution to the RHS, force and temporal terms
         call ProcPointer_nsm_elmrhu(e,dvol,testf,elext,eltemp,elrhs)
         
         !Residual Projection contributions to RHS : Block U
         call nsm_elmrhu_oss(e,tidiv,dvol,testf,gprep,elrhs)
         
         !Hook
         call ProcHook_InGaussElmats ! Terms involving pressure are inside this hook now
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer_PostGaussElmats
      
      !Matrix composition
      ! Assembly wrmat1  in elmuv
      forall (idime = 1:e%ndime)
         elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      end forall   
      
      !PreDirichlet
      call ProcHook_PreDirichlet
      
      !Dirichlet Boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      call php_elmdir(a,e,e%ndime,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
      
      do idime=1,e%ndime
         do inode=1,e%pnode
         if(abs(elmat(idime,inode,idime,inode)) <1e-14 .and. abs(elmat(idime,inode,idime,inode)) > 0.0_rp) then
            ipoin=e%lnods(inode)
            write(*,*) ielem,inode,ipoin,'1st','pivot_zero'
         endif
         end do
      end do
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   !Hook
   call ProcHook_Finalizations
   
   !Finalize Stats
   call a%FinalizeStats
   
   !Matrix Deallocations
   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsf_elmope_1st')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhs,'elrhs','nsf_elmope_1st')
   call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_1st')
   
   !Arrays Deallocations
   call DeallocateBaseElmopeArrays
   
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsf_elmope_1st')
   
   !if (a%istep == 2) then
   !   deb_PostprocessMatrix = 1
   !endif
   
end subroutine

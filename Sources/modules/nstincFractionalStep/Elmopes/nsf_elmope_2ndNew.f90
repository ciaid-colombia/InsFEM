module Mod_nsf_elmope_2ndnew
   use Mod_Mesh
   use Mod_Memor
   use Mod_nsFractionalStep
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_nsm_elmdir
   use Mod_NSF_Element
   
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
   use Mod_nsm_FreeSurfaceMatrices2nd   
   use Mod_nsm_HangingNodes
   use Mod_nsf_FreeSurfaceStabilization_elmope2nd
   implicit none   
   
   class(nsFractionalStepProblem), pointer :: b => NULL()
   
   real(rp) :: gplapu(3)      !Velocity laplacian
   real(rp), allocatable :: elveloc(:,:)   
   
   real(rp) :: divvel2
   
   logical :: IsCrankNicolson
   
   character(5) :: kfl_tsche_1st_current,kfl_tsche_2nd_current

   integer(ip) :: istep,elemStatus,inode,ipoin
   real(rp), allocatable :: divelvel(:,:)   
   
   
contains

   
   !---------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
   
      !External Procedures
      procedure() :: NULLSUB
      
      integer(ip) :: kfl_nonlinear, nelty
      
      call ResetProcedureComposition
      
      !------------------------------------------------
      !Defaults
      !Pointers
      ProcPointer_PostGaussElmats => PostGaussElmats
      
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB

      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersNonLinearDerivatives%Initialize
      call SetPointersExternalForces%Initialize 
      call SetPointersAdvectionVelocity%Initialize
      call SetPointersComputeTaus%Initialize
      call SetPointersPhysicalProperties%Initialize 
      call SetPointersTurbulenceModel%Initialize 
      call SetPointersDynSubsElmope2nd(0)
      call SetPointersInterpolateResidualProjection%Initialize
      call SetPointersInterpolateGradients%Initialize
      call SetPointersSwitchOff%Initialize
      call SetPointersLevelSetCoupling%Initialize
      call SetPointersFreeSurfaceMatrices2nd(0)     
      call SetPointersHangingNodes%Initialize
      call SetPointersFSS_Elmope2nd%Initialize
      
      !Now we set the pointers
      !------------------------------------------------------------------
      call SetPointersAdvectionVelocity%Set
      call SetPointersExternalForces%Set 
      call SetPointersComputeTaus%Set
      call SetPointersTurbulenceModel%Set 
      !We always need the gradients here
      call SetPointersInterpolateGradients%Set
      
      
      if (a%kfl_repro == 1) then
         call SetPointersInterpolateResidualProjection%Set
         call ConcatenateProcedures(ProcHook_InGaussElmats, InGaussElmrh2ndStep)
      elseif (a%kfl_repro == 2) then
         call SetPointersInterpolateResidualProjection%Set
         call ConcatenateProcedures(ProcHook_InGaussElmats, InGaussElmrh2ndStepSplit)
      else
         call runend('nsf_elmope2nd_onlyready for residual projection')
      endif
      
      !Repro Skip_FE type
      !-----------------------------------------------------------
      if(a%kfl_repro_SkipFE == 1)then
      
      elseif(a%kfl_repro_SkipFE == 2)then
         !The force is already in the residual projection, put it also in the regular residual
         call ConcatenateProcedures(ProcHook_InGaussElmats,ForceToelrhp)
      end if  
      !-----------------------------------------------------------    
      
      
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         if ( a%fvins > zensi ) then
            !call runend('Only Laplacian form ready for nonlinear elements')
         endif
         call SetPointersNonLinearDerivatives%Set
         call ConcatenateProcedures(ProcHook_Interpolates,InterpolateLapU)
         call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsNonLinear)
      endif
      
      if (kfl_nonlinear == 1 .or. a%kfl_colev==1) then
         call ConcatenateProcedures(ProcHook_InGauss,InGaussVolumesNonLinear)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer_PostGaussElmats)
         ProcPointer_PostGaussElmats => NULLSUB
      endif
      
      call SetPointersLevelSetCoupling%Set
      call SetPointersFreeSurfaceMatrices2nd(1)      
      call SetPointersPhysicalProperties%SetTask('Elmope') 
      call SetPointersPhysicalProperties%Set
      call SetPointersFSS_Elmope2nd%Set
      
      !Dynamic or non-linear subgrid scales
      if (a%kfl_tacsg /= 0) call SetPointersDynSubsElmope2nd(1)
      
      if (a%kfl_penal == 1) call ConcatenateProcedures(ProcHook_InGaussElmats,PenaltyPressure)
      
      !Hanging Nodes
      call SetPointersHangingNodes%Set
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange

      !-----------------------------------------------------------------------
      !Disconnection of elements
      if (a%kfl_SwitchOff == 1) call SetPointersSwitchOff%Set

      !-------------------------------------------------------------------
      !Finalize procedure flags so that we can use them again
      call SetPointersNonLinearDerivatives%Finalize
      call SetPointersAdvectionVelocity%Finalize
      call SetPointersExternalForces%Finalize 
      call SetPointersComputeTaus%Finalize
      call SetPointersPhysicalProperties%Finalize
      call SetPointersTurbulenceModel%Finalize 
      call SetPointersDynSubsElmope2nd(100)
      call SetPointersInterpolateResidualProjection%Finalize
      call SetPointersInterpolateGradients%Finalize
      call SetPointersSwitchOff%Finalize
      call SetPointersLevelSetCoupling%Finalize 
      call SetPointersFreeSurfaceMatrices2nd(100)      
      call SetPointersHangingNodes%Finalize
      call SetPointersFSS_Elmope2nd%Finalize
      
   end subroutine
   
   !----------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none
      
      ! Viscosity terms : we only consider mu*(grad q, grad p)
      call elmvis(e,dvolt0,1.0_rp,elmat)
   end subroutine
     
   !-------------------------------------------------------
   ! InGauss RHS contribution
   subroutine InGaussElmrh2ndStep
      implicit none
      call nsm_elmrhf_pre(e,dvol,acden,timom,LHSDtinv,divvel2,gprep,grpre,grvel,gpadv,elrhs)
   end subroutine
   
   subroutine InGaussElmrh2ndStepSplit
      implicit none
      call nsm_elmrhf_pre_split(e,dvol,acden,timom,LHSdtinv,divvel2,gprep(e%ndime+2),grpre,elrhs)
   end subroutine
   
   !----------------------------------------------------------
   !NonLinear Elements
   subroutine InGaussVolumesNonLinear
      implicit none
      !DvolsToZero
      dvolt0=0.0_rp
   end subroutine
   
   subroutine InterpolateLapU
      implicit none
      real(rp) :: aux1(e%pnode)
      integer(ip) :: inode
      integer(ip) :: idime
      
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo  
      do idime = 1,e%ndime                              ! Contribution from the laplacian term
         gplapu(idime) = dot_product(aux1(1:e%pnode),elvel(idime,1:e%pnode,1))
      end do
   end subroutine
   
   subroutine InGaussElmatsNonLinear
      implicit none
      
      call nsf_elm2nd_reslap(e,dvol,timom,acvis,gplapu,elrhs)
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
    end subroutine
   
    !For Skip_FE 2, the force term is put in the residual and its projection
    subroutine ForceToelrhp
       implicit none
       real(rp) :: eltemp(e%ndime)
       
       !Compute Elext
       elext=0.0_rp
       !Compute vector of external forces without rhs-temporal  
       call ProcPointer_ExternalForces      
       
       eltemp = 0.0_rp
       call nsm_elmrhp(e,timom,dvol,elext,eltemp,elrhs)
    end subroutine
   
   !PenaltyTerm
   subroutine PenaltyPressure
      implicit none
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elmat(1,inode,1,inode) = elmat(1,inode,1,inode) + a%penal*dvol
      enddo
   end subroutine
  
end module

subroutine nsf_elmope_2ndNew(NSProblem)
   use typre
   use Mod_nsf_elmope_2ndnew
   use Mod_nsFractionalStep
   use Mod_NSF_Element
   implicit none
   class(NSFractionalStepProblem), target :: NSProblem
   
   !real(rp) :: divvel2
   
   a => NSProblem%NavierStokesProblem
   b => NSProblem

   !Set Pointers for Execution
   ielem = 1
   call SetPointers
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsm_elmope')

   !Matrix Allocate
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmat,'elmat','nsf_elmope_2nd')
   call a%Memor%alloc(1,e%mnode,elrhs,'elrhs','nsf_elmope_2nd')
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','nsf_elmope_2nd')
   call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','nsf_elmope_2nd')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_2nd')
   call a%Memor%alloc(e%ndime,e%mnode,divelvel,'divelvel','nsf_elmope2nd')
   
   !Arrays Allocate
   call AllocateBaseElmopeArrays

   !Initializations
   call ProcHook_Initializations
   
   if (a%kfl_tsche_1st_current == 'CN   ') then
      IsCrankNicolson = .true.
      
      !The second step of the fractional step is solved as BDF1
      kfl_tsche_1st_current = 'BDF1 '
      call Integrator%Init(kfl_tsche_1st_current)
      call Integrator%GetLHSDtinv(a%dtinv,LHSdtinv)
      call Integrator%GetNumberOfTimeSteps(nsteps)
   else
      IsCrankNicolson = .false.
   endif
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)
     
      !Hook
      call ProcHook_OnIeltyChange
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      !Element length at center of gravity
      call e%elmlen
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp    

      !Gathers
      if (IsCrankNicolson) then
        call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,2)) !int.u_n+1/2
      else
        call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1)) !int.u_n+1
      end if
      call e%gather(1_ip,elpre(1:e%pnode,1),a%press(1:e%pnode,1))         ! p_n
      do istep=2,nsteps !bdf2 and others
         call e%gather(1_ip,elpre(:,istep),a%press(:,istep+1))
      enddo      
      call e%gather(e%ndime,divelvel,a%veloc(:,:,1))      
      !Hook
      call ProcHook_Gathers

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)
      
      !DvolsToZero
      dvolt0=0.0_rp

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

         call e%divergence(elvel(:,:,1),divvel)
         call e%divergence(divelvel,divvel2)
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
         !write(*,*) 'nsf_elmope_2ndNew: ','gpadv=',gpadv
         !Compute the stability parameters 
         !Hook
         call ProcHook_ComputeTaus
         
         !Compute Dvols
         dvolt0 = dvolt0 + (timom + 1.0_rp/(acden*LHSdtinv))*dvol
         
         !Hook
         call ProcHook_InGaussElmats
         
      enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer_PostGaussElmats

      !PreDirichlet
      call ProcHook_PreDirichlet
      
      
      !Dirichlet Boundary conditions
      call nsf_elmdir_press(b,e,elmat,elrhs)
      
      do inode=1,e%pnode
         if(abs(elmat(1,inode,1,inode)) <1e-14 .and. abs(elmat(1,inode,1,inode)) > 0.0_rp) then
            ipoin=e%lnods(inode)
            write(*,*) ielem,inode,ipoin,'2nd','pivot_zero'
         endif
      end do
      
      !Assembly
      call b%LinearSystemP%Assembly(e,elmat,elrhs)
   
   enddo elements
   
   !Initializations
   call ProcHook_Finalizations
   
   !Matrix Deallocate
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmat,'elmat','nsf_elmope_2nd')
   call a%Memor%dealloc(1,e%mnode,elrhs,'elrhs','nsf_elmope_2nd')
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','nsf_elmope_2nd')
   call a%Memor%dealloc(1,e%mnode,elrhp,'elrhp','nsf_elmope_2nd')
   call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_2nd')
   call a%Memor%dealloc(e%ndime,e%mnode,divelvel,'divelvel','nsf_elmope2nd')
   
   !Arrays DeAllocate
   call DeallocateBaseElmopeArrays

   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsm_elmope')
end subroutine

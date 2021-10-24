module Mod_nsf_HydrostaticPressureNew
   use typre
   use Mod_NSFractionalStep
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_NSF_Element
   use Mod_ConvectiveElement
   
   !For setting pointers
   use Mod_nsm_BaseElmope
   use Mod_nsm_TemperatureCoupling
   use Mod_nsm_NonLinearDerivatives
   use Mod_nsm_ExternalForces
   use Mod_nsm_LevelSetCoupling   
   use Mod_nsm_FreeSurfaceMatricesHydro 
   use Mod_nsm_HangingNodes
   
   implicit none
   
   class(NSFractionalStepProblem), pointer :: b => NULL()

contains
   
   !--------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      integer(ip) :: kfl_nonlinear, nelty

      
      !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
   
      !------------------------------------------------------------
      !Defaults
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !Pointers
      ProcPointer_PostGaussElmats => PostGaussElmats
      
      !-------------------------------------------------------------------
      !Initialize the procedure flags
      call SetPointersNonLinearDerivatives%Initialize
      call SetPointersExternalForces%Initialize
      call SetPointersTemperatureCoupling%Initialize
      call SetPointersLevelSetCoupling%Initialize 
      call SetPointersFreeSurfaceMatricesHydro(0)      
      call SetPointersHangingNodes%Initialize
   
      !Now we set the pointers
      call SetPointersExternalForces%Set
      call SetPointersTemperatureCoupling%Set
      call SetPointersLevelSetCoupling%Set
      call SetPointersFreeSurfaceMatricesHydro(1)  
      call SetPointersHangingNodes%Set
      
   
      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) call SetPointersNonLinearDerivatives%Set
      if (kfl_nonlinear == 1 .or. a%kfl_colev==1) then
         call ConcatenateProcedures(ProcHook_InGauss,InGaussVolumesNonLinear)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer_PostGaussElmats)
         ProcPointer_PostGaussElmats => NULLSUB
      endif       
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
      
      call SetPointersNonLinearDerivatives%Finalize
      call SetPointersExternalForces%Finalize
      call SetPointersTemperatureCoupling%Finalize 
      call SetPointersLevelSetCoupling%Finalize 
      call SetPointersFreeSurfaceMatricesHydro(100)    
      call SetPointersHangingNodes%Finalize
   
   end subroutine

   !---------------------------------------------------------------------
   !PostGauss Elmats
   subroutine PostGaussElmats
      implicit none
      
      !Viscosity terms : we only consider mu*(grad v, grad u)
      call elmvis(e,dvolt0,1.0_rp,elmat)
   end subroutine

   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none

      !DvolsToZero
      dvolt0=0.0_rp
   end subroutine   

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine


end module



subroutine nsf_HydrostaticPressureNew(NSFProblem)
   !This subroutine solves (grad q, grad p) = (grad q, f)
   !for Fractional Step in NavierStokes
   !This is important so that fractional step starts smoothly and does not 
   !go crazy
   use typre
   use Mod_NSFractionalStep
   use Mod_nsf_HydrostaticPressureNew
   implicit none
   class(NSFractionalStepProblem), target :: NSFProblem
   integer(ip) :: icomp, kfl_perio, kfl_HangingNodes 
   a => NSFProblem%NavierStokesProblem
   b => NSFProblem

   
   !--------------------------------------------------------
   !First we set the linear system to zero
   call b%LinearSystemP%ToZero
   
   
   
   !--------------------------------------------------------
   !Secondly we assembly the matrix
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers

   !Allocations
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsf_elmope_1st')
   
   !Matrix Allocations
   call a%Memor%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsf_hydro')
   call a%Memor%alloc(1_ip,e%mnode,elrhs,'elrhs','nsf_hydro')

   !Arrays Allocate
   call AllocateBaseElmopeArrays
   
   !Hook
   call ProcHook_Initializations
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook_OnIeltyChange

      !Hook
      call ProcHook_PreGauss
   
      !Initializations
      elmat=0.0_rp
      elrhs=0.0_rp
      
      !Hook
      call ProcHook_Gathers
      
      !Derivatives and Jacobian at the center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen
      
      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)      
      
      !DvolsToZero
      dvolt0=0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook_InGauss
         
         dvol = e%weigp(e%igaus)*e%detjm
         dvolt0=dvol + dvolt0                ! w(gp)*detjm
         
         !Hook
         call ProcHook_Interpolates
         
         !Physical Properties
         call a%GetPhysicalParameters(imat,acden,acvis)
         !Hook
         call ProcHook_PhysicalProp   
         
         !Compute g, Boussinesq, ...
         elext=0.0_rp
         !Compute vector of external forces
         call ProcPointer_ExternalForces   
         
         !Compute system rhs (grad q, grad f)
         call nsf_hydro_rhs(e,dvol,elext,elrhs)
         
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

      !Assembly
      call b%LinearSystemP%Assembly(e,elmat,elrhs)
   
   enddo elements
   
   !Deallocations
   
   !Hook
   call ProcHook_Finalizations
   
   !Matrix Allocations
   call a%Memor%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsf_hydro')
   call a%Memor%dealloc(1_ip,e%mnode,elrhs,'elrhs','nsf_hydro')
   
   !Arrays DeAllocate
   call DeAllocateBaseElmopeArrays
   
   !Element Deallocation
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsf_elmope_1st')
   
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(1_ip,b%LinearSystemP,a%Memor)
   
   !Hanging nodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(1_ip,b%LinearSystemP,a%Memor)
   
   !-----------------------------------------------------
   !Finally we solve the system
   if (a%MPIrank == a%MPIroot) write(a%lun_solve,101) 
   
   call b%LinearSystemP%Solve(b%unknoP)
   
   call a%Mesh%ArrayCommunicator%GhostCommunicate(1,b%unknoP)
   
   !HangingNodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(1_ip,b%unknoP)
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%MasterToSlave(1_ip,b%unknoP)
      
   do icomp = 1,size(a%press,2)
      a%press(:,icomp) = b%unknoP(1,:)
   enddo

   
   !Formats
   101 format('INITIAL HYDROSTATIC PRESSURE SOLVE: ')
   

   
end subroutine



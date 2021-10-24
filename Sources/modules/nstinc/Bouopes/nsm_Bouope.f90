module Mod_nsm_Bouope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope
   use Mod_nsm_NeumannBC
   !use Mod_nsm_ViscosityBouope
   use Mod_nsm_PhysicalProperties
   use Mod_nsm_DomainDecomposition
   use Mod_nsm_FSI
   use Mod_nsm_BoundarySGS
   implicit none

contains
   subroutine SetPointers
      implicit none
      call ResetProcedureComposition
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      call SetPointersNeumannBC%Initialize
      !call SetPointersViscosity%Initialize
      call SetPointersPhysicalProperties%Initialize 
      call SetPointersDomainD%Initialize
      call SetPointersFSI%Initialize
      call SetPointersBoundarySGS%Initialize
      call SetPointersNeumannBC%Set
      !call SetPointersViscosity%Set
      call SetPointersPhysicalProperties%SetTask('Bouope')
      call SetPointersPhysicalProperties%Set
      call SetPointersDomainD%Set
      call SetPointersFSI%Set
      call SetPointersBoundarySGS%Set
      call SetPointersNeumannBC%Finalize
      !call SetPointersViscosity%Finalize
      call SetPointersPhysicalProperties%Finalize 
      call SetPointersDomainD%Finalize
      call SetPointersFSI%Finalize
      call SetPointersBoundarySGS%Finalize
   end subroutine
end module

subroutine nsm_Bouope(NSProblem)
!-----------------------------------------------------------------------
!****f* Nstinc/nsm_bouope
! NAME 
!    nsm_bouope
! DESCRIPTION
!    Navier-Stokes boundary elemental operations:
!    1. Compute elemental matrix and RHS 
!    2. Impose Dirichlet boundary conditions
!    3. Assemble them
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   use Mod_nsm_elmdir
   use Mod_nsm_HangingNodes
   use Mod_nsm_Bouope
   implicit none
   class(NavierStokesProblem), target :: NSProblem
   integer(ip) :: kfl_HangingNodes
      integer(ip) :: inodb,inode
   
   !SetPointers for NavierStokesProblem
   a=>NSProblem
   
   !Pointers
   call SetPointers

   !Initializations
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetHanging(kfl_HangingNodes)
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   call AllocateBaseBouopeMatrices
   call AllocateBaseBouopeArrays
   call ProcHook_Initializations
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      call a%Mesh%BoundaryLoad(iboun,e) ! Load Element

      call BouopeElMatsToZero 
      call ProcHook_ElmatsToZero

      call BoundaryGathers
      call ProcHook_Gathers
      
      !Cartesian derivatives
      call e%elmdel
      !Element length
      call e%elmlen
      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,bovel,chale,a%kfl_hdifumin)
      dsurf0 = 0.0_rp

      call ProcHook_PreGauss
      
      !Gauss-Point Loop
      GaussPointLoop: do igaub=1,e%pgaub
         e%igaub = igaub
      
         !Derivatives at the boundary
         call e%elmderb
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         wmatr = 0.0_rp
         tract = 0.0_rp
        
         call BoundaryInterpolates 
         call ProcHook_Interpolates

         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)
         call ProcHook_PhysicalProp

         call ProcHook_InGauss

         !Velocity norm
         call vecnor(gpbve(:,1),ndime,gpvno,2)
         call ProcHook_InGaussElmats
         call ProcHook_InGaussElmatsAssembly
      end do GaussPointLoop
      call ProcHook_PreDirichlet

      if (kfl_HangingNodes == 1) call ModifyMatricesHanging
      
      !Assembly
      call nsm_elmdir(a,e,elmat,elrhs)
      call ProcHook_Assembly
   end do boundaries

   call ProcHook_Finalizations
   
   call DeallocateBaseBouopeArrays
   call DeallocateBaseBouopeMatrices
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
end subroutine nsm_bouope

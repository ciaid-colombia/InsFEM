module Mod_nsf_FreeSurfaceStabilization_elmope2nd
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeAdvectionVelocity
   use typre
   implicit none
   private
   public SetPointersFSS_Elmope2nd
   
   !For Free surface stabilization
   logical :: isFirst2Layers
   
   real(rp), allocatable :: projected_elGrpreAndForce(:,:), projected_GrPreAndForce(:)
   real(rp), allocatable :: elgrvel(:,:,:), elgrpre(:,:),elmass(:)
   real(rp), allocatable :: vmass(:)
   
   !For elmopes
   type, extends(PointerSetter) :: SPFreeSurfaceStabilization
contains
      procedure :: SpecificSet => SpecificSetFSS_Elmope2nd
   end type
   type(SPFreeSurfaceStabilization) :: SetPointersFSS_Elmope2nd
   
  

contains
   subroutine SpecificSetFSS_Elmope2nd(d)
      implicit none
      class(SPFreeSurfaceStabilization) :: d
      
      if (a%kfl_StabilizeFreeSurface == 1) then
         call ConcatenateProcedures(ProcHook_Initializations,Allocations)
         call ConcatenateProcedures(ProcHook_Finalizations,Deallocations)
         !Stabilization matrices are computed at the center of gravity
         !We prepend it to make sure it is done before the gauss points are changed by the cut elements
         call PrependProcedure(ProcHook_PreGauss,StabilizationMatrices_nsf2nd)
         call PrependProcedure(ProcHook_PreGauss,Gathers)
         call PrependProcedure(ProcHook_PreGauss,CheckIsFirst2Layers)
      endif
   end subroutine
   
   
   

   !------------------------------------------------------------------------
   !Free surface stabilization
   subroutine Allocations
      implicit none
      call a%Mesh%GetNdime(ndime)
      
!       call a%Memor%alloc(ndime,ndime,e%mnode,projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')
      call a%Memor%alloc(ndime,e%mnode,projected_elGrpreAndForce,'projected_elGrpreAndForce','alloc_freeSurfStabilization')
!       
!       call a%Memor%alloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')
      call a%Memor%alloc(ndime,projected_GrPreAndForce,'projected_GrPreAndForce','alloc_freeSurfStabilization')
      
      
      
   end subroutine
   
   subroutine DeAllocations
      implicit none
!       call a%Memor%dealloc(ndime,ndime,e%mnode,projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,size(projected_elGrpreAndForce,2),projected_elGrpreAndForce,'projected_elGrpreAndForce','alloc_freeSurfStabilization')
!       call a%Memor%dealloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,projected_GrPreAndForce,'projected_GrPreAndForce','alloc_freeSurfStabilization')
      
   end subroutine
   
   !call SetPointersInterpolateGradients(1)
   
   subroutine CheckIsFirst2Layers
      implicit none
      integer(ip) :: inode,ipoin,poinStatus
   
      isFirst2Layers = .false.
      
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         
         call a%CutMesh%GetPointType(ipoin,poinStatus) 
         if (poinStatus == 1) then
            !This means the element belongs to either the cut layer or the adjacent one inside
            !I want to add stabilization in these layers
            
            isFirst2Layers = .true.
            exit
         endif
      enddo
      
   end subroutine
   
   subroutine Gathers  
      implicit none
      if (isFirst2Layers) then
         call e%gather(ndime,projected_elGrpreAndForce,a%StabilizeFreeSurface_PressureGradientAndForce)
      endif
   
   end subroutine
   
   
   subroutine StabilizationMatrices_nsf2nd
      implicit none
      real(rp) :: velgraddiff(ndime,ndime), pregraddiff(ndime)
      integer(ip) :: idime,inode,jnode
      real(rp)  :: ghosttimom
      
      
      if (isFirst2Layers) then
      
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
            !Hook
            call ProcHook_PhysicalProp    
            
            !Advection velocity      
            call ProcPointer_ComputeAdvectionVelocity
            gpadv = gpvel(:,1)
            
            !Advection velocity norm
            call vecnor(gpadv,e%ndime,gpvno,2)
            
            !Compute a·grad(V)
            call ComputeAGradV(e,gpadv,AGradV)
            
            !Compute the stability parameters 
            call ProcHook_ComputeTaus
            
            call e%interpg(ndime,projected_elGrpreAndForce,projected_GrPreAndForce)
            
            !Compute Elext, Temporal Derivatives
            elext=0.0_rp
            !Compute vector of external forces
            call ProcPointer_ExternalForces  
         
            !first the right hand side
            do inode = 1,e%pnode
            
               !Stabilize
               elrhs(1,inode) = elrhs(1,inode) + a%StabilizeFreeSurface_Param*timom*dot_product(e%cartd(:,inode),projected_GrPreAndForce+elext)*e%weigp(e%igaus)*e%detjm 
            enddo
            
            !Second the left hand side
            do inode = 1,e%pnode
               do jnode = 1,e%pnode
                  !Stabilize the velocity array
                  elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + a%StabilizeFreeSurface_Param*timom*dot_product(e%cartd(:,inode),e%cartd(:,jnode))*e%weigp(e%igaus)*e%detjm 
               enddo
            enddo
            
         enddo gauss_points
            
       
      endif
      
   
   end subroutine
   

   
end module

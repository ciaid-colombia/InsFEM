module Mod_nsm_FreeSurfaceStabilization
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeAdvectionVelocity
   implicit none
   private
   public SetPointersFSS, SetPointersFSS_Endite
   
   logical :: isFirst2Layers
   real(rp), allocatable :: projected_elGrpreAndForce(:,:), projected_elgrvel(:,:,:), projected_GrPreAndForce(:), projected_grvel(:,:)
   real(rp), allocatable :: elgrvel(:,:,:), elgrpre(:,:),elmass(:)
   real(rp), allocatable :: vmass(:)
   integer(ip) :: ielem2 
   
   !For elmopes
   type, extends(PointerSetter) :: SPFreeSurfaceStabilization
contains
      procedure :: SpecificSet => SpecificSetFSS
   end type
   type(SPFreeSurfaceStabilization) :: SetPointersFSS
   
   !For endite
   type, extends(PointerSetter) :: SPFSS_Endite
contains
      procedure :: SpecificSet => SpecificSetFSS_Endite
   end type
   
   type(SPFSS_Endite) :: SetPointersFSS_Endite

contains
   subroutine SpecificSetFSS(d)
      implicit none
      class(SPFreeSurfaceStabilization) :: d
      
      if (a%kfl_StabilizeFreeSurface == 1) then
         call ConcatenateProcedures(ProcHook_Initializations,Allocations)
         call ConcatenateProcedures(ProcHook_Finalizations,Deallocations)
         call PrependProcedure(ProcHook_PreGauss,StabilizationMatrices)
         call PrependProcedure(ProcHook_PreGauss,Gathers)
         call PrependProcedure(ProcHook_PreGauss,CheckIsFirst2Layers)
      endif
   end subroutine
   
   subroutine SpecificSetFSS_Endite(d)
      implicit none
      class(SPFSS_Endite) :: d
      
      if (a%kfl_StabilizeFreeSurface == 1) then
         
         !We will need the interpolated gradients
         call SetPointersInterpolateGradients%Set
         call ConcatenateProcedures(ProcHook_Initializations,Allocations_Endite)
         !We prepend so that we make sure it is done with the full gauss points
         call PrependProcedure(ProcHook_PreGauss,InGaussOperations_Endite)
         call PrependProcedure(ProcHook_PreGauss,ElmatsToZero_Endite)
         call PrependProcedure(ProcHook_PreGauss,CheckIsFirst2Layers)
         call ConcatenateProcedures(ProcHook_AssemblyEndite,Assembly_Endite)
         call ConcatenateProcedures(ProcHook_Finalizations,Finalizations_Endite)
         call ConcatenateProcedures(ProcHook_Finalizations,Deallocations_Endite)
      endif
   end subroutine

   !------------------------------------------------------------------------
   !Free surface stabilization
   subroutine Allocations
      implicit none

      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,ndime,e%mnode,projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')
      call a%Memor%alloc(ndime,e%mnode,projected_elGrpreAndForce,'projected_elGrpreAndForce','alloc_freeSurfStabilization')
      call a%Memor%alloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')
      call a%Memor%alloc(ndime,projected_GrPreAndForce,'projected_GrPreAndForce','alloc_freeSurfStabilization')
   end subroutine
   
   subroutine DeAllocations
      implicit none

      call a%Memor%dealloc(ndime,ndime,size(projected_elgrvel,3),projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,size(projected_elGrpreAndForce,2),projected_elGrpreAndForce,'projected_elGrpreAndForce','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,projected_GrPreAndForce,'projected_GrPreAndForce','alloc_freeSurfStabilization')
   end subroutine
   
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
         call e%gather(ndime*ndime,projected_elgrvel,a%StabilizeFreeSurface_VelocityGradient)
         call e%gather(ndime,projected_elGrpreAndForce,a%StabilizeFreeSurface_PressureGradientAndForce)
      endif
   end subroutine
   
   subroutine StabilizationMatrices
      implicit none
      real(rp) :: velgraddiff(ndime,ndime), pregraddiff(ndime)
      integer(ip) :: idime,inode,jnode
      real(rp)  :: ghosttimom
      logical :: dopressure
      
      if (isFirst2Layers) then
         dopressure = .false.
         if (size(elmat,1) > e%ndime) dopressure = .true.
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            dvol = e%weigp(e%igaus)*e%detjm

            !Interpolate         
            call InterpolateGpVelocities
            !Hook
            call ProcHook_Interpolates
            !Physical Properties
            call a%GetPhysicalParameters(imat,acden,acvis)
            !Hook
            call ProcHook_PhysicalProp    

            !Advection velocity      
            call ProcPointer_ComputeAdvectionVelocity
            if (a%kfl_advec == 1) gpadv = gpvel(:,1)
            if (a%kfl_advco == 1) gpadv(1:e%ndime) = a%advco(1:e%ndime)

            !Advection velocity norm
            call vecnor(gpadv,e%ndime,gpvno,2)
            !Compute a·grad(V)
            call ComputeAGradV(e,gpadv,AGradV)
            
            !Compute the stability parameters 
            call ProcHook_ComputeTaus
            call e%interpg(ndime*ndime,projected_elgrvel,projected_grvel)
            call e%interpg(ndime,projected_elGrpreAndForce,projected_GrPreAndForce)
            
            !Compute Elext
            elext=0.0_rp
            !Compute vector of external forces
            call ProcPointer_ExternalForces  
         
            !Stabilize the velocity array
            !first the right hand side
            do inode = 1,e%pnode
               !Stabilize velocity array
               do idime = 1,ndime
                  elrhs(idime,inode) = elrhs(idime,inode) + a%StabilizeFreeSurface_Param*tidiv*dot_product(e%cartd(:,inode),projected_grvel(:,idime))*e%weigp(e%igaus)*e%detjm 
               enddo
            enddo
            !Second the left hand side
            do inode = 1,e%pnode
               do jnode = 1,e%pnode
                  !Stabilize the velocity array
                  do idime = 1,ndime
                     elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + a%StabilizeFreeSurface_Param*tidiv*dot_product(e%cartd(:,inode),e%cartd(:,jnode))*e%weigp(e%igaus)*e%detjm 
                  enddo
               enddo
            enddo
            
            !Stabilize pressure array
            if (dopressure) then
               do inode = 1,e%pnode
                   elrhs(ndime+1,inode) = elrhs(ndime+1,inode) + a%StabilizeFreeSurface_Param*timom*dot_product(e%cartd(:,inode),(projected_GrPreAndForce(:)+elext(1:e%ndime)))*e%weigp(e%igaus)*e%detjm 
               enddo
               !Second the left hand side
               do inode = 1,e%pnode
                  do jnode = 1,e%pnode
                     elmat(ndime+1,inode,ndime+1,jnode) = elmat(ndime+1,inode,ndime+1,jnode) + a%StabilizeFreeSurface_Param*timom*dot_product(e%cartd(:,inode),e%cartd(:,jnode))*e%weigp(e%igaus)*e%detjm 
                  enddo
               enddo
            endif
         enddo gauss_points
      endif
   
   end subroutine
   
   !------------------------------------------------------------------
   !For endite
   subroutine Allocations_Endite
      implicit none
      integer(ip) :: npoin
      
      a%StabilizeFreeSurface_VelocityGradient = 0.0_rp
      a%StabilizeFreeSurface_PressureGradientAndForce = 0.0_rp
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,ndime,e%mnode,elgrvel,'elgrvel','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%alloc(ndime,e%mnode,elgrpre,'elgrpre','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%alloc(e%mnode,elmass,'elmass','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%alloc(npoin,vmass,'auxvmass','Mod_nsm_FreeSurfaceStabilization')

   end subroutine
   
   subroutine Deallocations_Endite
      implicit none
      integer(ip) :: npoin
      
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(ndime,ndime,e%mnode,elgrvel,'elgrvel','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%dealloc(ndime,e%mnode,elgrpre,'elgrpre','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%dealloc(e%mnode,elmass,'elmass','Mod_nsm_FreeSurfaceStabilization')
      call a%Memor%dealloc(npoin,vmass,'auxvmass','Mod_nsm_FreeSurfaceStabilization')
   
   end subroutine
   
   subroutine ElmatsToZero_Endite
      implicit none
      
      !isFirst2Layers = .true.
      ielem2 = ielem
      if (isFirst2Layers) then
         elgrvel = 0.0_rp
         elgrpre = 0.0_rp
         elmass = 0.0_rp
      endif

   end subroutine
   
   subroutine InGaussOperations_Endite
      implicit none
      integer(ip) :: inode
   
      if (isFirst2Layers) then
         
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            dvol = e%weigp(e%igaus)*e%detjm
            
            !Interpolate         
            call InterpolateGpVelocities
            !Hook
            call ProcHook_Interpolates
            !Compute Elext, Temporal Derivatives
            elext=0.0_rp
            !Compute vector of external forces
            call ProcPointer_ExternalForces  
            do inode = 1,e%pnode
               elgrvel(:,:,inode) = elgrvel(:,:,inode) + e%shape(inode,e%igaus)*grvel*dvol
               elgrpre(:,inode)   = elgrpre(:,inode)   + e%shape(inode,e%igaus)*(grpre(1,:)- elext(1:e%ndime))*dvol
               elmass(inode) = elmass(inode) + e%shape(inode,e%igaus)*dvol
            enddo
         enddo gauss_points
      endif

   end subroutine
   
   subroutine Assembly_Endite
      implicit none
      if (isFirst2Layers) then
         call a%Mesh%AssemblyToArray(e,ndime*ndime,elgrvel,a%StabilizeFreeSurface_VelocityGradient) 
         call a%Mesh%AssemblyToArray(e,ndime,elgrpre,a%StabilizeFreeSurface_PressureGradientAndForce) 
         call a%Mesh%AssemblyToArray(e,1_ip,elmass,vmass) 
      endif

   end subroutine
   
   subroutine Finalizations_Endite
      use Mod_int2str
      implicit none
      interface
         subroutine AfterProjectOperations(a,ndofn,array)
            use typre
            use Mod_Memor
            use Mod_Mesh
            implicit none
            class(FemMesh) :: a
            integer(ip) :: ndofn
            real(rp)    :: array(:,:)
         end subroutine
      end interface
      
      integer(ip) :: ipoin,npoin,idime
      real(rp), pointer :: PointerArray(:,:) => NULL()
      character(150) :: auxname
      
      auxname = 'GraPress'//trim(adjustl(int2str(a%itera)))
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         if (vmass(ipoin) /= 0.0_rp) then
            a%StabilizeFreeSurface_VelocityGradient(:,:,ipoin) = a%StabilizeFreeSurface_VelocityGradient(:,:,ipoin)/vmass(ipoin)
            a%StabilizeFreeSurface_PressureGradientAndForce(:,ipoin) = a%StabilizeFreeSurface_PressureGradientAndForce(:,ipoin)/vmass(ipoin)
         endif
      enddo
      
      PointerArray(1:ndime*ndime,1:npoin) => a%StabilizeFreeSurface_VelocityGradient
      call AfterProjectOperations(a%Mesh,ndime*ndime,PointerArray)
      call AfterProjectOperations(a%Mesh,ndime      ,a%StabilizeFreeSurface_PressureGradientAndForce)
      
   end subroutine
   
end module

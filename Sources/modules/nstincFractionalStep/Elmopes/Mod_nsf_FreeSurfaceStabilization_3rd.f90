module Mod_nsf_FreeSurfaceStabilization_3rd
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeAdvectionVelocity
   use typre
   implicit none
   private
   public SetPointersFSS_3rd
   
   !For Free surface stabilization
   logical :: isFirst2Layers
   
   real(rp), allocatable :: projected_elgrvel(:,:,:),  projected_grvel(:,:)
   real(rp), allocatable :: elgrvel(:,:,:), elgrpre(:,:),elmass(:)
   real(rp), allocatable :: vmass(:)
   
   !For elmopes
   type, extends(PointerSetter) :: SPFreeSurfaceStabilization
contains
      procedure :: SpecificSet => SpecificSetFSS_3rd
   end type
   type(SPFreeSurfaceStabilization) :: SetPointersFSS_3rd
   
   

contains
   subroutine SpecificSetFSS_3rd(d)
      implicit none
      class(SPFreeSurfaceStabilization) :: d
      
      if (a%kfl_StabilizeFreeSurface == 1) then
         call ConcatenateProcedures(ProcHook_Initializations,Allocations)
         call ConcatenateProcedures(ProcHook_Finalizations,Deallocations)
         !Stabilization matrices are computed at the center of gravity
         !We prepend it to make sure it is done before the gauss points are changed by the cut elements
         call PrependProcedure(ProcHook_PreGauss,StabilizationMatrices_3rd)
         call PrependProcedure(ProcHook_PreGauss,Gathers)
         call PrependProcedure(ProcHook_PreGauss,CheckIsFirst2Layers)
      endif
   end subroutine
   
   
  


   !------------------------------------------------------------------------
   !Free surface stabilization
   subroutine Allocations
      implicit none
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%alloc(ndime,ndime,e%mnode,projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')

      call a%Memor%alloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')

      
      
      
   end subroutine
   
   subroutine DeAllocations
      implicit none
      call a%Memor%dealloc(ndime,ndime,size(projected_elgrvel,3),projected_elgrvel,'projected_elgrvel','alloc_freeSurfStabilization')
      call a%Memor%dealloc(ndime,ndime,projected_grvel,'projected_grvel','alloc_freeSurfStabilization')
      
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
         call e%gather(ndime*ndime,projected_elgrvel,a%StabilizeFreeSurface_VelocityGradient)
        
      endif
   
   end subroutine
   
   
   subroutine StabilizationMatrices_3rd
      implicit none
      real(rp) :: velgraddiff(ndime,ndime), pregraddiff(ndime)
      integer(ip) :: idime,inode,jnode
      real(rp)  :: ghosttimom
      
      
      if (isFirst2Layers) then
      
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
      
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus

            dvol = e%weigp(e%igaus)*e%detjm
            
            call e%interpg(ndime*ndime,projected_elgrvel,projected_grvel)
         
         
            !first the right hand side
            do inode = 1,e%pnode
            
                !Stabilize velocity array
                do idime = 1,ndime
                   elrhs(idime,inode) = elrhs(idime,inode) + a%StabilizeFreeSurface_Param*LHSdtinv*acden*chale(1)*chale(1)*dot_product(e%cartd(:,inode),projected_grvel(:,idime))*e%weigp(e%igaus)*e%detjm 
                enddo
              
            enddo
            
            !Second the left hand side
            do inode = 1,e%pnode
               do jnode = 1,e%pnode
               
                  !Stabilize the velocity array
                  do idime = 1,ndime
                     elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + a%StabilizeFreeSurface_Param*LHSdtinv*acden*chale(1)*chale(1)*dot_product(e%cartd(:,inode),e%cartd(:,jnode))*e%weigp(e%igaus)*e%detjm 
                  enddo
               enddo
            enddo
            
         enddo gauss_points
            
       
      endif
      
   
   end subroutine
   
 
   
end module

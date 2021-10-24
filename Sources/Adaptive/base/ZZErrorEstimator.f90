module Mod_ZZErrorEstimator
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_Element
   implicit none
   private
   public ZZErrorEstimator, GradientErrorEstimator
   
contains 

   subroutine ZZErrorEstimator(Mesh,Memor,ndofn,unkno,error) 
      implicit none
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: ndofn
      real(rp) :: unkno(:,:), error(:)
      
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable,  target :: gradient(:,:,:), elunk(:,:),grunk(:,:),elgrunk(:,:,:),gpgrunk(:,:),GrunkDiff(:,:)
      real(rp) :: DiffNorm,GrNorm,DiffNorm2
      real(rp) :: dvol
      integer(ip) :: ielem,nelem,ndime,npoin,inode,igaus
      
      call Mesh%GetNelem(nelem)
      
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','ZZErrorEstimator')
      call Mesh%GetNdime(ndime)
      call Mesh%GetNpoin(npoin)
      call Memor%alloc(ndofn,ndime,npoin,gradient,'gradient','ZZErrorEstimator')
      call Memor%alloc(ndofn,ndime,grunk,'grunk','ZZErrorEstimator')
      call Memor%alloc(ndofn,ndime,e%mnode,elgrunk,'elgrunk','ZZErrorEstimator')
      call Memor%alloc(ndofn,ndime,gpgrunk,'gpgrunk','ZZErrorEstimator')
      call Memor%alloc(ndofn,ndime,GrunkDiff,'GrunkDiff','ZZErrorEstimator')
      call Memor%alloc(ndofn,e%mnode,elunk,'elunk','ZZErrorEstimator')
      
      !First loop: compute and project the gradient
      do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         dvol = e%detjm
         
         call e%gather(ndofn,elunk,unkno)
         call e%gradient(ndofn,elunk,grunk) 
         
         elgrunk = 0.0_rp
         
         !Gauss Point Loop
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            
            
            dvol = e%weigp(e%igaus)*e%detjm
               
            do inode = 1,e%pnode
               elgrunk(:,:,inode) = elgrunk(:,:,inode) + e%shape(inode,e%igaus)*grunk*dvol
            enddo
         enddo gauss_points
         
         call Mesh%AssemblyToArray(e,ndofn*e%ndime,elgrunk,gradient)
      enddo
      
      !Project onto the nodes
      call Mesh%Smooth(ndofn*e%ndime,gradient) 
            
      !Second loop: compute the difference between the gradient and the projected gradient
      error(1:nelem) = 0

      !First loop: compute and project the gradient
      do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         dvol = e%detjm
         
         call e%gather(ndofn,elunk(:,:),unkno)
         call e%gradient(ndofn,elunk,grunk) 
         
         call e%gather(e%ndime*ndofn,elgrunk,gradient)
         call e%interpc(e%ndime*ndofn,elgrunk,gpGrunk)
         
         GrunkDiff = grunk - gpGrunk
         
         !call vecnor(GrunkDiff,e%ndime*ndofn,DiffNorm,2)
         !call vecnor(gpGrunk,e%ndime*ndofn,GrNorm,2)
         !error(ielem) = DiffNorm/GrNorm
         
         !call vecnor(GrunkDiff,e%ndime*ndofn,DiffNorm,2)
         !error(ielem) = DiffNorm*dvol
         
         call dot(GrunkDiff,GrunkDiff,e%ndime*ndofn,DiffNorm2)
         error(ielem) = DiffNorm2*dvol

         
      enddo

      call Memor%dealloc(ndofn,ndime,npoin,gradient,'gradient','ZZErrorEstimator')
      call Memor%dealloc(ndofn,ndime,grunk,'grunk','ZZErrorEstimator')
      call Memor%dealloc(ndofn,ndime,e%mnode,elgrunk,'elgrunk','ZZErrorEstimator')
      call Memor%dealloc(ndofn,ndime,gpgrunk,'gpgrunk','ZZErrorEstimator')
      call Memor%dealloc(ndofn,ndime,GrunkDiff,'GrunkDiff','ZZErrorEstimator')
      call Memor%dealloc(ndofn,e%mnode,elunk,'elunk','ZZErrorEstimator')
      
      call Mesh%ElementDeAlloc(e,Memor,'ForceClosedRule','ZZErrorEstimator')
      
   end subroutine
   
   subroutine GradientErrorEstimator(Mesh,Memor,ndofn,unkno,error)
      use typre
      use Mod_Memor
      use Mod_Mesh
      use Mod_Element
      implicit none
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: ndofn
      real(rp) :: unkno(:,:), error(:)
      
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable,  target ::  elunk(:,:),grunk(:,:)
      real(rp) :: grunknorm
      real(rp) :: dvol
      integer(ip) :: ielem,nelem,ndime,npoin,inode,igaus

      call Mesh%GetNelem(nelem)
   
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','ZZErrorEstimator')
      call Mesh%GetNdime(ndime)
      call Mesh%GetNpoin(npoin)
      call Memor%alloc(ndofn,ndime,grunk,'grunk','ZZErrorEstimator')
      call Memor%alloc(ndofn,e%mnode,elunk,'elunk','ZZErrorEstimator')
   
      !First loop: compute and project the gradient
      do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         dvol = e%detjm
         
         call e%elmlen
         
         call e%gather(ndofn,elunk,unkno)
         call e%gradient(ndofn,elunk,grunk) 
         
         call vecnor(grunk,e%ndime*ndofn,grunknorm,2)
         
         error(ielem) = grunknorm*(e%hleng(1)**2)
      enddo
      
      
      call Mesh%GetNdime(ndime)
      call Mesh%GetNpoin(npoin)
      call Memor%dealloc(ndofn,ndime,grunk,'grunk','ZZErrorEstimator')
      call Memor%dealloc(ndofn,e%mnode,elunk,'elunk','ZZErrorEstimator')
   
      call Mesh%ElementDeAlloc(e,Memor,'ForceClosedRule','ZZErrorEstimator')
   
   
   end subroutine
end module

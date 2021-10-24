module Mod_plcd_BMatrix_InfinitesimalStrain_2D
   use typre
   use Mod_plcd_BMatrix
   use Mod_Element
   use Mod_Memor
   implicit none
   private
   public BMatrix_IS2D
   
   integer(ip), parameter :: ndofn = 2, nVoigtComponents = 3

   type, extends(BMatrix) :: BMatrix_IS2D
   
            
      real(rp), allocatable :: B(:,:,:)
      integer(ip) :: pnode
contains
      procedure :: Alloc => Allocate_IS2D
      procedure :: Dealloc => DeAllocate_IS2D
      procedure :: Setup => Setup_IS2D
      procedure :: B_Times_Vector => B_Times_Vector_IS2D
      procedure :: Vector_Times_B => Vector_Times_B_IS2D
      procedure :: Bt_Times_Vector => Bt_Times_Vector_IS2D
      procedure :: Bt_Times_Matrix_Times_B => Bt_Times_Matrix_Times_B_IS2D
      procedure :: Betat_Times_Vector_Times_Beta => Betat_Times_Vector_Times_Beta_IS2D
      procedure :: SetDisplacementGradient => NULLSUB
   
   end type
   
contains

   subroutine Allocate_IS2D(a,e,Memor)
      class(BMatrix_IS2D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      call Memor%alloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','Allocate_IS2D')
      !call Memor%alloc(ndofn,e%mnode,a%Beta,'Beta','Allocate_IS2D')

   end subroutine
   
   subroutine DeAllocate_IS2D(a,e,Memor)
      class(BMatrix_IS2D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      
      call Memor%dealloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','DeAllocate_IS2D')
      !call Memor%dealloc(ndofn,e%mnode,a%Beta,'Beta','DeAllocate_IS2D')

   end subroutine

   subroutine Setup_IS2D(a,e)
      class(BMatrix_IS2D) :: a
      class(FiniteElement) :: e
      
      integer(ip) :: inode

      a%pnode = e%pnode
      do inode = 1,e%pnode
         !B matrix for Material Term
         a%B(1,1,inode) =e%cartd(1,inode)
         a%B(2,2,inode) =e%cartd(2,inode)
         !xy
         a%B(3,1,inode) =e%cartd(2,inode)
         a%B(3,2,inode) =e%cartd(1,inode)
         
      enddo
      
   end subroutine
   
   subroutine Bt_Times_Matrix_Times_B_IS2D(a,input,output)
      class(BMatrix_IS2D), target :: a
      real(rp), intent(in) :: input(:,:)
      real(rp), intent(out), target :: output(*)
      
      real(rp), pointer :: B2D(:,:)
      real(rp), pointer :: Output2D(:,:)
      
      !input(nVoigtComponents,nVoigtComponents)
      !output(ndofn*e%pnode,ndofn*e%pnode)
      call Pointer2D3D(nVoigtComponents,ndofn,a%pnode,a%B,B2D)
      !B2D(1:nVoigtComponents,1:ndofn*a%pnode) => a%B(:,:,1:a%pnode)
      Output2D(1:ndofn*a%pnode,1:ndofn*a%pnode) => output(1:ndofn*ndofn*a%pnode*a%pnode)
      
      Output2D = matmul(transpose(B2D),matmul(input,B2D))
   end subroutine
   
   subroutine Betat_Times_Vector_Times_Beta_IS2D(a,input,output)
      class(BMatrix_IS2D), target :: a
      real(rp), intent(in) :: input(:)
      real(rp), intent(out), target :: output(:,:)
      
      real(rp) :: Beta(ndofn,a%pnode)
      
      real(rp) :: inputMatrixform(ndofn,ndofn)
      !Input is in fact the stress vector
      
      Beta(1,1:a%pnode) = a%B(1,1,1:a%pnode)
      Beta(2,1:a%pnode) = a%B(2,2,1:a%pnode)
      
      inputMatrixform = 0.0_rp
      
      inputMatrixform =  reshape([input(1),input(3),input(3),input(2)],[ndofn,ndofn])
      
      output = matmul(transpose(Beta),matmul(inputMatrixform,Beta))

   end subroutine
      
   subroutine B_Times_Vector_IS2D(a,input,output)
      class(BMatrix_IS2D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
      call runend('B_Times_Vector_IS2D not ready')
   end subroutine
   
   subroutine Vector_Times_B_IS2D(a,input,output)
      class(BMatrix_IS2D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
      call runend('Vector_Times_B_IS2D not ready')
   end subroutine
   
   subroutine Bt_Times_Vector_IS2D(a,input,output)
      class(BMatrix_IS2D), target :: a
      real(rp), intent(in) :: input(:)
      real(rp), intent(out), target :: output(*)
      
      real(rp), pointer :: B2D(:,:)
      real(rp), pointer :: Output2D(:)
      
      !B2D(1:nVoigtComponents,1:ndofn*a%pnode) => a%B(:,:,1:a%pnode)
      call Pointer2D3D(nVoigtComponents,ndofn,a%pnode,a%B,B2D)
      Output2D(1:ndofn*a%pnode) => output(1:ndofn*a%pnode)
      
      Output2D = matmul(transpose(B2D),input(1:nVoigtComponents))
      
   end subroutine
   
   subroutine Pointer2D3D(sizeunk1,sizeunk2,sizeunk3,unk,ppointer)
      use typre
      implicit none
      real(rp), target :: unk(*)
      real(rp), pointer :: ppointer(:,:)
      integer(ip) :: sizeunk1,sizeunk2,sizeunk3
      
      ppointer(1:sizeunk1,1:sizeunk2*sizeunk3) => unk(1:sizeunk1*sizeunk2*sizeunk3)
   end subroutine
   
   subroutine NULLSUB(a,DisplacementGradient)
      class(BMatrix_IS2D), target :: a
      real(rp), intent(in) :: DisplacementGradient(:,:)
   end subroutine

end module
   
   

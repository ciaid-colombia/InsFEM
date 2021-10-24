module Mod_plcd_BMatrix_TotalLargeStrains_3D
   use typre
   use Mod_plcd_BMatrix
   use Mod_Element
   use Mod_Memor
   implicit none
   private
   public BMatrix_LSIS3D
   
   integer(ip), parameter :: ndofn = 3, nVoigtComponents = 6

   type, extends(BMatrix) :: BMatrix_LSIS3D
   
            
      real(rp), allocatable :: B(:,:,:)
      real(rp), allocatable :: Beta(:,:)
      real(rp), allocatable :: DisplacementGradient(:,:)
      integer(ip) :: pnode
contains
      procedure :: Alloc => Allocate_LSIS3D
      procedure :: Dealloc => DeAllocate_LSIS3D
      procedure :: Setup => Setup_LSIS3D
      procedure :: B_Times_Vector => B_Times_Vector_LSIS3D
      procedure :: Vector_Times_B => Vector_Times_B_LSIS3D
      procedure :: Bt_Times_Vector => Bt_Times_Vector_LSIS3D
      procedure :: Bt_Times_Matrix_Times_B => Bt_Times_Matrix_Times_B_LSIS3D
      procedure :: Betat_Times_Vector_Times_Beta => Betat_Times_Vector_Times_Beta_LSIS3D
      procedure :: SetDisplacementGradient => SetDisplacementGradient_LSIS3D
   
   end type
   
contains

   subroutine Allocate_LSIS3D(a,e,Memor)
      class(BMatrix_LSIS3D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      call Memor%alloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','Allocate_LSIS3D')
      call Memor%alloc(ndofn,e%mnode,a%Beta,'Beta','Allocate_LSIS3D')
      call Memor%alloc(ndofn,ndofn,a%DisplacementGradient,'DisplacementGradient','Allocate_LSIS3D')

   end subroutine
   
   
   subroutine DeAllocate_LSIS3D(a,e,Memor)
      class(BMatrix_LSIS3D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      call Memor%dealloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','DeAllocate_LSIS3D')
      call Memor%dealloc(ndofn,e%mnode,a%Beta,'Beta','DeAllocate_LSIS3D')
      call Memor%dealloc(ndofn,ndofn,a%DisplacementGradient,'DisplacementGradient','DeAllocate_LSIS3D')

   end subroutine

   subroutine Setup_LSIS3D(a,e)
      class(BMatrix_LSIS3D) :: a
      class(FiniteElement) :: e
      
      real(rp) :: F_mat(ndofn,ndofn)
      
      integer(ip) :: inode, idime
      
      a%pnode = e%pnode
      
      F_mat = 0.0_rp
      F_mat = a%DisplacementGradient
      forall(idime = 1:ndofn) F_mat(idime,idime) = F_mat(idime,idime) + 1.0_rp 
      
      do inode = 1,e%pnode
         !B matrix for Material Term
         a%B(1,1,inode) =e%cartd(1,inode)*F_mat(1,1)
         a%B(1,2,inode) =e%cartd(1,inode)*F_mat(2,1)
         a%B(1,3,inode) =e%cartd(1,inode)*F_mat(3,1)
         a%B(2,1,inode) =e%cartd(2,inode)*F_mat(1,2)
         a%B(2,2,inode) =e%cartd(2,inode)*F_mat(2,2)
         a%B(2,3,inode) =e%cartd(2,inode)*F_mat(3,2)
         a%B(3,1,inode) =e%cartd(3,inode)*F_mat(1,3)
         a%B(3,2,inode) =e%cartd(3,inode)*F_mat(2,3)
         a%B(3,3,inode) =e%cartd(3,inode)*F_mat(3,3)
         !xy
         a%B(4,1,inode) =e%cartd(1,inode)*F_mat(1,2)+e%cartd(2,inode)*F_mat(1,1)
         a%B(4,2,inode) =e%cartd(1,inode)*F_mat(2,2)+e%cartd(2,inode)*F_mat(2,1)
         a%B(4,3,inode) =e%cartd(1,inode)*F_mat(3,2)+e%cartd(2,inode)*F_mat(3,1)
         !xz
         a%B(5,1,inode) =e%cartd(1,inode)*F_mat(1,3)+e%cartd(3,inode)*F_mat(1,1)
         a%B(5,2,inode) =e%cartd(1,inode)*F_mat(2,3)+e%cartd(3,inode)*F_mat(2,1)
         a%B(5,3,inode) =e%cartd(1,inode)*F_mat(3,3)+e%cartd(3,inode)*F_mat(3,1)
         !yz
         a%B(6,1,inode) =e%cartd(2,inode)*F_mat(1,3)+e%cartd(3,inode)*F_mat(1,2)  
         a%B(6,2,inode) =e%cartd(2,inode)*F_mat(2,3)+e%cartd(3,inode)*F_mat(2,2)
         a%B(6,3,inode) =e%cartd(2,inode)*F_mat(3,3)+e%cartd(3,inode)*F_mat(3,2)
         
         !Beta matrix 
         a%Beta(1,inode) = e%cartd(1,inode)
         a%Beta(2,inode) = e%cartd(2,inode)
         a%Beta(3,inode) = e%cartd(3,inode)
      
      enddo
      
      a%DisplacementGradient = 0.0_rp
      
   end subroutine
   
   subroutine Bt_Times_Matrix_Times_B_LSIS3D(a,input,output)
      class(BMatrix_LSIS3D), target :: a
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
   
   subroutine Betat_Times_Vector_Times_Beta_LSIS3D(a,input,output)
      class(BMatrix_LSIS3D), target :: a
      real(rp), intent(in) :: input(:)
      real(rp), intent(out), target :: output(:,:)
      
      real(rp) :: inputMatrixform(ndofn,ndofn)
      !Input is in fact the stress vector
      
      inputMatrixform = 0.0_rp
      
      inputMatrixform =  reshape([input(1),input(4),input(5),input(4),input(2),input(6),input(5),input(6),input(3)],&
      [ndofn,ndofn])
      
      output = matmul(transpose(a%Beta),matmul(inputMatrixform,a%Beta))
      
   end subroutine
      
   subroutine B_Times_Vector_LSIS3D(a,input,output)
      class(BMatrix_LSIS3D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
      call runend('B_Times_Vector_LSIS3D not ready')
   end subroutine
   
   subroutine Vector_Times_B_LSIS3D(a,input,output)
      class(BMatrix_LSIS3D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
       call runend('Vector_Times_B_LSIS3D not ready')
   end subroutine
   
   subroutine Bt_Times_Vector_LSIS3D(a,input,output)
      class(BMatrix_LSIS3D), target :: a
      real(rp), intent(in) :: input(:)
      real(rp), intent(out), target :: output(*)
      
      real(rp), pointer :: B2D(:,:)
      real(rp), pointer :: Output2D(:)
      
      call Pointer2D3D(nVoigtComponents,ndofn,a%pnode,a%B,B2D)
      !B2D(1:nVoigtComponents,1:ndofn*a%pnode) => a%B(:,:,1:a%pnode)
      Output2D(1:ndofn*a%pnode) => output(1:ndofn*a%pnode)
      
      Output2D = matmul(transpose(B2D),input)
      
   end subroutine
   
   subroutine Pointer2D3D(sizeunk1,sizeunk2,sizeunk3,unk,ppointer)
      use typre
      implicit none
      real(rp), target :: unk(*)
      real(rp), pointer :: ppointer(:,:)
      integer(ip) :: sizeunk1,sizeunk2,sizeunk3
      
      ppointer(1:sizeunk1,1:sizeunk2*sizeunk3) => unk(1:sizeunk1*sizeunk2*sizeunk3)
   end subroutine
   
   subroutine SetDisplacementGradient_LSIS3D(a,DisplacementGradient)
      class(BMatrix_LSIS3D), target :: a
      real(rp),intent(in) :: DisplacementGradient(:,:)
      
      a%DisplacementGradient = DisplacementGradient
   
   end subroutine

end module
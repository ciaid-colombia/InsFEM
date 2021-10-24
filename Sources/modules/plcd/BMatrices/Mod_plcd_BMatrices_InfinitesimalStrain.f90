module Mod_plcd_BMatrix_InfinitesimalStrain_3D
   use typre
   use Mod_plcd_BMatrix
   use Mod_Element
   use Mod_Memor
   implicit none
   private
   public BMatrix_IS3D
   
   integer(ip), parameter :: ndofn = 3, nVoigtComponents = 6

   type, extends(BMatrix) :: BMatrix_IS3D
   
            
      real(rp), allocatable :: B(:,:,:)
      integer(ip) :: pnode
contains
      procedure :: Alloc => Allocate_IS3D
      procedure :: Dealloc => DeAllocate_IS3D
      procedure :: Setup => Setup_IS3D
      procedure :: B_Times_Vector => B_Times_Vector_IS3D
      procedure :: Vector_Times_B => Vector_Times_B_IS3D
      procedure :: Bt_Times_Vector => Bt_Times_Vector_IS3D
      procedure :: Bt_Times_Matrix_Times_B => Bt_Times_Matrix_Times_B_IS3D
      procedure :: Betat_Times_Vector_Times_Beta => Betat_Times_Vector_Times_Beta_IS3D
      procedure :: SetDisplacementGradient => NULLSUB
   
   end type
   
contains

   subroutine Allocate_IS3D(a,e,Memor)
      class(BMatrix_IS3D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      call Memor%alloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','Allocate_IS3D')

   end subroutine
   
   
   subroutine DeAllocate_IS3D(a,e,Memor)
      class(BMatrix_IS3D) :: a
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      
      call Memor%dealloc(nVoigtComponents,ndofn,e%mnode,a%B,'B','DeAllocate_IS3D')

   end subroutine

   subroutine Setup_IS3D(a,e)
      class(BMatrix_IS3D) :: a
      class(FiniteElement) :: e
      
      integer(ip) :: inode

      a%pnode = e%pnode
      do inode = 1,e%pnode
         !B matrix for Material Term
         a%B(1,1,inode) =e%cartd(1,inode)
         a%B(2,2,inode) =e%cartd(2,inode)
         a%B(3,3,inode) =e%cartd(3,inode)
         !xy
         a%B(4,1,inode) =e%cartd(2,inode)
         a%B(4,2,inode) =e%cartd(1,inode)
         !xz
         a%B(5,1,inode) =e%cartd(3,inode)
         a%B(5,3,inode) =e%cartd(1,inode)
         !yz
         a%B(6,2,inode) =e%cartd(3,inode)
         a%B(6,3,inode) =e%cartd(2,inode)
      enddo
   end subroutine
   
   subroutine Bt_Times_Matrix_Times_B_IS3D(a,input,output)
      class(BMatrix_IS3D), target :: a
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
   
   subroutine Betat_Times_Vector_Times_Beta_IS3D(a,input,output)
      class(BMatrix_IS3D), target :: a
      real(rp), intent(in) :: input(:)
      real(rp), intent(out), target :: output(:,:)
      
      real(rp) :: Beta(ndofn,a%pnode)
      
      real(rp) :: inputMatrixform(ndofn,ndofn)
      !Input is in fact the stress vector
      
      Beta(1,1:a%pnode) = a%B(1,1,1:a%pnode)
      Beta(2,1:a%pnode) = a%B(2,2,1:a%pnode)
      Beta(3,1:a%pnode) = a%B(3,3,1:a%pnode)
      
      inputMatrixform = 0.0_rp
      
      inputMatrixform =  reshape([input(1),input(6),input(5),input(6),input(2),input(4),input(5),input(4),input(3)],&
      [ndofn,ndofn])
      
      output = matmul(transpose(Beta),matmul(inputMatrixform,Beta))
      
   end subroutine
      
   subroutine B_Times_Vector_IS3D(a,input,output)
      class(BMatrix_IS3D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
      call runend('B_Times_Vector_IS3D not ready')
   end subroutine
   
   subroutine Vector_Times_B_IS3D(a,input,output)
      class(BMatrix_IS3D) :: a
      real(rp), intent(in) :: input(:)
      real(rp) :: output(:)
      
       call runend('Vector_Times_B_IS3D not ready')
   end subroutine
   
   subroutine Bt_Times_Vector_IS3D(a,input,output)
      class(BMatrix_IS3D), target :: a
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
   
   subroutine NULLSUB(a,DisplacementGradient)
      class(BMatrix_IS3D), target :: a
      real(rp), intent(in) :: DisplacementGradient(:,:)
   end subroutine

end module
   
   

module Mod_plcd_BMatrix
   use typre
   use Mod_Element
   use Mod_Memor
   implicit none
   private
   public BMatrix

   type, abstract :: BMatrix
   
   
contains
      procedure(Alloc), deferred :: Alloc
      procedure(DeAlloc), deferred :: DeAlloc
      procedure(Setup), deferred :: Setup
      procedure(B_Times_Vector), deferred :: B_Times_Vector
      procedure(Vector_Times_B), deferred :: Vector_Times_B
      procedure(Bt_Times_Vector), deferred :: Bt_Times_Vector
      procedure(Bt_Times_Matrix_Times_B), deferred :: Bt_Times_Matrix_Times_B
      procedure(Betat_Times_Vector_Times_Beta), deferred :: Betat_Times_Vector_Times_Beta
      procedure(SetDisplacementGradient), deferred :: SetDisplacementGradient
      
   
   end type

   abstract interface
   
      subroutine Alloc(a,e,Memor)
         import
         class(BMatrix) :: a
         class(FiniteElement) :: e
         type(MemoryMan) :: Memor
      end subroutine
      
      subroutine DeAlloc(a,e,Memor)
         import
         class(BMatrix) :: a
         class(FiniteElement) :: e
         type(MemoryMan) :: Memor
      end subroutine
   
      subroutine Setup(a,e)
         import
         class(BMatrix) :: a
         class(FiniteElement) :: e
      
      end subroutine
   
      subroutine B_Times_Vector(a,input,output)
         import
         class(BMatrix) :: a
         real(rp), intent(in) :: input(:)
         real(rp) :: output(:)
      end subroutine
      
      subroutine Vector_Times_B(a,input,output)
         import
         class(BMatrix) :: a
         real(rp), intent(in) :: input(:)
         real(rp) :: output(:)
      end subroutine
      
      subroutine Bt_Times_Vector(a,input,output)
         import
         class(BMatrix), target :: a
         real(rp), intent(in) :: input(:)
         real(rp), intent(out), target :: output(*)
      end subroutine
      
      subroutine Bt_Times_Matrix_Times_B(a,input,output)
         import
         class(BMatrix), target :: a
         real(rp), intent(in) :: input(:,:)
         real(rp), intent(out), target :: output(*)
      end subroutine
      
      subroutine Betat_Times_Vector_Times_Beta(a,input,output)
         import
         class(BMatrix), target :: a
         real(rp), intent(in) :: input(:)
         real(rp), intent(out), target :: output(:,:)
      end subroutine
   
      subroutine SetDisplacementGradient(a, DisplacementGradient)
         import
         class(BMatrix), target :: a
         real(rp), intent(in) :: DisplacementGradient(:,:)
      end subroutine
   
   end interface
   
contains

end module
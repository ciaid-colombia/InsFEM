module Mod_ElementDataStructures
   use typre
   implicit none
   
   type :: ElementDataStructure
      !Element jacobians and derivatives at gauss points and center of gravity
      real(rp), allocatable ::&
         cartd(:,:,:),&                 ! Cartesian derivatives
         cartg(:,:),  &                 ! Cartesian derivatives at the center of gravity
         tragl(:,:),&                   ! Inverse of jacobian at c.o.g
         hleng(:)                       ! Element length at c.o.g
      real(rp), allocatable ::&
         hessi(:,:,:)                   ! Hessian
      real(rp), allocatable ::&
         elvec(:,:,:,:),&               ! Normal and tangents
         elnor(:,:)
      
      !Element values at Gauss points
      real(rp), allocatable ::&
         detjm(:)                     ! Jacobian determinant
      real(rp)  :: detjmg             ! Jacobian determinant at center of gravity
   end type
   
   type :: BoundaryDataStructure
      real(rp), allocatable :: &
         cartb(:,:,:)                 ! Cartesian derivatives in the boundary
         !Boundary jacobians and derivatives
      real(rp), allocatable ::&
         baloc(:,:,:),   &              ! Local base of boundary element
         eucta(:)
   end type
end module



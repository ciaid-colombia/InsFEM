module Mod_supm_MatrixVector
     use typre
     use Mod_Element
  contains
  
   subroutine ConvertMatrixToVector(e,Matrix,Vector)

      implicit none
      class(FiniteElement) :: e
      real(rp) :: Matrix(e%ndime,e%ndime), Vector(e%ndime*e%ndime)
      integer(ip) :: ndime
      
      
      if (e%ndime==2) then
      
         Vector(1)=Matrix(1,1)
         Vector(2)=Matrix(2,2)
         Vector(3)=Matrix(1,2)
         Vector(4)=Matrix(2,1)
      
      elseif (e%ndime==3) then
      
         Vector(1)=Matrix(1,1)
         Vector(2)=Matrix(2,2)
         Vector(3)=Matrix(3,3)
         Vector(4)=Matrix(1,2)
         Vector(5)=Matrix(2,1)
         Vector(6)=Matrix(2,3)
         Vector(7)=Matrix(3,2)
         Vector(8)=Matrix(3,1)
         Vector(9)=Matrix(1,3)
      
      end if
   
   end subroutine
  
   subroutine ConvertVectorToMatrix(e,Vector,Matrix)

      implicit none
      class(FiniteElement) :: e
      real(rp) :: Matrix(e%ndime,e%ndime), Vector(e%ndime*e%ndime)
      integer(ip) :: ndime
   
      if (e%ndime==2) then
      
         Matrix(1,1)=Vector(1)
         Matrix(2,2)=Vector(2)
         Matrix(1,2)=Vector(3)
         Matrix(2,1)=Vector(4)
      
      elseif (e%ndime==3) then
      
         Matrix(1,1)=Vector(1)
         Matrix(2,2)=Vector(2)
         Matrix(3,3)=Vector(3)
         Matrix(1,2)=Vector(4)
         Matrix(2,1)=Vector(5)
         Matrix(2,3)=Vector(6)
         Matrix(3,2)=Vector(7)
         Matrix(3,1)=Vector(8)
         Matrix(1,3)=Vector(9)
      
      end if
   
   end subroutine
   
   subroutine ConvertVectorToDiagMatrix(ndime,Vector,Matrix)
   
   implicit none
      real(rp) :: Matrix(ndime,ndime), Vector(ndime*ndime)
      integer(ip) :: i,ndime
      
      Matrix=0.0_rp 
      do i=1, ndime
         Matrix(i,i)=Vector(i)
      end do
   end subroutine
   
    subroutine ConvertSigmaMatrixToVector(ndime,tens,Matrix,Vector)

      implicit none
      integer(ip) :: ndime,tens
      real(rp) :: Matrix(ndime,ndime), Vector(tens)
      
      
      
      if (ndime==2) then
      
         Vector(1)=Matrix(1,1)
         Vector(2)=Matrix(2,2)
         Vector(3)=Matrix(1,2)
      
      elseif (ndime==3) then
      
         Vector(1)=Matrix(1,1)
         Vector(2)=Matrix(2,2)
         Vector(3)=Matrix(3,3)
         Vector(4)=Matrix(2,3)
         Vector(5)=Matrix(1,3)
         Vector(6)=Matrix(1,2)
      end if
   
   end subroutine
   
   
   end module
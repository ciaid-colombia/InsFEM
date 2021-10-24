module Mod_MatVec4D
   implicit none
contains
   subroutine MatVec4D(matrix,vector,unkno)
      use typre
      implicit none
      real(rp), CONTIGUOUS, target :: matrix(:,:,:,:), vector(:,:), unkno(:,:)
      
      real(rp), pointer :: auxmat(:,:) => NULL(), auxvec(:) => NULL(), auxunkno(:) => NULL()
      
      integer(ip) :: s1,s2,s3,s4,s5,s6
      
      s1 = size(matrix,1)*size(matrix,2)
      s2 = size(matrix,3)*size(matrix,4)
      
      s3 = size(vector,1)*size(vector,2)
      S4 =size(unkno,1)*size(unkno,2)  
      auxmat(1:s1,1:s2) => matrix
      auxvec(1:s3) => vector
      auxunkno(1:s4) => unkno(:,:)
      
      auxunkno = matmul(auxmat,auxvec)

   end subroutine 
end module
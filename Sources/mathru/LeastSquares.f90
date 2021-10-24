subroutine LeastSquaresBuildSystem(msamp,ndof,X,y,systemMat,systemRHS)
   use typre
   implicit none
   !See http://en.wikipedia.org/wiki/Linear_least_squares_(mathematics) for notation
   
   integer(ip) :: ndof,msamp
   
   real(rp) :: X(msamp,ndof), y(msamp)
   real(rp) :: systemMAT(ndof,ndof), systemRHS(ndof)
   
   write(*,*) 'before leastsquresbuild'
   systemMat = 0.0_rp
   systemRHS = 0.0_rp
   write(*,*) 'before leastsquresbuild2'
   call LeastSquaresAddToSystem(msamp,ndof,X,y,systemMat,systemRHS)
   write(*,*) 'after leastsquresbuild'
end subroutine

subroutine LeastSquaresAddToSystem(msamp,ndof,X,y,systemMat,systemRHS)
   use typre
   implicit none
   !See http://en.wikipedia.org/wiki/Linear_least_squares_(mathematics) for notation
   
   integer(ip) :: ndof,msamp
   
   real(rp) :: X(msamp,ndof), y(msamp)
   real(rp) :: systemMAT(ndof,ndof), systemRHS(ndof)
   
   real(rp) :: XT(ndof,msamp)
   
   XT = transpose(X)
   
   systemMat = matmul(XT,X)
   systemRHS = matmul(XT,y)
end subroutine
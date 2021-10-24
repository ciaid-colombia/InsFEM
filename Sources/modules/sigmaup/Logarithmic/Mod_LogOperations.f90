module Mod_LogOperations
   use typre
   use Mod_supm_MatrixVector
   implicit none
   
   !-----------------------------
   !This module contains operation used by de Logarithm Conformation Reformulation Model
   !-----------------------------
   
   
contains 

   subroutine sup_ComputeExponential(ndime,auxtens,Vector,Expo_Vect,Expo_Matrix)
   !Compute Exp(Psi), where Psi is a tensor in voigt notation (i.e. a vector)
      
      implicit none
      integer(ip) :: auxtens,ndime
      real(rp)  :: D(ndime), V(ndime,ndime), V_inv(ndime,ndime)
      real(rp)  :: Exp1(ndime,ndime),D_Mat(ndime,ndime)
      real(rp)  :: Expo_Matrix(ndime,ndime), Expo_Vect(auxtens)
      real(rp)  :: Vector(auxtens), Matrix(ndime,ndime)
      real(rp)  :: ExpoDiagMatrix(ndime,ndime), ExpoD(ndime)
      real(rp)  ::  Sol1(ndime,ndime), Sol(ndime,ndime)
      integer   :: info,lwork,i
      real(rp)  :: work(ndime*4)
      lwork = 4*ndime
      
      call sup_SigmaMatrix(ndime,auxtens,Vector,Matrix)
      
      V=Matrix
      CALL DSYEV( 'Vectors', 'Upper', ndime, V, ndime, D, work, lwork, info)
  
      V_inv=transpose(V) ! V^-1 = V^T because of orthogonality
      
      ExpoD=exp(D)     
      call ConvertVectorToDiagMatrix(ndime,ExpoD,ExpoDiagMatrix)
    
      Exp1=matmul(V,ExpoDiagMatrix)
      Expo_Matrix=matmul(Exp1,V_inv)
      
      call ConvertSigmaMatrixToVector(ndime,auxtens,Expo_Matrix,Expo_Vect)
      
   end subroutine 
   

   subroutine sup_LogarithmicTransformation(ndime,ntens,auxL,sigma,psi)
   !Compute Log(auxL*S+I)
   
      implicit none
      real(rp) :: sigma(ntens),psi(ntens),auxL_inv,auxL
      real(rp) :: Sigma_Matrix(ndime,ndime), V(ndime,ndime),D(ndime),V_inv(ndime,ndime)
      real(rp) :: logD(ndime),LogDMatrix(ndime,ndime), Log_Matrix(ndime,ndime),tau(ntens)
      integer(ip) :: ntens,i,j,ndime
      integer   :: info,lwork
      real(rp)  :: work(ndime*4)
      lwork = 4*ndime
      call sup_ComputeConformationTensor(ndime,ntens,auxL,sigma,tau)

      call sup_SigmaMatrix(ndime,ntens,tau,Sigma_Matrix)
      V=Sigma_Matrix
      CALL DSYEV( 'Vectors', 'Upper', ndime, V, ndime, D, work, lwork, info)
  
      V_inv=transpose(V) ! V^-1 = V^T because of orthogonality
      
      LogD=log(D)     
      call ConvertVectorToDiagMatrix(ndime,LogD,LogDMatrix)
      
      Log_Matrix=matmul(matmul(V,LogDMatrix),V_inv)
      
      call ConvertSigmaMatrixToVector(ndime,ntens,Log_Matrix,psi)

   end subroutine
   

   subroutine sup_ComputeConformationTensor(ndime,auxtens,auxL,sigma,tau)
      !Compute tau=auxL*S+I
   
      implicit none
      integer(ip) :: i,j, auxtens,ndime
      real(rp) :: sigma(auxtens), tau(auxtens),auxL,auxL_inv
      
      auxL_inv=1.0_rp/auxL
      tau=0.0_rp
      tau= auxL_inv*sigma
      tau(1:ndime)=1.0_rp + tau(1:ndime) 
       
   end subroutine
   
   
   subroutine sup_ComputeEigenVectors(ndime,Matrix,eigenvectors,eigenvalues)
      implicit none
      real(rp)             :: Matrix(ndime,ndime),eigenvectors(ndime,ndime), eigenvalues(ndime)
      real(rp)             :: R(ndime,ndime), D(ndime)
      integer              :: info,lwork,ndime
      real(rp)             :: work(ndime*4)
      
      lwork = 4*ndime
      R=Matrix
      CALL DSYEV( 'Vectors', 'Upper', ndime, R, ndime, D, work, lwork, info)
      Eigenvectors=R
      Eigenvalues=D
   end subroutine
   

end module   
   
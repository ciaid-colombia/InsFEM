subroutine MatrixEigenValues(N,C,Eigenvalues)
   use typre
   implicit none
   integer(ip) :: N
   real(rp) :: C(N,N), Eigenvalues(N)
   
   real(rp) :: Caux(N,N)
   
   INTEGER          INFO, LWORK
      integer(ip)      :: LDA
      integer(ip) :: LWMAX
      real(rp) ::    WORK(N*N*5)

      Caux = C
      
      LWMAX = N*N
      LWORK = -1
      CALL DSYEV( 'Vectors', 'Upper', N, Caux, N, Eigenvalues, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      !Solve eingenvalue problem
      CALL DSYEV( 'Vectors', 'Upper', N, Caux, N, EigenValues, WORK, LWORK, INFO )

!       IF( INFO.GT.0 ) THEN
!          call runend('Failed to compute eigenvalues')
!       END IF


end subroutine 

subroutine sld_rotmat(imodi,nnode,mnode,ndime,ndof,amatr,bvect,&
  &               rotma)

!------------------------------------------------------------------------
!
! This routine modifies a matrix the components of which are
! submatrices by rotating them using a matrix R ( = ROTMA) as follows:
!
!                |       A_11  ....      A_1i  R   ....      A_1n   |
!                |   ............................................   |
!     AMATR <--  |   R^t A_i1  ....  R^t A_ii  R   ....  R^t A_in   |
!                |   ............................................   |
!                |       A_n1  ....      A_ni  R   ....      A_nn   |
!
! where i ( = IMODI) is a given position. Also, the i-subvector of the
! given vector BVECT is multiplied by R^t. 
!
! When A is the lumped mass matrix Aij=0 when i/=j and then
!
!                |       A_11  ....       0        ....        0    |
!                |   ............................................   |
!     AMATR <--  |         0   ....  R^t A_ii  R   ....        0    |
!                |   ............................................   |
!                |         0   ....       0        ....      A_nn   |
!
! and A_11 = int N1 * id(ndime,ndime)
!
! When A is the consistent mass matrix
!
!                |       A_11  ....      A_1i  R   ....      A_1n   |
!                |   ............................................   |
!     AMATR <--  |   R^t A_i1  ....  R^t A_ii  R   ....  R^t A_in   |
!                |   ............................................   |
!                |       A_n1  ....      A_ni  R   ....      A_nn   |
!
! and A_ij = int Ni Nj * id(ndime,ndime)
!------------------------------------------------------------------------
  use def_parame
  implicit none
  integer(ip), intent(in)    :: nnode,mnode,ndime,ndof
  real(rp),    intent(in)    :: rotma(ndime,ndime)
!  real(rp),    intent(inout) :: amatr(mtott,mtott)
  real(rp),    intent(inout) :: amatr(ndof,mnode,ndof,mnode)
  real(rp),    intent(inout) :: bvect(ndof,mnode)
  integer(ip)                :: imodi,inode,jnode
  integer(ip)                :: idime,jdime,kdime,idof,jdof
!
! Modifies column number IMODI of AMATR ( A_j,imodi <-- A_j,imodi R )
!
forall(idof=1:ndof,inode=1:nnode)
        amatr(idof,inode,1:ndime,imodi) = matmul(amatr(idof,inode,1:ndime,imodi),rotma)
end forall

forall(jdof=1:ndof,jnode=1:nnode)
        amatr(1:ndime,imodi,jdof,jnode) = matmul(transpose(rotma),amatr(1:ndime,imodi,jdof,jnode))
end forall

bvect(1:ndime,imodi) = matmul(transpose(rotma),bvect(1:ndime,imodi))
  
end subroutine sld_rotmat

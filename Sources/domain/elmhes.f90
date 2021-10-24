subroutine elmhes(heslo,hessi,ndime,nnode,ntens,xjaci,deriv,elcod)
!-----------------------------------------------------------------------
!****f* Domain/elmhes
! NAME
!    elmhes
! DESCRIPTION
!    This routine computes the Hessian matrix in the Cartesian system
!    of coordinates according to the rule
!
!    d^2 N / d x_i d x_j
!       = (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j)
!       + (d N / d s_k) (d^2 s_k / d x_i d x_j) 
! USES
!    vetoma
!    btdbma
!    matove
! SOURCE
!***
!-----------------------------------------------------------------------
  use typre
  implicit none

  integer(ip), intent(in)  :: ndime,nnode,ntens
  real(rp),    intent(in)  :: xjaci(ndime,ndime)
  real(rp),    intent(in)  :: deriv(ndime,nnode),elcod(ndime,nnode)
  real(rp),    intent(in)  :: heslo(ntens,nnode)
  real(rp),    intent(out) :: hessi(ntens,nnode) 
  
  real(rp) :: wmat1(ndime,ndime,nnode)      !really only 
  real(rp) :: wmat2(ndime,ndime,nnode)      !work not even
  real(rp) :: d2sdx(ndime,ndime,ndime)      !intent out
  integer(ip)              :: inode,idime,jdime,kdime,ldime
!
! Transforms the array HESLO to a symmetric matrix WMAT1
!
  do inode=1,nnode
     call vetoma(heslo(1,inode),wmat1(1,1,inode),ndime,ntens)
  end do
!
! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for
! each node
!
  do inode=1,nnode
     call btdbma(wmat2(1,1,inode),wmat1(1,1,inode),xjaci,ndime,ndime)
  end do
!
! Obtains (d^2 s_k / d x_i d x_j) as the solution of the system
! (d x_l / d s_k) (d^2 s_k / d x_i d x_j) 
!     = - (d^2 x_l / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j), 
! for l,i,j = 1,...,NDIME
!
  do kdime=1,ndime
     do idime=1,ndime
        do jdime=1,ndime
           d2sdx(kdime,idime,jdime)=0.0_rp
           do ldime=1,ndime
              do inode=1,nnode
                 d2sdx(kdime,idime,jdime)=d2sdx(kdime,idime,jdime)&
                      -xjaci(kdime,ldime)*wmat2(idime,jdime,inode)&
                      *elcod(ldime,inode)
              end do
           end do
        end do
     end do
  end do
!
! Computes the second Cartesian derivatives of the shape functions
!
  do inode=1,nnode
     do idime=1,ndime
        do jdime=1,ndime
           do kdime=1,ndime
              wmat2(idime,jdime,inode)=wmat2(idime,jdime,inode)&
                   +deriv(kdime,inode)*d2sdx(kdime,idime,jdime)
           end do
        end do
     end do
  end do
!
! Writes the Hessian matrix as an array
!
  do inode=1,nnode
     call matove(wmat2(1,1,inode),hessi(1,inode),ndime,ntens)
  end do
  
end subroutine elmhes

subroutine cderda(a)
!-----------------------------------------------------------------------
!****f* Domain/cderda
! NAME
!    cderda
! DESCRIPTION
!    This routine defines the following derivated parameters:
!       llapl(ielty)             (defined in setdim)
!       ltopo(ielty)             (defined in setdim)
!       lrule(ielty)             (defined in setdim)
!       nnodb(ielty)             (defined in setdim)
!       ngaub(ielty)             (defined in setdim)
!       hnatu(ielty)             (defined in setdim)
!
!    and allocates real arrays used for integration
!-----------------------------------------------------------------------
  use def_parame
  use Mod_Mesh
  implicit none
  
  class(FemMesh) :: a
  
  integer(ip) :: ielty,istat

  a%ndimb=a%ndime-1
  a%ntens=3*a%ndime-3

  !Allocate arrays that store mesh info
  call a%Memor%alloc(a%nelty,a%llapl,'LLAPL','cderda')
  call a%Memor%alloc(a%nelty,a%ltopo,'LTOPO','cderda')
  call a%Memor%alloc(a%nelty,a%lrule,'LRULE','cderda')
  call a%Memor%alloc(a%nelty,a%nnodb,'NNODB','cderda')
  call a%Memor%alloc(a%nelty,a%ngaub,'NGAUB','cderda')
  call a%Memor%alloc(a%nelty,a%hnatu,'HNATU','cderda')
  
  !Compute nnodb,ngaub,llapl,ltopo,lrule,hnatu
  do ielty=1,a%nelty
     call setdim(&
          a%nnode(ielty),a%ngaus(ielty),a%lquad(ielty),a%ndime,&
          a%llapl(ielty),a%ltopo(ielty),a%lrule(ielty),&
          a%nnodb(ielty),a%ngaub(ielty),a%hnatu(ielty))     
  end do

  !Maximum values
  a%mnodb=maxval(a%nnodb)
  a%mgaus=maxval(a%ngaus)
  a%mgaub=maxval(a%ngaub)
  a%mlapl=maxval(a%llapl)
  
  call a%Memor%alloc(a%ndime,a%mgaus,a%nelty,a%posgp,'POSGP','cderda')   

end subroutine cderda

subroutine codfix(ndofbc,kfl_fixno)
!-----------------------------------------------------------------------
! NAME 
!    codfix
! DESCRIPTION
!    This routine codes fixno
! USES
! USED BY
!    
!***
!-----------------------------------------------------------------------
  use typre
  implicit none
  integer(ip), intent(in)    :: ndofbc
  integer(ip), intent(inout) :: kfl_fixno(ndofbc)
  integer(ip)                :: icode,idofbc
  real(rp)                   :: rcode

  icode = kfl_fixno(1)
  do idofbc = 1,ndofbc
     rcode = icode/10**(ndofbc-idofbc)
     kfl_fixno(idofbc) = int(rcode)
     icode = icode - kfl_fixno(idofbc)*10**(ndofbc-idofbc)
  end do

end subroutine codfix


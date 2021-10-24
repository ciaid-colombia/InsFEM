subroutine lmn_outerr(a)
!------------------------------------------------------------------------
! NAME 
!    lmn_outerr
! DESCRIPTION
!    This routine checks if there are errors and warnings
!------------------------------------------------------------------------
  use Mod_int2str
  use typre
  use Mod_LowMach
  implicit none
  class(LowMachProblem) :: a
  integer(ip) :: ierro=0
  
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

!Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine lmn_outerr
  


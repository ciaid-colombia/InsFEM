subroutine plcd_outerr(a)
!------------------------------------------------------------------------
! NAME 
!    plcd_outerr
! DESCRIPTION
!    This routine checks if there are errors and warnings
!------------------------------------------------------------------------
  use typre
  use Mod_PLCD
  use Mod_int2str
  implicit none
  integer(ip) :: ierro=0,iwarn=0

  class(PLCDProblem) :: a
  
  
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

!Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine plcd_outerr
  


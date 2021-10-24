subroutine sld_outerr(a)
!------------------------------------------------------------------------
! NAME 
!    sld_outerr
! DESCRIPTION
!    This routine checks if there are errors and warnings
!------------------------------------------------------------------------
  use typre
  use Mod_int2str
  use Mod_Solids
  implicit none
  integer(ip) :: ierro=0,iwarn=0

  class(SolidsProblem) :: a
  
  
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

!Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine sld_outerr
  


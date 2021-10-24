subroutine nsc_outerr(a)
  use typre
  use Mod_NSCompressible
  use Mod_int2str
  implicit none
  integer(ip) :: ierro=0,iwarn=0

  class(NSCompressibleProblem) :: a
  
  ! Write Errors and Warnings
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

!Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine nsc_outerr
  

subroutine lev_outerr(a)
   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_int2str
   implicit none
   class(LevelSetProblem) :: a
   
   integer(ip) :: ierro = 0


   !Write Errors and Warnings
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

! Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine

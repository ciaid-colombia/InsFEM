subroutine php_closedatafi(a)
   use MPI
   use typre
   use def_parame
   use Mod_iofile
   use Mod_PhysicalProblem   
   class(PhysicalProblem) :: a
   character(150) :: outstr
   
   if (a%MPIrank == a%MPIroot) then
       !Output
       outstr = adjustl(trim(a%exmod))//'_CloseDataFile'
       call iofile(two,a%lun_pdata,' ','outstr','old')
   endif

end subroutine
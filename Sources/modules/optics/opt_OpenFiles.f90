subroutine opt_OpenFiles(a)
   use typre
   use def_parame
   use Mod_Optics
   use Mod_int2str
   use Mod_iofile
   implicit none
   class(OpticsProblem) :: a
   
   character(150) :: fil_outpu
   
   interface
      subroutine php_openfi(a)
         use typre
         use Mod_PhysicalProblem
         class(PhysicalProblem) :: a
      end subroutine
   end interface   
   
   !Open the usual files for a physical problem
   call php_openfi(a)
   
   !Now Compute and write down the various optical parameters
   if (a%MPIrank == a%MPIroot) then
      fil_outpu = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.res'
      call iofile(zero,a%lun_outres,fil_outpu,adjustl(trim(a%exmod))//' OUTPUT')
   endif
   
end subroutine

subroutine opt_CloseFiles(a)
   use typre
   use Mod_Optics
   use def_parame
   use Mod_iofile
   implicit none
   class(OpticsProblem) :: a
   
   character(150) :: fil_outpu
   
   interface
      subroutine php_closefi(a)
         use typre
         use Mod_PhysicalProblem
         class(PhysicalProblem) :: a
      end subroutine
   end interface   
   
   !Open the usual files for a physical problem
   call php_closefi(a)
   
   !Now Compute and write down the various optical parameters
   if (a%MPIrank == a%MPIroot) then
       call iofile(two,a%lun_outres,fil_outpu,adjustl(trim(a%exmod))//' OUTPUT')
   endif

end subroutine
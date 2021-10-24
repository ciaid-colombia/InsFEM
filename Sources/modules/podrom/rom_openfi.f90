subroutine rom_openfi(a)
   use typre
   use Mod_iofile
   use Mod_int2str
   use def_parame
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   character(150) :: fil_pdata,fil_outro,fil_solve
   character(3) :: exmod

   exmod = a%exmod
   if (a%MPIrank == a%MPIroot) then
      fil_pdata = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.rom.dat'
      call iofile(zero,a%lun_pdata,fil_pdata,'ROM DATA','old')
      fil_outro = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.rom.log'
      call iofile(zero,a%lun_outro,fil_outro,'ROM OUTPUT')
      fil_solve = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.rom.sol'
      call iofile(zero,a%lun_solro,fil_solve,'ROM OUTPUT')
   endif

end subroutine 

subroutine rom_closefi(a)
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a

   if (a%MPIrank == a%MPIroot) then
      close(a%lun_pdata)
      close(a%lun_outro)
      close(a%lun_solro)
   endif   
end subroutine

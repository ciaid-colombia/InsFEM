subroutine php_openfi(a)
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   
   class(PhysicalProblem) :: a
   
   character(150) :: fil_pdata,fil_outpu,fil_cpcon,fil_conve,fil_nolin,fil_solve,fil_adapt,fil_error,fil_ersgs
   
   if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
      fil_pdata = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(adjustl(trim(a%exmod))))//'.dat'
      
      call iofile(zero,a%lun_pdata,fil_pdata,adjustl(trim(a%exmod))//' DATA','old')
   endif
   
   if (a%MPIrank == a%MPIroot) then
      fil_outpu = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.log'
      fil_conve = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.cvg'
      fil_nolin = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.nli'
      fil_solve = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.sol'
      fil_cpcon = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.cp.cvg'
      fil_adapt = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.adp'
      fil_error = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.err'
      fil_ersgs = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.sgs'

      
      call iofile(zero,a%lun_outph,fil_outpu,adjustl(trim(a%exmod))//' OUTPUT')
      call iofile(zero,a%lun_conve,fil_conve,adjustl(trim(a%exmod))//' CONVERGENCE')
      call iofile(zero,a%lun_nolin,fil_nolin,adjustl(trim(a%exmod))//'NON LINEAR INFO')
      call iofile(zero,a%lun_solve,fil_solve,adjustl(trim(a%exmod))//' SOLVER')
      if(a%kfl_docoupconv) call iofile(zero,a%lun_cpconve,fil_cpcon,adjustl(trim(a%exmod))//'CP CONVERGENCE') 
      if(a%kfl_adap) call iofile(zero,a%lun_adapt,fil_adapt,adjustl(trim(a%exmod))//'ADAPTIVE') 
      call iofile(zero,a%lun_error,fil_error,adjustl(trim(a%exmod))//'ERROR')
      call iofile(zero,a%lun_ersgs,fil_ersgs,adjustl(trim(a%exmod))//'SGS ERROR')
   endif

end subroutine 

subroutine php_closefi(a)
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a

   if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
      close(a%lun_pdata)
      close(a%lun_outph)
      close(a%lun_conve)
      close(a%lun_nolin)
      close(a%lun_solve)
      if(a%kfl_docoupconv) close(a%lun_cpconve)
      if(a%kfl_adap) close(a%lun_adapt)
      close(a%lun_error)
      close(a%lun_ersgs)
   endif   
   
end subroutine

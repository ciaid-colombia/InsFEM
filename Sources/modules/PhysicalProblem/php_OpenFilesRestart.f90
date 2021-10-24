subroutine OpenFilesRestart(a,itask)
   !------------------------------------------------------------------------
   !    This routine opens the restart files
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip), intent(in) :: itask
   
   character(150) :: fil_rstar
   
   call a%Timer%Total%Tic
   call a%Timer%Restar%Tic

   select case (itask)

   case (1) !Open for reading

      if (a%kfl_iofor == 0) then
         a%fil_rstar = trim(a%OldRestartFolder)//'/'//adjustl(trim(a%oldnamda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.rst'
         call iofile(zero,a%lun_rstar,a%fil_rstar,'RESTART','old','unformatted')
         call iofile(zero,a%lun_rsta2,adjustl(trim(a%fil_rstar))//'2','RESTAR2','old','formatted')
      elseif (a%kfl_iofor == 1) then
         a%fil_rstar = trim(a%OldRestartFolder)//'/'//adjustl(trim(a%oldnamda))//'.'//adjustl(trim(a%exmod))//'.rst'
         if (a%MPIrank == a%MPIroot) call iofile(zero,a%lun_rstar,a%fil_rstar,'RESTART','old','unformatted')
      end if

   case (2) !Open for writing

      if (a%kfl_iofor == 0) then
         a%fil_rstar = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.rst'
         call iofile(zero,a%lun_rstar,a%fil_rstar,'RESTART','replace','unformatted')
         call iofile(zero,a%lun_rsta2,adjustl(trim(a%fil_rstar))//'2','RESTAR2','replace','formatted')
      elseif (a%kfl_iofor == 1) then
         a%fil_rstar = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(a%exmod))//'.rst'
         if (a%MPIrank == a%MPIroot) call iofile(zero,a%lun_rstar,a%fil_rstar,'RESTART','replace','unformatted')
      end if

   end select

end subroutine OpenFilesRestart

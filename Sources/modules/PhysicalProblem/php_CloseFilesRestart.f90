subroutine CloseFilesRestart(a,itask)
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

   case (1)
  
   !Close for reading
      if (a%kfl_iofor == 0) then
         call iofile(two,a%lun_rstar,a%fil_rstar,'RESTART','old','unformatted')
         call iofile(two,a%lun_rsta2,adjustl(trim(a%fil_rstar))//'2','RESTAR2','old','formatted')
      elseif (a%kfl_iofor == 1) then
         if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_rstar,a%fil_rstar,'RESTART','old','unformatted')
      end if

   case (2)
  
   !Close for writing
      if (a%kfl_iofor == 0) then
         call iofile(two,a%lun_rstar,fil_rstar,'RESTART','replace','unformatted')
         call iofile(two,a%lun_rsta2,adjustl(trim(fil_rstar))//'2','RESTAR2','replace','formatted')
      elseif (a%kfl_iofor == 1) then
         if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_rstar,a%fil_rstar,'RESTART','replace','unformatted')
      end if

   end select

end subroutine CloseFilesRestart

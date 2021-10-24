subroutine lev_output(a,itask)
!-----------------------------------------------------------------------
!    
! End of a Levelset time step 
! 
! itask = 0  When timemarching is true. There is output or post-process
!            of results if required.
! itask = 1  When timemarching is false. Output and/or post-process of
!            results is forced if they have not been written previously.
! 
!-----------------------------------------------------------------------

   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_Int2str
   use Mod_Iofile
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask
  
   ! integer(ip) :: itime
   integer(ip), save :: dopost(25)
   character(150) :: fil_outpu
   integer(ip) :: igauge
   
   
   
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
   
   !Do the actual postprocess
   !Level
   if (dopost(1) == 1) then
         call a%FilePostpr%postpr(a%level(:,1),'LEVEL',a%istep,a%ctime,a%Mesh)
   endif
   
  
   if (itask ==0 ) then
      !HeightGauges
      if (a%nHeightGauges >= 1) then
         !Only root postprocesses
         if (a%MPIrank == a%MPIroot) then
            if (a%istep == 1) then
               fil_outpu = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.height'
               call iofile(zero,a%lun_outHeight,fil_outpu,adjustl(trim(a%exmod))//' Height')
               
               write(a%lun_outHeight,*) 'Height graph'
               write(a%lun_outHeight,*) 'Time, HeightGauge1, HeightGauge2, ...'
            endif
            
            write(a%lun_outHeight,*) a%ctime, (a%HeightGauges(igauge)%height,igauge = 1,a%nHeightGauges)
            
            !Flush if necessary
            if (a%kfl_flush == 1) call flush(a%lun_outHeight)
         endif
      endif
    endif



end subroutine

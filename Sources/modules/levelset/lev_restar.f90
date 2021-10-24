subroutine lev_restar(a,itask)
   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_int2str
   use Mod_postpr
   use Mod_iofile
   use Mod_TimeIntegrator
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: icomp             ! component of the unknown vector
   integer(ip) :: ioerr             ! input/output error
   real(rp)    :: rdummy
   
   integer(ip) :: kfl_twost_pre
   integer(ip) :: lun_rstar,lun_rsta2
   integer(ip) :: ndime,npoin,pgaus,idime,ielem,igaus,ipoin,nelem,pnode
   character(150) :: fil_rstar
   integer(ip) :: iaux,iaux2,iaux3
   
  type(TimeIntegratorDt1) :: Integrator
  integer(ip)             :: nsteps,kfl_twost
   
   call a%Timer%Total%Tic
   call a%Timer%Restar%Tic
   
   kfl_twost_pre = 0_ip
   
   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   kfl_twost = nsteps-2
   
   select case (itask)
   
   case (1)
      fil_rstar = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.rst'
      if (a%MPIrank == a%MPIroot) call iofile(zero,lun_rstar,fil_rstar,'RESTART','old','unformatted')
      if (a%MPIrank == a%MPIroot) call iofile(zero,lun_rsta2,adjustl(trim(fil_rstar))//'2','RESTAR2','old','formatted')

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)
   
      ! Read the previous time scheme used: one integer number
      read(lun_rsta2,*) kfl_twost_pre
      
      if (kfl_twost_pre==0) then
         if (kfl_twost==0) then
            a%neule = 0_ip
         elseif(kfl_twost==1) then
            a%neule = 1_ip
         else
            call runend('restar: kfl_twost_pre not defined')
         endif
      elseif (kfl_twost_pre==1) then
         a%neule = 0_ip
      else
         call runend('restar: kfl_twost_pre not defined')
      endif
      
      ! first at previous time step (always one)
      icomp = 3
      read(lun_rstar) (a%level(ipoin,icomp),ipoin=1,npoin)

      ! then at pre previous time step
      icomp = 4
      if(kfl_twost>=1 .and. kfl_twost_pre == 1) then
         read(lun_rstar) (a%level(ipoin,icomp),ipoin=1,npoin)
      elseif (kfl_twost_pre == 1) then
         read(lun_rstar) (rdummy,ipoin=1,npoin)
      endif
      ! then at pre pre previous time step
      icomp = 5
      if(kfl_twost>=2) then
         call runend('restar: kfl_twost not defined')
      endif

      if (a%MPIrank == a%MPIroot) call iofile(two,lun_rstar,fil_rstar,'RESTART')
      
      ! Postprocess a%velocity and a%pressure at t=-dt
      if(kfl_twost>=1 .and. a%neule == 0) then
         icomp = 4
         call a%FilePostpr%postpr(a%level(:,icomp),'LEVEL',a%istep,a%ctime-a%dtime,a%Mesh)
      endif
      ! Postprocess a%velocity and a%pressure at t=0
      icomp = 3
      call a%FilePostpr%postpr(a%level(:,icomp),'LEVEL',a%istep,a%ctime,a%Mesh)
      
      ! Actualize Velocity an a%pressure at previous iteration and previous time step (guesses)
      icomp = 3
      a%level(:,1) = a%level(:,icomp)
      a%level(:,2) = a%level(:,icomp)      
      
      if (a%MPIrank == a%MPIroot) call iofile(two,lun_rsta2,adjustl(trim(fil_rstar))//'2','RESTAR2')
      
   case (2)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)
   
      fil_rstar = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.'//adjustl(trim(a%exmod))//'.rst'
      if (a%MPIrank == a%MPIroot) call iofile(zero,lun_rstar,fil_rstar,'RESTART','replace','unformatted') ! open file
      if (a%MPIrank == a%MPIroot) call iofile(zero,lun_rsta2,adjustl(trim(fil_rstar))//'2','RESTAR2','replace','formatted')
      
      ! write current time scheme
      write(lun_rsta2,*) kfl_twost
      
      ! t = -a%dtime
      icomp = 3
      write(lun_rstar) (a%level(ipoin,icomp),ipoin=1,npoin)

      
      ! t = -2*a%dtime
      icomp = 4
      if(kfl_twost>=1) then
         write(lun_rstar) (a%level(ipoin,icomp),ipoin=1,npoin)

      endif
      
      ! t = -3*a%dtime
      icomp = 5
      if(kfl_twost>=2) then
         call runend('restar: kfl_twost not defined')
      endif     
      

      if (a%MPIrank == a%MPIroot) call iofile(two,lun_rstar,fil_rstar,'RESTAR') !close file
      if (a%MPIrank == a%MPIroot) call iofile(two,lun_rsta2,adjustl(trim(fil_rstar))//'2','RESTAR2')
   end select
   
   call a%Timer%Total%Toc
   call a%Timer%Restar%Toc

100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)
   
end subroutine

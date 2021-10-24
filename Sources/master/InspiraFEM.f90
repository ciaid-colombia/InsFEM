program InspiraFEM
!-----------------------------------------------------------------------!
!                                                                       !
!     InspiraFEM                                                        !
!                                                                       !
!-----------------------------------------------------------------------!

   !This is the main routine of the program. It controls the main flow 
   !by calling routines that drive specific tasks 
   !Each of these routines, in turn, call modules to perform the given task.
   use def_master
   implicit none 
   
   integer(ip) :: kfl_gotim
   
   !InitializeMPI
   call InitializeMPI

   !Global Timer Initialize
   call cpu_total%Tic()

   !Read the problem data
   call Reapro
 
   !Turnon problems
   call Turnon 
   
   !Match the channels between Cases, Drivers, Channels
   !Also match the interpolators
   call MatchChannels
   
   !Setup the interpolators
   call SetupInterpolators
   
   !Initially setup channels
   call UpdateChannels

   kfl_gotim = 1
   time: do while (kfl_gotim /= 0)
  
      call DoTimeStep(kfl_gotim)
    
   end do time

   call Turnof
   
   
   
end program InspiraFEM

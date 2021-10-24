subroutine tem_output(a,itask)
!-----------------------------------------------------------------------
!    
! End of a temperature time step 
! 
! itask = 0  When timemarching is true. There is output or post-process
!            of results if required.
! itask = 1  When timemarching is false. Output and/or post-process of
!            results is forced if they have not been written previously.
! 
!-----------------------------------------------------------------------

   use typre
   use Mod_Temperature
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: itask
   integer(ip), save :: dopost(25)
   integer(ip) :: itime
 
   select case(itask)
   !End of a time step.
   case(zero)
      
      a%pos_alrea=0

      if(a%kfl_exacs/=0) then
         call runend('tem_output: kfl_exacs not ready')
         !call tem_exaerr(one)
      endif      

      ! Tracking of points.
      if(a%nptra > 0 ) then
         call runend('tem_output: tracking of points not ready')
         !call tem_outtpo
      end if
   
   !End of the run.
   case(one)

   
   end select
   
   
   !Do the actual postprocess
   
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
   !Temperature
   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%tempe(:,1),'Temperature',a%istep,a%ctime,a%Mesh)
   endif
   
   
   !Dissipation
   if (dopost(4) == 1) then
      call a%FilePostpr%postpr(a%Dissipation,'DISSI_TEMPE',a%istep,a%ctime,a%Mesh)
   endif
   
   !Averaged temperature along 1 dimension
   if (dopost(5) == 1) then
        call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%tempe(:,1),'Averaged_Tempe',a%istep,a%ctime,a%Mesh,a%Memor)
   endif
   
   !Averaged Dissipation along 1 dimension
   if (dopost(6) == 1) then
       call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%dissipation(:),'Averaged_Dissipation_TEM',a%istep,a%ctime,a%Mesh,a%Memor)
   endif
   
   !Temperature Gradient
   if (dopost(7) == 1) then
         call a%FilePostpr%postpr(a%gradient,'GRADIENT',a%istep,a%ctime,a%Mesh)
         if (allocated(a%GradientGaussPoints)) then
            call a%FilePostpr%pgpr2p(a%GradientGaussPoints,'GRADIENT_GAUSS',a%istep,a%ctime,a%Mesh)
         endif  
   endif
   
   !ShockCapturingViscosity
   if (dopost(8) == 1) then
       call a%FilePostpr%pgpr1p(a%ShockCapturingViscosity,'SHOCK_VISCO',a%istep,a%ctime,a%Mesh)
   endif
   
  !Sigma Term   sigma:GradSym(u)
   if (dopost(9) == 1) then
       if (a%kfl_CouplingThreeField==1) call a%FilePostpr%postgp(a%sigmatermarray(:),'Sigma_term',a%istep,a%ctime,a%Mesh)
       !if (a%kfl_CouplingThreeField==1) call a%FilePostpr%postpr(a%ViscousDissipation(:),'Viscous_Dis',a%istep,a%ctime,a%Mesh)
   endif 

   if (dopost(10) == 1) then
      call a%FilePostpr%postgp(a%tesgs,'Temperature SGS',a%istep,a%ctime,a%Mesh,'SCALAR')
   end if
end subroutine

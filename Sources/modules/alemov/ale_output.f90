subroutine ale_output(a,itask)
   !-----------------------------------------------------------------------
   !> End of an ALE time step 
   !! itask = 0  When timemarching is true. There is output or post-process
   !!            of results if required.
   !! itask = 1  When timemarching is false. Output and/or post-process of
   !!            results is forced if they have not been written previously.
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_Alemov
   use Mod_sldAlemov
   implicit none
   class(AlemovProblem) :: a
   integer(ip) :: itask
  
   integer(ip) :: itime
   integer(ip), save  :: dopost(25)
  
   select case(itask)

   !End of a time step.
   case(zero)
      
   !End of the run.
      !Tracking of points.
      if(a%nptra > 0) then
         call a%ale_outtpo()
      end if
   case(one)

      !Tracking of points.
      if(a%nptra > 0) then
         call a%ale_outtpo()
      end if

   end select

 
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)

   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%Displacement(:,:,1),'ALEdispl',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%Velocity(:,:,1),'ALEveloc',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(2) == 1) then
       select type (a)
       type is (sldAlemovProblem)
           if(a%kfl_printPrincipalStresses) then
               call a%FilePostpr%postgp(a%solid%printPStrain,'PStrain',a%istep,a%ctime,a%Mesh,'voigtT')
               call a%FilePostpr%postgp(a%solid%stress_g,'Stress',a%istep,a%ctime,a%Mesh,'voigtT')
               call a%FilePostpr%postgp(a%solid%strain_g,'Strain',a%istep,a%ctime,a%Mesh,'voigtT')
           end if
       end select
   end if

end subroutine

subroutine sld_output(a,itask)
    !-----------------------------------------------------------------------
    !****f* SOLIDS/sld_output
    ! NAME 
    !    sld_output
    ! DESCRIPTION
    !    End of a SOLIDS time step 
    !    itask = 0  When timemarching is true. There is output or post-process
    !               of results if required.
    !    itask = 1  When timemarching is false. Output and/or post-process of
    !               results is forced if they have not been written previously.
    ! USES
    !    output
    !    a%FilePostpr%postpr
    ! USED BY
    !    sld_endste (itask=0)
    !    sld_turnof (itask=1)
    !***
    !-----------------------------------------------------------------------
    use typre
    use Mod_Postpr
    use Mod_Solids
    implicit none
    class(SolidsProblem) :: a
    integer(ip)                :: itask
    integer(ip), save       :: dopost(49)
    integer(ip)             :: itime,istat,ifiel,ndime,aux4
    real(rp)                :: dummr
    
    select case(itask)

    !End of a time step.
   case(0)
      a%pos_alrea=0
     !------------------------------------------------------------------------- 
     ! Postprocess at a given time
   
      !Tracking of points.
      if(a%nptra > 0) then
         call a%PointTracking
      end if

   case(1)

   end select

   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
   !Do the actual postprocess
   !Displacements
   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%disp(:,:,1),'Disp',a%istep,a%ctime,a%Mesh)
   end if

   !Stress tensor
   if (dopost(2) == 1) then
       if(a%kfl_printStress .eqv. .true.) then
           call a%FilePostpr%postgp(a%stress_g,'Gauss Stress',a%istep,a%ctime,a%Mesh,'voigtT')
       end if
   end if

  !Strain
  if (dopost(21) == 1) then
      if(a%kfl_NodalStress.eqv. .true.) then
          if(a%kfl_printStress.eqv. .true.) then
              call a%FilePostpr%postpr(a%stress(:,:),'Nodal Stress',a%istep,a%ctime,a%Mesh,'voigtT')
          endif
      endif
  endif

   !Strain tensor
   if (dopost(3) == 1) then
       if(a%kfl_printStrain .eqv. .true.) then
           call a%FilePostpr%postgp(a%strain_g,'Gauss Strain',a%istep,a%ctime,a%Mesh,'voigtT')
       end if
   end if

  !Strain
  if (dopost(22) == 1) then
      if(a%kfl_NodalStrain .eqv. .true.) then
          if(a%kfl_printStrain.eqv. .true.) then
              call a%FilePostpr%postpr(a%strain(:,:),'Nodal Strain',a%istep,a%ctime,a%Mesh,'voigtT')
          endif
      endif
  endif

   !Internal tractions
   if (dopost(4) == 1) then
      call a%FilePostpr%postpr(a%btraction_nodal,'Sld_traction',a%istep,a%ctime,a%Mesh)
   end if

   !Fluid traction - nodal
   if (dopost(5) == 1) then
         if(associated(a%trac)) then
             call a%FilePostpr%postpr(-a%trac,'Fld_traction',a%istep,a%ctime,a%Mesh)
         end if
   endif

   !Total external traction - Gaussian
   if (dopost(6) == 1) then
       call a%FilePostpr%postgp(a%btraction,'Btraction',a%istep,a%ctime,a%Mesh)
   end if

   !Velocity
   if (dopost(7) == 1) then
       if(a%kfl_timei==1) then
           call a%FilePostpr%postpr(a%veloc(:,:,1),'Velocity',a%istep,a%ctime,a%Mesh)
       end if
   end if

   !Principal Strain tensor - generally auto used by sldALE to deform mesh
   if (dopost(8) == 1) then
       if(a%kfl_printPrincipalStresses) then
           call a%FilePostpr%postgp(a%printPStrain,'PpalStrain',a%istep,a%ctime,a%Mesh,'voigtT')
       end if
   endif

   !Acceleration
   if (dopost(9) == 1) then
       if(a%kfl_timei==1) then
           call a%FilePostpr%postpr(a%accel(:,:,1),'Acceleration',a%istep,a%ctime,a%Mesh)
       end if
   end if

   !Sigma
   if (dopost(11) == 1) then
       call a%FilePostpr%postpr(a%sigma(:,:,1),'Sigma',a%istep,a%ctime,a%Mesh,'voigtT')
   endif

   !Pressure
   if (dopost(12) == 1) then
       call a%FilePostpr%postpr(a%press(:,1),'Press',a%istep,a%ctime,a%Mesh)
   end if

   !Stress tensor
   if (dopost(20) == 1) then
        if(a%kfl_printGaussSigma) then
           call a%FilePostpr%postgp(a%sigma_g,'Gauss Sigma',a%istep,a%ctime,a%Mesh,'voigtT')
       endif
   end if

   if (dopost(23) == 1) then
        if(a%kfl_printGaussPress) then
           call a%FilePostpr%postgp(a%press_g,'Gauss Press',a%istep,a%ctime,a%Mesh)
       endif
   end if


   call a%SolidSpecificOutput(itask)

end subroutine sld_output

subroutine sldup_output(a,itask)
    !-----------------------------------------------------------------------
    !****f* SOLIDS/sldup_output
    ! NAME 
    !    sldup_output
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
    !    sldsup_endste (itask=0)
    !    sldsup_turnof (itask=1)
    !***
    !-----------------------------------------------------------------------
    use typre
    use Mod_Postpr
    use Mod_UPSolids
    implicit none
    class(UPSolidsProblem):: a
    integer(ip)            :: itask
    integer(ip), save      :: dopost(49)
    
    select case(itask)

    !End of a time step.
   case(0)
      a%pos_alrea=0

   case(1)

   end select

  !Decide which postprocesses need to be done
  call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)

  !Stress tensor
  if (dopost(20) == 1) then
       if(a%kfl_printGaussSigma) then
          call a%FilePostpr%postgp(a%sigma_g,'Gauss Sigma',a%istep,a%ctime,a%Mesh,'voigtT')
      endif
  end if


  !J2 stresses
  !if (dopost(15) == 1) then
  !    call a%FilePostpr%postpr(a%J2(:),'J2 Stresses',a%istep,a%ctime,a%Mesh)
  !endif

  !Residual
  if (dopost(16) == 1) then
      call a%FilePostpr%postgp(a%residualU,'ResidualU',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(a%residualP,'ResidualP',a%istep,a%ctime,a%Mesh)
  end if    

  call a%OutputExtend(itask)

end subroutine sldup_output

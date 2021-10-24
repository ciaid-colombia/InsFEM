subroutine sld_begste(a)
    ! DESCRIPTION
    !    This routine prepares for a new time step of the elastic solid prob
    !-----------------------------------------------------------------------
    use typre
    !use Mod_Postpr
    use Mod_Solids
    implicit none
    class(SolidsProblem) :: a
    integer(ip) :: ndime
    real(rp)    :: dtinv

    call a%Mesh%GetNdime(ndime)

    !call a%FilePostpr%postpr(a%disp(:,:,1),'Disp_ref',a%istep,a%ctime,a%Mesh)
    !if (a%kfl_timei == 1) then
    !    call a%FilePostpr%postpr(a%veloc(:,:,1),'Velocity_ref',a%istep,a%ctime,a%Mesh)
    !    call a%FilePostpr%postpr(a%accel(:,:,1),'Acceleration_ref',a%istep,a%ctime,a%Mesh)
    !endif

    a%force_factor = a%delta_force  

    !if( a%istep == 1) then
    !    !Compute initial accel a_0
    !    call a%Bouope
    !    call a%preElmope
    !end if

    !If dynamic
    if(a%kfl_timei == 1) then

        a%accel(:,:,2) = a%accel(:,:,1)
        a%veloc(:,:,2) = a%veloc(:,:,1)

        if(a%kfl_FSI_restart .eqv. .false.) then

            call a%sld_calcVelandAccel(a%dtinv)
        end if

    end if

    !This flag is only needed for the first step
    !Flag set in sld_restart, otherwise always false
    a%kfl_FSI_restart = .false.

end subroutine sld_begste

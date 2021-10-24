subroutine sld_endste(a,itask)
   !This routine ends a time step of the incompressible Solids equations.
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip) :: itask,icomp
   real(rp) :: dtinv,dtinv2
   
    !Pre-CrankNicolson
    if (itask == 1) then
        !Calculate stresses and strains
        call a%EndElmope('Endste')

        !Post-CrankNicolson
    elseif (itask == 2) then

        call a%UpdateDynamicComponents

        !Boundary operation at the end
        if(a%kfl_isSlave .eqv. .false.) then
            call a%EndBouope
        end if

        if(a%sld_type== 'NONLI' ) then
            !-------------------------------------------------------------------
            !ONLY MOVE mesh after all calculation involving derivatives have been
            !made
            ! !Dealloc and recompute ExtnorLpoty and Vmass
            call a%Mesh%DeallocExnorLpoty
            call a%Mesh%DeallocVmass

            !Recompute
            call a%Mesh%ComputeVmass
            call a%Mesh%ExtnorLpoty
        endif

   end if

end subroutine sld_endste

subroutine sld_updateDynamicComp(a)
    use typre
    use Mod_Solids
    implicit none
    class(SolidsProblem) :: a
    integer(ip)          :: icomp

    !Update a%displacements
    !Higher order components
    if (a%ncomp > 3) then
        do icomp = a%ncomp, 4,-1
            a%disp(:,:,icomp) = a%disp(:,:,icomp-1)
        enddo
    endif

    !Previous time step component
    a%disp(:,:,3) = a%disp(:,:,1)  ! Dn = Dn+1

end subroutine

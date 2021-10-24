subroutine sldup_updateDynamicComp(a)
    use typre
    use Mod_UPSolids
    implicit none
    class(UPSolidsProblem) :: a
    integer(ip) :: icomp

    if (a%ncomp > 3) then
        do icomp = a%ncomp, 4,-1
            a%disp(:,:,icomp)  = a%disp(:,:,icomp-1)
            a%press(:,icomp)   = a%press(:,icomp-1)
        enddo
    endif

    !Previous time step component
    a%disp(:,:,3)  = a%disp(:,:,1)  ! Vn = Vn+1
    a%press(:,3)   = a%press(:,1)  

end subroutine sldup_updateDynamicComp

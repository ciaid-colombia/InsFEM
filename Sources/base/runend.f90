subroutine runend_init(loutput)
    use typre
    use Mod_runend
    implicit none
    integer(ip) :: loutput

    lun_outpu = loutput

end subroutine

!subroutine runend(message,var)
subroutine runend(message)
    use typre
    use MPI
    use Mod_runend
    !-----------------------------------------------------------------------
    !
    ! This routine stops the run and writes the summary of CPU time.
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*) :: message
    logical(lg)  :: lopen, lexist
    integer :: STATUS = 1
    integer(ip)  :: ierr,status_mpi(MPI_STATUS_SIZE)

    !   !Write the table with the CPU times.
    !   call cputab

    !Write message and stop the run.
    if(message/='O.K.!') then
    inquire(unit=lun_outpu,opened=lopen,exist=lexist)
    if(lopen) then
        write(lun_outpu,900)
        write(lun_outpu,901) message
    end if
    !write(*,901) 'AN ERROR HAS BEEN DETECTED: '//message//var
    write(*,901) 'AN ERROR HAS BEEN DETECTED: '//adjustl(trim(message))

    !Flush everything to disk
    call flush

    !call MPI_Finalize(ierr)
    call MPI_ABORT(MPI_COMM_WORLD,status_mpi,ierr)
    !call exit(STATUS)
else
    !Flush everything to disk
    call flush
    inquire(unit=lun_outpu,opened=lopen,exist=lexist)
    if(lopen) then
        write(lun_outpu,'(//,a,/)') '     * * * END OF ANALYSIS * * *'
    end if
    write(*,'(//,a,/)')         '     * * * END OF ANALYSIS * * *'
endif

!
! Formats.
!
900 format(//,5x,34('* '),//,25x,'AN ERROR HAS BEEN DETECTED:',/)
901 format(/,5x,a,/)

end subroutine runend



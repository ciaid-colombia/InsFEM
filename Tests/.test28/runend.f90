

subroutine runend(message)
   use typre
!-----------------------------------------------------------------------
!
! This routine stops the run and writes the summary of CPU time.
!
!-----------------------------------------------------------------------
   implicit none
   character(*) :: message
   logical(lg)  :: lopen
   
!    !Write the table with the CPU times.
!    call cputab
   integer :: STATUS = 1
   
   !Write message and stop the run.
   if(message/='O.K.!') then
	 write(*,901) 'AN ERROR HAS BEEN DETECTED: '//message
   else
      write(*,'(//,a,/)')         '     * * * END OF ANALYSIS * * *'
   end if
 
   !
   ! Formats.
   !
   900 format(//,5x,34('* '),//,25x,'AN ERROR HAS BEEN DETECTED:',/)
   901 format(/,5x,a,/)

end subroutine runend



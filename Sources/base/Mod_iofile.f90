module Mod_iofile
  use typre
  use def_parame
  
  !Counter for file units
  integer(ip) :: funit = 10
  
contains 
  
  subroutine iofile(itask,nunit,files,messa,stato,formo,posio,ioerr2)
  !------------------------------------------------------------------------
  !
  ! This routine opens/closes logical units and treat errors
  !
  !------------------------------------------------------------------------
    implicit none
    integer(ip),   intent(in)    :: itask 
    integer(ip),   intent(inout) :: nunit
    character*(*), intent(in)    :: messa,files
    character*(*), optional      :: stato,formo,posio
    integer(ip), optional        :: ioerr2
    integer                      :: ioerr
    character(7)                 :: statu
    character(6)                 :: posit
    character(11)                :: forma

    select case (itask)
!
! Open unit
!
    case(zero)

       ! Optional arguments
       if(present(stato)) then
          statu=stato
       else
          statu='unknown'
       end if
       if(present(formo)) then
          forma=formo
       else
          forma='formatted'
       end if
       if(present(posio)) then
          posit=posio
       else
          posit='asis'
       end if

       funit = funit+1
       nunit = funit
       
       !open(nunit,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr,position=posit,access='STREAM')
       open(nunit,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr,position=posit,access='SEQUENTIAL')
       if(ioerr/=0) then
         if (present(ioerr2)) then
            !error flag
            ioerr2 = -1
         else
            write(*,*) files
            write(*,*)'Error :', ioerr
            call runend('ERROR WHEN OPENING THE '//trim(messa)//' FILE')
         endif
      endif
!
! Close unit
!
    case(two)
       close(nunit,iostat=ioerr)
       if(ioerr/=0) then
             call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE')
       end if
    end select
!
! Format
!
101 format(5x,'WARNING: ',a)
    
  end subroutine iofile

end module Mod_iofile

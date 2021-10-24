subroutine OpenFileDom(DataFolder,namda,lun_pdata_dom)
  use   typre
  use   Mod_iofile
  implicit none

  character(150), intent(in) :: DataFolder,namda
  integer(ip), intent(inout)    :: lun_pdata_dom
  character(150) :: fil_pdata_dom, fil_outpu

  fil_pdata_dom = trim(DataFolder)//'/'//adjustl(trim(namda))//'.dom.dat'
  call iofile(zero,lun_pdata_dom,fil_pdata_dom,fil_pdata_dom,'old')
  
end subroutine OpenFileDom

subroutine CloseFileDom(lun_pdata_dom)
  use typre
  use Mod_iofile
  implicit none
  
  integer(ip), intent(inout)    :: lun_pdata_dom
  call iofile(two,lun_pdata_dom,'','DOMAIN','old')

end subroutine CloseFileDom

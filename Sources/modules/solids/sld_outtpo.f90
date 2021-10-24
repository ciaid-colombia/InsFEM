subroutine sld_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the SLD problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_Solids
   implicit none

   class(SolidsProblem) :: a

   integer(ip) :: ndime,npoin
   real(rp)    :: disp(3_ip,a%nptra)
   real(rp)    :: veloc(3_ip,a%nptra)
   real(rp)    :: accel(3_ip,a%nptra)

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   if(a%nptra>0) then
      call a%TrackingInterpolator%Interpolate(ndime,a%disp(:,:,1),disp(1:ndime,:))
      !If dynamic
      if(a%kfl_timei == 1) then
          call a%TrackingInterpolator%Interpolate(ndime,a%veloc(:,:,1),veloc(1:ndime,:))
          call a%TrackingInterpolator%Interpolate(ndime,a%accel(:,:,1),accel(1:ndime,:))
      endif

      if (a%MPIrank == a%MPIroot) then
          !If dynamic
          if(a%kfl_timei == 1) then
              write(a%lun_trap,20)(a%ctime),(disp(1:ndime,1:a%nptra)),(veloc(1:ndime,1:a%nptra)),(accel(1:ndime,1:a%nptra))
          else
              write(a%lun_trap,20)(a%ctime),(disp(1:ndime,1:a%nptra))
          endif
          if (a%kfl_flush == 1) call flush(a%lun_trap)
      endif  
  endif

  
   10 format(1x,'Nodal Point number:  ',10(1x,i12  ))
   11 format(1x,'X-Coordinate:        ',10(1x,e14.7))
   12 format(1x,'Y-Coordinate:        ',10(1x,e14.7))
   13 format(1x,'Z-Coordinate:        ',10(1x,e14.7))
   15 format(1x,/, &
             1x,'Time  Point1  Point2 ...',/, &
             1x,'------------------------')
   20 format(*(1x,e14.7))
  return

end subroutine sld_outtpo


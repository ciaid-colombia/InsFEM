subroutine ale_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the SLD problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_Alemov
   implicit none

   class(AlemovProblem) :: a
   real(rp)    :: disp(3_ip,a%nptra)
   real(rp)    :: veloc(3_ip,a%nptra)
   integer(ip) :: ndime,npoin

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   if(a%nptra>0) then

      call a%TrackingInterpolator%Interpolate(ndime,a%Displacement(:,:,1),disp(1:ndime,:))
      call a%TrackingInterpolator%Interpolate(ndime,a%Velocity(:,:,1),veloc(1:ndime,:))

      if (a%MPIrank == a%MPIroot) then

          write(a%lun_trap,20)(a%ctime),(disp(1:ndime,1:a%nptra)),(veloc(1:ndime,1:a%nptra))
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

end subroutine ale_outtpo


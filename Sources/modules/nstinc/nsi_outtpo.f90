subroutine nsi_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the NSI problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   implicit none

   class(NavierStokesProblem) :: a

   integer(ip) :: ndime,npoin
   real(rp) :: veloc(3_ip,a%nptra), press(1_ip,a%nptra)
   real(rp), allocatable :: NSpress(:,:)

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   call a%Memor%alloc(1,npoin,NSpress,'NSpress','nsi_outtpo')
   NSpress(1,:) = a%press(:,1)

   if(a%nptra>0) then
      call a%TrackingInterpolator%Interpolate(ndime,a%veloc(:,:,1),veloc(1:ndime,:))
      call a%TrackingInterpolator%Interpolate(1,NSpress,press)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_trap,20)(a%ctime),(veloc(1:ndime,1:a%nptra)),(press(1_ip,1:a%nptra))
         if (a%kfl_flush == 1) call flush(a%lun_trap)
      endif  
   endif

   call a%Memor%dealloc(1,npoin,NSpress,'NSpress','nsi_outtpo')
  
   10 format(1x,'Nodal Point number:  ',10(1x,i12  ))
   11 format(1x,'X-Coordinate:        ',10(1x,e14.7))
   12 format(1x,'Y-Coordinate:        ',10(1x,e14.7))
   13 format(1x,'Z-Coordinate:        ',10(1x,e14.7))
   15 format(1x,/, &
             1x,'Time  Point1  Point2 ...',/, &
             1x,'------------------------')
   20 format(*(1x,e14.7))
  return

end subroutine nsi_outtpo


subroutine sup_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the SLD problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_ThreeField
   implicit none

   class(ThreeFieldNSProblem) :: a

   integer(ip) :: ndime,npoin,tn
   real(rp) :: veloc(3_ip,a%nptra)
   real(rp) :: press(1_ip,a%nptra)
   real(rp),allocatable :: sigma(:,:)

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   tn = (ndime*(ndime+1))/2

   call a%Memor%alloc(tn,a%nptra,sigma,'sigma'   ,'sup_outtpo')

   if(a%nptra>0) then
      call a%TrackingInterpolator%Interpolate(tn   ,a%sigma(:,:,1),sigma(1:tn,:))
      call a%TrackingInterpolator%Interpolate(ndime,a%veloc(:,:,1),veloc(1:ndime,:))
      call a%TrackingInterpolator%Interpolate(1,a%press(:,1),press)

      if (a%MPIrank == a%MPIroot) then
          !If dynamic
              write(a%lun_trap,20)a%ctime,veloc(1:ndime,1:a%nptra),press(1,1:a%nptra),sigma(1:tn,1:a%nptra)
          if (a%kfl_flush == 1) call flush(a%lun_trap)
      endif  
  endif

  call a%Memor%dealloc(tn,a%nptra,sigma,'sigma'   ,'sup_outtpo')
  
   10 format(1x,'Nodal Point number:  ',10(1x,i12  ))
   11 format(1x,'X-Coordinate:        ',10(1x,e14.7))
   12 format(1x,'Y-Coordinate:        ',10(1x,e14.7))
   13 format(1x,'Z-Coordinate:        ',10(1x,e14.7))
   15 format(1x,/, &
             1x,'Time  Point1  Point2 ...',/, &
             1x,'------------------------')
   20 format(*(1x,e14.7))
  return

end subroutine sup_outtpo


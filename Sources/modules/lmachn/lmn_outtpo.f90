subroutine lmn_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the LowMach problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_LowMach
   implicit none

   class(LowMachProblem) :: a

   integer(ip) :: ndime,npoin
   real(rp) :: tempe(1_ip,a%nptra), veloc(3_ip,a%nptra), press(1_ip,a%nptra)

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)


   if(a%nptra>0) then
      call a%TrackingInterpolator%Interpolate(1,a%tempe(1:npoin,1),tempe)
      call a%TrackingInterpolator%Interpolate(ndime,a%veloc(:,:,1),veloc(1:ndime,:))
      call a%TrackingInterpolator%Interpolate(1,a%press(1:npoin,1),press)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_trap,20)(a%ctime),(tempe(1_ip,1:a%nptra)),(press(1_ip,1:a%nptra)),(veloc(1:ndime,1:a%nptra))
         if (a%kfl_flush == 1) call flush(a%lun_trap)
      endif
  endif
  
   20 format(21(1x,e14.7))
  return

end subroutine lmn_outtpo

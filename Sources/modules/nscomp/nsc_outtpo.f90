subroutine nsc_outtpo(a)
!-----------------------------------------------------------------------
!
! This routine tracks points for the NSC problem.
!
!-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_PhysicalProblem
   use Mod_NSCompressible
   use Mod_MeshInterpolator
   implicit none

   class(NSCompressibleProblem) :: a

   integer(ip) :: ndime,npoin
   real(rp) :: densf(1_ip,a%nptra), veloc(3_ip,a%nptra), press(1_ip,a%nptra)

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)


   if(a%nptra>0) then
      call a%TrackingInterpolator%Interpolate(1,a%densf(1:npoin,1),densf)
      call a%TrackingInterpolator%Interpolate(ndime,a%veloc(:,:,1),veloc(1:ndime,:))
      call a%TrackingInterpolator%Interpolate(1,a%press(1:npoin,1),press)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_trap,20)(a%ctime),(densf(1_ip,1:a%nptra)),(press(1_ip,1:a%nptra)),(veloc(1:ndime,1:a%nptra))
         if (a%kfl_flush == 1) call flush(a%lun_trap)
      endif
  endif
  
   20 format(21(1x,e14.7))
  return

end subroutine nsc_outtpo


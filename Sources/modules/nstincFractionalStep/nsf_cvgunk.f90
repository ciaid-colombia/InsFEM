subroutine nsf_InnerResiduals(a,rinsi,rprnsi)
   !This subroutine performs convergence checks for the fractional step NavierStokes
   use typre
   use Mod_NSFractionalStep
   use Mod_NavierStokes, only : zensi
   use MPI
   use def_parame
   implicit none
   class(NSFractionalStepProblem) :: a
   
    real(rp) :: rinsi,rprnsi
   
   integer(ip) :: ndime, npoinLocal
   
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,a%unkno,a%veloc(1,1,1),rinsi,zensi,a%MPIcomm,a%MPIroot)
   rprnsi = 0.0_rp
   
   
end subroutine

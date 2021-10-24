subroutine Smooth(a,ndofn,array)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   real(rp)    :: array(ndofn,*)
   
   integer(ip) :: ipoin
   
   do ipoin = 1,a%npoinLocal
      array(:,ipoin) = array(:,ipoin)*a%vmass(ipoin)
   enddo
   
   call a%ArrayCommunicator%GhostCommunicate(ndofn,array(:,1:a%npoin))
   
   !Hanging nodes
   if (a%kfl_HangingNodes .eqv. .true.) call a%InterpolateHangingValues(ndofn,array)
   
   if (a%kfl_perio == 1) call a%MasterToSlave(ndofn,array(:,1:a%npoin))
   
end subroutine   
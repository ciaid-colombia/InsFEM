subroutine rom_SolutionToPhysicalProblem(a)
   use typre
   use Mod_PodRom
   use Mod_PhysicalProblem
   implicit none
   class(PodRomProblem) :: a
   integer(ip)          :: npoinLocal,idofn,idofr,ipoin
   real(rp)             :: Solution(a%ndofr)
   integer(ip)          :: kfl_HangingNodes,kfl_perio

   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%EigenSystem%GetSolution(Solution)
   
   a%FOMSolution = 0.0_rp
   do idofr = 1, a%ndofr
      do ipoin=1,npoinLocal
         do idofn =1,a%ndofn
            a%FOMSolution(idofn,ipoin) = a%FOMSolution(idofn,ipoin) + a%Basis(idofn,ipoin,idofr)*Solution(idofr)
         end do
      end do
   end do
   a%FOMSolution = a%FOMSOlution + a%SnapMean
      
   !Ghostcommunicate
   call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%FOMSolution)

   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%MasterToSlave(a%ndofn,a%FOMSolution)

   !HangingNodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(a%ndofn,a%FOMSolution)
   
   call a%Problem%SetUnkno(a%FOMSolution)

end subroutine
      

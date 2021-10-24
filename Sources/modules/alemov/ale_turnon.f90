subroutine ale_turnon(a)
   !-----------------------------------------------------------------------
   !> This routine sets the mesh to ALE. The displacement and velocity pointers
   !! of the mesh are also initialized. 
   !-----------------------------------------------------------------------
   use typre  
   use Mod_Alemov
   implicit none

   class(AlemovProblem)  :: a
   integer(ip):: npoin,ipoin,idime,ndime

   call a%Mesh%GetNpoin(npoin)
  
   call a%Mesh%SetALE(1_ip)
   call a%Mesh%SetDisplacements(a%Displacement)
   call a%Mesh%SetVelocities(a%Velocity)

   a%kfl_fixno0=a%kfl_fixno
   call a%deferredTurnon

end subroutine ale_turnon

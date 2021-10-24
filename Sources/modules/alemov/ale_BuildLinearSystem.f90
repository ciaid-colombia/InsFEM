subroutine ale_BuildLinearSystem(a,currentbvess)
   !-----------------------------------------------------------------------
   !> This routine computes the elemental matrix and RHS for the mesh displacements.
   !-----------------------------------------------------------------------
   use typre
   use Mod_Alemov   
   implicit none
   class(AlemovProblem) :: a
   integer(ip) :: kfl_perio, kfl_HangingNodes
   integer(ip) :: currentbvess
   
   interface
      subroutine ale_elmope (a,currentbvess)
         use typre
         import AlemovProblem
         implicit none
         integer(ip) :: currentbvess
	      class(AlemovProblem) :: a
      end subroutine
   end interface

   ! Initializations
   call a%LinearSystem%ToZero
   
   call ale_elmope(a,currentbvess)
   
   call a%Bouope
   
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(a%ndofn,a%LinearSystem,a%Memor)
   
   !Hanging nodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(a%ndofn,a%LinearSystem,a%Memor)
   
end subroutine ale_BuildLinearSystem

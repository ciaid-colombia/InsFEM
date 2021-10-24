subroutine php_BuildLinearSystem(a)
   !-----------------------------------------------------------------------
   !    This routine computes elemental matrix and RHS.
   !-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: kfl_perio, kfl_HangingNodes

   ! Initializations
   if (a%kfl_elmat_current==1) then
      call a%LinearSystem%ToZero
   else
      call a%LinearSystem%ToZeroRhs
   end if
   
   call a%Elmope
   
   call a%Bouope
   
   call a%PoinOpe
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(a%ndofn,a%LinearSystem,a%Memor)

   !Hanging nodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(a%ndofn,a%LinearSystem,a%Memor)
   
end subroutine php_BuildLinearSystem

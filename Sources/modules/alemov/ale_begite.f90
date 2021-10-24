subroutine ale_begite(a)
   use typre
   use Mod_Alemov
   implicit none
   class(AlemovProblem) :: a
   
   !We do it always because of adaptive mesh refinement requires it to be done here
   !Not so expensive anyway
   
   !The global (case) iteration counter should come here
   !if (a%OutIiter > 1 .or. .true.) call a%Updbcs
   if (a%OutIiter > 1) call a%Updbcs
   
end subroutine

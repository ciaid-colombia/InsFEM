subroutine nsf_restar(a,itask)
   use typre
   use Mod_NSFractionalStep
   implicit none
   class(NSFractionalStepProblem) :: a
   integer(ip), intent(in) :: itask
   
   interface
      subroutine nsi_restar(a,itask)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine
   end interface
   
   call nsi_restar(a,itask)
end subroutine

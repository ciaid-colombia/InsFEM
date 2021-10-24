subroutine sldup_begite(SldProblem)
!-----------------------------------------------------------------------
! NAME 
!    sldsup_begite
! DESCRIPTION
!    This routine starts an internal iteration for the elastic solid prob. 
!-----------------------------------------------------------------------
   use typre
   use Mod_sld_begiteElmope
   implicit none
   class(UPSolidsProblem_NH), target :: SldProblem

   a      => SldProblem
   up     => SldProblem
   up_NH  => SldProblem

   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   up_NH%disp(:,:,1)  = up_NH%disp(:,:,2)
   up_NH%press(:,1)   = up_NH%press(:,2)

end subroutine sldup_begite

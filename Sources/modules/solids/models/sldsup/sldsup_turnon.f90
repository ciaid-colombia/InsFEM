subroutine sldsup_turnon(a)
   use typre
   use Mod_SUPSolids
   implicit none
   class(SUPSolidsProblem) :: a
   integer :: ielem,nelem
   
   call a%Mesh%GetNelem(nelem)

   a%kfl_printNodalSigma = .true.
   a%kfl_printGaussSigma = .false.
   call a%SUPModelTurnon

end subroutine sldsup_turnon

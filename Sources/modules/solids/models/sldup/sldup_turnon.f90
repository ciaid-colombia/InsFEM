subroutine sldup_turnon(a)
   use typre
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   integer :: ielem,nelem
   
   call a%Mesh%GetNelem(nelem)

   call a%UPModelTurnon

end subroutine sldup_turnon

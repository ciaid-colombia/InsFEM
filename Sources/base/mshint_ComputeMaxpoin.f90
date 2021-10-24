subroutine mshint_ComputeMaxpoin(Mesh,maxpoin)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: Mesh
   integer(ip) :: maxpoin
   
   integer(ip) :: ndime
   
   call Mesh%GetNdime(ndime)
   if (ndime == 2) then
      maxpoin = 20
   elseif (ndime == 3) then
      maxpoin = 20
   endif
   
   
   
end subroutine
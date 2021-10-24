subroutine sldsup_SetExmod(a)
   use typre
   use Mod_SUPSolids
   implicit none
   
   class(SUPSolidsProblem) :: a
   
   a%exmod = 'sld'
   a%namod = 'SUPSOLIDS'
   
end subroutine

subroutine sldsup_SetNdofn(a)
   use typre
   use Mod_SUPSolids
   implicit none
   
   class(SUPSolidsProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)

   !             sigma           + disp  + press
   a%ndofn = (ndime*(ndime+1))/2 + ndime + 1
   a%udofn = ndime
   
end subroutine

subroutine sldsup_SetNdofbc(a)
   use typre
   use Mod_SUPSolids
   implicit none
   
   class(SUPSolidsProblem) :: a
   integer(ip) :: ndime,tn
   integer(ip) :: u1,uf,s1,sf,p1,bc
   
   call a%Mesh%GetNdime(ndime)
   call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

   a%ndofbc      = ndime
   a%ndofbcstart = bc
   
end subroutine



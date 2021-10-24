subroutine sldup_SetExmod(a)
   use typre
   use Mod_UPSolids
   implicit none
   
   class(UPSolidsProblem) :: a
   
   a%exmod = 'sld'
   a%namod = 'UPSOLIDS'
   
end subroutine

subroutine sldup_SetNdofn(a)
   use typre
   use Mod_UPSolids
   implicit none
   
   class(UPSolidsProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)

   !         disp  + press
   a%ndofn = ndime + 1
   a%udofn = ndime
   
end subroutine

subroutine sldup_SetNdofbc(a)
   use typre
   use Mod_UPSolids
   implicit none
   
   class(UPSolidsProblem) :: a
   integer(ip) :: ndime
   integer(ip) :: u1,uf,s1,sf,p1,bc
   
   call a%Mesh%GetNdime(ndime)
   call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

   a%ndofbc      = ndime
   a%ndofbcstart = bc
   
end subroutine



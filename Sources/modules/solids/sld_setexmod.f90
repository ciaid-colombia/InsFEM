subroutine sld_SetExmod(a)
   use typre
   use Mod_Solids
   implicit none
   
   class(SolidsProblem) :: a
   
   a%exmod = 'sld'
   a%namod = 'SOLIDS'
   
end subroutine


subroutine sld_SetNdofn(a)
   use typre
   use Mod_Solids
   implicit none
   
   class(SolidsProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = ndime
   a%udofn = ndime
   
end subroutine

subroutine sld_SetNdofbc(a)
   use typre
   use Mod_Solids
   implicit none
   
   class(SolidsProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofbc=ndime
   a%ndofbcstart=0
   
end subroutine



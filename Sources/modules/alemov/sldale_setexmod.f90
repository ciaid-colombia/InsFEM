subroutine sldale_SetNdofn(a)
   use typre
   use Mod_sldAlemov
   implicit none
   
   class(sldAlemovProblem) :: a
   integer(ip) :: ndime
   call a%Mesh%GetNdime(ndime)

   a%ndofn = ndime
end subroutine


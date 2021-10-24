subroutine ale_SetExmod(a)
   use typre
   use Mod_Alemov
   implicit none
   
   class(AlemovProblem) :: a
   
   a%exmod = 'ale'
   a%namod = 'Alemov'
end subroutine


subroutine ale_SetNdofn(a)
   use typre
   use Mod_Alemov
   implicit none
   
   class(AlemovProblem) :: a

   a%ndofn = 1
end subroutine

subroutine ale_SetNdofbc(a)
   use typre
   use Mod_Alemov
   implicit none
   
   class(AlemovProblem) :: a
   integer(ip) :: ndime
   call a%Mesh%GetNdime(ndime)

   a%ndofbc = ndime
   a%ndofbcstart = 0
end subroutine

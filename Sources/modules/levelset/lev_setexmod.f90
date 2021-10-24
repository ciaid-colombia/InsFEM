subroutine lev_SetExmod(a)
   use typre
   use Mod_LevelSet
   implicit none
   
   class(LevelSetProblem) :: a
   
   a%exmod = 'lev'
   a%namod = 'LevSet'
end subroutine


subroutine lev_SetNdofn(a)
   use typre
   use Mod_LevelSet
   implicit none
   
   class(LevelSetProblem) :: a
   
   a%ndofn = 1
end subroutine

subroutine lev_SetNdofbc(a)
   use typre
   use Mod_LevelSet
   implicit none
   
   class(LevelSetProblem) :: a
   integer(ip) :: ndime
   
   a%ndofbc = 1
   a%ndofbcstart = 0
end subroutine

subroutine plcd_SetExmod(a)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a
   
   a%exmod = 'plcd'
   a%namod = 'PLCD '
   
end subroutine


subroutine plcd_SetNdofn(a)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = ndime
      
end subroutine

subroutine plcd_SetNdofbc(a)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofbc=ndime
      
end subroutine



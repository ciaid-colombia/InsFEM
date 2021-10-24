subroutine nsi_SetExmod(a)
   use typre
   use Mod_NavierStokes
   implicit none
   
   class(NavierStokesProblem) :: a
   
   a%exmod = 'nsi'
   a%namod = 'NSTINC'
   
end subroutine


subroutine nsi_SetNdofn(a)
   use typre
   use Mod_NavierStokes
   implicit none
   
   class(NavierStokesProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = ndime+1
   
   
end subroutine

subroutine nsi_SetNdofbc(a)
   use typre
   use Mod_NavierStokes
   implicit none
   
   class(NavierStokesProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(a%ndofbc)
   a%ndofbcstart=0
   
   
end subroutine



subroutine nsc_SetExmod(a)
   use typre
   use Mod_NSCompressible
   implicit none
   
   class(NSCompressibleProblem) :: a
   
   a%exmod = 'nsc'
   a%namod = 'NSCOMP'
   
end subroutine


subroutine nsc_SetNdofn(a)
   use typre
   use Mod_NSCompressible
   implicit none
   
   class(NSCompressibleProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = ndime+2
   
   
end subroutine

subroutine nsc_SetNdofbc(a)
   use typre
   use Mod_NSCompressible
   implicit none
   
   class(NSCompressibleProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofbc=ndime+2
   
end subroutine



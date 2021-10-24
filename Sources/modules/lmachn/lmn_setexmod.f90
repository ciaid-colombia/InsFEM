subroutine lmn_SetExmod(a)
   use typre
   use Mod_LowMach
   implicit none
   
   class(LowMachProblem) :: a
   
   a%exmod = 'lmn'
   a%namod = 'LMACHN'
   
end subroutine


subroutine lmn_SetNdofn(a)
   use typre
   use Mod_LowMach
   implicit none
   
   class(LowMachProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = ndime+2
   
   
end subroutine

subroutine lmn_SetNdofbc(a)
   use typre
   use Mod_LowMach
   implicit none
   
   class(LowMachProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofbc=ndime+1
   a%ndofbcstart=0
   
   
end subroutine



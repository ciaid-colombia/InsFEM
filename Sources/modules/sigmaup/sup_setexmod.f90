subroutine sup_SetExmod(a)
   use typre
   use Mod_ThreeField
   implicit none
   
   class(ThreeFieldNSProblem) :: a
   
   a%exmod = 'nsi'
   a%namod = 'sigmaup'
   
end subroutine


subroutine sup_SetNdofn(a)
   use typre
   use Mod_ThreeField

   implicit none
   
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   a%ndofn = (ndime-1)*(ndime-1)+2+ndime+1
   
   
end subroutine

subroutine sup_SetNdofbc(a)
   use typre
   use Mod_ThreeField
   implicit none
   integer(ip) :: aux,aux1,aux2
   
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   
!    Viscous case
    a%ndofbc = ndime
    a%ndofbcstart=(ndime-1)*(ndime-1)+2    
   
   
end subroutine



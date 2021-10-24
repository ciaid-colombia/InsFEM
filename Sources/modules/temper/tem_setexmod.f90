subroutine tem_SetExmod(a)
   use typre
   use Mod_Temperature
   implicit none
   
   class(TemperatureProblem) :: a
   
   a%exmod = 'tem'
   a%namod = 'TEMPER'
end subroutine


subroutine tem_SetNdofn(a)
   use typre
   use Mod_Temperature
   implicit none
   
   class(TemperatureProblem) :: a
   
   a%ndofn = 1
end subroutine

subroutine tem_SetNdofbc(a)
   use typre
   use Mod_Temperature
   implicit none
   
   class(TemperatureProblem) :: a
   integer(ip) :: ndime
   
   a%ndofbc = 1
   a%ndofbcstart = 0
end subroutine

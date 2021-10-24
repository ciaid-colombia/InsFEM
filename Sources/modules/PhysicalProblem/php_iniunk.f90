subroutine php_iniunk(a)
!   This routine sets up the initial conditions for the a%velocity.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   real(rp)                :: venew,veold,dummr
   real(rp), external      :: funcre
   integer(ip)             :: ibopo,idime,ipoin,ndime,npoin
   
   !Update boundary conditions
   call a%Updbcs
   
   !Initialize a%kfl_stead
   a%kfl_stead = 0

   !SpecificIniunk
   call a%SpecificIniunk
   
   
end subroutine php_iniunk


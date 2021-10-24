subroutine nsc_ManufacturedBoundaryConditions(a)
   use typre
   use Mod_NSCompressible  
   implicit none
   class(NSCompressibleProblem) :: a   
   integer(ip)  :: ipoin, npoin, ndime, ibopo
   real(rp), pointer  :: exnor(:,:) => NULL()
   
   
   if (a%ManufacturedBoundaryCondition == 1) then

   endif


end subroutine

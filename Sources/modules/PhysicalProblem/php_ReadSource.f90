subroutine php_ReadSource(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: ipoin
   
   !Read point number
   ipoin                   = int(a%Listener%param(1)) 
   
   !To local numbering
   call a%Mesh%Global2Local(ipoin,ipoin)  

   !Specific on nodes
   a%gipoin = ipoin
   call a%SpecificReadSource
   
   
   
end subroutine

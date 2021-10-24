subroutine php_ReadOnNodes(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: ipoin
   
   !Read point number
   ipoin                   = int(a%Listener%param(1)) 
   
   !To local numbering
   call a%Mesh%Global2Local(ipoin,ipoin)  
   
   !Read nodal Fixity conditions
   a%kfl_fixno(1,ipoin)   = int(a%Listener%param(2))
   a%bvess(1:a%ndofbc,ipoin,1)=a%Listener%param(3:2+a%ndofbc)
   call codfix(a%ndofbc,a%kfl_fixno(1,ipoin))
   
   !Non-constant boundary conditions
   if (a%kfl_conbc/=1) then   
      a%kfl_funno(ipoin)  = int(a%Listener%param(3+a%ndofbc))
      !if( (maxval(a%kfl_fixno(:,ipoin)) > 0).and.a%kfl_funno(ipoin)/=0) then
         a%bvess(1:a%ndofbc,ipoin,2)=a%bvess(1:a%ndofbc,ipoin,1)
      !end if
   endif   
   
   !Specific on nodes
   a%gipoin = ipoin
   call a%SpecificReadOnNodes
   
end subroutine

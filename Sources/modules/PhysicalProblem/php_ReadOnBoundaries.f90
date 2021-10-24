subroutine php_InitializeOnBoundariesConditions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a

   a%bvnat_coun = 0                                   ! Dimension of the bvnat array
end subroutine
   
   
subroutine php_ReadOnBoundaries(a,knodb)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: knodb(:)
   
   integer(ip) :: pnodb,aux_fixbo
   integer(ip) :: iboun,idime,inodb,ipsta,ndofi
   
   !Read boundary nodes
   pnodb=int(a%Listener%param(2))
   knodb(1:pnodb)=int(a%Listener%param(3:2+pnodb))
   
   !To local numbering
   call a%Mesh%Global2Local(pnodb,knodb(1:pnodb),knodb(1:pnodb))
   
   !Find which is the boundary
   call a%Mesh%Finbou(pnodb,knodb(1:pnodb),iboun)    
   if(iboun==0) call runend('php_ReadOnBoundaries: Boundary not found')
   
   !Check that the boundary wasn't already assigned a condition
   !if (a%kfl_fixbo(iboun) /= 0) call runend('php_ReadOnBoundaries: tried to assign condition to boundary, but boundary was already assigned a condition')
   
   ipsta=3+pnodb
   aux_fixbo = int(a%Listener%param(ipsta))
   !write(*,*) 'aux_fixbo         =',aux_fixbo
   !write(*,*) 'a%kfl_fixbo(iboun)=',a%kfl_fixbo(iboun)
   if ((a%kfl_fixbo(iboun)>=0).AND.(a%kfl_fixbo(iboun) /= aux_fixbo)) call runend('php_ReadOnBoundaries: tried to re-assign a different condition to boundary')
   a%kfl_fixbo(iboun) = aux_fixbo
   
   ! Dirichlet
   if(a%kfl_fixbo(iboun)==1) then
      if (associated(a%bvnat(iboun)%a)) then
         a%bvnat_coun = a%bvnat_coun - size(a%bvnat(iboun)%a)
         deallocate(a%bvnat(iboun)%a)
      endif
      
      allocate(a%bvnat(iboun)%a(a%ndofbc*pnodb))
      ndofi=0
      do inodb=1,pnodb
         do idime=1,a%ndofbc
            ipsta=ipsta+1
            ndofi=ndofi+1
            a%bvnat(iboun)%a(ndofi)=a%Listener%param(ipsta)
         end do
      end do
      if(a%kfl_conbc==1) then
         
      else
         a%kfl_funbo(iboun) = int(a%Listener%param(ipsta+1))
      end if
      
   !Neumann Boundary condition
   else if(a%kfl_fixbo(iboun)==2) then
      if (associated(a%bvnat(iboun)%a)) then
         a%bvnat_coun = a%bvnat_coun - size(a%bvnat(iboun)%a)
         deallocate(a%bvnat(iboun)%a)
      endif
   
      allocate(a%bvnat(iboun)%a(1))
      a%bvnat(iboun)%a(1)=a%Listener%param(ipsta+1)
      if (a%kfl_conbc /= 1) a%kfl_funbo(iboun) =int(a%Listener%param(ipsta+2))
   end if
   
   !On boundaries specific
   a%giboun = iboun
   a%gipsta = ipsta
   call a%SpecificReadOnBoundaries
   
   if (associated(a%bvnat(iboun)%a)) then
      a%bvnat_coun = a%bvnat_coun + size(a%bvnat(iboun)%a)
   endif
   
end subroutine

subroutine php_FinalizeOnBoundaries(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: iboun,nboun
   
   call a%Memor%allocObj(0,'bvnat%a','php_FinalizeOnBoundaries',a%bvnat_coun*rp)
   
   call a%Mesh%GetNboun(nboun)
   !if boundary was not read, mark it as free
   do iboun=1,nboun
      if (a%kfl_fixbo(iboun)<0) then
         a%kfl_fixbo(iboun) = 0
      end if
   end do
   
end subroutine

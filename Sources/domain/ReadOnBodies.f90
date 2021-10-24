subroutine ReadOnBodies(a,knodb)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: knodb(*)
   
   integer(ip) :: ibody,auxnbound
   
   integer(ip) :: pnodb,nbody,ierr
   integer(ip) :: iboun,idime,inodb,ipsta,ndofi
   
   !testing
   auxnbound=a%nboun
   
   !body number initialization
   ibody=0
   
   !Read boundary nodes
   pnodb=int(a%Listener%param(2))
   knodb(1:pnodb)=int(a%Listener%param(3:2+pnodb))   
   ibody = int(a%Listener%param(2+pnodb+1))
   
   !To local numbering
   call a%Global2Local(pnodb,knodb(1:pnodb),knodb(1:pnodb))  
   
   !Find which is the boundary
   call a%Finbou(pnodb,knodb(1:pnodb),iboun)    
   if (iboun==0) call runend('ReadOnBodies: Boundary not found')   
   
   a%lbody(iboun) = ibody 
   
   !Body numbers
   nbody=a%nbody
   a%nbody=max(ibody,nbody)
   
end subroutine   



subroutine nsi_RVE(a,itask)
   !This routine ends a time step of the incoma%pressible NS equations.
   use MPI
   use typre
   use Mod_Element
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable    :: elvel(:,:),gpvel(:),auxvel(:),auxvelg(:)
   real(rp)    :: fact,dvol,vnor
   integer(ip) :: itask,ierr
   integer(ip) :: ielem,nelem,ndime,nboun,iboun,npoinlocal,idime,ipoin,pnodb
   integer(ip) :: igaus,inodb,nodecount,inode,jnode,jnodeS,jnodeM,ibopo
   real(rp), pointer :: exnor(:,:) => NULL()
   integer(ip), pointer :: lnodb(:) => NULL()
   logical  :: cycleflag

   if (itask == 1) then
      a%RVETensor = reshape([ 0.0_rp, 1.0_rp, 0.0_rp, &
                              1.0_rp, 0.0_rp, 0.0_rp, &
                              0.0_rp, 0.0_rp, 0.0_rp],[3,3])

   !Post-CrankNicolson
   elseif (itask == 2) then
   
      !Element Allocation
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_RVE')
   
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoinLocal(npoinLocal)      
      
      boundaries: do iboun=1,nboun
         !Load Element
         call a%Mesh%BoundaryLoad(iboun,e)
         call e%elmdel

         cycleflag = .false.
         do inodb = 1,e%pnodb
            ipoin = e%lnodb(inodb)
            call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
            if (ibopo == 0) then
               cycleflag =  .true.
            else
               call vecnor(exnor(:,1),e%ndime,vnor,2)
               if (vnor == 0.0_rp) cycleflag =  .true. 
            endif
         enddo

         if (cycleflag) cycle
            
         call e%bounor
         call a%Mesh%GetLnodb(iboun,pnodb,lnodb)
         do jnode = 1,e%pnodb
            jnodeS = lnodb(jnode)
            call a%Mesh%GetMaster(jnodeS,jnodeM)
            if (jnodeS /= jnodeM) a%veloc(1:ndime,jnodeS,3) = a%veloc(1:ndime,jnodeS,3) + matmul(a%RVETensor(1:e%ndime,1:e%ndime),e%baloc(:,ndime))
         enddo

      enddo boundaries
   
      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','nsi_RVE')
      call a%Memor%alloc(ndime,gpvel,'gpvel','nsi_RVE')
      call a%Memor%alloc(ndime,auxvel,'auxvel','nsi_RVE')
      call a%Memor%alloc(ndime,auxvelg,'auxvelg','nsi_RVE')
      
      call a%Mesh%ElementSetPointers(e)
      elements : do ielem = 1,nelem   
         call a%Mesh%ElementLoad(ielem,e)  
      
         nodecount = 0
         do inode = 1,e%pnode
            if (e%lnods(inode) <= npoinLocal) then
               nodecount = nodecount+1
            endif
         enddo
         fact = real(nodecount)/real(e%pnode)
      
         call e%gather(ndime,elvel(:,1:e%pnode),a%veloc(:,:,3))

         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            !Cartesian derivatives and Jacobian at center of gravity
            call e%elmdcg
            dvol = e%weigp(e%igaus)*e%detjm
            call e%interpg(ndime,elvel,gpvel)
            auxvel = auxvel + fact*dvol*gpvel
         enddo gauss_points
      enddo elements
      
      CALL MPI_ALLREDUCE(auxvel,auxvelg,ndime,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
      do idime=1,ndime
         a%veloc(idime,:,3) = a%veloc(idime,:,3) - auxvelg(idime)
      enddo
      
      call a%Memor%dealloc(ndime,auxvelg,'auxvelg','nsi_RVE')
      call a%Memor%dealloc(ndime,auxvel,'auxvel','nsi_RVE')
      call a%Memor%dealloc(ndime,gpvel,'gpvel','nsi_RVE')
      call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','nsi_RVE')
      !Element Deallocation
      call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_RVE')
      
   endif

end subroutine nsi_RVE

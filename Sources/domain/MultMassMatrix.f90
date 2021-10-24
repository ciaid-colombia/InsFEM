subroutine MultMassMatrix(a,array,idofn,ndofn)
   use typre
   use Mod_Memor
   use Mod_ParallelSystemInterface
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(ParallelSystemInterface),pointer :: LinearSystem => NULL()
   class(FemMesh) :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: inode,jnode,igaus,ielem,idofn,ipoin,npoin,nelem,npoinLocal,ndofn
   real(rp)    :: array(ndofn,a%npoin),adiag,dvol,vecout(a%npoin),vecin(a%npoin)
   real(rp), allocatable :: elmat(:,:),elrhs(:)
   character(150) :: auxstring, auxstring2

   call a%ElementAlloc(e,a%Memor,'DefaultRule','L2Projector')
   call a%Memor%Alloc(e%mnode,e%mnode,elmat,'elmat','L2Projector')
   call a%Memor%Alloc(e%mnode,elrhs,'elrhs','L2Projector')
   call a%GetNelem(nelem)
   call a%GetNpoinLocal(npoinLocal)

   auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.dom.sol '
   auxstring2 = 'dom '
  
   call a%ParallelLibrary%CreateSystem(LinearSystem, a%Memor)
   call a%InitializeSystem(1_ip,0_ip,LinearSystem,a%Memor,auxstring,auxstring2,a%lun_solve_dom)

   call LinearSystem%ToZero
   vecin(:) = array(idofn,:)

   !Element loop
   elements : do ielem = 1,nelem
      !Load Element
      call a%ElementLoad(ielem,e)   
      call e%elmdel
      
      elmat = 0.0_rp
     
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         do inode = 1,e%pnode
            do jnode = 1,e%pnode
               elmat(inode,jnode) = elmat(inode,jnode) + dvol*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            enddo
         enddo
       enddo gauss_points

       !Assembly
       call LinearSystem%Assembly(e,elmat,elrhs)
   enddo elements  

   call LinearSystem%Matvec(vecin,vecout)
   array(idofn,:) = vecout 

   !Communicate Ghosts
   call a%ArrayCommunicator%GhostCommunicate(1,array(idofn,:)) 

   call a%Memor%Dealloc(e%mnode,e%mnode,elmat,'elmat','L2Projector')
   call a%Memor%Dealloc(e%mnode,elrhs,'elrhs','L2Projector')
   call a%ElementDealloc(e,a%Memor,'DefaultRule','Projector')
   call LinearSystem%Deallocate
   call a%ParallelLibrary%DeallocateSystem(LinearSystem,a%Memor)
end subroutine 

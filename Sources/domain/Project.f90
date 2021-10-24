subroutine Project(a,ndofn,array)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   real(rp)    :: array(ndofn,a%npoin)
   integer(ip) :: kfl_perio,idofn,ipoin

   interface
      subroutine AfterProjectOperations(a,ndofn,array)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
         integer(ip) :: ndofn
         real(rp)    :: array(:,:)
      end subroutine
   end interface  
 
   !If the mass matrix has not been computed yet, compute it
   if (a%kfl_ProjectorReady == 0) then
      a%kfl_ProjectorReady = 1
      
      call a%ComputeMassMatrix
   endif
   
   !Just one degree of freedom
   if (ndofn == 1 .or. a%kfl_ProjectionFull == 1) then
      call a%L2ProjectorSystem%CopyToRHS(array)
      call a%L2ProjectorSystem%Solve(array)
      
   !More than one degree of freedom   
   else
      !copy each dof to auxunkno, solve, bring it back
      do idofn = 1,ndofn
         !Solve
         call a%L2ProjectorSystem%CopyToRHS(array(idofn,:))
         call a%L2ProjectorSystem%Solve(array(idofn,:))
      enddo   
      
   endif
   
   call AfterProjectOperations(a,ndofn,array)
   
end subroutine

subroutine AfterProjectOperations(a,ndofn,array)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   real(rp)    :: array(:,:)
   
   integer(ip) :: kfl_perio,idofn,ipoin
   
   !Hanging nodes
   if (a%kfl_HangingNodes .eqv. .true.) call a%InterpolateHangingValues(ndofn,array)
   
   !Communicate Ghosts
   call a%ArrayCommunicator%GhostCommunicate(ndofn,array)
   
   !Periodic boundary conditions
   call a%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%MasterToSlave(ndofn,array)
   
end subroutine   

subroutine ComputeMassMatrix(a)
   use typre
   use def_parame
   use Mod_int2str
   use Mod_Iofile
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   character(150) :: auxstring, auxstring2
   character(150):: fil_solve
   
   real(rp), allocatable :: elmat(:,:,:,:), elrhs(:,:)
   class(FiniteElement), pointer :: e => NULL()
   
   real(rp) :: dvol
   integer(ip) :: inode,jnode,ielem,igaus
   
   !Open file for postprocessing the results of the system solve
   if (a%MPIrank == a%MPIroot) then
      fil_solve = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.dom.sol'
      call iofile(zero,a%lun_solve_dom,fil_solve,'MESH SOLVE')
   endif
   
   !Create the LinearSystem
   auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.dom.sol '
   auxstring2 = 'dom '
   call a%ParallelLibrary%CreateSystem(a%L2ProjectorSystem, a%Memor) 
   call a%L2ProjectorSystem%SetFlush(a%kfl_flush)
   call a%L2ProjectorSystem%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call a%InitializeSystem(1_ip,0_ip,a%L2ProjectorSystem,a%Memor,auxstring,auxstring2,a%lun_solve_dom)
   
   
   !Now we compute the mass matrix
   call a%L2ProjectorSystem%ToZero
   
   !Allocations
   call a%ElementAlloc(e,a%Memor,'DefaultRule','L2Projector')
   
   call a%Memor%Alloc(1,e%mnode,1,e%mnode,elmat,'elmat','L2Projector')
   call a%Memor%Alloc(1,e%mnode,elrhs,'elrhs','L2Projector')
   
   !Element loop
   elements : do ielem = 1,a%nelem
      !Load Element
      call a%ElementLoad(ielem,e)   
      call e%elmdel
      
      elmat = 0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         do inode = 1,e%pnode
            do jnode = 1,e%pnode
               elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) +dvol*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            enddo
         enddo
      enddo gauss_points
      
      !Assembly
      call a%L2ProjectorSystem%Assembly(e,elmat,elrhs)
      
   enddo elements  
   
   if (a%kfl_HangingNodes .eqv. .true.) call a%AssemblyHangingNodesDiag(1_ip,a%L2ProjectorSystem,a%Memor)
   if (a%kfl_perio == 1) call a%AssemblyPeriodicBC(1_ip,a%L2ProjectorSystem,a%Memor)
   
   !Deallocations
   call a%Memor%Dealloc(1,e%mnode,1,e%mnode,elmat,'elmat','L2Projector')
   call a%Memor%Dealloc(1,e%mnode,elrhs,'elrhs','L2Projector')
   
   call a%ElementDealloc(e,a%Memor,'DefaultRule','Projector')
  
end subroutine

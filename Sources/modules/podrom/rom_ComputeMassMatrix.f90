subroutine rom_ComputeMassMatrix(a)
   use typre
   use def_parame
   use Mod_int2str
   use Mod_Iofile
   use Mod_Element
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable :: elmat(:,:,:,:), elrhs(:,:)
   character(150) :: auxstring, auxstring2
   character(150) :: fil_solve
   character(3)   :: exmod
   real(rp)       :: dvol
   integer(ip)    :: inode,jnode,ielem,igaus,npoin,nelem
   integer(ip)    :: auxkfl_alemov = 0
   integer(ip), save :: ipass=0
  
   exmod = a%exmod
   a%Mesh%kfl_ProjectorReady = 1
   a%Mesh%kfl_ProjectionFull = 1
   !Open file for postprocessing the results of the system solve
   if (a%Mesh%MPIrank == a%Mesh%MPIroot .and. ipass==0) then
      ipass=1
      fil_solve = trim(a%Mesh%OutputFolder)//'/'//adjustl(trim(a%Mesh%namda))//adjustl(trim(int2str(a%Mesh%MPIrank)))//'.'//adjustl(trim(exmod))//'.dom.sol'
      call iofile(zero,a%Mesh%lun_solve_dom,fil_solve,'MESH SOLVE')
   endif
   
   !Create the LinearSystem
   auxstring = trim(a%Mesh%InputFolder)//'/'//adjustl(trim(a%Mesh%namda))//'.dom.sol '
   auxstring2 = 'dom '
   call a%EigenLibrary%CreateSystemProj(a%Mesh%L2ProjectorSystem,a%Mesh%Memor) 
   call a%Mesh%L2ProjectorSystem%SetFlush(a%Mesh%kfl_flush)
   call a%Mesh%L2ProjectorSystem%SetMPI(a%Mesh%MPIcomm,a%Mesh%MPIsize,a%Mesh%MPIroot,a%Mesh%MPIrank)
   call a%Mesh%InitializeSystem(a%ndofn,0_ip,a%Mesh%L2ProjectorSystem,a%Mesh%Memor,auxstring,auxstring2,a%Mesh%lun_solve_dom)
   
   
   !Now we compute the mass matrix
   call a%Mesh%L2ProjectorSystem%ToZero
   
   !Allocations
   call a%Mesh%ElementAlloc(e,a%Mesh%Memor,'DefaultRule','L2Projector')
   
   call a%Mesh%Memor%Alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','L2Projector')
   call a%Mesh%Memor%Alloc(a%ndofn,e%mnode,elrhs,'elrhs','L2Projector')
   
   !Element loop
   elements : do ielem = 1,a%Mesh%nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)   
      call e%elmdel
      
      elmat = 0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         do inode = 1,e%pnode
            do jnode = 1,e%pnode
               elmat(1:a%ndofn,inode,1:a%ndofn,jnode) = elmat(1:a%ndofn,inode,1:a%ndofn,jnode) +dvol*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            enddo
         enddo
      enddo gauss_points
      
      !Assembly
      call a%Mesh%L2ProjectorSystem%Assembly(e,elmat,elrhs)
      
   enddo elements  
   
   if (a%Mesh%kfl_HangingNodes .eqv. .true.) call a%Mesh%AssemblyHangingNodesDiagToZero(a%ndofn,a%Mesh%L2ProjectorSystem,a%Mesh%Memor)
   if (a%Mesh%kfl_perio == 1) call a%Mesh%AssemblyPeriodicBCToZero(a%ndofn,a%Mesh%L2ProjectorSystem,a%Mesh%Memor)
   
   !Deallocations
   call a%Mesh%Memor%Dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','L2Projector')
   call a%Mesh%Memor%Dealloc(a%ndofn,e%mnode,elrhs,'elrhs','L2Projector')
   
   call a%Mesh%ElementDealloc(e,a%Mesh%Memor,'DefaultRule','Projector')
 
   !Set extra arrays 
   call a%EigenSystem%SetBasisToLS(a%Mesh%L2ProjectorSystem,a%Basis)
   call a%Mesh%L2ProjectorSystem%Init2(a%kfl_massMatrix,a%kfl_Precondition)

end subroutine

subroutine rom_ComputeLumpedMassMatrix(a)
   use typre
   use def_parame
   use Mod_int2str
   use Mod_Iofile
   use Mod_Element
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable :: elmat(:,:,:,:), elrhs(:,:),vmass(:),auxshape(:)
   character(150) :: auxstring, auxstring2
   character(150) :: fil_solve
   character(3)   :: exmod
   real(rp)       :: dvol
   integer(ip)    :: inode,jnode,ielem,igaus,npoin,nelem
   integer(ip)    :: auxkfl_alemov = 0
   integer(ip), save :: ipass=0

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   !For the first pass, a%Mesh%displ has not been set yet
   if ((a%Mesh%kfl_alemov == 1)  ) then
      if (.not. (associated(a%Mesh%displ)) ) then
         auxkfl_alemov = 1
         a%Mesh%kfl_alemov = 0
      elseif (size(a%Mesh%displ,2) /= npoin) then
         !If adaptive this might not be the same size
         auxkfl_alemov = 1
         a%Mesh%kfl_alemov = 0
      endif
   endif
   
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
   
   call a%Mesh%ElementAlloc(e,a%Memor,'ForceClosedRule','vmass')
   !Compute the lumped mass vector
   call a%Memor%alloc(npoin,vmass,'VMASS','exnor')
   call a%Memor%alloc(e%mnode,auxshape,'auxshape','rom_lumped')
   !We force a closed integration rule for computing vmass
   !Loop over elements and compute vmass
   
   do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      call e%elmdel
      
      !Gauss Point Loop
      do igaus=1,e%pgaus
        e%igaus = igaus
        call e%elmder
        dvol = e%weigp(e%igaus)*e%detjm
        auxshape(1:e%pnode) = dvol*e%shape(1:e%pnode,e%igaus)
        call a%Mesh%AssemblyToArray(e,1_ip,auxshape,vmass)
      enddo
   enddo 
   !call a%Mesh%EndAssemblyToArray(1_ip,vmass)
   
   call a%Mesh%ElementDealloc(e,a%Memor,'ForceClosedRule','vmass')

   
   !Communicate between subdomains
   call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,vmass)
   
   !Periodic Boundary conditions
   if (a%Mesh%kfl_perio == 1) call a%Mesh%MasterToSlave(1_ip,vmass)
   
   !Assign mass vector to mass matrix
   call a%Mesh%L2ProjectorSystem%ToZero
   
   !Allocations
   call a%Mesh%ElementAlloc(e,a%Mesh%Memor,'DefaultRule','L2Projector')
   
   call a%Mesh%Memor%Alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','L2Projector')
   call a%Mesh%Memor%Alloc(a%ndofn,e%mnode,elrhs,'elrhs','L2Projector')
   
   !Element loop
   elements : do ielem = 1,nelem
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
            elmat(1:a%ndofn,inode,1:a%ndofn,inode) = vmass(inode)
         enddo
      enddo gauss_points
      
      !Assembly
      call a%Mesh%L2ProjectorSystem%Assembly(e,elmat,elrhs)
      
   enddo elements  
   
   if (a%Mesh%kfl_HangingNodes .eqv. .true.) call a%Mesh%AssemblyHangingNodesDiagToZero(a%ndofn,a%Mesh%L2ProjectorSystem,a%Mesh%Memor)
   if (a%Mesh%kfl_perio == 1) call a%Mesh%AssemblyPeriodicBCToZero(a%ndofn,a%Mesh%L2ProjectorSystem,a%Mesh%Memor)
   
   !Deallocations
   call a%Memor%dealloc(npoin,vmass,'VMASS','exnor')
   call a%Memor%dealloc(e%mnode,auxshape,'auxshape','rom_lumped')
   call a%Mesh%Memor%Dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','L2Projector')
   call a%Mesh%Memor%Dealloc(a%ndofn,e%mnode,elrhs,'elrhs','L2Projector')
   
   call a%Mesh%ElementDealloc(e,a%Mesh%Memor,'DefaultRule','Projector')
 
   !Set extra arrays 
   call a%EigenSystem%SetBasisToLS(a%Mesh%L2ProjectorSystem,a%Basis)
   call a%Mesh%L2ProjectorSystem%Init2(a%kfl_massMatrix,a%kfl_Precondition)

end subroutine

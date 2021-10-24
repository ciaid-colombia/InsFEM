subroutine MeshInitializeAdaptiveRefinement(a,Refiner)
   use typre
   use Mod_Mesh
   use Mod_AdaptiveInterface
   implicit none
   class(FemMesh), target :: a
   class(AdaptiveRefinerInterface) :: Refiner

   call a%SetAdaptiveRefiner(Refiner)
   
  
end subroutine

subroutine MeshMatchOldBoundariesToRefiner(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh), target :: a
   integer(ip), pointer :: pboel(:) => NULL()
   integer(ip) :: iboun 
   
   call a%Timer%Refine%Tic
   
   !For boundary info trespassing
   if (allocated(a%auxBoundaryMatch1)) then
      call a%Memor%realloc(a%nboun,a%auxBoundaryMatch1,'auxBoundaryMatch1','MeshInitializeAdaptiveRefinement')
   else
      call a%Memor%alloc(a%nboun,a%auxBoundaryMatch1,'auxBoundaryMatch1','MeshInitializeAdaptiveRefinement')
   endif
   
   !aux pboel if required
   if (a%nelty == 1) then
      call a%Memor%palloc(a%nboun+1,pboel,'pboel','MeshInitializeAdaptiveRefinement')
      pboel(1) = 1
      do iboun = 1,a%nboun
         pboel(iboun+1) = pboel(iboun) + a%pboel(1) + 1
      enddo
   else
      pboel => a%pboel
   endif
   call a%Refiner%MatchOldLboelToOldExternalBoundaries(a%nboun,pboel,a%lboel,a%auxBoundaryMatch1,a%auxOldNboun)
   
   !deallocate auxpboel if required
   if (a%nelty == 1) then
      call a%Memor%pdealloc(a%nboun+1,pboel,'pboel','MeshInitializeAdaptiveRefinement')
   endif
   
   call a%Timer%Refine%Toc
end subroutine
   

subroutine MeshRefine(a, itask)
   use MPI
   use typre
   use def_parame
   use Mod_int2str
   use Mod_iofile
   use Mod_Mesh
   use Mod_Partitionate
   use Mod_HangingNodes
   implicit none
   class(FemMesh), target :: a
   character(6) :: itask
   
   integer(ip) :: newnelem,newnpoin,newnpoinLocal,newgnpoin,newlnodsSize,hanginglistsize,newnboun
   integer(ip), allocatable :: pnods(:),lnods(:)
   integer(ip), allocatable :: LocalToGlobal(:),ProcessorList(:)
   real(rp), pointer :: coord(:,:) => NULL()
   real(rp), allocatable :: newcoord(:,:)
   real(rp), allocatable :: newblock(:)
   integer(ip), allocatable :: pHangingList(:),lHangingList(:)
   real(rp), allocatable :: rHangingList(:)
   
   integer(ip), allocatable :: pelpo(:),lelpo(:),ia(:),ja(:)
   
   logical, allocatable :: isExnor(:)
   type(r1p), allocatable :: ExternalNormal(:)
   
   integer(ip) :: icount,ipoin
   real(rp), allocatable :: oldExwor(:,:), newExwor(:,:)
   integer(ip) :: ndime,oldnboun,iboun,iboundarymatchcoun,ielem
   
   integer(ip), allocatable :: lbody(:),auxlbody(:)
   integer(ip), pointer :: pboel(:) => NULL()
   
   integer(ip) :: iskfl_alemov, ierr, maxsizePHangingList
   character(150) :: fil_solve

   integer(ip), pointer :: iauxBoundaryMatch1(:)
   type(i1p), pointer   :: iauxBoundaryMatch2(:)

   interface
      subroutine ComputeRanges(a)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine PeriodicBCModify(a)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
   end interface
   
   call a%Timer%Refine%Tic
   
    !Deallocate projector if necessary
    if (a%kfl_ProjectorReady == 1) then
      call a%L2ProjectorSystem%Deallocate
      call a%ParallelLibrary%DeallocateSystem(a%L2ProjectorSystem, a%Memor) 
      if (a%MPIrank == a%MPIroot) then
         fil_solve = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.dom.sol'
         call iofile(two,a%lun_solve_dom,fil_solve,'MESH SOLVE')
      endif   
      a%kfl_ProjectorReady = 0
   endif
   
   !Cases for which mesh refinement is not ready
   if (count(a%isExnor) /= 0) then
      call runend('Mesh refinement not ready for isExnor')
   endif
   oldnboun = a%nboun
   
   !New Dimensions
   call a%Refiner%GetPointDimensions(newnpoinLocal,newnpoin,newgnpoin)
   call a%Refiner%GetElementDimensions(newnelem,newLnodsSize)
   
   !New Pnods and Lnods
   call a%Memor%alloc(newnelem+1,pnods,'pnods','MeshRefine')
   call a%Memor%alloc(newLnodsSize,lnods,'lnods','MeshRefine')
   call a%Refiner%GetLnods(pnods,lnods)
   
   !New Local Ordering and ProcessorList
   call a%Memor%alloc(newnpoin,LocalToGlobal,'LocalToGlobal','MeshRefine')
   call a%Memor%alloc(newnpoin,ProcessorList,'ProcessorList','MeshRefine')
   call a%Refiner%GetLocalOrdering(LocalToGlobal,ProcessorList)

   !New coord array
   call a%Memor%alloc(a%ndime,newnpoin,newcoord,'coord','MeshRefine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndime,a%coord,newcoord)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndime,a%coord,newcoord)
   endif

   !New block array
   call a%Memor%alloc(newnpoin,newblock,'geoblock','MeshRefine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%geoblock,newblock)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%geoblock,newblock)
   endif
   
   !New HangingList info
   call a%Refiner%GetHangingListDimensions(HangingListSize)
   call a%Memor%alloc(newnpoin+1,pHangingList,'pHangingList','MeshRefine')
   call a%Memor%alloc(HangingListSize,lHangingList,'lHangingList','MeshRefine')
   call a%Memor%alloc(HangingListSize,rHangingList,'rHangingList','MeshRefine')
   call a%Refiner%GetHangingList(pHangingList,lHangingList,rHangingList)

   !Externally fixed normal
   call a%Memor%alloc(newnpoin,isExnor,'isExnor','MeshRefine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%isExnor,isExnor,'minim')
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%isExnor,isExnor)
   endif
   !Communicate the external normal
   call a%Memor%alloc(a%ndime,a%npoin,oldExwor,'oldExwor','MeshRefine')
   call a%Memor%alloc(a%ndime,newnpoin,newExwor,'newExwor','MeshRefine')
   do ipoin = 1,a%npoin
      if (a%isExnor(ipoin) .eqv. .true.) then
         oldExwor(:,ipoin) = a%ExternalNormal(ipoin)%a
      endif
   enddo
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndime,oldExwor,NewExwor)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndime,oldExwor,NewExwor)
   endif
   call a%Memor%alloc(newnpoin,ExternalNormal,'ExternalNormal','MeshRefine')
   icount = 0
   do ipoin = 1,newnpoin
      if (isExnor(ipoin) .eqv. .true.) then
         allocate(ExternalNormal(ipoin)%a(a%ndime))
         icount = icount + a%ndime
         ExternalNormal(ipoin)%a = NewExwor(:,ipoin)
      endif
   enddo
   call a%Memor%AllocObj(0_ip,'ExternalNormal','MeshRefine',icount*rp)
   call a%Memor%dealloc(a%ndime,a%npoin,oldExwor,'oldExwor','MeshRefine')
   call a%Memor%dealloc(a%ndime,newnpoin,newExwor,'newExwor','MeshRefine')
   !------------------------------------------------------------------------------------

   !Copy lbody for later use
   call a%Memor%alloc(oldnboun,lbody,'lbody','MeshRefine')
   lbody = a%lbody
   
!--------------------------------------------------------------------------------------------------------
!If Alemov, non of these operations can use displacements or velocities, since they are in the modules
!and still in the previous mesh
!so we temporarily deactivate kfl_alemov
iskfl_alemov = 0
if (a%kfl_alemov == 1) then
   iskfl_alemov = 1
   a%kfl_alemov = 0
endif
!-------------------------------------------------------------------------------------------------------   
   
   
   
   !We need to reallocate everything
   call a%DeallocateLocalArrays
   call a%BuildLnodsFromArray(newnelem,pnods,lnods)
   call StatisticsAndElementPostprocessing(a)
   call a%SetNpoinAndBuildOrderingFromArray(newgnpoin,newnpoin,newnpoinLocal,LocalToGlobal,ProcessorList,1)
   call a%BuildCoordinatesFromArray(newcoord)
   call a%BuildBlocksFromArray(newblock)
   call a%BuildExnorFromArray(isExnor,ExternalNormal)
   
   
   call MPI_REDUCE( size(pHangingList,1), maxsizePHangingList, 1, MPI_INTEGER4, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
   call MPI_BCAST(maxsizePHangingList, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   
   if (maxsizePHangingList > 1) then
      call a%BuildHangingNodesFromArray(pHangingList,lHangingList,rHangingList)
      !Compute Mnode for meshes with hanging nodes
      call ComputeMnodeForHangingNodesMeshes(a)
   else
      call a%SetHanging(0)
   endif
   
   !Auxiliary arrays Deallocations
   call a%Memor%dealloc(size(pnods),pnods,'pnods','MeshRefine')
   call a%Memor%dealloc(size(lnods),lnods,'lnods','MeshRefine')
   call a%Memor%dealloc(newnpoin,LocalToGlobal,'LocalToGlobal','MeshRefine')
   call a%Memor%dealloc(newnpoin,ProcessorList,'ProcessorList','MeshRefine')
   call a%Memor%dealloc(a%ndime,newnpoin,newcoord,'coord','MeshRefine')
   call a%Memor%dealloc(newnpoin,newblock,'geoblock','MeshRefine')
   call a%Memor%dealloc(newnpoin+1,pHangingList,'pHangingList','MeshRefine')
   call a%Memor%dealloc(HangingListSize,lHangingList,'lHangingList','MeshRefine')
   call a%Memor%dealloc(HangingListSize,rHangingList,'rHangingList','MeshRefine')
   icount = 0
   do ipoin = 1,newnpoin
      if (isExnor(ipoin) .eqv. .true.) then
         deallocate(ExternalNormal(ipoin)%a)
         icount = icount + a%ndime
      endif
   enddo
   call a%Memor%DeAllocObj(0_ip,'ExternalNormal','MeshRefine',icount*rp)
   call a%Memor%dealloc(newnpoin,ExternalNormal,'ExternalNormal','MeshRefine')
   call a%Memor%dealloc(newnpoin,isExnor,'isExnor','MeshRefine')

   
   !Mesh Operations
   call a%LocalGraph
   call a%BuildArrayCommunicator
   !Compute the coordinates Ranges
   if (itask == 'Rebala') call ComputeRanges(a)
   if (a%kfl_perio == 1) then
      call PeriodicBCModify(a)
         
      !Build local arrays from Scattered Data
      call a%GhostCoord
   endif
   call a%Dombou     
   call a%ComputeVmass
   call a%ExtnorLpoty
   call a%Blbopo
   
   !aux pboel if required
   if (a%nelty == 1) then
      call a%Memor%palloc(a%nboun+1,pboel,'pboel','MeshInitializeAdaptiveRefinement')
      pboel(1) = 1
      do iboun = 1,a%nboun
         pboel(iboun+1) = pboel(iboun) + a%pboel(1) + 1
      enddo
   else
      pboel => a%pboel
   endif
   
   !For communicating boundary data
   call a%GetNboun(newnboun)
   
   !New BoundaryMatch
   if (allocated(a%auxBoundaryMatch2)) then
      call a%Memor%dealloc(size(a%auxBoundaryMatch2,1),size(a%auxBoundaryMatch2,2),a%auxBoundaryMatch2,'auxBoundaryMatch2','MeshRefine')
   endif   
   call a%Memor%alloc(2,newnboun,a%auxBoundaryMatch2,'auxBoundaryMatch2','MeshRefine')
   call a%Refiner%MatchNewLboelToNewElementsAndFaces(newnboun,pboel,a%lboel,a%auxBoundaryMatch2)
   
   if (itask == 'Refine') then
      if (allocated(a%BoundaryNewToOld)) then
         call a%Memor%dealloc(size(a%BoundaryNewToOld),a%BoundaryNewToOld,'BoundaryNewToOld','MeshRefine')
      endif
      call a%Memor%alloc(newnboun,a%BoundaryNewToOld,'BoundaryNewToOld','MeshRefine')
      call a%Refiner%MatchOldBoundariesToNewBoundaries(oldnboun,a%auxBoundaryMatch1,newnboun,a%auxBoundaryMatch2,a%BoundaryNewToOld)
   elseif (itask == 'Rebala') then
      if (allocated(a%iauxBoundaryMatch1)) then
         call a%Memor%dealloc(size(a%iauxBoundaryMatch1),a%iauxBoundaryMatch1,'iauxBoundaryMatch1','MeshRefine')
         iboundarymatchcoun = 0
         do ielem = 1,size(a%iauxBoundaryMatch2)
            iboundarymatchcoun = iboundarymatchcoun + size(a%iauxBoundaryMatch2(ielem)%l)
            deallocate(a%iauxBoundaryMatch2(ielem)%l)
         enddo
         call a%Memor%deallocObj(0,'iBoundaryMatch','MeshTurnof',iboundarymatchcoun*ip)
         call a%Memor%dealloc(size(a%iauxBoundaryMatch2),a%iauxBoundaryMatch2,'iBoundaryMatch','MeshTurnof')
      endif
      
      call a%Memor%alloc(a%auxOldNboun,a%iauxBoundaryMatch1,'iauxBoundaryMatch1','MeshRefine')
      !iauxBoundaryMatch2 is allocated inside the routine
      !call a%Memor%alloc(newnboun,a%iauxBoundaryMatch2,'iauxBoundaryMatch2','MeshRefine') 
      
      !Inverse matching
      call a%Refiner%GetiBoundaryMatch(oldnboun,a%auxBoundaryMatch1,a%iauxBoundaryMatch1)
      call a%Refiner%GetiBoundaryToElementsAndFacesMatch(newnboun,a%auxBoundaryMatch2,a%iauxBoundaryMatch2,a%Memor)
   endif
   !deallocate auxpboel if required
   if (a%nelty == 1) then
      call a%Memor%pdealloc(a%nboun+1,pboel,'pboel','MeshInitializeAdaptiveRefinement')
   endif
   
   !Update Lbody
   call a%Memor%alloc(newnboun,auxlbody,'auxlbody', 'MeshRefine')
   if (itask == 'Refine') then
      do iboun = 1,a%nboun
         if (a%BoundaryNewToOld(iboun) > 0) then
            auxlbody(iboun) = lbody(a%BoundaryNewToOld(iboun))
         else
            auxlbody(iboun) = 0
         endif
      enddo
      call move_alloc(auxlbody,lbody)
      call a%Memor%deallocObj(0,'auxlbody','MeshRefine',ip*oldnboun)
   elseif (itask == 'Rebala') then
      call a%GetRebalanceInverseBoundaryMatches(iauxBoundaryMatch1,iauxBoundaryMatch2)
      call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,1,lbody,auxlbody)
      call move_alloc(auxlbody,lbody)
      call a%Memor%deallocObj(0,'auxlbody','MeshRefine',ip*oldnboun)
   endif
   call a%BuildListOfBodiesFromArray(lbody)
   
   !Auxiliary array deallocations
   call a%Memor%dealloc(newnboun, lbody,'lbody', 'MeshRefine')
 
  
 
!--------------------------------------------------------------------
!For ALE, end of temporary deactivation of alemov  
if (iskfl_alemov == 1) then
   a%kfl_alemov = 1
endif   
!---------------------------------------------------------------------
   
   !Postprocess mesh info
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outpu_dom,*) 'New global number of points: ', newgnpoin
      call flush(a%lun_outpu_dom)
   endif
   
   call a%Timer%Refine%Toc
   
end subroutine


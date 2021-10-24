subroutine MeshTurnof(a)
   use typre
   use Mod_iofile
   use def_parame
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip) :: ielem,bvnatcoun
   character(150) :: fil_outpu
   
   call a%DeallocateInitial
   call a%DeallocateMeshInfo
   call a%DeallocateLocalArrays
   
   !This is to be deallocated only at the very end (cannot be done in local arrays
   if (allocated(a%auxBoundaryMatch1)) then
      call a%Memor%dealloc(size(a%auxBoundaryMatch1),a%auxBoundaryMatch1,'auxBoundaryMatch1','MeshTurnof')
   endif
   if (allocated(a%iauxBoundaryMatch1)) then
      call a%Memor%dealloc(size(a%iauxBoundaryMatch1),a%iauxBoundaryMatch1,'iauxBoundaryMatch1','MeshTurnof')
   endif
   if (allocated(a%iauxBoundaryMatch2)) then
      bvnatcoun = 0
      do ielem = 1,size(a%iauxBoundaryMatch2)
         bvnatcoun = bvnatcoun + size(a%iauxBoundaryMatch2(ielem)%l)
         deallocate(a%iauxBoundaryMatch2(ielem)%l)
      enddo
      call a%Memor%deallocObj(0,'iBoundaryMatch','MeshTurnof',bvnatcoun*ip)
      call a%Memor%dealloc(size(a%iauxBoundaryMatch2),a%iauxBoundaryMatch2,'iBoundaryMatch','MeshTurnof')
   endif
   if (allocated(a%auxBoundaryMatch2)) then
      call a%Memor%dealloc(size(a%auxBoundaryMatch2,1),size(a%auxBoundaryMatch2,2),a%auxBoundaryMatch2,'auxBoundaryMatch2','MeshRefine')
   endif   
   
   if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_outpu_dom,fil_outpu,'MESH OUTPUT')
   
end subroutine

subroutine MeshDeallocateMeshInfo(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip) :: truemnode
         
   !DeAllocate arrays that store mesh info
   
   truemnode = size(a%ltype)
   
   !Element data info      
   call a%Memor%dealloc(a%nelty,a%nnode,'NNODE','ReadatMPI')
   call a%Memor%dealloc(a%nelty,a%lquad,'LQUAD','ReadatMPI')
   call a%Memor%dealloc(a%nelty,a%ngaus,'NGAUS','ReadatMPI')
   call a%Memor%dealloc(a%nelty,a%nface,'NFACE','ReadatMPI')
   call a%Memor%dealloc(truemnode,a%ltype,'LTYPE','ReadatMPI')   

   call a%Memor%dealloc(a%nelty,a%llapl,'LLAPL','cderda')
   call a%Memor%dealloc(a%nelty,a%ltopo,'LTOPO','cderda')
   call a%Memor%dealloc(a%nelty,a%lrule,'LRULE','cderda')
   call a%Memor%dealloc(a%nelty,a%nnodb,'NNODB','cderda')
   call a%Memor%dealloc(a%nelty,a%ngaub,'NGAUB','cderda')
   call a%Memor%dealloc(a%nelty,a%hnatu,'HNATU','cderda')
   call a%Memor%dealloc(a%ndime,a%mgaus,a%nelty,a%posgp,'POSGP','cderda')   
   
   !Element a%shape function and a%derivatives (a%shape,a%deriv,a%heslo,a%weigp) 
   call a%Memor%dealloc(truemnode,a%mgaus,a%nelty,a%shape,'shape','cshder')
   call a%Memor%dealloc(a%ndime,truemnode,a%mgaus,a%nelty,a%deriv,'deriv','cshder')
   call a%Memor%dealloc(a%ntens,truemnode,a%mgaus,a%nelty,a%heslo,'heslo','cshder')
   call a%Memor%dealloc(a%mgaus,a%nelty,a%weigp,'weigp','cshder') 
   
   !Element Center of gravity a%shape function and a%derivatives (a%shacg,a%dercg,a%weicg) 
   call a%Memor%dealloc(truemnode,a%nelty,a%shacg,'shacg','cshder')
   call a%Memor%dealloc(a%ndime,truemnode,a%nelty,a%dercg,'dercg','cshder')
   call a%Memor%dealloc(a%nelty,a%weicg,'weicg','cshder')  
   
   !Closed Rule integration data
   call a%Memor%dealloc(truemnode,truemnode,a%nelty,a%shape_clo,'shape_clo','cshder')
   call a%Memor%dealloc(a%ndime,truemnode,truemnode,a%nelty,a%deriv_clo,'deriv_clo','cshder')
   call a%Memor%dealloc(a%ntens,truemnode,truemnode,a%nelty,a%heslo_clo,'heslo_clo','cshder')
   call a%Memor%dealloc(truemnode,a%nelty,a%weigp_clo,'weigp_clo','cshder') 
   
   call a%Memor%dealloc(a%mnodb,a%mgaub,a%nelty,a%shapb,'shapb','cshder')
   call a%Memor%dealloc(a%ndimb,a%mnodb,a%mgaub,a%nelty,a%derib,'derib','cshder')        
   call a%Memor%dealloc(a%mgaub,a%nelty,a%weigb,'weigb','cshder')
   
   call a%Memor%dealloc(a%mgaus,truemnode,a%nelty,a%shaga,'shaga','cshder')
   call a%Memor%dealloc(a%mgaub,a%mnodb,a%nelty,a%shagb,'shagb','cshder')
   
   call a%Memor%dealloc(a%mnodb,a%ltypb,'ltypb','dombou')
   
end subroutine
   
subroutine MeshDeallocateLocalArrays(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   interface
      subroutine msh_FinalizeElementDataStructures(a)
         import
         implicit none
         class(FemMesh) :: a
      end subroutine
   end interface
   
   !blbopo
   call a%Memor%dealloc(a%pbopo(a%npoin+1),a%lbopo,'lbopo','blbopo')
   call a%Memor%dealloc(a%npoin+1,a%pbopo,'pbopo','blbopo')
   
   !List of body deallocation
   call a%Memor%dealloc(a%nboun,a%lbody,'lbody','lbody')       
   
      if (a%nelty > 1) then
      call a%Memor%dealloc(a%pboel(a%nboun+1)-1,a%lboel,'lboel','dombou')
      call a%Memor%dealloc(a%nboun+1,a%pboel,'pboel','dombou')
   else
      call a%Memor%dealloc(a%nboun*(a%mnodb+1),a%lboel,'lboel','dombou')
      call a%Memor%dealloc(1,a%pboel,'pboel','dombou')
   endif
   
   !FOR EXNOR and LPOTY 
   call a%Memor%dealloc(a%npoin,a%lpoty,'lpoty','exnor')
   call a%Memor%dealloc(a%ndime,a%ndime,a%nbopo,a%exnor,'exnor','exnor')
   
   !Lelpo and Pelpo
   call a%Memor%dealloc(size(a%lelpo),a%lelpo,'lelpo','BuildElpo')
   call a%Memor%dealloc(size(a%pelpo),a%pelpo,'pelpo','BuildElpo')
   
   !Graph
   call a%Memor%dealloc(a%ia(a%npoin+1)-1,a%ja,'ja','BuildGraph')
   call a%Memor%dealloc(a%npoin+1,a%ia,'ia','BuildGraph')
   
   !Lnods and Coordinates
   if (a%nelty > 1) then 
      call a%Memor%dealloc(a%pnods(a%nelem+1)-1,a%lnods,'lnods','DeallocateLocals')
      call a%Memor%dealloc(a%nelem+1,a%pnods,'pnods','DeallocateLocals')
   else
      call a%Memor%dealloc(size(a%lnods),a%lnods,'lnods','DeallocateLocals')
      call a%Memor%dealloc(1,a%pnods,'pnods','DeallocateLocals')
   endif
   call a%Memor%dealloc(a%ndime,a%npoin,a%coord,'coord','DeallocateLocals')
   call a%Memor%dealloc(a%mnodb,a%mface,a%nelty,a%cfael,'cfael','DeallocateLocals')
   call a%Memor%dealloc(a%npoin,a%geoblock,'geoblock','DeallocateLocals')
   !if(a%kfl_ReadBlock)call a%Memor%dealloc(a%numGeoBlocks,a%blockName,'BlockName','DeallocateLocals')

   !FOR VMASS
   call a%Memor%dealloc(a%npoin,a%vmass,'VMASS','exnor')
   
   !ParallelOrdering
   call a%LocalOrdering%Deallocate(a%Memor)
   
   !ParallelCommunicator 
   call a%ArrayCommunicator%Deallocate
   
   call a%ParallelLibrary%DeallocateCommunicator(a%ArrayCommunicator,a%Memor)
   call a%ParallelLibrary%DeallocateOrdering(a%LocalOrdering,a%Memor)
   
   !Periodic boundary conditions
   if (allocated(a%MasterSlave)) call a%Memor%dealloc(a%npoin,a%MasterSlave,'MasterSlave','DeallocateLocals')
   if (allocated(a%SlaveList))   call a%Memor%dealloc(a%nslave,a%SlaveList,'SlaveList','DeallocateLocals')
   
   
   !Deallocate Projector system if it is allocated
   if (a%kfl_ProjectorReady == 1) then
      call a%L2ProjectorSystem%Deallocate
      call a%ParallelLibrary%DeallocateSystem(a%L2ProjectorSystem, a%Memor) 
      a%kfl_ProjectorReady = 0
   endif
   
   !Read Exnor
   !Allocate the arrays for storing the Exnor systems
   call a%Memor%dealloc(a%npoin,a%IsExnor,'isExnor','ReaExnor')
   call a%Memor%dealloc(a%npoin,a%ExternalNormal,'ExternalNormal','ReaExnor')
   
   if (a%kfl_HangingNodes .eqv. .true.) then
      call a%Memor%dealloc(size(a%pHangingList),a%pHangingList,'pHangingList','DeallocateLocals')
      call a%Memor%dealloc(size(a%lHangingList),a%lHangingList,'lHangingList','DeallocateLocals')
      call a%Memor%dealloc(size(a%rHangingList),a%rHangingList,'rHangingList','DeallocateLocals')
      
      call a%Memor%dealloc(size(a%iwaHanging),a%iwaHanging,'iwaHanging','DeallocateLocals')
      call a%Memor%dealloc(size(a%seenHanging),a%seenHanging,'seenHanging','DeallocateLocals')
   endif
   
   if (allocated(a%BoundaryNewToOld)) then
      call a%Memor%dealloc(size(a%BoundaryNewToOld),a%BoundaryNewToOld,'BoundaryNewToOld','DeallocateLocals')
   endif
   
   if (a%kfl_UseElementDataStructures) then
      call msh_FinalizeElementDataStructures(a)
   endif
   
end subroutine

subroutine MeshComputeRebalancingNumbering(a,PointGlobNumber,PointProcNumber)
   use typre
   use Mod_ParallelPartitionerInterface
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: PointGlobNumber(:), PointProcNumber(:)

   class(ParallelPartitionerInterface), pointer :: Partitioner => NULL()
   integer(ip), allocatable :: auxja(:)

   integer(ip), allocatable :: localid(:), globalid(:)
   character(6) :: GraphPartType
   
   integer(ip) :: ipoin
   
   !-----------------------------------------------------------------------------------------------------------
   !Partitionate
   call a%ParallelLibrary%CreateRebalancePartitioner(Partitioner,a%Memor)
   call Partitioner%Initialize(a%Memor)
   
   call Partitioner%GetGraphPartType(graphPartType)
   
   if (GraphPartType .eq. 'Graph') then
      call a%Memor%alloc(size(a%ja),auxja,'auxja','MeshComputeRebalancingNumbering')
      call a%LocalOrdering%Local2Global(size(a%ja),a%ja,auxja)
      
      call Partitioner%GraphPart(a%MPIcomm,a%npoinLocal,a%gnpoin,a%ia,auxja,PointProcNumber,PointGlobNumber)
      call a%Memor%dealloc(size(a%ja),auxja,'auxja','MeshComputeRebalancingNumbering')
   
   elseif (GraphPartType .eq. 'Geom') then
      call a%Memor%alloc(a%npoinlocal,localid,'localid','MeshComputeRebalancingNumbering')
      call a%Memor%alloc(a%npoinlocal,globalid,'globalid','MeshComputeRebalancingNumbering')
      
      do ipoin = 1,a%npoinLocal
         localid(ipoin) = ipoin
      enddo
      call a%LocalOrdering%Local2Global(size(localid),localid,globalid)
  
      call Partitioner%GraphPartGeom(a%MPIcomm,a%ndime,a%npoinlocal,a%gnpoin,localid,globalid,a%coord,pointProcNumber(1:a%npoinLocal),pointGlobNumber(1:a%npoinLocal))
   
      call a%Memor%dealloc(a%npoinlocal,localid,'localid','MeshComputeRebalancingNumbering')
      call a%Memor%dealloc(a%npoinlocal,globalid,'globalid','MeshComputeRebalancingNumbering')
   
   endif
   
   call a%ArrayCommunicator%GhostCommunicate(1_ip,PointProcNumber(1:a%npoin))
   call a%ArrayCommunicator%GhostCommunicate(1_ip,PointGlobNumber(1:a%npoin))
   
   call a%ParallelLibrary%DeallocatePartitioner(Partitioner,a%Memor)
   
   !-----------------------------------------------------------------------------------------------------------
end subroutine   

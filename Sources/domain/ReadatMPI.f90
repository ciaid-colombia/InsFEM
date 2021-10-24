subroutine ReadatMPI(a)
   use MPI
   use typre
   use Mod_Mesh
   implicit none
   
   class(FemMesh), intent(inout) :: a
   
   integer(ip) :: ierr,j
   
   !Communicate dimensions
   CALL MPI_BCAST(a%numGeoBlocks, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   !do j=1,a%numGeoBlocks
   !    CALL MPI_BCAST(a%blockName(j), 10, MPI_CHAR, a%MPIroot, a%MPIcomm, ierr)
   !enddo
   CALL MPI_BCAST(a%gnpoin, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%gnelem, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%gnskew, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%ndime, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%nelty, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%mitsm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_perio, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_ForcePeriodicMesh, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%ForcePeriodicMeshTol, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_ReadBlock, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%PerDime, 3, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_exnpo, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%mnode, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   CALL MPI_BCAST(a%tolsm, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%scale, 3, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_TestCommunications, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_MeshInputType, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_MeshInputType == 1) then
      CALL MPI_BCAST(a%ManufacturedNpoinX, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%ManufacturedType, 6, MPI_CHAR, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%kfl_ManuMeshType, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      if (a%kfl_ManuMeshType == 1) then
         CALL MPI_BCAST(a%ManufacturedInternalRadius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedExternalRadius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedDivisionsInAngle, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      elseif (a%kfl_ManuMeshType == 2) then
         CALL MPI_BCAST(a%ManufacturedInternalRadius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedExternalRadius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedAxialLength, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedDivisionsInAngle, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%ManufacturedNpoinZ, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      endif
   endif
   CALL MPI_BCAST(a%kfl_UseElementDataStructures, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
   
   !Allocate
   if (a%MPIrank /= a%MPIroot) then
      call a%Memor%alloc(a%nelty,a%nnode,'NNODE','ReadatMPI')
      call a%Memor%alloc(a%nelty,a%nface,'NFACE','ReadatMPI')
      call a%Memor%alloc(a%nelty,a%lquad,'LQUAD','ReadatMPI')
      call a%Memor%alloc(a%nelty,a%ngaus,'NGAUS','ReadatMPI')
      call a%Memor%alloc(a%mnode,a%ltype,'LTYPE','ReadatMPI')
   endif
   
   !Communicate arrays
   call MPI_BCAST(a%nnode,a%nelty,MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%nface,a%nelty,MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%lquad,a%nelty,MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%ngaus,a%nelty,MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   call MPI_BCAST(a%ltype,size(a%ltype),MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
end subroutine

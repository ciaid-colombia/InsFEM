!   subroutine NewDomain(itask,MeshName,Mesh1)
!   !-----------------------------------------------------------------------
!   !    This routine loads a mesh
!   !-----------------------------------------------------------------------
!      use typre
!      use def_master
!      use Mod_Mesh
!      use Mod_int2str
!      implicit none
!      class(FemMesh) :: Mesh1

!      integer(ip) :: itask
!      character(150) :: MeshName

!      !Initialize Mesh
!      if (itask == 1) then

!   
!         !Mesh Input data
!         call Mesh1%SetMPI(MPIsize,MPIroot,MPIrank)
!         call Mesh1%SetReadType(ReadTypeString)
!         call Mesh1%SetInputFolder(BaseDataFolder)
!         call Mesh1%SetInputFile(NameMesh)
!         call Mesh1%SetOutputFolder(ResultsFolder)
!         call Mesh1%SetOutputFiles(lun_memor,lun_outpu)
!         call Mesh1%SetParallelLibrary(ParallelLibrary)
!         call Mesh1%SetFlush(kfl_flush)

!         !Read data
!         call Mesh1%Readat
!      
!         !Initial operations (Graph and Scatter)
!         call Mesh1%Initialize

!         !Deallocate Read Globals
!         elseif (itask == 2) then

!         call Mesh1%DeallocateReadGlobals
!      
!         !Deallocate Locals
!         elseif (itask == 3) then
!         call Mesh1%Turnof
!      endif
!   end subroutine NewDomain

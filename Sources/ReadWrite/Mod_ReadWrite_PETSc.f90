module Mod_ReadWrite_PETSc
   use typre
   use Mod_ReadWrite
#include <petscversion.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>
   use petscvec
   implicit none
   private
   public Reader_PETSc, Writer_PETSc, Reader_PETSc_Const, Writer_PETSc_Const

   type, extends(Writer) :: Writer_PETSc
      Vec              :: auxVec

   contains

#ifdef HDF5
      procedure :: WriteAttributeReal
      procedure :: WriteAttributeInteger
      procedure :: WriteDataReal1
      procedure :: WriteDataReal2
      procedure :: WriteArray1
      procedure :: WriteArray2
      procedure :: WriteArray3
      procedure :: InitWriter
      procedure :: EndWriter
#endif

   end type

   type, extends(Reader) :: Reader_PETSc
      Vec              :: auxVec

   contains

#ifdef HDF5
      procedure :: ReadAttributeReal
      procedure :: ReadAttributeInteger
      procedure :: ReadDataReal1
      procedure :: ReadDataReal2
      procedure :: ReadArray1
      procedure :: ReadArray2
      procedure :: ReadArray3
      procedure :: InitReader
      procedure :: EndReader
#endif

   end type

   interface Writer_PETSc_Const
      procedure constructorW
   end interface Writer_PETSc_Const

   interface Reader_PETSc_Const
      procedure constructorR
   end interface Reader_PETSc_Const

contains

   function constructorW()
        class(Writer_PETSc), pointer :: constructorW

        allocate(constructorW)

    end function constructorW

   function constructorR()
        class(Reader_PETSc), pointer :: constructorR

        allocate(constructorR)

    end function constructorR

#ifdef HDF5
   subroutine InitWriter(a)
      implicit none
      class(Writer_PETSc) :: a
      PetscErrorCode   :: ierr

      call VecCreate(a%MPIcomm,a%auxVec,ierr)
      call VecSetFromOptions(a%auxVec,ierr)
      call VecSetSizes(a%auxVec,a%sizeLocal,a%sizeGlobal,ierr)
   end subroutine
   
   subroutine EndWriter(a)
      implicit none
      class(Writer_PETSc) :: a
      PetscErrorCode   :: ierr

      call VecDestroy(a%auxVec,ierr)
   end subroutine
   
   subroutine WriteAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      integer(ip)      :: Array(:)
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call PetscViewerHDF5WriteAttribute(viewer,PETSC_NULL_CHARACTER,nameV,PETSC_INT,Array,ierr)

      call PetscViewerDestroy(viewer,ierr)
   end subroutine
   
   subroutine WriteAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call PetscViewerHDF5WriteAttribute(viewer,PETSC_NULL_CHARACTER,nameV,PETSC_REAL,Array,ierr)

      call PetscViewerDestroy(viewer,ierr)
   end subroutine

   subroutine WriteDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,j,sg
      real(rp)         :: Array(:)
      Vec              :: auxVec
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      sg = size(Array,1)
      call PetscViewerCreate(MPI_COMM_SELF,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5Open(MPI_COMM_SELF,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call VecCreate(MPI_COMM_SELF,auxVec,ierr)
      call VecSetFromOptions(auxVec,ierr)
      call VecSetSizes(auxVec,sg,sg,ierr)

      call VecSetValues(auxVec,sg,(/(j,j=0,sg-1)/),Array,INSERT_VALUES,ierr)
      call VecAssemblyBegin(auxVec,ierr)
      call VecAssemblyEnd(auxVec,ierr)

      call PetscObjectSetName(auxVec,adjustl(trim(nameV)),ierr)
      call VecView(auxVec,viewer,ierr)

      call PetscViewerDestroy(viewer,ierr)
      call VecDestroy(auxVec,ierr)
   end subroutine

   subroutine WriteDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,j,sg(2),i
      real(rp)         :: Array(:,:)
      Vec              :: auxVec
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)
      call PetscViewerCreate(MPI_COMM_SELF,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5Open(MPI_COMM_SELF,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call VecCreate(MPI_COMM_SELF,auxVec,ierr)
      call VecSetFromOptions(auxVec,ierr)
      call VecSetSizes(auxVec,sg(1)*sg(2),sg(1)*sg(2),ierr)

      do i =1,sg(2)
         call VecSetValues(auxVec,sg(1),(/(j,j=(i-1)*sg(1),i*sg(1)-1)/),Array(:,i),INSERT_VALUES,ierr)
      end do
      call VecAssemblyBegin(auxVec,ierr)
      call VecAssemblyEnd(auxVec,ierr)

      call PetscObjectSetName(auxVec,adjustl(trim(nameV)),ierr)
      call VecView(auxVec,viewer,ierr)

      call PetscViewerDestroy(viewer,ierr)
      call VecDestroy(auxVec,ierr)
   end subroutine

   subroutine WriteArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,j,ind(a%sizeLocal)
      real(rp)         :: Array(:)
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr
      integer(ip), optional :: ncol

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)
      if (present(ncol)) call PetscViewerHDF5SetTimestep(viewer,ncol,ierr)

      ind = a%LocalToGlobal
      call VecSetValues(a%auxVec,a%sizeLocal,ind,Array,INSERT_VALUES,ierr)
      call VecAssemblyBegin(a%auxVec,ierr)
      call VecAssemblyEnd(a%auxVec,ierr)

      call PetscObjectSetName(a%auxVec,adjustl(trim(nameV)),ierr)
      call VecView(a%auxVec,viewer,ierr)

      call PetscViewerDestroy(viewer,ierr)
   end subroutine

   subroutine WriteArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG,string,nameV2
      integer(ip)      :: fil_unit,sg(2),i
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)

      if (sg(1) == a%sizeLocal) then
         do i =1,sg(2)
            write(string,"(I3)") i
            nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
            call a%WriteArray(fil_name,fil_unit,Array(:,i),nameV2,nameG,ncol)
         end do
      elseif (sg(2) == a%sizeLocal) then
         do i =1,sg(1)
            write(string,"(I3)") i
            nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
            call a%WriteArray(fil_name,fil_unit,Array(i,:),nameV2,nameG,ncol)
         end do
      endif
   end subroutine

   subroutine WriteArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG,string,nameV2
      integer(ip)      :: fil_unit,sg(3),icount,i,j
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)
      sg(3) = size(Array,3)

      if (sg(1) == a%sizeLocal) then
         do i =1,sg(2)
            do j =1,sg(3)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%WriteArray(fil_name,fil_unit,Array(:,i,j),nameV2,nameG,ncol)
            end do
         end do
      elseif (sg(2) == a%sizeLocal) then
         do i =1,sg(1)
            do j =1,sg(3)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%WriteArray(fil_name,fil_unit,Array(i,:,j),nameV2,nameG,ncol)
            end do
         end do
      elseif (sg(3) == a%sizeLocal) then
         do i =1,sg(1)
            do j =1,sg(2)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%WriteArray(fil_name,fil_unit,Array(i,j,:),nameV2,nameG,ncol)
            end do
         end do
      endif
   end subroutine

   subroutine InitReader(a)
      implicit none
      class(Reader_PETSc) :: a
      PetscErrorCode   :: ierr

      call VecCreate(a%MPIcomm,a%auxVec,ierr)
      call VecSetFromOptions(a%auxVec,ierr)
      call VecSetSizes(a%auxVec,a%sizeLocal,a%sizeGlobal,ierr)
      call VecAssemblyBegin(a%auxVec,ierr)
      call VecAssemblyEnd(a%auxVec,ierr)
   end subroutine
   
   subroutine EndReader(a)
      implicit none
      class(Reader_PETSc) :: a
      PetscErrorCode   :: ierr

      call VecDestroy(a%auxVec,ierr)
   end subroutine
   
   subroutine ReadAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      integer(ip)      :: Array(:)
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call PetscViewerHDF5ReadAttribute(viewer,PETSC_NULL_CHARACTER,nameV,PETSC_INT,Array,ierr)

      call PetscViewerDestroy(viewer,ierr)
   end subroutine
   
   subroutine ReadAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_APPEND,ierr)
      call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_APPEND,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call PetscViewerHDF5ReadAttribute(viewer,PETSC_NULL_CHARACTER,nameV,PETSC_REAL,Array(:),ierr)

      call PetscViewerDestroy(viewer,ierr)
   end subroutine

   subroutine ReadDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,j,sg(2),i
      real(rp)         :: Array(:,:)
      Vec              :: auxVec
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)
      call PetscViewerCreate(MPI_COMM_SELF,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_READ,ierr)
      call PetscViewerHDF5Open(MPI_COMM_SELF,fil_name,FILE_MODE_READ,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call VecCreate(MPI_COMM_SELF,auxVec,ierr)
      call VecSetFromOptions(auxVec,ierr)
      call VecSetSizes(auxVec,sg(1)*sg(2),sg(1)*sg(2),ierr)
      call VecAssemblyBegin(auxVec,ierr)
      call VecAssemblyEnd(auxVec,ierr)

      call PetscObjectSetName(auxVec,adjustl(trim(nameV)),ierr)
      call VecLoad(auxVec,viewer,ierr)
      do i =1,sg(2)
         call VecGetValues(auxVec,sg(1),(/(j,j=(i-1)*sg(1),i*sg(1)-1)/),Array(:,i),ierr)
      end do

      call PetscViewerDestroy(viewer,ierr)
      call VecDestroy(auxVec,ierr)
   end subroutine
   
   subroutine ReadDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,j,sg
      real(rp)         :: Array(:)
      Vec              :: auxVec
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      sg = size(Array,1)
      call PetscViewerCreate(MPI_COMM_SELF,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_READ,ierr)
      call PetscViewerHDF5Open(MPI_COMM_SELF,fil_name,FILE_MODE_READ,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)

      call VecCreate(MPI_COMM_SELF,auxVec,ierr)
      call VecSetFromOptions(auxVec,ierr)
      call VecSetSizes(auxVec,PETSC_DECIDE,sg,ierr)
      call VecAssemblyBegin(auxVec,ierr)
      call VecAssemblyEnd(auxVec,ierr)

      call PetscObjectSetName(auxVec,adjustl(trim(nameV)),ierr)
      call VecLoad(auxVec,viewer,ierr)
      call VecGetValues(auxVec,sg,(/(j,j=0,sg-1)/),Array,ierr)

      call PetscViewerDestroy(viewer,ierr)
      call VecDestroy(auxVec,ierr)
   end subroutine 

   subroutine ReadArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      integer(ip)      :: ipoin,irank,icount
      real(rp), pointer        :: arraypointer(:) => NULL()
      type(r1p), allocatable   :: ArrayINeed(:),ArrayOthersNeed(:)
      integer(ip), allocatable :: countPoint(:)
      integer, allocatable     :: irequest01(:),irequest02(:)
      integer, parameter       :: mtag1 = 1
      integer(ip), optional :: ncol
      PetscViewer      :: viewer
      PetscErrorCode   :: ierr

      call PetscViewerCreate(a%MPIcomm,viewer,ierr)
      call PetscViewerFileSetMode(viewer,FILE_MODE_READ,ierr)
      !call PetscViewerHDF5SetCollective(viewer,PETSC_TRUE,ierr)
      call PetscViewerHDF5Open(a%MPIcomm,fil_name,FILE_MODE_READ,viewer,ierr)
      call PetscViewerHDF5PushGroup(viewer,nameG,ierr)
      if (present(ncol)) call PetscViewerHDF5SetTimestep(viewer,ncol,ierr)

      call PetscObjectSetName(a%auxVec,adjustl(trim(nameV)),ierr)
      call VecLoad(a%auxVec,viewer,ierr)

      call VecGetArrayF90(a%auxVec,arraypointer,ierr)

      call a%Memor%alloc(a%MPIsize,ArrayOthersNeed,'ArrayOthersNeed','Distribute')
      call a%Memor%alloc(a%MPIsize,ArrayINeed,'ArrayINeed','Distribute')
      
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NpoinOthersNeed(irank+1),ArrayOthersNeed(irank+1)%a,'ArrayOthersNeed','ReadArray')
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NpoinINeed(irank+1),ArrayINeed(irank+1)%a,'ArrayINeed','ReadArray')
         endif
      enddo

      call a%Memor%alloc(a%MPIsize,countPoint,'countPoint','ReadArray')
      countPoint = 0
      !Read coordinates and fill the communications structure
      do ipoin = 1,a%sizeLocal
         countPoint(a%ProcNumber(ipoin)+1) = countPoint(a%ProcNumber(ipoin)+1)+1
         ArrayOthersNeed(a%ProcNumber(ipoin)+1)%a(countPoint(a%ProcNumber(ipoin)+1)) = arraypointer(ipoin)
      enddo
      call a%Memor%dealloc(a%MPIsize,countPoint,'countPoint','ReadArray')
      
      call a%Memor%alloc(a%MPIsize,irequest01,'irequest01','ReadArray')
      call a%Memor%alloc(a%MPIsize,irequest02,'irequest02','ReadArray')
      
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(ArrayOthersNeed(irank+1)%a,a%NpoinOthersNeed(irank+1),MPI_REAL8,irank,mtag1,a%MPIcomm,irequest01(irank+1),ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(ArrayINeed(irank+1)%a,a%NpoinINeed(irank+1),MPI_REAL8,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01,MPI_STATUSES_IGNORE,ierr)
      call MPI_WAITALL(a%MPIsize,irequest02,MPI_STATUSES_IGNORE,ierr)
         
      call a%Memor%dealloc(a%MPIsize,irequest01,'irequest01','ReadArray')
      call a%Memor%dealloc(a%MPIsize,irequest02,'irequest02','ReadArray')

      icount = 1
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            do ipoin = 1,a%NpoinINeed(irank+1)
               Array(a%GlobalToLocal(icount)) = ArrayINeed(irank+1)%a(ipoin)
               icount = icount +1
            enddo
         endif
      enddo
      
      !Deallocate buffers
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NpoinOthersNeed(irank+1),ArrayOthersNeed(irank+1)%a,'ArrayOthersNeed','ReadArray')
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NpoinINeed(irank+1),ArrayINeed(irank+1)%a,'ArrayINeed','ReadArray')
         endif
      enddo

      call a%Memor%dealloc(a%MPIsize,ArrayOthersNeed,'ArrayOthersNeed','Distribute')
      call a%Memor%dealloc(a%MPIsize,ArrayINeed,'ArrayINeed','Distribute')

      call VecRestoreArrayF90(a%auxVec,arraypointer,ierr)
      call PetscViewerDestroy(viewer,ierr)
   end subroutine

   subroutine ReadArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG,string,nameV2
      integer(ip)      :: fil_unit,sg(2),i
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)

      if (sg(1) == a%sizeLocal) then
         do i =1,sg(2)
            write(string,"(I3)") i
            nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
            call a%ReadArray(fil_name,fil_unit,Array(:,i),nameV2,nameG,ncol)
         end do
      elseif (sg(2) == a%sizeLocal) then
         do i =1,sg(1)
            write(string,"(I3)") i
            nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
            call a%ReadArray(fil_name,fil_unit,Array(i,:),nameV2,nameG,ncol)
         end do
      endif

   end subroutine

   subroutine ReadArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader_PETSc) :: a
      character(150)   :: fil_name,nameV,nameG,string,nameV2
      integer(ip)      :: fil_unit,sg(3),i,j
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol

      sg(1) = size(Array,1)
      sg(2) = size(Array,2)
      sg(3) = size(Array,3)

      if (sg(1) == a%sizeLocal) then
         do i =1,sg(2)
            do j =1,sg(3)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%ReadArray(fil_name,fil_unit,Array(:,i,j),nameV2,nameG,ncol)
            end do
         end do
      elseif (sg(2) == a%sizeLocal) then
         do i =1,sg(1)
            do j =1,sg(3)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%ReadArray(fil_name,fil_unit,Array(i,:,j),nameV2,nameG,ncol)
            end do
         end do
      elseif (sg(3) == a%sizeLocal) then
         do i =1,sg(1)
            do j =1,sg(2)
               write(string,"(I3,I3)") i,j
               nameV2 = adjustl(trim(nameV))//" "//adjustl(trim(string))
               call a%ReadArray(fil_name,fil_unit,Array(i,j,:),nameV2,nameG,ncol)
            end do
         end do
      endif
   end subroutine
#endif

end module Mod_ReadWrite_PETSc

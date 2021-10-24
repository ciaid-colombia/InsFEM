module Mod_ReadWrite
   use typre
   use Mod_Memor
#include <petscversion.h>
#if PETSC_VERSION_MINOR >= 8
#include <petsc/finclude/petsc.h>
   use MPI
   implicit none
#else
   implicit none
#include <petsc/finclude/petsc.h90>
#endif
   private
   public Reader, Writer

   type, abstract :: Writer
      type(MemoryMan), pointer   :: Memor => NULL()
      integer(ip) :: sizeLocal,sizeGlobal
      integer(ip), allocatable :: LocalToGlobal(:)
      integer(ip) :: MPIsize,MPIrank,MPIroot,MPIcomm

   contains
      procedure :: SetWriterMPI
      procedure :: SetWriter
      procedure :: DeallocateWriter
      procedure :: WriteAttributeReal
      procedure :: WriteAttributeInteger
      procedure :: WriteDataReal1
      procedure :: WriteDataReal2
      procedure :: WriteArray1
      procedure :: WriteArray2
      procedure :: WriteArray3
      procedure :: InitWriter
      procedure :: EndWriter

      generic :: WriteArray => WriteArray1,WriteArray2,WriteArray3
      generic :: WriteAttribute => WriteAttributeInteger,WriteAttributeReal
      generic :: WriteData => WriteDataReal1,WriteDataReal2

   end type

   type, abstract :: Reader
      type(MemoryMan), pointer   :: Memor => NULL()
      integer(ip) :: sizeLocal,sizeGlobal
      integer(ip), allocatable :: NpoinINeed(:),NpoinOthersNeed(:)
      integer(ip), allocatable :: GlobalToLocal(:),InitialToGlobal(:),ProcNumber(:)
      integer(ip) :: MPIsize,MPIrank,MPIroot,MPIcomm

   contains
      procedure :: SetReaderMPI
      procedure :: SetReader
      procedure :: SetG2L
      procedure :: DeallocateReader
      procedure :: ReadAttributeReal
      procedure :: ReadAttributeInteger
      procedure :: ReadDataReal1
      procedure :: ReadDataReal2
      procedure :: ReadArray1
      procedure :: ReadArray2
      procedure :: ReadArray3
      procedure :: InitReader
      procedure :: EndReader

      generic :: ReadArray => ReadArray1,ReadArray2,ReadArray3
      generic :: ReadAttribute => ReadAttributeInteger,ReadAttributeReal
      generic :: ReadData => ReadDataReal1,ReadDataReal2

   end type

contains

   subroutine SetWriterMPI(a,MPIcomm,MPIsize,MPIroot,MPIrank)
      implicit none
      class(Writer) :: a
      integer(ip) :: MPIsize,MPIrank,MPIroot,MPIcomm

      a%MPIsize = MPIsize
      a%MPIrank = MPIrank
      a%MPIroot = MPIroot
      a%MPIcomm = MPIcomm
   end subroutine

   subroutine SetWriter(a,sizeLocal,sizeGlobal,LocalToGlobal,Memor)
      implicit none
      class(Writer) :: a
      type(MemoryMan), target   :: Memor
      integer(ip)   :: sizeLocal,sizeGlobal,sg
      integer(ip)   :: LocalToGlobal(:)

      sg = size(LocalToGlobal,1)
      a%Memor => Memor
      a%sizeLocal = sizeLocal
      a%sizeGlobal = sizeGlobal

      call a%Memor%alloc(sg,a%LocalToGlobal,'LocalToGlobal','ReadWrite')
      a%LocalToGlobal = LocalToGlobal-1

      call a%InitWriter
   end subroutine

   subroutine DeallocateWriter(a)
      implicit none
      class(Writer) :: a
      integer(ip)   :: sg

      call a%EndWriter
      sg = size(a%LocalToGlobal,1)
      call a%Memor%dealloc(sg,a%LocalToGlobal,'LocalToGlobal','ReadWrite')
   end subroutine

   subroutine WriteAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      integer(ip)      :: Array(:)

      call runend('WriteData1 not defined')
   end subroutine

   subroutine WriteAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)

      call runend('WriteData0 not defined')
   end subroutine

   subroutine WriteDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:)

      call runend('WriteData1 not defined')
   end subroutine

   subroutine WriteDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)

      call runend('WriteData1 not defined')
   end subroutine

   subroutine WriteArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      integer(ip), optional :: ncol

      call runend('WriteArray1 not defined')
   end subroutine

   subroutine WriteArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol
      
      call runend('WriteArray2 not defined')
   end subroutine

   subroutine WriteArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Writer) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol
      
      call runend('WriteArray3 not defined')
   end subroutine

   subroutine InitWriter(a)
      implicit none
      class(Writer) :: a
      
      call runend('InitWriter not defined')
   end subroutine

   subroutine EndWriter(a)
      implicit none
      class(Writer) :: a
      
      call runend('EndWriter not defined')
   end subroutine

   subroutine SetReaderMPI(a,MPIcomm,MPIsize,MPIroot,MPIrank)
      implicit none
      class(Reader) :: a
      integer(ip) :: MPIsize,MPIrank,MPIroot,MPIcomm

      a%MPIcomm = MPIcomm
      a%MPIsize = MPIsize
      a%MPIrank = MPIrank
      a%MPIroot = MPIroot
   end subroutine

   subroutine SetReader(a,sizeLocal,sizeGlobal,GlobalToLocal,InitialToGlobal,ProcNumber,Memor)
      implicit none
      class(Reader) :: a
      type(MemoryMan), target  :: Memor
      integer(ip), allocatable :: countPoint(:)
      type(i1p), allocatable   :: PointsINeed(:),PointsOthersNeed(:)
      integer, allocatable     :: irequest01(:),irequest02(:)
      integer, parameter       :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer(ip)   :: sizeLocal,sizeGlobal,sg
      integer(ip)   :: GlobalToLocal(:),InitialToGlobal(:),ProcNumber(:)
      integer(ip) :: irank,ipoin,ierr,icount

      a%Memor => Memor
      a%sizeLocal = sizeLocal
      a%sizeGlobal = sizeGlobal
      sg = size(GlobalToLocal,1)

      call a%Memor%alloc(sg,a%GlobalToLocal,'GlobalToLocal','ReadWrite')
      
      call a%Memor%alloc(sg,a%InitialToGlobal,'InitialToGlobal','ReadWrite')
      a%InitialToGlobal = InitialToGlobal
      
      call a%Memor%alloc(sg,a%ProcNumber,'ProcNumber','ReadWrite')
      a%ProcNumber = ProcNumber
      
      call a%Memor%alloc(a%MPIsize,a%NpoinINeed,'NpoinINeed','ReadWrite')
      call a%Memor%alloc(a%MPIsize,a%NpoinOthersNeed,'NpoinOthersNeed','ReadWrite')
      
      !Count the number of points others need
      do ipoin = 1,a%sizeLocal
         a%NpoinOthersNeed(a%ProcNumber(ipoin)+1) = a%NpoinOthersNeed(a%ProcNumber(ipoin)+1)+1
      enddo
      call MPI_AllToAll(a%NpoinOthersNeed,1,MPI_INTEGER4,a%NpoinINeed,1,MPI_INTEGER4,a%MPIcomm,ierr)

      !Allocate communications structure
      call a%Memor%alloc(a%MPIsize,PointsOthersNeed,'PointsOthersNeed','Distribute')
      call a%Memor%alloc(a%MPIsize,PointsINeed,'PointsINeed','Distribute')

      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NpoinOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','Distribute')
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NpoinINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','Distribute')
         endif
      enddo

      call a%Memor%alloc(a%MPIsize,countPoint,'countPoint','Distribute')
      countPoint = 0
      !Read coordinates and fill the communications structure
      do ipoin = 1,a%sizeLocal
         countPoint(a%ProcNumber(ipoin)+1) = countPoint(a%ProcNumber(ipoin)+1)+1
         PointsOthersNeed(a%ProcNumber(ipoin)+1)%l(countPoint(a%ProcNumber(ipoin)+1)) = a%InitialToGlobal(ipoin)
      enddo
      call a%Memor%dealloc(a%MPIsize,countPoint,'countPoint','Distribute')
      
      call a%Memor%alloc(a%MPIsize,irequest01,'irequest01','Distribute')
      call a%Memor%alloc(a%MPIsize,irequest02,'irequest02','Distribute')
      
      irequest01 = MPI_REQUEST_NULL
      irequest02 = MPI_REQUEST_NULL
      !Communicate the coordinates
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(PointsOthersNeed(irank+1)%l,a%NpoinOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest01(irank+1),ierr)
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call MPI_IRECV(PointsINeed(irank+1)%l,a%NpoinINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest02(irank+1),ierr)
         endif
      enddo
      call MPI_WAITALL(a%MPIsize,irequest01,MPI_STATUSES_IGNORE,ierr)
      call MPI_WAITALL(a%MPIsize,irequest02,MPI_STATUSES_IGNORE,ierr)
         
      call a%Memor%dealloc(a%MPIsize,irequest01,'irequest01','Distribute')
      call a%Memor%dealloc(a%MPIsize,irequest02,'irequest02','Distribute')

      icount = 1
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            do ipoin = 1,a%NpoinINeed(irank+1)
               GlobalToLocal(icount) = PointsINeed(irank+1)%l(ipoin)
               icount = icount +1
            enddo
         endif
      enddo
      
      !Deallocate buffers
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NpoinOthersNeed(irank+1),PointsOthersNeed(irank+1)%l,'PointsOthersNeed','Distribute')
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NpoinINeed(irank+1),PointsINeed(irank+1)%l,'PointsINeed','Distribute')
         endif
      enddo
      call a%Memor%dealloc(a%MPIsize,PointsOthersNeed,'PointsOthersNeed','Distribute')
      call a%Memor%dealloc(a%MPIsize,PointsINeed,'PointsINeed','Distribute')

      call a%InitReader

   end subroutine

   subroutine SetG2L(a,GlobalToLocal)
      implicit none
      class(Reader) :: a
      integer(ip)   :: GlobalToLocal(:)

      a%GlobalToLocal = GlobalToLocal
      
   end subroutine

   subroutine DeallocateReader(a)
      implicit none
      class(Reader) :: a
      integer(ip)   :: sg

      call a%EndReader
      sg = size(a%GlobalToLocal,1)
      call a%Memor%dealloc(sg,a%GlobalToLocal,'GlobalToLocal','ReadWrite')
      call a%Memor%dealloc(sg,a%InitialToGlobal,'InitialToGlobal','ReadWrite')
      call a%Memor%dealloc(sg,a%ProcNumber,'ProcNumber','ReadWrite')
      call a%Memor%dealloc(a%MPIsize,a%NpoinINeed,'NpoinINeed','ReadWrite')
      call a%Memor%dealloc(a%MPIsize,a%NpoinOthersNeed,'NpoinOthersNeed','ReadWrite')
   end subroutine

   subroutine ReadAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      
      call runend('ReadData1 not defined')
   end subroutine

   subroutine ReadAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      integer(ip)      :: Array(:)
      
      call runend('ReadData0 not defined')
   end subroutine

   subroutine ReadDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:)
      
      call runend('ReadData1 not defined')
   end subroutine

   subroutine ReadDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      
      call runend('ReadData1 not defined')
   end subroutine

   subroutine ReadArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:)
      integer(ip), optional :: ncol
      
      call runend('ReadArray1 not defined')
   end subroutine

   subroutine ReadArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol
      
      call runend('ReadArray2 not defined')
   end subroutine

   subroutine ReadArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      implicit none
      class(Reader) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol
      
      call runend('ReadArray3 not defined')
   end subroutine

   subroutine InitReader(a)
      implicit none
      class(Reader) :: a
      
      call runend('InitReader not defined')
   end subroutine

   subroutine EndReader(a)
      implicit none
      class(Reader) :: a
      
      call runend('EndReader not defined')
   end subroutine

end module Mod_ReadWrite

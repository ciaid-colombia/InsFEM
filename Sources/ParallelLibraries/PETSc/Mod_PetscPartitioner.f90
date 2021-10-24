module Mod_PetscPartitioner
   use typre
   use Mod_ParallelPartitionerInterface
   implicit none
   private
   public PetscPartitioner1 
   
   type, extends(ParallelPartitionerInterface) :: PetscPartitioner1
      integer(ip) :: MPIcomm
contains  
   
      procedure :: GraphPart => GraphPartPETSc
      procedure :: GetGraphPartType => GetGraphPartTypePETSc
   end type

contains

   subroutine GetGraphPartTypePETSc(a,type)
      use typre
      implicit none
      class(PetscPartitioner1) :: a
      character(6) :: type
      
      !Type will be either: "Graph" or "Geom"
      
      type = "Graph"
   end subroutine

   subroutine GraphPartPETSc(a,MPIcomm,npoin,gnpoin,ia,ja,pointProcNumber,pointGlobNumber)
      use typre
#include <petsc/finclude/petsc.h>
      use petsc
      implicit none
      class(PetscPartitioner1) :: a
      integer(ip) :: npoin,gnpoin,MPIcomm
      integer(ip) :: ia(*)
      integer(ip) :: ja(*)
      integer(ip) :: pointProcNumber(*),pointGlobNumber(*)
      
      !Petsc
      Mat              ::  Adjacency
      MatPartitioning  ::  Partition
      IS               ::  isPointProcNumber, isPointGlobNumber
      MatPartitioningType :: parttype
      integer(ip), pointer :: p_is(:) => NULL()
      integer :: ierr

      integer(ip) :: inipoin,jasize
      
      a%MPIcomm = MPIcomm
      
      !Adjust for ParMetis
      inipoin = ia(1)
      jasize  = ia(npoin+1)-ia(1)
      
      ia(1:npoin+1) = ia(1:npoin+1) - inipoin
      ja(1:jasize) = ja(1:jasize)-1
      
      !Partitionate
      call MatCreateMPIAdj(a%MPIcomm,npoin,gnpoin,ia,ja,PETSC_NULL_INTEGER,Adjacency,ierr)
      call MatPartitioningCreate(a%MPIcomm,Partition,ierr)
      call MatPartitioningSetAdjacency(Partition,Adjacency,ierr);
      call MatPartitioningSetType(Partition,"parmetis",ierr)
      
      if (ierr /= 0) call runend('Partitionate: Partitioning library not found')
!      call MPI_BARRIER(a%MPIcomm,ierr)
      
      call MatPartitioningApply(Partition,isPointProcNumber,ierr)
      call MatPartitioningDestroy(Partition,ierr)
      call MatDestroy(Adjacency,ierr)
      call ISPartitioningToNumbering(isPointProcNumber,isPointGlobNumber,ierr);
      
      !extract the data
      call ISGetIndicesF90(isPointProcNumber,p_is,ierr)
      PointProcNumber(1:npoin) = p_is(1:npoin)
      call ISRestoreIndicesF90(isPointProcNumber,p_is,ierr)
    
      call ISGetIndicesF90(isPointGlobNumber,p_is,ierr)
      PointGlobNumber(1:npoin) = p_is(1:npoin)+1  !Adjust numbering to Fortran
      call ISRestoreIndicesF90(isPointGlobNumber,p_is,ierr)
   
      call ISDestroy(isPointProcNumber,ierr)
      call ISDestroy(isPointGlobNumber,ierr)
      
      ia(1:npoin+1) = ia(1:npoin+1) + inipoin
      ja(1:jasize) = ja(1:jasize) + 1
      
   end subroutine

end module

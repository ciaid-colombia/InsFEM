module Mod_ParallelPartitionerInterface
   use Mod_Memor
   implicit none
   private
   public ParallelPartitionerInterface

   type ParallelPartitionerInterface


      type(MemoryMan), pointer :: Memor => NULL()

   contains   
      procedure :: Initialize
      procedure :: GraphPart                              
      procedure :: GraphPartGeom

      procedure :: GetGraphPartType

   end type

contains

   subroutine GetGraphPartType(a,type)
      use typre
      implicit none
      class(ParallelPartitionerInterface) :: a
      character(6) :: type

      !Type will be either: "Graph" or "Geom"

      call runend('ParallelPartitionerInterface_ GetGraphPartType not defined')
   end subroutine

   subroutine GraphPart(a,MPIcomm,npoin,gnpoin,ia,ja,pointProcNumber,pointGlobNumber)
      use typre
      implicit none

      class(ParallelPartitionerInterface) :: a
      integer(ip) :: npoin,gnpoin,MPIcomm
      integer(ip) :: ia(*)
      integer(ip) :: ja(*)
      integer(ip) :: pointProcNumber(*),pointGlobNumber(*)

      call runend('ParallelPartitionerInterface: GraphPart not ready')
   end subroutine

   subroutine GraphPartGeom(a,MPIcomm,ndime,npoinlocal,gnpoin,localid,globalid,coord,pointProcNumber,pointGlobNumber)
      use typre
      implicit none

      class(ParallelPartitionerInterface) :: a
      integer(ip) :: ndime,npoinlocal,gnpoin,MPIcomm
      integer(ip), target :: localid(:),globalid(:),pointProcNumber(:),pointGlobNumber(:)
      real(rp), target :: coord(:,:)

      call runend('ParallelPartitionerInterface: GraphPartGeom not ready')
   end subroutine

   subroutine Initialize(a,Memor)
      implicit none
      class(ParallelPartitionerInterface) :: a
      type(MemoryMan), target :: Memor

      a%Memor => Memor
   end subroutine
end module

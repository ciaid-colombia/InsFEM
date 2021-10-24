module Mod_RebalanceCommunicatorInterface
   use typre
   use Mod_MPIObject
   implicit none
   private
   public RebalanceCommunicatorInterface

   type, abstract, extends(MPIObject) :: RebalanceCommunicatorInterface
contains
      procedure(Initialize), deferred :: Initialize
      procedure(CommunicateReal), deferred :: CommunicateReal
      procedure(Dealloc), deferred :: Deallocate
   end type

   abstract interface
      subroutine Initialize(a,OldNpoinLocal,NewNpoinLocal,RebalancePointNumbering,RebalanceProcessorList)
         use typre
         use Mod_LinkedList
         use MPI
         import RebalanceCommunicatorInterface
         implicit none
         class(RebalanceCommunicatorInterface) :: a
         integer(ip) :: OldNPoinLocal, NewNpoinLocal,RebalancePointNumbering(*),RebalanceProcessorList(*)
      end subroutine
      
      subroutine CommunicateReal(a,ndofn,OldArray,NewArray)
         use typre
         use MPI 
         import RebalanceCommunicatorInterface
         implicit none
         class(RebalanceCommunicatorInterface) :: a
         integer(ip) :: ndofn
         real(rp) :: OldArray(ndofn,*), NewArray(ndofn,*)
      end subroutine
      
      subroutine Dealloc(a)
         use typre
         import RebalanceCommunicatorInterface
         implicit none
         class(RebalanceCommunicatorInterface) :: a
      end subroutine
   end interface


end module

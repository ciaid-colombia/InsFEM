module Mod_ParallelCommunicatorInterface
   implicit none
   private
   public ParallelCommunicatorInterface
   
   type ParallelCommunicatorInterface
   private 
   
   
contains   
   
      procedure :: Init
      procedure :: GhostCommunicateInteger
      procedure :: GhostCommunicateReal
      procedure :: GhostCommunicateInteger1
      procedure :: GhostCommunicateReal1
      procedure :: GhostCommunicateLogical
      procedure :: GhostCommunicateLogical1
      generic :: GhostCommunicate => GhostCommunicateReal, &
                                     GhostCommunicateInteger, &
                                     GhostCommunicateReal1, &
                                     GhostCommunicateInteger1, &
                                     GhostCommunicateLogical, &
                                     GhostCommunicateLogical1 
      procedure :: Deallocate
      
!       procedure :: AssemblyGhostContributionsInteger
!       procedure :: AssemblyGhostContributionsReal
!       procedure :: AssemblyGhostContributionsInteger1
!       procedure :: AssemblyGhostContributionsReal1
!       generic :: AssemblyGhostContributions => AssemblyGhostContributionsReal, &
!                                      AssemblyGhostContributionsInteger, &
!                                      AssemblyGhostContributionsReal1, &
!                                      AssemblyGhostContributionsInteger1 
                                     
      
   
   end type

contains

   subroutine Init(a,MPIcomm,MPIsize,MPIrank,npoinLocal,npoinGhost,gnpoin,ParallelOrdering,Memor)
      use typre
      use Mod_Memor
      use Mod_ParallelOrderingInterface
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: MPIsize,MPIrank,MPIcomm
      class(ParallelOrderingInterface)     :: ParallelOrdering
      integer(ip) :: npoinLocal,npoinGhost,gnpoin
      type(MemoryMan), target :: Memor
   
      call runend('ParallelCommunicator, Init procedure not implemented')
   end subroutine
   
   subroutine GhostCommunicateReal(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      real(rp)    :: array(ndime,*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine GhostCommunicateReal1(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      real(rp),target    :: array(*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine GhostCommunicateInteger(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      integer(ip) :: array(ndime,*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine GhostCommunicateInteger1(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      integer(ip),target :: array(*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine GhostCommunicateLogical(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      logical :: array(ndime,*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine GhostCommunicateLogical1(a,ndime,array)
      use typre
      implicit none
      
      class(ParallelCommunicatorInterface) :: a
      integer(ip) :: ndime
      logical, target :: array(*)
      integer(ip) :: ierr
      
      call runend('ParallelCommunicator, Communicate procedure not implemented')
   
   end subroutine
   
   subroutine Deallocate(a)
      use typre
      implicit none
      class(ParallelCommunicatorInterface) :: a
      
      call runend('ParallelCommunicator, Deallocate procedure not implemented')
   
   end subroutine
   
   
!       subroutine AssemblyGhostContributionsReal(a,ndime,array)
!       use typre
!       implicit none
!       
!       class(ParallelCommunicatorInterface) :: a
!       integer(ip) :: ndime
!       real(rp)    :: array(ndime,*)
!       integer(ip) :: ierr
!       
!       call runend('ParallelCommunicator, Communicate procedure not implemented')
!    
!    end subroutine
!    
!    subroutine AssemblyGhostContributionsReal1(a,ndime,array)
!       use typre
!       implicit none
!       
!       class(ParallelCommunicatorInterface) :: a
!       integer(ip) :: ndime
!       real(rp),target    :: array(*)
!       integer(ip) :: ierr
!       
!       call runend('ParallelCommunicator, Communicate procedure not implemented')
!    
!    end subroutine
!    
!    subroutine AssemblyGhostContributionsInteger(a,ndime,array)
!       use typre
!       implicit none
!       
!       class(ParallelCommunicatorInterface) :: a
!       integer(ip) :: ndime
!       integer(ip) :: array(ndime,*)
!       integer(ip) :: ierr
!       
!       call runend('ParallelCommunicator, Communicate procedure not implemented')
!    
!    end subroutine
!    
!    subroutine AssemblyGhostContributionsInteger1(a,ndime,array)
!       use typre
!       implicit none
!       
!       class(ParallelCommunicatorInterface) :: a
!       integer(ip) :: ndime
!       integer(ip),target :: array(*)
!       integer(ip) :: ierr
!       
!       call runend('ParallelCommunicator, Communicate procedure not implemented')
!    
!    end subroutine


end module
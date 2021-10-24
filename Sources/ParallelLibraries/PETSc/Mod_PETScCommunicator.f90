module Mod_PETScCommunicator
   use typre
   use Mod_ParallelCommunicatorInterface
   use Mod_Memor
#include <petsc/finclude/petscvec.h>
   use petscvec
   implicit none
   private
   public PETScCommunicator
   
   type, extends(ParallelCommunicatorInterface) ::  PETScCommunicator
      
      integer(ip) :: npoinLocal, npoinGhost, gnpoin,MPIcomm
      integer(ip), allocatable :: GhostLocal2Global(:) 
      
      Vec       :: VectorCommunicator(20)
      logical   :: kfl_VectorInitialized(20) = .false.   
      logical   :: kfl_VectorInitializedInteger(20) = .false.
      real(rp)  :: dummyr
      
      integer(ip) ::               npoinLocalBoundary     !Number of local nodes which belong to a boundary element
      integer(ip), allocatable ::  lnodsLocalBoundary(:)  !List of local nodes which belong to a boundary element    
      type(r2p) :: Integer2RealArray(20)
      
      type(MemoryMan), pointer :: Memor => NULL()
   
contains   
   
      procedure :: Init => PETScInitCommunicator
      procedure :: GhostCommunicateReal => PETScGhostCommunicateReal
      procedure :: GhostCommunicateInteger => PETScGhostCommunicateInteger
      procedure :: GhostCommunicateReal1 => PETScGhostCommunicateReal1
      procedure :: GhostCommunicateInteger1 => PETScGhostCommunicateInteger1
      procedure :: GhostCommunicateLogical => PETScGhostCommunicateLogical
      procedure :: GhostCommunicateLogical1 => PETScGhostCommunicateLogical1
      procedure :: Deallocate => PETScDeallocate
   
   end type

contains

   subroutine PETScInitCommunicator(a,MPIcomm,MPIsize,MPIrank,npoinLocal,npoinGhost,gnpoin,ParallelOrdering,Memor)
      use typre
      use Mod_Memor
      use Mod_ParallelOrderingInterface
      use Mod_PETScOrdering
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: MPIrank,MPIsize,MPIcomm
      class(ParallelOrderingInterface)  :: ParallelOrdering
      integer(ip) :: npoinLocal,npoinGhost,gnpoin
      type(MemoryMan),target :: Memor
      
      integer(ip) :: ipoin,jpoin,poinj
      
      integer(ip), allocatable :: iwa(:)
      logical,     allocatable :: lwa(:)

      Vec :: UNKNO
      real(rp), allocatable :: realUNKNO(:)
      integer(ip), allocatable :: auxLocal2Global(:)
      integer(ip) :: ierr,iaux,ndofn,npoin
      
      a%Memor => Memor
      a%npoinLocalBoundary = 0
      a%npoinLocal = npoinLocal
      a%npoinGhost = npoinGhost
      a%gnpoin     = gnpoin
      npoin = npoinLocal+npoinGhost
      a%MPIcomm = MPIcomm
      
      !Build the list of the global numbering of the Ghost points
      call a%Memor%alloc(npoinGhost,a%GhostLocal2Global,'GhostLocal2Global','PETScCommunicator_Init')
      do ipoin = 1,a%npoinGhost
         a%GhostLocal2Global(ipoin) = a%npoinLocal + ipoin
      enddo
      call ParallelOrdering%Local2Global(a%npoinGhost,a%GhostLocal2Global,a%GhostLocal2Global)
      a%GhostLocal2Global = a%GhostLocal2Global - 1
      
      call a%Memor%alloc(npoin,auxLocal2Global,'auxLocal2Global','PetscCommunicator_Init')
      do ipoin = 1,npoin
         auxLocal2Global(ipoin) = ipoin
      enddo
      call ParallelOrdering%Local2Global(npoin,auxLocal2Global,auxLocal2Global)
      auxLocal2Global = auxLocal2Global-1
      call a%Memor%alloc(npoin,realUNKNO,'realUnkno','PETScCommunicator_Init')
      
      !Build the list of the local poins which need to be sent to other processes
      !Unkno PETSC vector
      ndofn = 1_ip
      call VecCreateMPIWithArray(MPIcomm,ndofn,npoinLocal*ndofn,gnpoin*ndofn,realUNKNO,UNKNO,ierr)
      call VecSetBlockSize(UNKNO,ndofn,ierr);
      iaux = 1
      call Memor%allocObj(ierr,'UNKNO','PETScSystem_Init',iaux)
      
      realUnkno(1:npoinLocal) = 0.0_rp
      realUnkno(npoinLocal+1:npoin) = 1.0_rp
      
      !call VecSetValues(UNKNO,npoinLocal,auxLocal2Global(1:npoinLocal),realUnkno(1:npoinLocal),ADD_VALUES,ierr)
      if (a%npoinGhost > 0) then
         call VecSetValues(UNKNO,npoinGhost,auxLocal2Global(npoinLocal+1),realUnkno(npoinLocal+1),ADD_VALUES,ierr)
      endif
      

      call VecAssemblyBegin(UNKNO,ierr);
      call VecAssemblyEnd(UNKNO,ierr);

      call a%Memor%alloc(npoinLocal,iwa,'iwa','PETScInit')
      call a%Memor%alloc(npoinLocal+npoinGhost,lwa,'lwa','PETScInit')
      
      do ipoin = 1,npoinLocal
         if (realUNKNO(ipoin) /= 0.0_rp) then
            a%npoinLocalBoundary = a%npoinLocalBoundary+1
            iwa(a%npoinLocalBoundary) = ipoin
            lwa(ipoin) = .true.
         endif 
      enddo
         
      call VecDestroy(UNKNO,ierr);
      call a%Memor%dealloc(npoin,auxLocal2Global,'auxLocal2Global','PetscCommunicator_Init')
      call Memor%deallocObj(ierr,'UNKNO','PETScSystem_Init',iaux)
      call a%Memor%dealloc(npoin,realUNKNO,'realUnkno','PETScCommunicator_Init')
      
      call a%Memor%alloc(a%npoinLocalBoundary,a%lnodsLocalBoundary,'lnodsLocalBoundary','PETScCommunicator_Init')
      a%lnodsLocalBoundary = iwa(1:a%npoinLocalBoundary)
      
      call a%Memor%dealloc(npoinLocal,iwa,'iwa','PETScInit')
      call a%Memor%dealloc(npoinLocal+npoinGhost,lwa,'lwa','PETScInit')
   end subroutine
   
   subroutine PETScGhostCommunicateReal1(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      real(rp), target    :: array(*)
      
      real(rp), pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call PETScGhostCommunicateReal(a,ndime,aux_array)
   end subroutine
   
   subroutine PETScGhostCommunicateInteger1(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      integer(ip),target :: array(*)
      
      integer(ip), pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call PETScGhostCommunicateInteger(a,ndime,aux_array)
   end subroutine
   
   subroutine PETScGhostCommunicateLogical1(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      logical,target :: array(*)
      
      logical, pointer :: aux_array(:,:) => NULL()
      
      aux_array(1:1,1:1) => array(1:1)
      call PETScGhostCommunicateLogical(a,ndime,aux_array)
   end subroutine
   
   
   subroutine PETScGhostCommunicateReal(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      real(rp)    :: array(ndime,*)
      integer(ip) :: ierr
      
      if (a%kfl_VectorInitialized(ndime) .eqv. .false.) then
         call VecCreateGhostBlockWithArray(a%MPIcomm,ndime,ndime*a%npoinLocal,ndime*a%gnpoin,a%npoinGhost,a%GhostLocal2Global,a%dummyr,a%VectorCommunicator(ndime),ierr)
         a%kfl_VectorInitialized(ndime) = .true.
      endif
     
      call VecPlaceArray(a%VectorCommunicator(ndime),array,ierr)
      
      call VecGhostUpdateBegin(a%VectorCommunicator(ndime),INSERT_VALUES,SCATTER_FORWARD,ierr)
      call VecGhostUpdateEnd  (a%VectorCommunicator(ndime),INSERT_VALUES,SCATTER_FORWARD,ierr)
      
      call VecResetArray(a%VectorCommunicator(ndime),ierr)
   
   end subroutine
   
   subroutine PETScGhostCommunicateInteger(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      integer(ip) :: array(ndime,*)
      
      !Allocate if it's the first time
      if (a%kfl_VectorInitializedInteger(ndime) .eqv. .false.) then
         call a%Memor%palloc(ndime,a%npoinLocal+a%npoinGhost,a%Integer2RealArray(ndime)%a,'Integer2RealArray','PETScGhostCommunicateInteger')
         a%kfl_VectorInitializedInteger(ndime) = .true.
      endif
      
      !Integer 2 Real, only local nodes which are to be communicated
      a%Integer2RealArray(ndime)%a(:,a%lnodsLocalBoundary) = array(:,a%lnodsLocalBoundary)
      
      !Communicate reals
      call PETScGhostCommunicateReal(a,ndime,a%Integer2RealArray(ndime)%a)
      
      !Update ghost points in the integer array
      array(:,a%npoinLocal+1:a%npoinLocal+a%npoinGhost) = nint(a%Integer2RealArray(ndime)%a(:,a%npoinLocal+1:a%npoinLocal+a%npoinGhost))
      
   end subroutine 
   
   subroutine PETScGhostCommunicateLogical(a,ndime,array)
      use typre
      implicit none
      
      class(PETScCommunicator) :: a
      integer(ip) :: ndime
      logical :: array(ndime,*)
      
      integer(ip) :: ipoin,poini,idime
      
      !Allocate if it's the first time
      if (a%kfl_VectorInitializedInteger(ndime) .eqv. .false.) then
         call a%Memor%palloc(ndime,a%npoinLocal+a%npoinGhost,a%Integer2RealArray(ndime)%a,'Integer2RealArray','PETScGhostCommunicateInteger')
         a%kfl_VectorInitializedInteger(ndime) = .true.
      endif
      
     !Integer 2 Real, only local nodes which are to be communicated
      do poini = 1,a%npoinLocalBoundary
         ipoin = a%lnodsLocalBoundary(poini)
         do idime = 1,ndime
            if (array(idime,ipoin) .eqv. .true.) then
               a%Integer2RealArray(ndime)%a(idime,a%lnodsLocalBoundary) = 1.0_rp
            else
               a%Integer2RealArray(ndime)%a(idime,a%lnodsLocalBoundary) = 0.0_rp
            endif
         enddo
      enddo
      
      !Communicate reals
      call PETScGhostCommunicateReal(a,ndime,a%Integer2RealArray(ndime)%a)
      
      !Update ghost points in the integer array
      do ipoin = a%npoinLocal+1,a%npoinLocal+a%npoinGhost
         do idime = 1,ndime
            if (a%Integer2RealArray(ndime)%a(idime,ipoin) == 1.0_rp) then
               array(idime,ipoin) = .true.
            else
               array(idime,ipoin) = .false.
            endif
         enddo
      enddo
   end subroutine 
      
   subroutine PETScDeallocate(a)
      implicit none
      class(PETScCommunicator) :: a
      integer(ip) :: ierr,ndime
      
      call a%Memor%dealloc(a%npoinGhost,a%GhostLocal2Global,'GhostLocal2Global','PETScCommunicator_Deallocate')
      call a%Memor%dealloc(a%npoinLocalBoundary,a%lnodsLocalBoundary,'lnodsLocalBoundary','PETScCommunicator_Deallocate')
      
      do ndime = 1,20
         if (a%kfl_VectorInitialized(ndime) .eqv. .true.) then
            call VecDestroy(a%VectorCommunicator(ndime),ierr)
            a%kfl_VectorInitialized(ndime) = .false.
         endif
         if (a%kfl_VectorInitializedInteger(ndime) .eqv. .true.) then
            call a%Memor%pdealloc(ndime,a%npoinLocal+a%npoinGhost,a%Integer2RealArray(ndime)%a,'Integer2RealArray','PETScGhostCommunicateInteger')
            a%kfl_VectorInitializedInteger(ndime) = .false.
         endif
      enddo   
      
      a%npoinLocalBoundary = 0
      a%npoinLocal = 0
      a%npoinGhost = 0
      a%gnpoin = 0
   
   end subroutine
   
end module

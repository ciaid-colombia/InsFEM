module Mod_PointToProc
   use typre
   use Mod_Memor
   implicit none
   private
   public PointToProc
   
   type :: PointToProc
      !MPI
      integer(ip) :: MPIsize, MPIrank, MPIroot,MPIcomm
      !Memor
      type(MemoryMan), pointer :: Memor => NULL()
      integer(ip), allocatable :: poin0(:)

contains   
   
      procedure :: SetMPI
      procedure :: SetMemor
      procedure :: Initialize
      procedure :: InitializeFromList
      procedure :: GetProc
      procedure :: Finalize
   end type

contains
   
   subroutine SetMPI(a,comm,size,root,irank)
      implicit none
      class(PointToProc) :: a
      integer(ip) :: irank,size,root,comm
      
      a%MPIrank = irank
      a%MPIsize = size
      a%MPIroot = root
      a%MPIcomm = comm
   end subroutine
   
   subroutine SetMemor(a,Memor)
      implicit none
      class(PointToProc) :: a
      type(MemoryMan), target :: Memor
      
      a%Memor => Memor
   end subroutine
   
   subroutine Initialize(a,poin0)
      use MPI
      implicit none
      class(PointToProc) :: a
      integer(ip) :: poin0
      
      integer(ip), allocatable :: auxpoin0(:)
      integer(ip) :: ierr
      
      call a%Memor%alloc(a%MPIsize,a%poin0,'poin0','InitializePointToProc')
      !Allocate working
      call a%Memor%alloc(a%MPIsize,auxpoin0,'poin0','InitializePointToProc')
      
      !Communicate the poin0 of each processor
      auxpoin0 = poin0
      call MPI_AllToAll(auxpoin0, 1, MPI_INTEGER4, a%poin0, 1, MPI_INTEGER4, a%MPIcomm,ierr);
      
      !Deallocate working
      call a%Memor%dealloc(a%MPIsize,auxpoin0,'poin0','InitializePointToProc')

   end subroutine
   
   subroutine InitializeFromList(a,ListNpoinLocal)
      implicit none
      class(PointToProc) :: a
      integer(ip) :: ListNpoinLocal(*)
      
      integer(ip) :: irank
      
      if (.not. allocated(a%poin0)) then
         call a%Memor%alloc(a%MPIsize,a%poin0,'poin0','InitializePointToProc')
      endif
      
      !Allocate working
      a%poin0(1) = 0
      do irank = 0,a%MPIsize-2
         a%poin0(irank+2) = a%poin0(irank+1)+ListNpoinLocal(irank+1)
      enddo
   end subroutine
      
   
   subroutine GetProc(a,gipoin,irank)
      implicit none
      class(PointToProc) :: a
      integer(ip) :: gipoin
      integer(ip) :: irank
      
      integer(ip) :: rankStep,StepDirection
      
      integer(ip) :: ipass
      
      irank = -1
      StepDirection = 1
      ipass = 0
      
      !This is an order(log(MPIsize)) algorithm for finding the processor to which a node belongs
      !We start in the mid processor
      do while (StepDirection /= 0)
         ipass = ipass+1
         rankStep = a%MPIsize/2**ipass+1
         
         irank = irank + rankStep*StepDirection
         
         if (irank < 0) irank = 0
         if (irank > a%MPIsize -1) irank = a%MPIsize-1
         
         call CheckFound(a,irank,gipoin,StepDirection)
      end do
      
   end subroutine
   
   subroutine CheckFound(a,irank,gipoin,StepDirection)
      implicit none
      class(PointToProc) :: a
      integer(ip) :: gipoin
      integer(ip) :: irank,StepDirection
      
      if (irank < a%MPIsize-1) then
         if (gipoin > a%poin0(irank+1) .and. gipoin <= a%poin0(irank+2)) then
            StepDirection = 0
         elseif (gipoin <= a%poin0(irank+1)) then
            StepDirection = -1
         elseif (gipoin > a%poin0(irank+2)) then
            StepDirection = 1
         endif
      elseif (irank == a%MPIsize -1) then
         if (gipoin > a%poin0(irank+1)) then
            StepDirection = 0
         else
            StepDirection = -1
         endif
      endif
   end subroutine

   subroutine Finalize(a)
      implicit none
      class(PointToProc) :: a

      !Deallocate
      call a%Memor%dealloc(a%MPIsize,a%poin0,'poin0','InitializePointToProc')
   end subroutine
end module

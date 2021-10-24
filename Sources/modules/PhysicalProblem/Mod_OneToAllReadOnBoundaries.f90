module Mod_OneToAllReadOnBoundaries
   use typre
   use Mod_Listen
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop
   implicit none
   private
   public OneToAllReadOnBoundaries
   
   character(150), parameter :: outstr = 'OneToAllReadOnBoundaries'

   type, extends(OneToAllLoop) :: OneToAllReadOnBoundaries
      
      class(PhysicalProblem), pointer :: Problem => NULL()
      class(FemMesh), pointer         :: Mesh => NULL()
      type(ListenFile), pointer       :: Listener => NULL()
      
      !Number of boundaries read
      integer(ip) :: cboun
      integer(ip), allocatable :: iBounFlag(:)
      
      !Working arrays
      integer(ip), allocatable :: lnode(:), aux_lnode(:)
      
      !Periodic boundary conditions
      integer(ip) :: kfl_perio
      
   
contains
      procedure :: Initialize
      procedure :: Finalize

      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheck 
      procedure :: SpecificRootAddToBuffer => RootAddToBuffer
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperateOnBoundaries
   
   end type
   
contains
   
   subroutine Initialize(a,Problem)
      implicit none
      class(OneToAllReadOnBoundaries) :: a
      class(PhysicalProblem),  target :: Problem
      
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      integer(ip) :: mnode
      
      a%Problem  => Problem
      a%Mesh     => Problem%Mesh
      a%Listener => Problem%Listener
      
      !Set MPI and Memor
      call a%SetMPI(Problem%MPIcomm,Problem%MPIsize,Problem%MPIroot,Problem%MPIrank)
      call a%SetMemor(Problem%Memor)
      
      !SetBufferDimensions
      !500 reads per processor or 25 Mb buffer for the listener%param
      BufMaxNsends = 1000
      BufMaxSize = 25e6/a%MPIsize/rp
      BufSafetySend = 5e2
      call a%SetBufferSizes(BufMaxNsends,BufMaxSize,BufSafetySend)
      
      !Allocation of working arrays
      call a%Mesh%GetMnode(mnode)
      !Allocation of working arrays
      call a%Memor%alloc(mnode,a%lnode,'lnode',outstr)
      call a%Memor%alloc(mnode,a%aux_lnode,'aux_lnode',outstr)
      call a%Memor%alloc(a%MPIsize,a%iBounFlag,'iBounFlag',outstr)
      a%iBounFlag = -1
      a%cboun = 0
      
      !Periodic boundary conditions
      call a%Mesh%GetPerio(a%kfl_perio)
   end subroutine
   
   subroutine Finalize(a)
      implicit none
      class(OneToAllReadOnBoundaries) :: a
      
      integer(ip) :: mnode
      
      !Allocation of working arrays
      call a%Mesh%GetMnode(mnode)
      !Allocation of working arrays
      call a%Memor%dealloc(mnode,a%lnode,'lnode',outstr)
      call a%Memor%dealloc(mnode,a%aux_lnode,'aux_lnode',outstr)
      call a%Memor%dealloc(a%MPIsize,a%iBounFlag,'iBounFlag',outstr)
   end subroutine
   
   subroutine RootGetDataAndCompletionCheck(a)
      implicit none
      class(OneToAllReadOnBoundaries) :: a
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'ENDON') a%kfl_CompletionCheck = .true.
   end subroutine
   
   subroutine RootAddToBuffer(a)
      implicit none
      class(OneToAllReadOnBoundaries) :: a
      
      integer(ip) :: cproc    !Number of processors to which this boundary belongs
      
      integer(ip) :: inode,pnode, ipoin, irank
      
      a%cboun = a%cboun +1 !Total global Number of boundaries read
      cproc = 0        !Number of processes to which this boundary belongs   
      !Read boundary nodes
      pnode=int(a%Listener%param(2))
      a%lnode(1:pnode)=int(a%Listener%param(3:2+pnode))

      !Decide to which processes the boundary belongs and build processor list
      do inode = 1,pnode
         ipoin = a%lnode(inode)
         
         !Periodic boundary conditions
         if (a%kfl_perio == 1) then
            !Send the boundary condition to the domain
            !of my master if I am slave
            call a%Mesh%Initial2MasterInitial(ipoin,ipoin)
         endif   
         
         call a%Mesh%Initial2ProcNumber(ipoin,irank)
         if (a%iBounFlag(irank+1) /= a%cboun) then
            a%iBounFlag(irank+1) = a%cboun
            cproc = cproc+1

            !To parallel numbering for sending
            call a%Mesh%Initial2Global(pnode,a%lnode,a%aux_lnode)  
            a%Listener%param(3:2+pnode) = a%aux_lnode(1:pnode)

            !We add Listener to the corresponding buffer
            call a%Buffer%AddToBuffer(irank,a%Listener%nnpar,a%Listener%param)
         endif
      enddo
   end subroutine
   
   subroutine GetDataFromBufferAndOperateOnBoundaries(a)
      implicit none
      class(OneToAllReadOnBoundaries) :: a    
      
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on boundaries information
      !call php_ReadOnBoundaries(a%Problem,a%lnode)
      call a%Problem%ReadOnBoundaries(a%lnode)
   end subroutine
   
end module

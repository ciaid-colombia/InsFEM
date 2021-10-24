module Mod_OneToAllReadOnNodes
   use typre
   use Mod_Listen 
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop
   implicit none
   private
   public OneToAllReadOnNodes
   
   character(150), parameter :: outstr = 'OneToAllReadOnNodes'

   type, extends(OneToAllLoop) :: OneToAllReadOnNodes
      
      class(PhysicalProblem), pointer :: Problem => NULL()
      class(FemMesh), pointer         :: Mesh => NULL()
      type(ListenFile), pointer       :: Listener => NULL()

   
contains
      procedure :: InitializeOnNodes
      procedure :: FinalizeOnNodes

      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheckOnNodes
      procedure :: SpecificRootAddToBuffer => RootAddToBufferOnNodes
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperateOnNodes
   
   end type
   
contains
   
   subroutine InitializeOnNodes(a,Problem)
      implicit none
      class(OneToAllReadOnNodes) :: a
      class(PhysicalProblem),  target :: Problem
      
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      integer(ip) :: mnodb
      
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
      
      
   end subroutine
   
   subroutine FinalizeOnNodes(a)
      implicit none
      class(OneToAllReadOnNodes) :: a
      
      integer(ip) :: mnodb
      
   end subroutine
   
   subroutine RootGetDataAndCompletionCheckOnNodes(a)
      implicit none
      class(OneToAllReadOnNodes) :: a
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'ENDON') a%kfl_CompletionCheck = .true.
   end subroutine
   
   subroutine RootAddToBufferOnNodes(a)
      implicit none
      class(OneToAllReadOnNodes) :: a      
      
      integer(ip) :: irank,cpoin,wpoin
      
      !readpoin number
      cpoin=int(a%Listener%param(1))
      !Decide to which processes the boundary belongs and build processor list          
      call a%Mesh%Initial2ProcNumber(cpoin,irank)
      !To parallel numbering for sending
      call a%Mesh%Initial2Global(cpoin,wpoin)    
      a%Listener%param(1) = wpoin   
      !We add Listener to the corresponding buffer
      call a%Buffer%AddToBuffer(irank,a%Listener%nnpar,a%Listener%param)

   end subroutine
   
   subroutine GetDataFromBufferAndOperateOnNodes(a)
      implicit none
      class(OneToAllReadOnNodes) :: a      
   
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on nodes information
      call a%Problem%ReadOnNodes
   end subroutine
   
end module

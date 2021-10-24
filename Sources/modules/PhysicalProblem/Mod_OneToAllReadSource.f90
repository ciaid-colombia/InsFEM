module Mod_OneToAllReadSource
   use typre
   use Mod_Listen 
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop
   implicit none
   private
   public OneToAllReadSource
   
   character(150), parameter :: outstr = 'OneToAllReadSource'

   type, extends(OneToAllLoop) :: OneToAllReadSource
      
      class(PhysicalProblem), pointer :: Problem => NULL()
      class(FemMesh), pointer         :: Mesh => NULL()
      type(ListenFile), pointer       :: Listener => NULL()

   
contains
      procedure :: InitializeSource
      procedure :: FinalizeSource

      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheckSource
      procedure :: SpecificRootAddToBuffer => RootAddToBufferSource
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperateSource
   
   end type
   
contains
   
   subroutine InitializeSource(a,Problem)
      implicit none
      class(OneToAllReadSource) :: a
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
   
   subroutine FinalizeSource(a)
      implicit none
      class(OneToAllReadSource) :: a
      
      integer(ip) :: mnodb
      
   end subroutine
   
   subroutine RootGetDataAndCompletionCheckSource(a)
      implicit none
      class(OneToAllReadSource) :: a
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'ENDSO') a%kfl_CompletionCheck = .true.
   end subroutine
   
   subroutine RootAddToBufferSource(a)
      implicit none
      class(OneToAllReadSource) :: a      
      
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
   
   subroutine GetDataFromBufferAndOperateSource(a)
      implicit none
      class(OneToAllReadSource) :: a      
   
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on nodes information
      call a%Problem%ReadSource
   end subroutine
   
end module

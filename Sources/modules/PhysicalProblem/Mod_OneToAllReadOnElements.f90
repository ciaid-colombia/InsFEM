module Mod_OneToAllReadOnElements
   use typre
   use Mod_Listen
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop
   implicit none
   private
   public OneToAllReadOnElements
   
   character(150), parameter :: outstr = 'OneToAllReadOnElements'

  type, extends(OneToAllLoop) :: OneToAllReadOnElements
      
      class(PhysicalProblem), pointer :: Problem => NULL()
      class(FemMesh), pointer         :: Mesh => NULL()
      type(ListenFile), pointer       :: Listener => NULL()
      
      
   
contains
      procedure :: Initialize
      procedure :: Finalize

      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheck 
      procedure :: SpecificRootAddToBuffer => RootAddToBuffer
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperateOnElements
   
   end type
   
contains

   subroutine Initialize(a,Problem)
      implicit none
      class(OneToAllReadOnElements) :: a
      class(PhysicalProblem),  target :: Problem
      
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      
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
   
   subroutine Finalize(a)
      implicit none
      class(OneToAllReadOnElements) :: a
      
      integer(ip) :: mnode
   end subroutine

   subroutine RootGetDataAndCompletionCheck(a)
      implicit none
      class(OneToAllReadOnElements) :: a
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'ENDON') a%kfl_CompletionCheck = .true.
   end subroutine

   subroutine RootAddToBuffer(a)
      implicit none
      class(OneToAllReadOnElements) :: a
      
      integer(ip) ::  irank, ranki, InitialElement,InitialElementSend
      
     
      InitialElement = int(a%Listener%param(1))
      
      if (a%Problem%kfl_Readtype == 0) then
         InitialElementSend = InitialElement
      elseif (a%Problem%kfl_Readtype == 1) then
         call runend('Read on elements implemented but not tested for kfl_readtype = 1, comment this line and check')
         call a%Mesh%ElementNaive2Initial%Global2Local(InitialElement,InitialElementSend)
      endif   
      
      if (associated(a%Mesh%gElementInitialSendToProcessor(InitialElementSend)%l)) then
         do ranki = 1,size(a%Mesh%gElementInitialSendToProcessor(InitialElementSend)%l)
            irank = a%Mesh%gElementInitialSendToProcessor(InitialElementSend)%l(ranki)
            
            call a%Buffer%AddToBuffer(irank,a%Listener%nnpar,a%Listener%param)
         enddo   
      endif
   end subroutine
   
  
   
   subroutine GetDataFromBufferAndOperateOnElements(a)
      implicit none
      class(OneToAllReadOnElements) :: a    
      
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on Elements information
      call a%Problem%ReadOnElements
   end subroutine
   
end module

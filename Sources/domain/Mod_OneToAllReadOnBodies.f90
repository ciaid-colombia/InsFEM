module Mod_OneToAllReadOnBodies
   use typre
   use Mod_Mesh
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop
   use Mod_Listen
   use MPI
   implicit none
   private
   public OneToAllReadOnBodies
   
   character(150), parameter :: outstr = 'OneToAllReadOnBodies'

   type, extends(OneToAllLoop) :: OneToAllReadOnBodies
      
      class(FemMesh), pointer         :: Mesh => NULL()
      type(ListenFile), pointer       :: Listener => NULL()
      
      !Number of boundaries read
      integer(ip) :: cboun
      integer(ip), allocatable :: iBounFlag(:)
      
      !Working arrays
      integer(ip), allocatable :: knodbo(:), aux_knodbo(:)
      
      
   
contains
      procedure :: Initialize
      procedure :: Finalize

      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheckBody
      procedure :: SpecificRootAddToBuffer => RootAddToBufferBody
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperateBody
   
   end type
   
contains
   
   subroutine Initialize(a,Mesh)
      implicit none
      class(OneToAllReadOnBodies) :: a
      class(FemMesh),  target :: Mesh      
      
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      integer(ip) :: mnodb
      
          
      a%Mesh  => Mesh
      a%Listener => Mesh%Listener      

      Mesh%nbody=0
      
      !Set MPI and Memor
      call a%SetMPI(Mesh%MPIcomm,Mesh%MPIsize,Mesh%MPIroot,Mesh%MPIrank)
      call a%SetMemor(Mesh%Memor)
      
      !SetBufferDimensions
      !500 reads per processor or 25 Mb buffer for the listener%param
      BufMaxNsends = 1000
      BufMaxSize = 25e6/a%MPIsize/rp
      BufSafetySend = 5e2
      call a%SetBufferSizes(BufMaxNsends,BufMaxSize,BufSafetySend)
      
      !Allocation of working arrays
      !Now we have Mesh=FemMesh
      call a%Memor%alloc(Mesh%mnodb,a%knodbo,'knodbo',outstr)
      call a%Memor%alloc(Mesh%mnodb,a%aux_knodbo,'aux_knodbo',outstr)
      call a%Memor%alloc(a%MPIsize,a%iBounFlag,'iBounFlag',outstr)
      a%iBounFlag = -1
      a%cboun = 0      

   end subroutine
   
   subroutine Finalize(a)
      implicit none
      class(OneToAllReadOnBodies) :: a
      integer(ip) :: ierr
      
      CALL MPI_BCAST(a%Mesh%nbody, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

      !Allocation of working arrays
      call a%Memor%dealloc(a%mesh%mnodb,a%knodbo,'knodbo',outstr)
      call a%Memor%dealloc(a%mesh%mnodb,a%aux_knodbo,'aux_knodbo',outstr)
      call a%Memor%dealloc(a%MPIsize,a%iBounFlag,'iBounFlag',outstr)
        
      
   end subroutine
   
   subroutine RootGetDataAndCompletionCheckbody(a)
      implicit none
      class(OneToAllReadOnBodies) :: a
      integer(ip) :: nbody,ierr      
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'BODYE') a%kfl_CompletionCheck = .true.
   end subroutine
   
   subroutine RootAddToBufferbody(a)
      implicit none
      class(OneToAllReadOnBodies) :: a
      
      integer(ip) :: cproc    !Number of processors to which this boundary belongs
      
      integer(ip) :: inodb,pnodb,ipoin,irank,ibody
      
      a%cboun = a%cboun +1 !Total global Number of boundaries read
      cproc = 0        !Number of processes to which this boundary belongs   
      !Read boundary nodes
      pnodb=int(a%Listener%param(2))
      ibody=int(a%Listener%param(2+pnodb+1))
      a%knodbo(1:pnodb)=int(a%Listener%param(3:2+pnodb))

      !Decide to which processes the boundary belongs and build processor list
      do inodb = 1,pnodb
         ipoin = a%knodbo(inodb)
  
         
         call a%Mesh%Initial2ProcNumber(ipoin,irank)
         if (a%iBounFlag(irank+1) /= a%cboun) then
            a%iBounFlag(irank+1) = a%cboun
            cproc = cproc+1
            
            !Definition of nbody
            ibody = a%Listener%param(2+pnodb+1)            
            a%Mesh%nbody=max(a%Mesh%nbody,ibody)

            !To parallel numbering for sending
            call a%Mesh%Initial2Global(pnodb,a%knodbo,a%aux_knodbo)  
            a%Listener%param(3:2+pnodb) = a%aux_knodbo
            
            !We add Listener to the corresponding buffer
            call a%Buffer%AddToBuffer(irank,a%Listener%nnpar,a%Listener%param)
         endif
      enddo
   end subroutine
   
   subroutine GetDataFromBufferAndOperateBody(a)
      implicit none
      class(OneToAllReadOnBodies) :: a
      
      interface
      
         subroutine ReadOnBodies(a,knodbo)
            use typre
            use Mod_Mesh
            use Mod_Listen
            implicit none
            class(FemMesh) :: a
            !type(ListenFile), pointer       :: Listener => NULL() 
            integer(ip) :: knodbo(*)            

         end subroutine      
         
      end interface        
      
      
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on nodes information
      call ReadOnBodies(a%Mesh,a%knodbo)
               
     
   end subroutine
   
end module

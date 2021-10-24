module Mod_OneToAllReadExnor
   use typre
   use Mod_Listen
   use Mod_Mesh
   use Mod_OneToAllBuffer
   use Mod_OneToAllLoop
   implicit none
   private
   public OneToAllReadExnor

   type, extends(OneToAllLoop)  :: OneToAllReadExnor
      class(FemMesh), pointer   :: Mesh => NULL()
      type(ListenFile), pointer :: Listener => NULL()
      
      integer(ip) :: ndime
      
      integer(ip) :: AllocCount
      
contains
      procedure :: Initialize
      procedure :: Finalize 
      
      procedure :: SpecificRootGetDataAndCompletionCheck => RootGetDataAndCompletionCheck
      procedure :: SpecificRootAddToBuffer => RootAddToBuffer
      procedure :: SpecificGetDataFromBufferAndOperate => GetDataFromBufferAndOperate
   end type
   
   character(150) :: outstr = 'OneToAllReadExnor'

contains

   subroutine Initialize(a,Mesh)
      implicit none
      class(OneToAllReadExnor) :: a
      class(FemMesh), target :: Mesh
      
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      
      a%Mesh => Mesh
      a%Listener => Mesh%Listener
      call a%Mesh%GetNdime(a%ndime)
      
      !Set MPI and Memor
      call a%SetMPI(Mesh%MPIcomm,Mesh%MPIsize,Mesh%MPIroot,Mesh%MPIrank)
      call a%SetMemor(Mesh%Memor)
      
      !SetBufferDimensions
      !500 reads per processor or 25 Mb buffer for the listener%param
      BufMaxNsends = 1000
      BufMaxSize = 25e6/a%MPIsize/rp
      BufSafetySend = 5e2
      call a%SetBufferSizes(BufMaxNsends,BufMaxSize,BufSafetySend)    
      
      a%AllocCount = 0_ip
   end subroutine
      
   subroutine Finalize(a)
      implicit none
      class(OneToAllReadExnor) :: a
      
      !Memory Counter
      call a%Memor%allocObj(0_ip,'ExnorSystems',outstr,a%AllocCount)
      
   end subroutine
   
   subroutine RootGetDataAndCompletionCheck(a)
      implicit none
      class(OneToAllReadExnor) :: a
      
      !Reads a line from the file
      call a%Listener%listen(outstr)
      
      !Checks wether we are done reading on boundaries boundary conditions
      if(a%Listener%words(1) == 'ENDEX') a%kfl_CompletionCheck = .true.
   end subroutine
   
   subroutine RootAddToBuffer(a)
      implicit none
      class(OneToAllReadExnor) :: a      
      
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
   
   subroutine GetDataFromBufferAndOperate(a)
      implicit none
      class(OneToAllReadExnor) :: a   
      
      integer(ip) :: ipoin
      
 
   
      !Get data from Buffer
      call a%Buffer%GetFromBuffer(a%isend,a%Listener%nnpar,a%Listener%param)
      
      !Do what needs to be done with on nodes information
      !Read point number
      ipoin                   = int(a%Listener%param(1)) 
   
      !To local numbering
      call a%Mesh%Global2Local(ipoin,ipoin)  
      
      a%Mesh%IsExnor(ipoin) = .true.
      allocate(a%Mesh%ExternalNormal(ipoin)%a(a%ndime))
      a%AllocCount = a%AllocCount + 1
      
      a%Mesh%ExternalNormal(ipoin)%a = a%Listener%param(2:a%ndime+1)
      a%AllocCount = a%AllocCount + a%ndime*rp
      
   end subroutine
end module

subroutine reaExnor(a)
   use typre
   use Mod_Mesh
   use Mod_OneToAllReadExnor
   use Mod_iofile
   use MPI
   implicit none
   class(FemMesh) :: a
   character(150) :: fil_dom_fix 
   
   type(OneToAllReadExnor) :: OTAReadExnor
   
   character(150) :: outstr = 'reaExnor'
   integer(ip) :: icount,counter
   
   !MPI
   integer :: ierr
   
   !Allocate the arrays for storing the Exnor systems
   call a%Memor%alloc(a%npoin,a%IsExnor,'isExnor','ReaExnor')
   call a%Memor%alloc(a%npoin,a%ExternalNormal,'ExternalNormal','ReaExnor')
   
   !We will need a buffer structure, etc
   
   !Reach the Exnor Systems Section
   if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
      !Reach the ExnorSystems section.
      call a%Listener%listen(outstr)
      
      !We count how many lines we read so that we can rewind them if there were no EXNOR 
      !(only for old cases, new ones always have EXNOR)
      !(in old cases it only reads 4 lines so it is not computationally costly)
      counter = 0
      do while(a%Listener%words(1)/='EXTER' .and. a%Listener%words(1)/='ENDGE')
         call a%Listener%listen(outstr)
         counter = counter +1
      end do
   endif
   
   if (a%kfl_ReadType == 0) &
      CALL MPI_BCAST(a%Listener%words(1), 5, MPI_CHARACTER, a%MPIroot, a%MPIcomm, ierr)
   
   !If I arrived to the end of the geometry, the rewind just the counted ones
   if (a%Listener%words(1) == 'ENDGE') then
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         do icount = 1,counter
            backspace(a%Listener%nunit)
         enddo
      endif
      
   !Else do what we need to do   
   elseif (a%Listener%words(1)=='EXTER') then 
           
      !ExnorSystems
      call OTAReadExnor%SetType(a%kfl_ReadType)
      call OTAReadExnor%Initialize(a)
      call OTAReadExnor%Loop
      call OTAReadExnor%Finalize
   endif
   
!    if (a%kfl_inter==1) then
!       fil_dom_fix = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.dom.fix' 
!       call iofile(two,a%Listener%nunit,fil_dom_fix,'DOMAIN','old','formatted')
!    endif

end subroutine


Module Mod_OneToAllBuffer
   use typre
   use Mod_Memor
   use MPI
   private
   public OneToAllBuffer
   
   !MPI parameters
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   
   !> Object for sending data from one processor to the others <BR>
   !> Processor "root" Adds data to the buffer <BR>
   !> All the processors receive data from the buffer <BR>
   !> See an example of use in type(OneToAllLoop)
   type :: OneToAllBuffer
      integer(ip) :: kfl_type=-1  !buffer type (0: OneToAll. 1:AllToAll
      
      !MPI
      integer(ip)    :: MPIrank
      integer(ip)    :: MPIroot
      integer(ip)    :: MPIsize
      integer(ip)    :: MPIComm
      type(MemoryMan), pointer :: Memor => NULL()
      
      integer(ip), allocatable :: irequest01(:), irequest02(:), irequest03(:)
      
      !Buffers for sending (CSR format)
      integer(ip), allocatable :: Snsend(:)
      type(i1p), allocatable ::  Sia(:)
      type(r1p), allocatable ::    Sja(:)
      logical , allocatable ::    kfl_Sending(:)
      
      
      !Buffers for receiving (CSR format)
      integer(ip) :: Rnsend
      integer(ip), allocatable :: Ria(:)
      real(rp), allocatable :: Rja(:)
      
      integer(ip) :: MaxNSends,MaxSize, SafetySend
      
      integer(ip) :: kfl_doreceiver
      
      !Procedure pointers, depends on type of buffer
      procedure(InitializeOneToAll), pointer :: Initialize => NULL()
      procedure(AddToBufferOneToAll), pointer :: AddToBuffer => NULL()
      procedure(SendBufferOneToAll), pointer :: SendBuffer => NULL()
      procedure(FinalSendOneToAll), pointer :: FinalSendOnlyRoot => NULL()
      procedure(FinalizeOneToAll), pointer :: Finalize => NULL()
   
contains

      procedure :: SetMPI
      procedure :: SetBufferTypeString
      procedure :: SetBufferTypeInteger
      
      procedure :: WaitBuffer
      procedure :: GetFromBuffer
      
      procedure :: ReceiveBufferOneToAll
      procedure :: ReceiveBufferAllToAll
      
      generic :: ReceiveBuffer => ReceiveBufferOneToAll, ReceiveBufferAllToAll
      generic :: SetBufferType => SetBufferTypeInteger, SetBufferTypeString

   end type

      

contains

!Base subroutines
   subroutine SetBufferTypeInteger(a,TypeInteger)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: TypeInteger
      
      a%kfl_type = TypeInteger
      call subSetBufferType(a)
   end subroutine
   
   subroutine SetBufferTypeString(a,TypeString)
      implicit none
      class(OneToAllBuffer) :: a
      character(8) :: TypeString
      
      if (trim(TypeString) == 'OneToAll') then
         a%kfl_type = 0
      elseif (trim(TypeString) == 'AllToAll') then
         a%kfl_type = 1
      endif
      call subSetBufferType(a)
   end subroutine   
   
   subroutine subSetBufferType(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      if (a%kfl_type == 0) then
         a%Initialize => InitializeOneToAll
         a%AddToBuffer => AddToBufferOneToAll
         a%SendBuffer => SendBufferOneToAll
         a%FinalSendOnlyRoot => FinalSendOneToAll
         a%Finalize => FinalizeOneToAll
         
      elseif (a%kfl_type == 1) then
         a%Initialize => InitializeAllToAll
         a%AddToBuffer => AddToBufferAllToAll
         a%SendBuffer => SendBufferAllToAll
         a%FinalSendOnlyRoot => FinalSendAllToAll
         a%Finalize => FinalizeAllToAll
      endif
      
   end subroutine   
   
   

   subroutine SetMPI(a,MPIcomm,MPIsize,MPIroot,MPIrank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: MPIcomm,MPIrank,MPIsize,MPIroot
      
      a%MPIrank = MPIrank
      a%MPIsize = MPIsize
      a%MPIroot = MPIroot
      a%MPIcomm = MPIcomm
   end subroutine
      
   subroutine AllocateSendStructure(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      integer(ip) :: irank
      
      call a%Memor%alloc(a%MPIsize,a%Snsend,'Snsend','OneToAllBufferInitialize')
      call a%Memor%alloc(a%MPIsize,a%Sia,'Sia','OneToAllBufferInitialize')
      do irank = 0,a%MPIsize-1
         call a%Memor%palloc(a%MaxNSends+1,a%Sia(irank+1)%l,'Sia','OneToAllBufferInitialize')
      enddo
      call a%Memor%alloc(a%MPIsize,a%Sja,'Sja','OneToAllBufferInitialize')
      do irank = 0,a%MPIsize-1
         call a%Memor%palloc(a%MaxSize,a%Sja(irank+1)%a,'Sja','OneToAllBufferInitialize')
      enddo
      call a%Memor%alloc(a%MPIsize,a%kfl_Sending,'kfl_Sending','OneToAllBufferInitialize')
      
      call a%Memor%alloc(a%MPIsize,a%irequest01,'irequest01','OneToAllBufferInitialize')
      call a%Memor%alloc(a%MPIsize,a%irequest02,'irequest02','OneToAllBufferInitialize')
      call a%Memor%alloc(a%MPIsize,a%irequest03,'irequest03','OneToAllBufferInitialize')
      a%irequest01 = MPI_REQUEST_NULL
      a%irequest02 = MPI_REQUEST_NULL
      a%irequest03 = MPI_REQUEST_NULL
      
      !First value for CSR
      do irank = 0,a%MPIsize -1
         a%Sia(irank+1)%l(1) = 1
      enddo
   end subroutine
   
   subroutine AllocateReceiveStructure(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      call a%Memor%alloc(a%MaxNSends+1,a%Ria,'Ria','OneToAllBufferInitialize')
      call a%Memor%alloc(a%MaxSize,a%Rja,'Rja','OneToAllBufferInitialize')
      a%Rnsend = 0
   end subroutine
    
   subroutine DeallocateSendStructure(a)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
      
      call a%Memor%dealloc(a%MPIsize,a%Snsend,'Snsend','OneToAllBufferFinalize')
      
      do irank = 0,a%MPIsize-1
         call a%Memor%pdealloc(size(a%Sia(irank+1)%l),a%Sia(irank+1)%l,'Sia','OneToAllBufferInitialize')
      enddo
      call a%Memor%dealloc(a%MPIsize,a%Sia,'Sia','OneToAllBufferInitialize')
      do irank = 0,a%MPIsize-1
         call a%Memor%pdealloc(size(a%Sja(irank+1)%a),a%Sja(irank+1)%a,'Sja','OneToAllBufferInitialize')
      enddo
      call a%Memor%dealloc(a%MPIsize,a%Sja,'Sja','OneToAllBufferInitialize')
      
      call a%Memor%dealloc(a%MPIsize,a%kfl_Sending,'kfl_Sending','OneToAllBufferFinalize')
      call a%Memor%dealloc(a%MPIsize,a%irequest01,'irequest01','OneToAllBufferFinalize')
      call a%Memor%dealloc(a%MPIsize,a%irequest02,'irequest02','OneToAllBufferFinalize')
      call a%Memor%dealloc(a%MPIsize,a%irequest03,'irequest03','OneToAllBufferFinalize')
   end subroutine   
   
   subroutine DeallocateReceiveStructure(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      call a%Memor%dealloc(size(a%Ria),a%Ria,'Ria','OneToAllBufferFinalize')
      call a%Memor%dealloc(size(a%Rja),a%Rja,'Rja','OneToAllBufferFinalize')
   end subroutine
    
   subroutine AddToBuffer(a,irank,nsize,array)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank,nsize
      real(rp) :: array(nsize)
      
      integer(ip) :: ispos, endpos
      
      !If we just sent the buffer, wait for sending completion
      if (a%kfl_Sending(irank+1) .eqv. .true.) then
         call WaitBuffer(a,irank)
      endif   
      
      !Add everything to the Buffer
      a%Snsend(irank+1) = a%Snsend(irank+1)+1
      a%Sia(irank+1)%l(a%Snsend(irank+1)+1) = a%Sia(irank+1)%l(a%Snsend(irank+1)) + nsize
      ispos = a%Sia(irank+1)%l(a%Snsend(irank+1))
      endpos = ispos + nsize -1
      a%Sja(irank+1)%a(ispos:endpos) = array
   end subroutine
   
   subroutine ReallocateBuffer(a,irank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
      
      integer(ip) :: MaxNSends, MaxSize
      
      if (  (a%Snsend(irank+1) == size(a%Sia(irank+1)%l)-1)) then
         MaxNSends = size(a%Sia(irank+1)%l)*2
         call a%Memor%prealloc(MaxNSends+1,a%Sia(irank+1)%l,'Sia','OneToAllBufferInitialize')
      endif
      
      if ( (a%Sia(irank+1)%l(a%Snsend(irank+1)+1) > size(a%Sja(irank+1)%a) - a%SafetySend)) then
         MaxSize = size(a%Sja(irank+1)%a)*2
         call a%Memor%prealloc(MaxSize,a%Sja(irank+1)%a,'Sja','OneToAllBufferInitialize')
      endif   
      
   end subroutine
   
   subroutine WaitBuffer(a,irank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      !Wait for the sending of the buffer to complete
      call MPI_WAIT(a%irequest01(irank+1), status, ierr)
      call MPI_WAIT(a%irequest02(irank+1), status, ierr)
      call MPI_WAIT(a%irequest03(irank+1), status, ierr)
      !Reset the buffer to empty
      a%Snsend(irank+1) = 0
      
      !We have completed send
      a%kfl_Sending(irank+1) = .false.
   end subroutine
   
   subroutine SendBuffer(a,irank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
      
      integer(ip) :: endpos
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      if (a%Snsend(irank+1) > 0) then
         endpos = a%Sia(irank+1)%l(a%Snsend(irank+1)+1)-1
         call MPI_ISEND(a%Snsend(irank+1),1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,a%irequest01(irank+1), ierr)
         call MPI_ISEND(a%Sia(irank+1)%l(1),a%Snsend(irank+1)+1,MPI_INTEGER4, irank, mtag2, a%MPIcomm,a%irequest02(irank+1), ierr)
         call MPI_ISEND(a%Sja(irank+1)%a(1),endpos,MPI_REAL8, irank, mtag3, a%MPIcomm,a%irequest03(irank+1), ierr)
         a%kfl_Sending(irank+1) = .true.
      endif
   end subroutine
   
   subroutine GetFromBuffer(a,isend,nsize,array)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: isend,nsize
      real(rp)    :: array(*)
      
      integer(ip) :: ispos,endpos
      
      nsize = a%Ria(isend+1) - a%Ria(isend)
      ispos = a%Ria(isend)
      endpos = a%Ria(isend+1)-1
      array(1:nsize) = a%Rja(ispos:endpos)
   end subroutine
   
   
   

!One to all subroutines
   
   subroutine InitializeOneToAll(a,MaxNSends,MaxSize,SafetySend,Memor)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: MaxNSends, MaxSize, SafetySend
      type(MemoryMan), target :: Memor
      
      a%MaxNSends = MaxNSends
      a%MaxSize   = MaxSize
      a%SafetySend= SafetySend
      a%Memor => Memor
      
      !Default is do receive for everyone except for root
      a%kfl_doreceiver = 1
      
      !Memory allocation
      if (a%MPIrank == a%MPIroot) then
         call AllocateSendStructure(a)
         !Default is do not receive (for root)
         a%kfl_doreceiver = 0
      endif
      call AllocateReceiveStructure(a)
   end subroutine 
    
   subroutine AddToBufferOneToAll(a, irank, nsize, array)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank, nsize
      real(rp) :: array(nsize)
      
      integer(ip) :: ispos, endpos
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      if (a%MPIrank /= a%MPIroot) then
         call runend('OneToAllBuffer: only root can add to buffer')
      endif
      
      !After checks, add the data to the buffer
      call AddToBuffer(a,irank,nsize,array)
      
      !If the buffer is full enough, send
      if (a%Snsend(irank+1) == a%MaxNSends .or. a%Sia(irank+1)%l(a%Snsend(irank+1)+1) > a%MaxSize - a%SafetySend) then
         call SendBufferOneToAll(a,irank)
      endif
   end subroutine
      
   subroutine SendBufferOneToAll(a,irank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
   
      if (a%MPIrank /= a%MPIroot) then
         call runend('OneToAllBuffer: only root can send buffer')
      endif
      
      call SendBuffer(a,irank)
      
      !If I sent to root, I need to receive
      if (irank == a%MPIroot) then
         if (a%Snsend(irank+1) > 0) a%kfl_doreceiver = 1
      endif
      
   end subroutine
      
   subroutine ReceiveBufferOneToAll(a,nsend)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: nsend
      
      integer(ip) :: endpos
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      nsend = 0_ip
      if (a%kfl_doreceiver == 1) then
         call MPI_RECV(nsend,1_ip, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,status, ierr) 
         
         !We are done with receiving
         if (nsend == -1) then
         
         !Receive the buffer
         elseif (nsend /= 0) then
            a%Rnsend = nsend
            
            call MPI_RECV(a%Ria,a%Rnsend+1, MPI_INTEGER4,a%MPIroot, mtag2, a%MPIcomm,status, ierr) 
            endpos = a%Ria(a%Rnsend+1)-1
            call MPI_RECV(a%Rja,endpos, MPI_REAL8, a%MPIroot, mtag3, a%MPIcomm,status, ierr) 
         endif
      endif
      
      !Default for root is do not receive
      if (a%MPIrank == a%MPIroot) then
         a%kfl_doreceiver = 0_ip
      endif
         
   end subroutine
   
   subroutine FinalSendOneToAll(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      integer(ip) :: irank
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      if (a%MPIrank /= a%MPIroot) then
         call runend('OneToAllBuffer: only root can FinalSendOnlyRoot')
      endif
      
      do irank = 0,a%MPIsize-1
         !Wait for everyone to complete send
         call WaitBuffer(a,irank)
         !Tell everybody we are done
         
         call MPI_ISEND(-1_ip,1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,a%irequest01(irank+1), ierr)
      enddo
      
      a%kfl_doreceiver = 1
   end subroutine
      
   subroutine FinalizeOneToAll(a)
      implicit none
      class(OneToAllBuffer) :: a
      type(MemoryMan) :: Memor
      
      integer(ip) :: irank
      
      !Memory deallocation
      if (a%MPIrank == a%MPIroot) then
         !Wait for buffers before deallocating
         do irank = 0,a%MPIsize-1
            call WaitBuffer(a,irank)
         enddo
      
         call DeallocateSendStructure(a)
      endif
      call DeallocateReceiveStructure(a)
   end subroutine   
   
!All to all subroutines
   subroutine InitializeAllToAll(a,MaxNSends,MaxSize,SafetySend,Memor)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: MaxNSends, MaxSize, SafetySend
      type(MemoryMan), target :: Memor
      
      a%MaxNSends = MaxNSends
      a%MaxSize   = MaxSize
      a%SafetySend= SafetySend
      a%Memor => Memor

      call AllocateReceiveStructure(a)
      call AllocateSendStructure(a)
   end subroutine
   
   subroutine AddToBufferAllToAll(a, irank, nsize, array)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank, nsize
      real(rp) :: array(nsize)
      
      integer(ip) :: ispos, endpos
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      
      !After checks, add the data to the buffer
      call AddToBuffer(a,irank,nsize,array)
      
      !If the buffer is full, reallocate it
      if (  (a%Snsend(irank+1) == size(a%Sia(irank+1)%l)-1) .or. (a%Sia(irank+1)%l(a%Snsend(irank+1)+1) > size(a%Sja(irank+1)%a) - a%SafetySend)  ) then
         call ReallocateBuffer(a,irank)
      endif
   end subroutine
   
   subroutine SendBufferAllToAll(a,irank)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: irank
      
      integer(ip) :: endpos
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      if (a%Snsend(irank+1) > 0) then
         call SendBuffer(a,irank)
      !just tell them we are done
      else
         call WaitBuffer(a,irank)
         !Tell everybody we are done
         call MPI_ISEND(-1_ip,1,MPI_INTEGER4, irank, mtag1, a%MPIcomm,a%irequest01(irank+1), ierr)
      endif
   end subroutine
   
   subroutine ReceiveBufferAllToAll(a,irank,nsend)
      implicit none
      class(OneToAllBuffer) :: a
      integer(ip) :: nsend,irank
      
      integer(ip) :: endpos
      
      !MPI
      integer            :: ierr
      integer            :: status(MPI_STATUS_SIZE)
      
      nsend = 0_ip
      
      call MPI_RECV(nsend,1_ip, MPI_INTEGER4, irank, mtag1, a%MPIcomm,status, ierr) 
      
      !We are done with receiving
      if (nsend == -1) then
      
      !Receive the buffer
      elseif (nsend /= 0) then
         a%Rnsend = nsend
         !Reallocate if necessary
         if (a%Rnsend+1 > size(a%Ria)) then
            call a%Memor%realloc(a%Rnsend+1,a%Ria,'Ria','ReceiveBufferAllToAll')
         endif
         call MPI_RECV(a%Ria,a%Rnsend+1, MPI_INTEGER4,irank, mtag2, a%MPIcomm,status, ierr) 
         endpos = a%Ria(a%Rnsend+1)-1
         !Reallocate if necessary
         if (endpos > size(a%Rja) ) then
            call a%Memor%realloc(endpos,a%Rja,'Rja','ReceiveBufferAllToAll')
         endif
         call MPI_RECV(a%Rja,endpos, MPI_REAL8, irank, mtag3, a%MPIcomm,status, ierr) 
      endif
   end subroutine
   
   subroutine FinalSendAllToAll(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      call runend('Do not call final send if all to all (not needed, a single send to each processor is used')
      
   end subroutine
   
   subroutine FinalizeAllToAll(a)
      implicit none
      class(OneToAllBuffer) :: a
      
      integer(ip) :: irank
      
      !Wait for buffers before deallocating
      do irank = 0,a%MPIsize-1
         call WaitBuffer(a,irank)
      enddo
      
      call DeallocateSendStructure(a)
         
      call DeallocateReceiveStructure(a)
   end subroutine
   
end module

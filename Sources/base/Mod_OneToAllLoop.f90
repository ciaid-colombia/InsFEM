module Mod_OneToAllLoop
   use typre
   use Mod_Memor
   use Mod_OneToAllBuffer
   private
   public OneToAllLoop
   
   type, abstract :: OneToAllLoop
      
      integer(ip) :: kfl_type = -1
      
      !Memor
      type(MemoryMan), pointer :: Memor => NULL()
      
      !Buffer
      type(OneToAllBuffer) :: Buffer
      integer(ip) :: BufMaxNsends, BufMaxSize, BufSafetySend
      
      !MPI
      integer(ip)    :: MPIrank
      integer(ip)    :: MPIroot
      integer(ip)    :: MPIsize
      integer(ip)    :: MPIcomm
      
      logical :: kfl_CompletionCheck
      
      !Counter for looping trough buffer when receiving
      integer(ip) :: isend
       
      procedure(LoopOneToAll), pointer :: Loop => NULL()
   
contains
      
      procedure :: SetMPI
      procedure :: SetLoopTypeString
      procedure :: SetLoopTypeInteger
      procedure :: SetMemor
      procedure :: SetBufferSizes
      procedure(SpecificRootGetDataAndCompletionCheck), deferred  :: SpecificRootGetDataAndCompletionCheck
      procedure(SpecificRootAddToBuffer), deferred                :: SpecificRootAddToBuffer
      procedure(SpecificGetDataFromBufferAndOperate), deferred    :: SpecificGetDataFromBufferAndOperate
      
      generic :: SetType => SetLoopTypeInteger, SetLoopTypeString
      
   end type
   
   abstract interface
      subroutine SpecificRootGetDataAndCompletionCheck(a)
         import OneToAllLoop
         implicit none
         class(OneToAllLoop) :: a
      end subroutine
      
      subroutine SpecificRootAddToBuffer(a)
         import OneToAllLoop
         implicit none
         class(OneToAllLoop) :: a
      end subroutine
      
      subroutine SpecificGetDataFromBufferAndOperate(a)
         import OneToAllLoop
         implicit none
         class(OneToAllLoop) :: a
      end subroutine
   end interface

contains

   subroutine SetMPI(a,MPIcomm,MPIsize,MPIroot,MPIrank)
      implicit none
      class(OneToAllLoop) :: a
      integer(ip) :: MPIrank,MPIsize,MPIroot,MPIcomm
      
      a%MPIcomm = MPIcomm
      a%MPIrank = MPIrank
      a%MPIsize = MPIsize
      a%MPIroot = MPIroot
   end subroutine

   subroutine SetMemor(a,Memor)
      implicit none
      class(OneToAllLoop) :: a
      type(MemoryMan), target :: Memor
      
      a%Memor => Memor
   end subroutine
   
   subroutine SetLoopTypeString(a,TypeString) 
      implicit none
      class(OneToAllLoop) :: a
      character(8) :: TypeString
      
      if (trim(TypeString) == 'OneToAll') then
         a%kfl_type = 0
         
      elseif (trim(TypeString) == 'AllToAll') then
         a%kfl_type = 1   
      endif
      call subSetLoopType(a)
   end subroutine
   
   subroutine SetLoopTypeInteger(a,TypeInteger) 
      implicit none
      class(OneToAllLoop) :: a
      integer(ip) :: TypeInteger
      
      a%kfl_type = TypeInteger
      call subSetLoopType(a)
   end subroutine
   
   subroutine subSetLoopType(a)
      implicit none
      class(OneToAllLoop) :: a
      
      if (a%kfl_type == 0) then
         a%Loop => LoopOneToAll
         
      elseif (a%kfl_type == 1) then
         a%Loop => LoopAllToAll
         
      endif
   end subroutine
   
   subroutine SetBufferSizes(a,MaxNSends,MaxSize,SafetySend)
      implicit none
      class(OneToAllLoop) :: a
      integer(ip) :: MaxNSends, MaxSize, SafetySend
      
      a%BufMaxNsends = MaxNSends
      a%BufMaxSize = MaxSize
      a%BufSafetySend = SafetySend
   end subroutine
   
   subroutine LoopInitializations(a)
      implicit none
      class(OneToAllLoop) :: a
      
      !Initializations
      a%kfl_CompletionCheck = .false.
      !Initialize Buffer
      call a%Buffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Buffer%SetBufferType(a%kfl_type)
      call a%Buffer%Initialize(a%BufMaxNsends,a%BufMaxSize,a%BufSafetySend,a%Memor)
   end subroutine
      
   
   subroutine LoopOneToAll(a)
      implicit none
      class(OneToAllLoop) :: a
      
      integer(ip) :: kfl_Continue, kfl_DoRoot
      integer(ip) :: isend,nsend,irank
      
      if (a%kfl_type == -1) then
         call runend('OneToAllLoop: loop type not set')
      endif
      
      !Initializations
      call LoopInitializations(a)
      
      
      kfl_Continue = 1
      kfl_DoRoot = 1
      do while (kfl_Continue == 1)
         !To do by root
         if (a%MPIrank == a%MPIroot) then
      
            !We are done with boundaries, second pass
            if (kfl_DoRoot == 0) then
               call a%Buffer%FinalSendOnlyRoot
            
            !Listen
            elseif (kfl_DoRoot == 1) then
               
               call a%SpecificRootGetDataAndCompletionCheck
               !call a%Listener%listen(outstr)
               
               !We are done with boundaries, first pass
               !if(a%Listener%words(1) == 'ENDON') then
               if (a%kfl_CompletionCheck .eqv. .true.) then
                  
                  !If something still needs to be sent, we send it
                  do irank = 0,a%MPIsize-1
                     call a%Buffer%SendBuffer(irank)
                  enddo
                  !next one is the last one, and I do not want to read
                  kfl_DoRoot = 0
               
               !OnBoundaries information
               else
                  call a%SpecificRootAddToBuffer

               endif
            endif
         endif    
         
         !----------------
         !To do by everybody
         !Receive Buffer
         call a%Buffer%ReceiveBuffer(nsend)
         if (nsend == -1) then
            kfl_Continue = 0 !We are done with boundaries
         elseif (nsend /= 0) then
            !Receive the buffer
            !Now get data from buffer and do what needs to be done
            do isend = 1,nsend
               a%isend = isend
               call a%SpecificGetDataFromBufferAndOperate

            enddo
         endif
      enddo
      
      !Finalize Buffer
      call a%Buffer%Finalize
   end subroutine
   
   subroutine LoopAllToAll(a)
      use MPI
      implicit none
      class(OneToAllLoop) :: a
      
      integer(ip) :: kfl_Continue, kfl_DoRoot
      integer(ip) :: isend,nsend,irank
      
      integer numtasks, rank, ierr, rc, len, i
      character*(MPI_MAX_PROCESSOR_NAME) name
      
      if (a%kfl_type == -1) then
         call runend('OneToAllLoop: loop type not set')
      endif
      
      call MPI_GET_PROCESSOR_NAME(name, len, ierr)
      
      !Initializations
      call LoopInitializations(a)
      
      a%kfl_CompletionCheck = .false.
      call a%SpecificRootGetDataAndCompletionCheck
      do while (a%kfl_CompletionCheck .eqv. .false.) 
         call a%SpecificRootAddToBuffer
         call a%SpecificRootGetDataAndCompletionCheck
      enddo
      
      !Send Buffer
       do irank = 0,a%MPIsize-1
          call a%Buffer%SendBuffer(irank)
       enddo
      
      !Receive all buffers
      do irank = 0,a%MPIsize-1
         call a%Buffer%ReceiveBuffer(irank,nsend)
         if (nsend == -1) then
            !We are done with boundaries
         elseif (nsend /= 0) then
            !Receive the buffer
            !Now get data from buffer and do what needs to be done
            do isend = 1,nsend
               a%isend = isend
               call a%SpecificGetDataFromBufferAndOperate
            enddo
         endif
      enddo
      !Finalize Buffer
      call a%Buffer%Finalize
      
   end subroutine   

end module

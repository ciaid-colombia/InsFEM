module Mod_ReaLnods
   use typre
   implicit none

contains   
   subroutine reaLnodsSerial(a)
      use MPI
      use typre
      use Mod_Int2Str
      use Mod_Mesh
      use Mod_NaivePartition
      implicit none
      class(FemMesh) :: a
      
      integer(ip) :: kfl_ContinueReading,kfl_newformat,kfl_doreceiver
      
      integer(ip) :: nelem, maxnelem
      integer(ip) :: lnods(a%mnode)
      integer(ip) :: ielem,gielem,inode,ipoin,ispos,nnode
      integer(ip), allocatable :: ielemFlag(:)
      
      integer(ip) :: iaux,newmaxnelem,inipoin,endpoin
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer            :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize)
      integer            :: status(MPI_STATUS_SIZE)
      integer(ip) :: irank
      
      integer(ip) :: myIOStat1
      
      integer(ip), allocatable :: SBuffer_Nelem(:)
      integer(ip), allocatable :: SBuffer_Pnods(:,:,:)  !1: pnods 2: element number
      integer(ip), allocatable :: SBuffer_Lnods(:,:)
      
      integer(ip)              :: RBuffer_Nelem
      integer(ip), allocatable :: RBuffer_Pnods(:,:)    !1: pnods 2: element number
      integer(ip), allocatable :: RBuffer_Lnods(:)
      
      integer(ip) :: SBBytesNelem,SBBytesProc,SBNelemProc
      
      integer(ip) :: SBielem, SBispos
      integer(ip) :: lnods_size
      
      integer(ip) :: NelemCounter(a%MPIsize)
         
      type(NaivePartition) :: NaivePartitioner
      
      integer(ip) :: aux
      integer(ip), allocatable :: gLocalNaiveToInitialElementNumbering(:)
      
      
      !Initialize NaivePartitioner
      call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
      
      !Auxiliar value to know which is the processor for each node
      call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
      a%npoinLocal = endpoin - inipoin + 1
      
      !Approximate number of elements per node
      maxnelem = a%gnelem/a%MPIsize*1.5_rp
      
      call a%Memor%alloc(maxnelem+1,a%gpnods,'gpnods','reageo')
      call a%Memor%alloc(maxnelem,gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
      call a%Memor%alloc(maxnelem*a%mnode,a%glnods,'glnods','reageo')
      
      
      !Allocate a 5 MBytes Buffer for avoiding too many send and receive operations
      SBBytesNelem = ip*(1+a%mnode)
      SBBytesProc  = SBBytesNelem*a%MPIsize
      SBNelemProc = min(50e6/SBBytesProc,1.0_rp*a%gnelem)
      
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%alloc(a%MPIsize,SBuffer_Nelem,'SBuffer_Nelem','reageo')
         call a%Memor%alloc(2,SBNelemProc+1,a%MPIsize,SBuffer_Pnods,'SBuffer_Lnods','reageo')
         call a%Memor%alloc(SBNelemProc*a%mnode,a%MPIsize,SBuffer_Lnods,'SBuffer_Lnods','reageo')
         
         !We initialize SBuffer_Pnods
         SBuffer_Pnods(1,1,:) = 1
         
      endif
      
      call a%Memor%alloc(2,SBNelemProc+1,RBuffer_Pnods,'RBuffer_Pnods','reageo')
      call a%Memor%alloc(SBNelemProc*a%mnode,RBuffer_Lnods,'RBuffer_Lnods','reageo')
      
      !Initialize the element counter for each processor
      NelemCounter = 0
      
      !Root Reach the elements section
      if (a%MPIrank == a%MPIroot) then
         call ReachElementsSection(a,kfl_newformat)
      endif
      
      !Allocate ielemFlag
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%alloc(a%MPIsize,ielemFlag,'ielemFlag','reabcs')
         ielemFlag = -1
         
      endif
      
      a%gpnods(1) = 1
      
      !Still not done flag
      kfl_ContinueReading = 0
      
      !Number of read elements
      if (a%MPIrank == a%MPIroot) gielem = 0
      
      !Number of elements per processor
      a%nelem = 0
      
      irequest1 = MPI_REQUEST_NULL
      irequest2 = MPI_REQUEST_NULL
      irequest3 = MPI_REQUEST_NULL
      
      do while (kfl_ContinueReading /= -1)
      
         !Default is do what receivers need to do
         kfl_doreceiver = 1
                  
         !To be done by root
         if (a%MPIrank == a%MPIroot) then
            
            !For MPIrank, default is do not do receive
            kfl_doreceiver = 0
         
            gielem = gielem +1
            
            !If we are done, first pass we send all that still is in the buffer
            if (gielem == a%gnelem +1 ) then
               !If something is left to be sent, we send it
               do irank = 0,a%MPIsize-1
                  if (SBuffer_Nelem(irank+1) > 0 .and. SBuffer_Nelem(irank+1) /= SBNelemProc) then
                     call realnods_sendbuffers(a%MPIcomm,irank,SBuffer_Nelem(irank+1),SBuffer_Pnods(1,1,irank+1),SBuffer_Lnods(1,irank+1),irequest1(irank+1),irequest2(irank+1),irequest3(irank+1))
                     !If I send it to root, then root must receive
                     if (irank == a%MPIroot) kfl_doreceiver = 1
                  endif
               enddo
         
            !Second pass after we are done, we tell everybody that we are done
            elseif (gielem == a%gnelem+2) then
               
               !We tell everybody we are done
               do irank = 0,a%MPIsize-1
                  call MPI_WAIT(irequest1(irank+1), status, ierr)
                  call MPI_ISEND(-1_ip,1, MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest1(irank+1), ierr) 
               enddo
               !If we are done, root must do receive
               kfl_doreceiver = 1
               
            
            !Continue sending
            else   
               call ReadElement(a,kfl_newformat,ielem,nnode,lnods)
               
               !Send to the corresponding processor depending on the node number
               do inode = 1,nnode
                  ipoin = lnods(inode)
                  call NaivePartitioner%GetPointRank(ipoin,irank)
                  
                  if (ielemFlag(irank+1) /= gielem) then
                     ielemFlag(irank+1) = gielem
                     
                     !Update the counter
                     NelemCounter(irank+1) = NelemCounter(irank+1)+1
                     
                     !If we just sent the buffer, we wait for the send to be complete
                     if (SBuffer_Nelem(irank+1) == SBNelemProc) then
                        call MPI_WAIT(irequest1(irank+1), status, ierr)
                        call MPI_WAIT(irequest2(irank+1), status, ierr)
                        call MPI_WAIT(irequest3(irank+1), status, ierr)
                        
                        !We start the buffer again
                        SBuffer_Nelem(irank+1) = 0
                     endif
                     
                     !First we place the element in the buffer
                     SBuffer_Nelem(irank+1) = SBuffer_Nelem(irank+1)+1
                     
                     SBielem = SBuffer_Nelem(irank+1)
                     SBuffer_Pnods(2,SBielem,irank+1) = ielem
                     SBuffer_Pnods(1,SBielem+1,irank+1) =  SBuffer_Pnods(1,SBielem,irank+1) + nnode
                     SBispos = SBuffer_Pnods(1,SBielem,irank+1)-1
                     SBuffer_Lnods(SBispos+1:SBispos+nnode,irank+1) = lnods(1:nnode)
                     
                     !If the Buffer is full, we send it
                     if (SBuffer_Nelem(irank+1) == SBNelemProc) then
                        call realnods_sendbuffers(a%MPIcomm,irank,SBuffer_Nelem(irank+1),SBuffer_Pnods(1,1,irank+1),SBuffer_Lnods(1,irank+1),irequest1(irank+1),irequest2(irank+1),irequest3(irank+1))
                        !If I sent it to root, then root must receive
                        if (irank == a%MPIroot) kfl_doreceiver = 1
                     endif
                  endif
               enddo
            endif
               
         endif
         
         if (kfl_doreceiver == 1) then
            !To be done by everybody
            call MPI_RECV(kfl_ContinueReading,1_ip, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,status, ierr) 
            if (kfl_continueReading > 0_ip) then
               
               RBuffer_Nelem = kfl_continueReading
               
               !Reallocate arrays if the actual number of elements is too large
               if (RBuffer_Nelem + a%nelem > maxnelem) then
                  newmaxnelem = max(maxnelem*1.5_rp,1.0_rp*(RBuffer_Nelem + a%nelem))
                  
                  call a%Memor%realloc(newmaxnelem,gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
                  call a%Memor%realloc(newmaxnelem*a%mnode,a%glnods,'glnods','reageo')
                  call a%Memor%realloc(newmaxnelem+1,a%gpnods,'gpnods','reageo')
                  maxnelem = newmaxnelem
               endif
                        
               call MPI_RECV(RBuffer_Pnods,2*(RBuffer_Nelem+1), MPI_INTEGER4, a%MPIroot, mtag2, a%MPIcomm,status, ierr)
               
               !Update Pnods
               a%gpnods(a%nelem+2:a%nelem+RBuffer_Nelem+1) = a%gpnods(a%nelem+1)+RBuffer_Pnods(1,2:RBuffer_Nelem+1)-1
               
               !Update Lnods
               lnods_size = RBuffer_Pnods(1,RBuffer_Nelem+1)-1
               gLocalNaiveToInitialElementNumbering(a%nelem+1:a%nelem+RBuffer_Nelem) = RBuffer_Pnods(2,1:RBuffer_Nelem)
               ispos = a%gpnods(a%nelem+1)
               
               call MPI_RECV(a%glnods(ispos),lnods_size, MPI_INTEGER4, a%MPIroot, mtag3, a%MPIcomm,status, ierr)
               
               !New number of elements
               a%nelem = a%nelem + RBuffer_Nelem
            endif  
         endif 
         
      enddo
      
      !Wait for all communications to be done
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest1(irank+1), status, ierr)
            call MPI_WAIT(irequest2(irank+1), status, ierr)
            call MPI_WAIT(irequest3(irank+1), status, ierr)
         enddo
      endif
      
      !Send the elem0 value to each processor
      if (a%MPIrank == a%MPIroot) then
         aux = 0
         do irank = 0,a%MPIsize-1
            call MPI_ISEND(aux,1, MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest1(irank+1), ierr) 
            aux = aux+NelemCounter(irank+1)
         enddo      
      endif
      call MPI_RECV(a%elem0,1, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,status, ierr)
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest1(irank+1), status, ierr)
         enddo
      endif
      
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%dealloc(a%MPIsize,ielemFlag,'ielemFlag','reabcs')
      endif
      
      call CreateInitialElementNumberingOrdering(a,a%nelem,gLocalNaiveToInitialElementNumbering)
      
      call a%Memor%dealloc(size(gLocalNaiveToInitialElementNumbering),gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
      call a%Memor%realloc(a%gpnods(a%nelem+1)-1,a%glnods,'glnods','reageo')
      call a%Memor%realloc(a%nelem+1,a%gpnods,'gpnods','reageo')
      
      !Deallocate the 5 MBytes Buffer for avoiding too many send and receive operations
      if (a%MPIrank == a%MPIroot) then
         call a%Memor%dealloc(a%MPIsize,SBuffer_Nelem,'SBuffer_Nelem','reageo')
         call a%Memor%dealloc(2,SBNelemProc+1,a%MPIsize,SBuffer_Pnods,'SBuffer_Lnods','reageo')
         call a%Memor%dealloc(SBNelemProc*a%mnode,a%MPIsize,SBuffer_Lnods,'SBuffer_Lnods','reageo')
      endif
      call a%Memor%dealloc(2,SBNelemProc+1,RBuffer_Pnods,'RBuffer_Pnods','reageo')
      call a%Memor%dealloc(SBNelemProc*a%mnode,RBuffer_Lnods,'RBuffer_Lnods','reageo')

      
   end subroutine   

   subroutine reaLnodsPartitioned(a)
      use MPI
      use typre
      use Mod_Int2Str
      use Mod_Mesh
      use Mod_NaivePartition
      implicit none
      class(FemMesh) :: a
      
      integer(ip) :: inipoin,endpoin
      integer(ip) :: kfl_newformat
      
      type(NaivePartition) :: NaivePartitioner
      integer(ip) :: ispos,nnode,lnods(a%mnode)
      integer(ip) :: ielem
      integer(ip) :: maxnelem,newmaxnelem
      
      integer(ip), allocatable :: gLocalNaiveToInitialElementNumbering(:)
      
      !Initialize NaivePartitioner
      call NaivePartitioner%Initialize(a%gnpoin,a%MPIsize)
      
      !Auxiliar value to know which is the processor for each node
      call NaivePartitioner%GetIniEndPoint(a%MPIrank,inipoin,endpoin)
      a%npoinLocal = endpoin - inipoin + 1
      
      !Approximate number of elements per node
      maxnelem = a%gnelem/a%MPIsize*1.5_rp
      
      call a%Memor%alloc(maxnelem+1,a%gpnods,'gpnods','reageo')
      call a%Memor%alloc(maxnelem*a%mnode,a%glnods,'glnods','reageo')
      call a%Memor%alloc(maxnelem,gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
      
      call ReachElementsSection(a,kfl_newformat)
      
      ielem = 0
      a%nelem = 0
      a%gpnods(1) = 1
      do while (ielem /= -1)
         call ReadElement(a,kfl_newformat,ielem,nnode,lnods)
         
         if (ielem /= -1) then
            !reallocate if necessary
            if (a%nelem+1 > maxnelem) then
               newmaxnelem = maxnelem*1.5_rp
               
               call a%Memor%realloc(newmaxnelem,gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
               call a%Memor%realloc(newmaxnelem*a%mnode,a%glnods,'glnods','reageo')
               call a%Memor%realloc(newmaxnelem+1,a%gpnods,'gpnods','reageo')
                        
               maxnelem = newmaxnelem
            endif
            
            a%nelem = a%nelem+1
            a%gpnods(a%nelem+1) = a%gpnods(a%nelem)+nnode
            ispos = a%gpnods(a%nelem)-1
            a%glnods(ispos+1:ispos+nnode) = lnods(1:nnode)
            gLocalNaiveToInitialElementNumbering(a%nelem) = ielem
         endif
      enddo
      
      call CreateInitialElementNumberingOrdering(a,a%nelem,gLocalNaiveToInitialElementNumbering)
      
      call a%Memor%dealloc(size(gLocalNaiveToInitialElementNumbering),gLocalNaiveToInitialElementNumbering,'gLocalNaiveToInitialElementNumbering','reageo')
      call a%Memor%realloc(a%gpnods(a%nelem+1)-1,a%glnods,'glnods','reageo')
      call a%Memor%realloc(a%nelem+1,a%gpnods,'gpnods','reageo')
   end subroutine
   
   subroutine CreateInitialElementNumberingOrdering(a,nelemNaive,ElementNaive2Initial)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: nelemNaive,ElementNaive2Initial(nelemNaive)
      
      integer(ip) :: ProcList(nelemNaive)
      
      ProcList = 0_ip
      !Create and Initialize the LocalOrdering
      call a%ParallelLibrary%CreateOrdering(a%ElementNaive2Initial,a%Memor)  !FACTORY
      call a%ElementNaive2Initial%Init(a%MPIcomm,nelemNaive,ElementNaive2Initial,ProcList,a%Memor)
   
   
   end subroutine

   subroutine ReachElementsSection(a,kfl_newformat)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_newformat
      
      !Read vectors.
      call a%Listener%rewind()
      call a%Listener%listen('REAGEO')
      do while(a%Listener%words(1)/='GEOME')
         call a%Listener%listen('REAGEO')
      end do
      
      !Reach elements section
      do while(a%Listener%words(1)/='ELEME')
         call a%Listener%listen('reageo')
      enddo
      
      !New format reads faster
      if (a%Listener%words(2)=='NEWFO') then
         kfl_newformat = 1
         !ielem nnode lnods
      else
         kfl_newformat = 0
         !ielem lnods
      endif
   end subroutine
   
   subroutine ReadElement(a,kfl_newformat,ielem,nnode,lnods)
      use typre
      use Mod_Int2Str
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: kfl_newformat,ielem, nnode, lnods(*)

      
!       integer, parameter :: line_buf_len= 1024*4
!       character(LEN=line_buf_len) :: InS
!       integer(ip) :: size
      
      integer(ip) :: MyIostat1
   
      if (kfl_newformat == 0) then
         call a%Listener%listen('reageo')
         ielem=int(a%Listener%param(1))
         nnode = a%Listener%nnpar-1
         lnods(1:nnode) = int(a%Listener%param(2:nnode+1))
      else
         !Read everything at once
         read(a%Listener%nunit,*,IOStat=myIoStat1) ielem,nnode, lnods(1:nnode)
         if (myIOStat1/=0) call runend('Reageo: Error reading element '//adjustl(trim(int2str(a%nelem))))
      endif
      
      if(ielem<-1.or.ielem>a%gnelem) &
            call runend('REAGEO: WRONG NUMBER OF ELEMENTS') 
      if (nnode /= 0) then 
         if((nnode > a%mnode) .or. (a%ltype(nnode) == 0))&
            call runend('REAGEO: WRONG TYPE OF ELEMENT '//int2str(ielem))
      endif
   end subroutine
      
   subroutine realnods_sendbuffers(MPIcomm,irank,SBuffer_Nelem,SBuffer_Pnods,SBuffer_Lnods,irequest1,irequest2,irequest3)
      use MPI
      use typre
      implicit none
      
      integer(ip), intent(in) :: SBuffer_Nelem,irank,MPIcomm
      integer(ip), intent(in) :: SBuffer_Pnods(2,*),SBuffer_Lnods(*)
      integer(ip), intent(in) :: irequest1,irequest2,irequest3
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
      integer     :: ierr   
      
      integer(ip) :: lnods_size
      
      call MPI_ISEND(SBuffer_Nelem,1, MPI_INTEGER4, irank, mtag1, MPIcomm,irequest1, ierr)
      !Pnods
      call MPI_ISEND(SBuffer_Pnods(1,1),2*(SBuffer_Nelem+1), MPI_INTEGER4, irank, mtag2, MPIcomm,irequest2, ierr) 
      !Lnods
      lnods_size = SBuffer_Pnods(1,SBuffer_Nelem+1)-1
      call MPI_ISEND(SBuffer_Lnods(1),lnods_size, MPI_INTEGER4, irank, mtag3, MPIcomm,irequest3, ierr) 
   end subroutine
   
end module

subroutine reaLnods(a)
   use MPI
   use typre
   use Mod_Int2Str
   use Mod_Mesh
   use Mod_NaivePartition
   use Mod_ReaLnods
   implicit none
   class(FemMesh) :: a

   call a%Timer%reaLnods%Tic
   
   !Serial
   if (a%kfl_ReadType == 0) then
      call reaLnodsSerial(a)
      
   !Partitioned
   elseif (a%kfl_ReadType == 1) then
      call reaLnodsPartitioned(a)
   endif
   
      
   call a%Timer%reaLnods%Toc
   
end subroutine

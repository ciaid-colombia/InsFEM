module Mod_MeshInterpolator
   use typre
   use MPI
   use Mod_Memor
   use Mod_Octree
   use Mod_MpiObject
   use Mod_Mesh
   use Mod_Element
   use Mod_postpr
   implicit none
   private
   public Interpolator
   public SubdomainRangeGiver
   public ElemRangeGiver

   type, extends(MpiObject)   :: Interpolator
      class(FemMesh), pointer :: Mesh => NULL()
      class(FiniteElement), pointer :: e => NULL()
      type(i1p), allocatable :: PointToProcOthersNeed(:),IsFoundINeed(:),IsFoundOthersNeed(:)
      type(i1p), allocatable :: PointsINeed(:),PointsOthersNeed(:),NodesINeed(:),NodesOthersNeed(:)
      type(r2p), allocatable :: CoordsINeed(:), CoordsOthersNeed(:),ShafunOthersNeed(:)
      type(r2p), allocatable :: ArrayINeed(:),ArrayOthersNeed(:)
      integer(ip), allocatable :: NpoinINeed(:),NpoinOthersNeed(:),NcoordINeed(:),NcoordOthersNeed(:)
      real(rp), allocatable :: range(:,:,:)
      real(rp), pointer :: coord(:,:) => NULL()
      real(rp) :: Tol = 0.0001_rp
      integer(ip) :: npoin1,npoin2,ndime
      integer(ip) :: kfl_Init = -1_ip

contains
      procedure :: Initialize_Interp
      procedure :: Interpolate
      procedure :: Finalize
      procedure :: IsInitialized
      procedure :: SetTol
      
      !This is implemented in a particular subroutine because a different extension is needed
      ! (for periodic boundary conditions)
      procedure :: ModifyShafun
   end type

   type, extends(RangeGiver)  :: ElemRangeGiver
      class(FemMesh), pointer :: Mesh => NULL()
      class(FiniteElement), pointer :: e => NULL()
      type(MemoryMan), pointer      :: Memor => NULL()
      integer(ip)             :: ndime,nelem
      real(rp), allocatable   :: range(:,:,:)
contains      
      procedure :: GiveMeRange => ElemRange
      procedure :: InitializeRange
      procedure :: FinalizeRange
   end type

   type, extends(RangeGiver) :: SubdomainRangeGiver
      integer(ip)            :: ndime
      real(rp), pointer      :: range(:,:,:) => NULL()
contains
      procedure :: GiveMeRange => SDGiveRange
      procedure :: SDInitialize
   
   end type

contains
      
   subroutine ElemRange(a,ielem,range)
      use typre
      implicit none
      class(ElemRangeGiver) :: a
      integer(ip) :: ielem
      real(rp)    :: range(*)

      integer(ip) :: idime,pnode
      

      do idime = 1,a%ndime
         range((idime-1)*2+1) = a%range(1,idime,ielem)
         range((idime-1)*2+2) = a%range(2,idime,ielem)
      end do
   end subroutine

   subroutine InitializeRange(a,Mesh,Memor)
      use Mod_Mesh
      use Mod_Element
      use Mod_Memor
      implicit none
      class(ElemRangeGiver) :: a
      class(FemMesh), target :: Mesh
      type(MemoryMan), target :: Memor
      integer(ip) :: idime,ielem,pnode,ndime,nelem
      
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
      a%ndime = ndime
      a%nelem = nelem
      a%Memor => Memor
      call a%Memor%alloc(2,a%ndime,a%nelem,a%range,'range','MeshInit')

      a%Mesh => Mesh
      call a%Mesh%ElementAlloc(a%e,a%Memor,'ForceClosedRule','MeshInit')   

      do ielem = 1,a%nelem
         call a%Mesh%ElementLoad(ielem,a%e)
         pnode = a%e%pnode

         do idime = 1,a%ndime
            a%range(1,idime,ielem) = minval(a%e%elcod(idime,1:pnode))
            a%range(2,idime,ielem) = maxval(a%e%elcod(idime,1:pnode))
         end do
      enddo
      call a%Mesh%ElementDeAlloc(a%e,a%Memor,'ForceClosedRule','MeshInit')   
   end subroutine
   
   subroutine FinalizeRange(a)
      implicit none
      class(ElemRangeGiver) :: a
      
      call a%Memor%dealloc(2,a%ndime,a%nelem,a%range,'range','MeshInit')
   end subroutine   
      
   subroutine SDGiveRange(a,ielem,range)
      implicit none
      class(SubdomainRangeGiver) :: a
      integer(ip) :: ielem
      real(rp) :: range(*)
      
      integer(ip) :: icount, idime, j
      
      icount = 0
      do idime = 1,a%ndime
         do j = 1,2
            icount = icount +1
            range(icount) = a%range(j,idime,ielem)
         enddo
      enddo
   end subroutine

   subroutine SDInitialize(a,nelem,ndime,range)
      implicit none
      class(SubdomainRangeGiver) :: a
      integer(ip) :: ndime
      integer(ip) :: nelem
      real(rp), target :: range(2,ndime,nelem)
      
      a%ndime = ndime
      a%range => range
   end subroutine

   subroutine Initialize_Interp(a,Mesh,coord,DoFlag)
      use typre
      use Mod_Memor
      use Mod_Octree
      use Mod_Mesh
      use Mod_MpiObject
      use MPI
      implicit none

      class(FemMesh), target :: Mesh
      real(rp), target :: coord(:,:)      
      logical, optional :: DoFlag(:)
      
      type(Octree) :: SubdomainTree,ElemTree
      class(Interpolator) :: a
      type(ElemRangeGiver) :: PowerRanger
      type(SubdomainRangeGiver) :: SDPowerRanger
      integer(ip) :: nelem,maxpoin,ndime,npoin1,npoin2,ierr
      real(rp) :: myrange(2,3)
      logical :: isALE
      
      interface 
         subroutine mshint_ComputeMaxpoin(Mesh,maxpoin)
            use typre
            use Mod_Mesh
            implicit none
            class(FemMesh) :: Mesh
            integer(ip) :: maxpoin
         end subroutine
      end interface

      call mshint_ComputeMaxpoin(Mesh,maxpoin)
     
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
      call Mesh%GetNpoin(npoin1)
      a%Mesh => Mesh
      a%coord => coord
      call a%Mesh%GetNdime(ndime)
      a%ndime = ndime
      a%npoin1 = npoin1
      a%npoin2 = size(coord,2)

      call a%SetLOutputFile(Mesh%lun_memo,Mesh%lun_outpu)

      call a%SetMPI(Mesh%MPIcomm,Mesh%MPIsize,Mesh%MPIroot,Mesh%MPIrank)

      call a%Mesh%ElementAlloc(a%e,a%Memor,'ForceClosedRule','MeshInit')

      !First compute my range
      call a%Mesh%GetALE(isALE)
      if (isALE) then
         call a%Mesh%GetDeformedRange(myrange)
      else
         call a%Mesh%GetRange(myrange)
      endif
   
      !Now send my range to everybody, and receive everybody's range
      !First allocate SubdomainRange 
      call a%Memor%alloc(2,a%ndime,a%MPIsize,a%range,'range','SubdomainRange')
      call MPI_AllGather(myrange, 2*a%ndime, MPI_REAL8, a%range, 2*a%ndime, MPI_REAL8, a%MPIcomm,ierr)

      !Prepare an Octree of Subdomains
      call SDPowerRanger%SDInitialize(a%MPIsize,ndime,a%range)
      call SubdomainTree%InitializeElementOctree(ndime,a%MPIsize,SDPowerRanger,maxpoin,a%Memor)

      !Prepare an Octree of Elements
      call PowerRanger%InitializeRange(Mesh,a%Memor)
      call ElemTree%InitializeElementOctree(ndime,nelem,PowerRanger,maxpoin,a%Memor)

      if (.not. present(DoFlag)) then
         call LinkNodesWithSubdomains(a,SubdomainTree,SDPowerRanger)
      else
         call LinkNodesWithSubdomains(a,SubdomainTree,SDPowerRanger,DoFlag)
      endif

      call Prepare_Interpolation(a,ElemTree,PowerRanger)
      call PowerRanger%FinalizeRange
      a%kfl_Init = 1_ip
      call a%Memor%dealloc(2,a%ndime,a%MPIsize,a%range,'range','SubdomainRange')
      call a%Mesh%ElementDeAlloc(a%e,a%Memor,'ForceClosedRule','MeshInit')
   end subroutine

   subroutine Finalize(a)
      use typre
      use Mod_Memor
      use Mod_Octree
      use Mod_Mesh
      use Mod_MpiObject
      use Mod_LinkedList
      implicit none
      class(Interpolator) :: a
      integer(ip) :: irank,ndofn

      do irank = 0,a%MPIsize-1
         if (a%NcoordOthersNeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NcoordOthersNeed(irank+1),a%IsFoundOthersNeed(irank+1)%l,'IsFoundOthersNeed','Finalize')
            call a%Memor%pdealloc(a%NcoordOthersNeed(irank+1),a%PointToProcOthersneed(irank+1)%l,'PointToProcOthersneed','Finalize')
            call a%Memor%pdealloc(a%NcoordOthersNeed(irank+1),a%NodesOthersneed(irank+1)%l,'PointToProcOthersneed','Finalize')
            call a%Memor%pdealloc(size(a%ShafunOthersNeed(irank+1)%a,1),size(a%ShafunOthersNeed(irank+1)%a,2),a%ShafunOthersNeed(irank+1)%a,'shafun','Finalize')
            call a%Memor%pdealloc(a%NcoordOthersNeed(irank+1),a%PointsOthersNeed(irank+1)%l,'PointsINeed','Finalize')
            call a%Memor%pdealloc(a%ndime,a%NcoordOthersNeed(irank+1),a%CoordsOthersNeed(irank+1)%a,'CoordsOthersNeed','Finalize')
         endif
         if (a%NcoordINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NcoordINeed(irank+1),a%IsFoundINeed(irank+1)%l,'IsFoundINeed','Finalize')
            call a%Memor%pdealloc(a%ndime,a%NcoordINeed(irank+1),a%CoordsINeed(irank+1)%a,'CoordsINeed','Finalize')
            call a%Memor%pdealloc(a%NcoordINeed(irank+1),a%PointsINeed(irank+1)%l,'PointsINeed','Finalize')
         endif
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(a%NpoinINeed(irank+1),a%NodesINeed(irank+1)%l,'NodesINeed','Finalize')
         endif
      enddo

      call a%Memor%dealloc(a%MPIsize,a%NcoordOthersNeed,'NcoordOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%NcoordINeed,'NcoordINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%PointToProcOthersNeed,'PoinToProcOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%NodesINeed,'NodesINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%IsFoundINeed,'IsFoundINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%NpoinINeed,'NpoinINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%NpoinOthersNeed,'NpoinOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%ShafunOthersNeed,'ShafunOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%CoordsINeed,'CoordsINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%CoordsOthersNeed,'CoordsOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%PointsINeed,'PointsINeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%PointsOthersNeed,'PointsOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%NodesOthersNeed,'NodesOthersNeed','Finalize')
      call a%Memor%dealloc(a%MPIsize,a%IsFoundOthersNeed,'IsFoundOthersNeed','Finalize')

      a%kfl_Init = -1

   end subroutine

   subroutine LinkNodesWithSubdomains(a,SubdomainTree,SDPowerRanger,DoFlag)
      use typre
      use Mod_Memor
      use Mod_Octree
      use Mod_Mesh
      use Mod_LinkedList
      use MPI
      implicit none
      class(Interpolator) :: a
      type(Octree) :: SubdomainTree
      type(SubdomainRangeGiver) :: SDPowerRanger
      logical, optional :: DoFlag(:)
      
      class(LinkedListHeader), pointer  :: PerHead => NULL()
      type(LinkedListHeader), allocatable  :: pointlistH(:)
      class(LinkedListStorage), pointer :: PerStore => NULL()
      type(LinkedListStorage)  :: pointlistS
      real, pointer :: coord(:) => NULL()
      integer(ip) :: ipoin,irank,idime,value,ipos,gpoin
      real(rp), allocatable :: coordlist(:,:,:)
      real(rp) :: Tol
      logical :: IsIn
      integer :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize),istatus,i
      integer, parameter :: mtag1 = 1
      integer status(MPI_STATUS_SIZE)

      Tol = 0.0001_rp

      call a%Memor%alloc(a%MPIsize,a%NcoordOthersNeed,'NcoordOthersNeed','LinkNodesWithSubdomains')
      call a%Memor%alloc(a%MPIsize,a%NcoordINeed,'NcoordINeed','LinkNodesWithSubdomains')

      allocate(pointlistH(a%MPIsize),stat=istatus)
      call a%Memor%allocObj(istatus,'pointlistH','LinkNodesWithSubdomains',a%MPIsize)
      do irank = 0,a%MPIsize-1
         call pointlistH(irank+1)%Initialize
      enddo

      call pointlistS%Initialize(a%npoin2,'NonRemov','Repeatable',a%Memor)

      do ipoin=1,a%npoin2
         !If the optional parameter is present, we will only interpolate the marked ones
         if (present(DoFlag)) then
            if (DoFlag(ipoin) .eqv. .false.) cycle
         endif
         
         call SubdomainTree%GetListForCoord(a%coord(:,ipoin),PerHead,PerStore)
         call PerStore%GoToFirstOfList(PerHead)
         call PerStore%GetNext(irank)
         do while (irank /= -1)
            !Check if I really need to test it
            IsIn = .true.
            do idime = 1,a%ndime
               if (a%coord(idime,ipoin) >= SDPowerRanger%range(1,idime,irank)-Tol .and. a%coord(idime,ipoin) <= SDPowerRanger%range(2,idime,irank)+Tol) then  
                  !The node is not out
               else
                  IsIN = .false.
               endif
            enddo
            !If it is in, then add the subdomain to the list of check subdomains
            if (IsIn .eqv. .true.) then
               a%NcoordINeed(irank) = a%NcoordINeed(irank) + 1
               call pointlistS%Add(pointlistH(irank),ipoin)
            endif
            call PerStore%GetNext(irank)
         enddo
      enddo

      call MPI_AllToAll(a%NcoordINeed, 1, MPI_INTEGER4, a%NcoordOthersNeed, 1, MPI_INTEGER4, a%MPIcomm,ierr)

      call a%Memor%alloc(a%MPIsize,a%CoordsINeed,'CoordsINeed','LinkNodesWithSubdomains')
      call a%Memor%alloc(a%MPIsize,a%CoordsOthersNeed,'CoordsOthersNeed','LinkNodesWithSubdomains')
      call a%Memor%alloc(a%MPIsize,a%PointsINeed,'PointsINeed','LinkNodesWithSubdomains')
      call a%Memor%alloc(a%MPIsize,a%PointsOthersNeed,'PointsOthersNeed','LinkNodesWithSubdomains')

      do irank = 0,a%MPIsize-1
         if (a%NcoordINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%ndime,a%NcoordINeed(irank+1),a%CoordsINeed(irank+1)%a,'CoordsINeed','LinkNodesWithSubdomains')
            call a%Memor%palloc(a%NcoordINeed(irank+1),a%PointsINeed(irank+1)%l,'PointsINeed','LinkNodesWithSubdomains')
         endif
      enddo
      
      do irank=0,a%MPIsize-1
         call pointlistS%GoToFirstOfList(pointlistH(irank+1))
         do ipoin=1,a%NcoordINeed(irank+1)
            call pointlistS%GetNext(value)
            a%PointsINeed(irank+1)%l(ipoin) = value
            a%CoordsINeed(irank+1)%a(:,ipoin) = a%coord(1:a%ndime,value)
         enddo
      enddo

      call pointlistS%Dealloc
      deallocate(pointlistH,stat=istatus)
      call a%Memor%deallocObj(istatus,'pointlistH','LinkNodesWithSubdomains',a%MPIsize)

      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%ndime,a%NcoordOthersNeed(irank+1),a%CoordsOthersNeed(irank+1)%a,'CoordsOthersNeed','LinkNodesWithSubdomains')
            call MPI_IRECV(a%CoordsOthersNeed(irank+1)%a,a%NcoordOthersNeed(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr) 
         endif                         
      enddo

      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordINeed(irank+1) /= 0) then
            call MPI_ISEND(a%CoordsINeed(irank+1)%a,a%NcoordINeed(irank+1)*a%ndime,MPI_REAL8,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr)
         endif 
      enddo    

      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1), status, ierr)
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo

      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NcoordOthersNeed(irank+1),a%PointsOthersNeed(irank+1)%l,'PointsINeed','LinkNodesWithSubdomains')
            call MPI_IRECV(a%PointsOthersNeed(irank+1)%l,a%NcoordOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr)  
         endif                         
      enddo

      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordINeed(irank+1) /= 0) then
            call MPI_ISEND(a%PointsINeed(irank+1)%l,a%NcoordINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr)
         endif 
      enddo    

      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1), status, ierr)
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo

      call SubdomainTree%Dealloc(a%Memor)
   end subroutine
   
   subroutine ModifyShafun(a,newshape)
      class(Interpolator) :: a
      real(rp) :: newshape(:)
      
      !This does nothing by default
      !Required for the round interpolator in PeriodicBoundaryConditions
       
   end subroutine

   subroutine Prepare_Interpolation(a,ElemTree,PowerRanger)
      use typre
      use Mod_Memor
      use Mod_Octree
      use Mod_LinkedList
      use Mod_Mesh
      use MPI
      implicit none
      class(Interpolator) :: a
      class(ElemRangeGiver) :: PowerRanger
      type(Octree) :: ElemTree
      type(MemoryMan) :: Memor
      class(LinkedListHeader), pointer :: Head => NULL()
      class(LinkedListStorage), pointer :: Store => NULL()
      real(rp), allocatable :: coord(:)
      real(rp) :: xloc(a%ndime),newshape(200),isorange(2),Tol
      integer(ip) :: irank,ipos,inode,ielem,foundelem,ipoin,poini,idime,flag
      integer :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize)
      integer, parameter :: mtag1 = 1, mtag2 = 2
      integer status(MPI_STATUS_SIZE)

      call a%Memor%alloc(a%MPIsize,a%PointToProcOthersNeed,'PoinToProcOthersNeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%NodesINeed,'NodesINeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%NodesOthersNeed,'NodesOthersNeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%IsFoundINeed,'IsFoundINeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%IsFoundOthersNeed,'IsFoundOthersNeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%NpoinINeed,'NpoinINeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%NpoinOthersNeed,'NpoinOthersNeed','Prepare_Interpolation')
      call a%Memor%alloc(a%MPIsize,a%ShafunOthersNeed,'ShafunOthersNeed','Prepare_Interpolation')
      call a%Memor%alloc(PowerRanger%ndime,coord,'coord','Prepare_Interpolation')
      call a%Mesh%ElementAlloc(PowerRanger%e,PowerRanger%Memor,'ForceClosedRule','Prepare_Interpolation')   

      do irank = 0,a%MPIsize-1
         if (a%NcoordOthersNeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NcoordOthersNeed(irank+1),a%IsFoundOthersNeed(irank+1)%l,'IsFoundOthersNeed','Prepare_Interpolation')
            call a%Memor%palloc(a%NcoordOthersNeed(irank+1),a%PointToProcOthersneed(irank+1)%l,'PointToProcOthersneed','Prepare_Interpolation')
            call a%Memor%palloc(a%NcoordOthersNeed(irank+1),a%NodesOthersneed(irank+1)%l,'PointToProcOthersneed','Prepare_Interpolation')
            call a%Memor%palloc(a%NcoordOthersNeed(irank+1),a%e%mnode,a%ShafunOthersNeed(irank+1)%a,'shafun','Prepare_Interpolation')
        
            do ipos = 1,a%NcoordOthersNeed(irank+1)
               coord(:) = a%CoordsOthersNeed(irank+1)%a(:,ipos)

               !Octree for elements
               call ElemTree%GetListForCoord(coord,Head,Store)
               call Store%GoToFirstOfList(Head)      
               call Store%GetNext(ielem)
               flag = 0_ip

               elemlist: do while (ielem /= -1 .and. flag==0_ip)
                  call PowerRanger%Mesh%ElementLoad(ielem,PowerRanger%e)
                  call PowerRanger%Mesh%GetIsoRange(PowerRanger%e%pnode,isorange)
                  call PowerRanger%e%isoparinv(coord,xloc) 
                  if (isorange(1)==0.0) then 
                     if (sum(xloc)>1.0) then
                        foundelem = -1
                     else
                        foundelem = 0
                     endif
                  else
                     foundelem = 0
                  endif
                  do idime = 1,PowerRanger%ndime
                     if (xloc(idime) >= isorange(1) .and. xloc(idime) <= isorange(2)) then
                        foundelem = foundelem + 1
                     endif
                  enddo

                  if (foundelem == PowerRanger%ndime) then
                     !Do not consider elements with negative jacobian (folded elements FMALE)
                     call PowerRanger%e%elmdcg
                     if (PowerRanger%e%detjm < 0.0_rp) then
                        !do nothing, the element is folded
                     else
                        flag = 1_ip
                        a%NpoinOthersNeed(irank+1) = a%NpoinOthersNeed(irank+1) + 1
                        a%IsFoundOthersNeed(irank+1)%l(ipos) = 1_ip
                        a%PointToProcOthersNeed(irank+1)%l(a%NpoinOthersNeed(irank+1)) = ielem
                        a%NodesOthersNeed(irank+1)%l(a%NpoinOthersNeed(irank+1)) = a%PointsOthersNeed(irank+1)%l(ipos)
                        
                        call PowerRanger%e%ComputeShafun(xloc,newshape)
                        call a%ModifyShafun(newshape)
                        do inode = 1,PowerRanger%e%pnode
                           a%ShafunOthersNeed(irank+1)%a(a%NpoinOthersNeed(irank+1),inode) = newshape(inode)
                        end do
                     endif
                  endif 
                  call Store%GetNext(ielem) 
               enddo elemlist

               if (flag==0_ip) then
                  call Store%GoToFirstOfList(Head)
                  call Store%GetNext(ielem)
                  do while (ielem /= -1 .and. flag==0_ip)
                     call PowerRanger%Mesh%ElementLoad(ielem,PowerRanger%e)
                     call PowerRanger%Mesh%GetIsoRange(PowerRanger%e%pnode,isorange)
                     call PowerRanger%e%isoparinv(coord,xloc) 
                     call PowerRanger%e%ComputeShafun(xloc,newshape)
                     
                     if (isorange(1)==0.0) then 
                        if (sum(xloc) > 1.0_rp+a%Tol) then
                           foundelem = -1
                        else
                           foundelem = 0
                        endif
                     else
                        foundelem = 0
                     endif
                     do idime = 1,PowerRanger%ndime
                        if (xloc(idime) >= (isorange(1) - a%Tol) .and. xloc(idime) <= (isorange(2) + a%Tol)) then
                           foundelem = foundelem + 1
                        endif
                     enddo
                  
                     if (foundelem == PowerRanger%ndime) then
                        flag = 1_ip
                        a%NpoinOthersNeed(irank+1) = a%NpoinOthersNeed(irank+1) + 1
                        a%IsFoundOthersNeed(irank+1)%l(ipos) = 1_ip
                        a%PointToProcOthersNeed(irank+1)%l(a%NpoinOthersNeed(irank+1)) = ielem
                        a%NodesOthersNeed(irank+1)%l(a%NpoinOthersNeed(irank+1)) = a%PointsOthersNeed(irank+1)%l(ipos)
                        call a%ModifyShafun(newshape)
                        do inode = 1,PowerRanger%e%pnode
                           a%ShafunOthersNeed(irank+1)%a(a%NpoinOthersNeed(irank+1),inode) = newshape(inode)
                        end do
                     endif
                     call Store%GetNext(ielem)
                  enddo  
               endif
            enddo
         endif
      enddo

      call MPI_AllToAll(a%NpoinOthersNeed, 1, MPI_INTEGER4, a%NpoinINeed, 1, MPI_INTEGER4, a%MPIcomm,ierr)

      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(a%IsFoundOthersNeed(irank+1)%l,a%NcoordOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm, &
            irequest1(irank+1),ierr)
         endif 
      enddo

      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NcoordINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NcoordINeed(irank+1),a%IsFoundINeed(irank+1)%l,'IsFoundINeed','Prepare_Interpolation')
            call MPI_IRECV(a%IsFoundINeed(irank+1)%l,a%NcoordINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr)
         endif
      enddo

      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1),status,ierr)
         call MPI_WAIT(irequest2(irank+1),status,ierr)
      enddo

      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(a%NodesOthersNeed(irank+1)%l(1:a%NpoinOthersNeed(irank+1)),a%NpoinOthersNeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm, &
            irequest1(irank+1),ierr)
         endif 
      enddo

      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%palloc(a%NpoinINeed(irank+1),a%NodesINeed(irank+1)%l,'NodesINeed','Prepare_Interpolation')
            call MPI_IRECV(a%NodesINeed(irank+1)%l,a%NpoinINeed(irank+1),MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr)
         endif
      enddo

      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1),status,ierr)
         call MPI_WAIT(irequest2(irank+1),status,ierr)
      enddo

      call a%Memor%dealloc(PowerRanger%ndime,coord,'coord','Prepare_Interpolation')
      call a%Mesh%ElementDeAlloc(PowerRanger%e,PowerRanger%Memor,'ForceClosedRule','Prepare_Interpolation')   
      call ElemTree%Dealloc(a%Memor)
   end subroutine

   subroutine Interpolate(a,ndofn,array1,array2)
      use typre
      use Mod_Mesh
      use Mod_Element
      use MPI
      implicit none

      class(Interpolator) :: a
      class(FiniteElement), pointer :: e => NULL()
      real(rp) :: array1(ndofn,a%npoin1),array2(ndofn,a%npoin2)
      integer(ip) :: irank,ipos,ipoin,jpoin,npoin,ielem,inode,poini,poinj,iposaux,gpoin,ndofn
      integer :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize)
      integer status(MPI_STATUS_SIZE)
      integer, parameter :: mtag1 = 1, mtag2 = 2

      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','Interpolate')
      call a%Memor%alloc(a%MPIsize,a%ArrayINeed,'ArrayINeed','Interpolate')
      call a%Memor%alloc(a%MPIsize,a%ArrayOthersNeed,'ArrayOthersNeed','Interpolate')

      do irank = 0,a%MPIsize-1
         if (a%NPoinOthersNeed(irank+1)/=0) then
            call a%Memor%palloc(ndofn,a%NPoinOthersNeed(irank+1),a%ArrayOthersNeed(irank+1)%a,'ArrayOthersNeed','Interpolate')
            do ipos = 1,a%NpoinOthersNeed(irank+1)
               ielem = a%PointToProcOthersNeed(irank+1)%l(ipos)
               call a%Mesh%ElementLoad(ielem,e)
               do inode = 1,e%pnode
                  jpoin = e%lnods(inode)
                  a%ArrayOthersNeed(irank+1)%a(1:ndofn,ipos) = a%ArrayOthersNeed(irank+1)%a(1:ndofn,ipos) + &
                  a%ShafunOthersNeed(irank+1)%a(ipos,inode)*array1(1:ndofn,jpoin)
               end do
            end do
         endif
      enddo

      irequest2 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%palloc(ndofn,a%NpoinINeed(irank+1),a%ArrayINeed(irank+1)%a,'ArrayINeed','Interpolate')
            call MPI_IRECV(a%ArrayINeed(irank+1)%a,ndofn*a%NpoinINeed(irank+1),MPI_REAL8,irank,mtag1,a%MPIcomm,irequest2(irank+1),ierr)
         endif
      enddo

      irequest1 = MPI_REQUEST_NULL
      do irank = 0,a%MPIsize-1
         if (a%NpoinOthersNeed(irank+1) /= 0) then
            call MPI_ISEND(a%ArrayOthersNeed(irank+1)%a,ndofn*a%NPoinOthersNeed(irank+1),MPI_REAL8,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr)
         endif
      enddo

      do irank = 0,a%MPIsize-1
         call MPI_WAIT(irequest1(irank+1), status, ierr)
         call MPI_WAIT(irequest2(irank+1), status, ierr)
      enddo

      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1)/=0) then
            iposaux = 0
            do ipos = 1,a%NcoordINeed(irank+1)
               if (a%IsFoundINeed(irank+1)%l(ipos)==1) then
                  iposaux = iposaux + 1
                  ipoin = a%NodesINeed(irank+1)%l(iposaux)
                  array2(1:ndofn,ipoin) = a%ArrayINeed(irank+1)%a(1:ndofn,iposaux)
               endif
            enddo
         endif
      enddo

      do irank = 0,a%MPIsize-1
         if (a%NpoinINeed(irank+1) /= 0) then
            call a%Memor%pdealloc(ndofn,a%NpoinINeed(irank+1),a%ArrayINeed(irank+1)%a,'ArrayINeed','Interpolate')
         endif
         if (a%NPoinOthersNeed(irank+1)/=0) then
            call a%Memor%pdealloc(ndofn,a%NPoinOthersNeed(irank+1),a%ArrayOthersNeed(irank+1)%a,'ArrayOthersNeed','Interpolate')
         endif
      enddo

      call a%Memor%dealloc(a%MPIsize,a%ArrayINeed,'ArrayINeed','Interpolate')
      call a%Memor%dealloc(a%MPIsize,a%ArrayOthersNeed,'ArrayOthersNeed','Interpolate')
      call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','Interpolate')

   end subroutine

   subroutine IsInitialized(a,kfl_Init)
      use typre
      implicit none
      class(Interpolator) :: a
      integer(ip) :: kfl_Init

      kfl_Init = a%kfl_Init
   end subroutine

   subroutine SetTol(a,Tol)
      use typre
      implicit none
      class(Interpolator) :: a
      real(rp) :: Tol

      a%Tol = Tol
   end subroutine

end module

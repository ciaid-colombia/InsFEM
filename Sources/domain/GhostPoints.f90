subroutine GhostPointsAndLocalOrdering(a)
   use MPI
   use typre
   use Mod_PointToProc
   use Mod_Mesh
   implicit none
   
   class(FemMesh), intent(inout) :: a

   logical, allocatable     :: lwa(:)
   integer(ip), allocatable :: iwa(:)
   
   integer(ip) :: minpo,maxpo,ielem,ipoin,pnode,inode,ispos,iwasize
   
   integer(ip), allocatable :: LocalToGlobal(:)
   type(PointToProc) :: Point2Proc
   integer(ip), allocatable :: ProcList(:)
   
   !MPI
   integer :: ierr
   
   call a%Timer%GhostPoints%Tic
   
   !Decide which are my ghost points
   minpo = a%poin0 +1
   maxpo = a%poin0 + a%npoinLocal
   
   call a%Memor%alloc(a%npoinLocal,iwa ,'iwa' ,'ScatterData')
   iwasize = a%npoinLocal
   call a%Memor%alloc(a%gnpoin,lwa,'lwa','ScatterData')
   
   a%npoinGhost = 0
   do ielem = 1,a%nelem
      call getlnods(ielem,pnode,ispos,a%pnods,a%lnods,a%nelty,a%nnode)
      do inode = 1,pnode
         if (a%lnods(ispos+inode) > maxpo .or. a%lnods(ispos+inode) < minpo) then !point is a ghost!
            if (lwa(a%lnods(ispos+inode)) .eqv. .false.) then
               lwa(a%lnods(ispos+inode)) = .true.
               a%npoinGhost = a%npoinGhost+1
               if (a%npoinGhost > iwasize) then
                  iwasize = iwasize*2
                  call a%Memor%Realloc(iwasize,iwa,'iwa','ScatterData')
               endif
               iwa(a%npoinGhost) = a%lnods(ispos+inode)
            endif
         endif
      enddo
   enddo
   a%npoin = a%npoinLocal + a%npoinGhost
   
   !Build  the LocalToGlobal Array
   call a%Memor%alloc(a%npoinLocal+a%npoinGhost,LocalToGlobal,'LocalToGlobal','ScatterData')
   do ipoin = 1,a%npoinLocal
      LocalToGlobal(ipoin) = ipoin + a%poin0
   enddo
   LocalToGlobal(a%npoinLocal+1:a%npoinLocal+a%npoinGhost) = iwa(1:a%npoinGhost)
   
      
   call a%Memor%dealloc(a%gnpoin,lwa,'lwa','ScatterData')
   call a%Memor%dealloc(size(iwa), iwa,'iwa','ScatterData')
   
   call MPI_REDUCE( a%npoinGhost, a%gnpoinGhost, 1, MPI_INTEGER4, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call MPI_REDUCE( a%npoinGhost, a%MaxNpoinGhostPerProc, 1, MPI_INTEGER4, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
   
   call a%Timer%GhostPoints%Toc
   
!----------------------------------------------------------------------------   
   
   call a%Timer%BuildLocalOrdering%Tic
   
   
   
   !I also want to store the processor number of each ghost point. At this point, I do not know it
   call a%Memor%alloc(a%npoin,ProcList,'ProcList','GhostPoints')
   
   call Point2Proc%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call Point2Proc%SetMemor(a%Memor)
   call Point2Proc%Initialize(a%poin0)
   ProcList(1:a%npoinLocal) = a%MPIrank
   do ipoin = 1,a%npoinGhost
      call Point2Proc%GetProc(LocalToGlobal(a%npoinLocal+ipoin),ProcList(a%npoinLocal+ipoin))
   enddo
   call Point2Proc%Finalize
   
   !Create and Initialize the LocalOrdering
   call a%ParallelLibrary%CreateOrdering(a%LocalOrdering,a%Memor)  !FACTORY
   call a%LocalOrdering%Init(a%MPIcomm,a%npoinLocal+a%npoinGhost,LocalToGlobal,ProcList,a%Memor)
   
   
   !Ghost Point List deallocate
   call a%Memor%dealloc(a%npoin,ProcList,'ProcList','BuildArrayCommunicator')
   call a%Memor%dealloc(a%npoinLocal+a%npoinGhost,LocalToGlobal,'LocalToGlobal','ScatterData')

   
   call a%Timer%BuildLocalOrdering%Toc
end subroutine  

subroutine BuildArrayCommunicator(a)
   use typre
   use Mod_Mesh
   implicit none
   
   class(FemMesh), intent(inout) :: a
   
   integer(ip), allocatable :: ProcList(:)
   
   call a%Timer%BuildArrayCommunicator%Tic
   
   call a%ParallelLibrary%CreateCommunicator(a%ArrayCommunicator,a%Memor)
   call a%ArrayCommunicator%Init(a%MPIcomm,a%MPIsize,a%MPIrank,a%npoinLocal,a%npoinGhost,a%gnpoin,a%LocalOrdering,a%Memor)
   
   call a%Timer%BuildArrayCommunicator%Toc
   
end subroutine


   
subroutine Ghostcoord(a)
   use MPI
   use typre
   use Mod_Mesh
   implicit none
   
   class(FemMesh), intent(inout) :: a
   
   integer(ip) :: ierr
   
   call a%Timer%GhostCoord%Tic

   call a%ArrayCommunicator%GhostCommunicate(a%ndime,a%coord)
   
   call a%Timer%GhostCoord%Toc
   
end subroutine   

subroutine LnodsToLocal(a)
   use typre
   use Mod_Mesh
   implicit none
   
   class(FemMesh), intent(inout) :: a   
   
   integer(ip) :: iaux,idummy,ierr
   
   call a%Timer%LnodsToLocal%Tic
   
   !Transform lnods to local numbering
   iaux = size(a%lnods)
   call a%LocalOrdering%Global2Local(iaux,a%lnods,a%lnods)   

   call a%Timer%LnodsToLocal%Tic
end subroutine
   

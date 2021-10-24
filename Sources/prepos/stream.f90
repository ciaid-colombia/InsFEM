module Mod_stream
contains

subroutine stream(vecto,ndime,istep,ttime,Mesh,Memor,FilePostpr,MPIcomm,MPIroot,MPIrank,MPISize)
   !This subroutine computes the stream function in parallel
   use typre
   use Mod_Mesh
   use Mod_Postpr
   use Mod_Element
   use Mod_memor
   use MPI
   implicit none
   type(MemoryMan) :: Memor
   class(FemMesh)    :: Mesh  
   class(FiniteElement), pointer :: e => NULL() 
   class(PostprFile) :: FilePostpr   
   
   real(rp), intent(in)     :: vecto(ndime,*),ttime
   integer(ip)              :: istat,kpoin,itouc,MPIrank,MPIsize,MPIroot,MPIcomm
   integer(ip) , intent(in)  :: istep  
   integer(ip)              :: metmp(2)=0
   real(rp),    allocatable :: strea(:),auxstrea(:)
   integer(ip), allocatable :: lpoin(:)
   integer(ip), allocatable :: nocha(:,:),noinv(:,:)
   real(rp),    allocatable :: westr(:,:),elvec(:,:)
   integer(ip)              :: ielty,npoin,ndime,ipoin,idime,nelem,ielem,nelty,&
                                 jpoin,jnode,inode,ispos,one,kelem
   integer(ip), pointer     :: nnode(:) => NULL(),ngaus(:) => NULL(),ltopo(:) => NULL(),lelpo(:) => NULL()
   integer(ip)              :: pelpo
   real(rp), pointer        :: coord(:)  => NULL()
   integer(ip)              :: auxpoin
   
   integer(ip)              :: IAmDone,WeAreDone
   integer                  :: ierr
   
   integer(ip)              :: IAmLeveled,WeAreLeveled
   integer(ip)              :: npoinLocal,npoinGhost
   integer(ip)              :: nleveledpoints
   real(rp)                 :: leveljump
   real(rp)                 :: nan
   integer(ip)              :: ielemcounter
   logical, allocatable     :: ElementDone(:)
   integer(ip), allocatable :: PointList(:)
   integer(ip)              :: countPointList,elemi,poini

   
   !Mesh parameters used
   call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
   call Mesh%GetNdime(ndime)
   call Mesh%GetNelem(nelem)
   call Mesh%GetNpoin(npoin) 
   !Element object Initialization
   call Mesh%ElementAlloc(e,Memor,'DefaultRule','stream')
   
   !Allocate memory
   call Memor%alloc(16,nelty,nocha,'nocha','stream') 
   call Memor%alloc(16,nelty,noinv,'noinv','stream')
   call Memor%alloc(16,nelty,westr,'westr','stream')    
   call Memor%alloc(npoin,strea,'strea','stream')
   call Memor%alloc(npoin,lpoin,'lpoin','stream')  
   call Memor%alloc(ndime,e%mnode,elvec,'elvec','stream')  
   call Memor%alloc(nelem,ElementDone,'ElementDone','stream')
   call Memor%alloc(npoin,PointList,'PointList','stream')
   
 
   !Parameters used in the generation of stream-functions  
   do ielty=1,nelty
      call chanum(&
          ltopo(ielty),nnode(ielty),nocha(1,ielty),&
          noinv(1,ielty),westr(1,ielty))
   end do

   !Compute stream function
   ielem=0
   kpoin=1
   
   ielty = 1   
   
   ielemcounter=0
   ElementDone = .false.
   
   strea(1)=0.0_rp
   lpoin(1)=1
   PointList(1) = 1
   countPointList = 1
   do poini = 1,npoin
      ipoin = PointList(poini)
      call Mesh%GetLelpo(ipoin,pelpo,lelpo)
      do elemi = 1,pelpo
         ielem = lelpo(elemi)
         if (ElementDone(ielem) .eqv. .false.) then
            ElementDone(ielem) = .true.
            ielemcounter=ielemcounter+1
            call Mesh%ElementLoad(ielem,e)
            
            !nodes in the element are done one more time
            do inode = 1,e%pnode
               ipoin = e%lnods(inode)
               if (lpoin(ipoin) == 0) then
                  !Add the point to the list
                  countPointList = countPointList+1
                  PointList(countPointList) = ipoin
               endif
            enddo   
            call e%gather(ndime,elvec,vecto) 
            
            call strloc(e%pnode,ndime,npoin,&
               e%lnods,lpoin,kpoin,e%elcod,elvec,strea,&
               nocha(1,e%ielty),noinv(1,e%ielty),westr(1,e%ielty))
         endif
      enddo
   enddo                                   

   do ipoin=1,npoin
      strea(ipoin)=strea(ipoin)/real(lpoin(ipoin))
   end do
   
   !-------------------------------------------------------------------
   !Levelling
   IAmLeveled=0
   WeAreLeveled=0
   call Memor%alloc(npoin,auxstrea,'auxstrea','stream')
   nan = 0.0_rp
   nan = nan/nan
   auxstrea = nan
   call Mesh%GetNpoinLocal(npoinLocal)
   npoinGhost = npoin-npoinLocal
   
   if (MPIrank == MPIroot) then
      auxstrea = strea
      IAmLeveled=1
   endif
   
   call MPI_REDUCE( IAmLeveled, WeAreLeveled, 1, MPI_INTEGER4, MPI_SUM, MPIroot,MPIcomm, ierr )
   call MPI_BCAST(WeAreLeveled, 1, MPI_INTEGER4, MPIroot, MPIcomm, ierr)
   
   do while (WeAreLeveled /= MPIsize) 
      call Mesh%ArrayCommunicator%GhostCommunicate(1_ip,auxstrea)
      if (IAmLeveled == 0) then
         nleveledpoints=0
         leveljump=0.0_rp
         do ipoin = npoinLocal+1,npoin
            if (isNaN(auxstrea(ipoin)) .eqv. .false.) then
               nleveledpoints = nleveledpoints+1
               leveljump=leveljump+auxstrea(ipoin)-strea(ipoin)
            endif
         enddo
         if (nleveledpoints > 0) then
            leveljump = leveljump/nleveledpoints
            strea = strea + leveljump
            auxstrea = strea
         
            IAmLeveled = 1
         endif
      endif
   
      call MPI_REDUCE( IAmLeveled, WeAreLeveled, 1, MPI_INTEGER4, MPI_SUM, MPIroot,MPIcomm, ierr )
      call MPI_BCAST(WeAreLeveled, 1, MPI_INTEGER4, MPIroot, MPIcomm, ierr)
   enddo   
   
   
   call Memor%dealloc(npoin,auxstrea,'auxstrea','stream')
   
   
   
   !Postprocess stream function
   call FilePostpr%postpr(strea,'stream',istep,ttime,Mesh)
   
   !Deallocate memory
   call Memor%dealloc(16,nelty,nocha,'nocha','stream') !antes nnode era 16 el máximo número de nodos por elemento 2d 
   call Memor%dealloc(16,nelty,noinv,'noinv','stream')
   call Memor%dealloc(16,nelty,westr,'westr','stream')    
   call Memor%dealloc(npoin,strea,'strea','stream')
   call Memor%dealloc(npoin,lpoin,'lpoin','stream')  
   call Memor%dealloc(ndime,e%mnode,elvec,'elvec','stream')  
   call Memor%Dealloc(nelem,ElementDone,'ElementDone','stream')
   call Memor%Dealloc(npoin,PointList,'PointList','stream')
   call Mesh%ElementDealloc(e,Memor,'DefaultRule','stream')
   
   

end subroutine stream

end module

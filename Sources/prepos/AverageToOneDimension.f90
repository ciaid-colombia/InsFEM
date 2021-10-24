subroutine AverageToOneDimension(a,AvgDime,array,wopos,istep,ttime,Mesh,Memor)
   use typre
   use Mod_Element
   use Mod_Memor
   use Mod_Mesh
   use MPI
   use Mod_Postpr
   implicit none
   class(PostprFile) :: a
   integer(ip) :: AvgDime
   real(rp) :: array(:,:)
   character(*), intent(in)  :: wopos
   integer(ip) :: istep
   real(rp) :: ttime
   class(FemMesh)   :: Mesh
   type(MemoryMan) :: Memor
   
   integer(ip) :: ndofn
   
   real(rp) :: range(2,3),subcoord
   integer(ip) :: nsubdivisions
   real(rp)    :: sublength
   
   integer(ip) :: nelem,ielem,isub,maxnode,minnode,maxpoin,minpoin,minsub,maxsub,insub,ipoin
   real(rp)    :: maxcoord,mincoord,dvol,alpha
   
   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable :: AuxSubArray(:,:)
   real(rp),   pointer   :: coord(:) => NULL()
   
   integer(ip) :: npoin,gnpoin
   
   integer(ip) :: gminmaxsub(2)
   integer(ip), allocatable :: rGminmaxsub(:,:)
   real(rp),    allocatable :: rbuffer(:,:)
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer status(MPI_STATUS_SIZE),irequest1,irequest2
   integer :: ierr
   integer(ip) :: irank
   
   integer(ip) :: npoinLocal,local_nodes,inode,ndime
   
   !Open the graf file if I did not do it
   call OpenGrafFileTrue(a)
   
   ndofn = size(array,1)
   
   !GetDimensions
   call Mesh%GetGRange(range)
   call Mesh%GetGNpoin(gnpoin)
   call Mesh%GetNpoinLocal(npoinLocal)
   call Mesh%GetNdime(ndime)
   call Mesh%GetNelem(nelem)
   
   !We divide the AvgDime in how many parts? Depending on the number of points
   !We want 10 subdivisions per point, just to make sure everything is fine
   nsubdivisions = 10.0_rp*(gnpoin**(1.0_rp/ndime))
   sublength = (range(2,Avgdime)-range(1,Avgdime))/nsubdivisions
   
   !Memor allocations
   call Memor%alloc(ndofn+1,nsubdivisions+1,AuxSubArray,'AuxSubArray','AverageToOneDimension')
   call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','vmass')
   if (Mesh%MPIrank == Mesh%MPIroot) then
      call Memor%alloc(2,Mesh%MPIsize,rGminmaxsub,'rGminmaxsub','AverageToOneDimension')
      call Memor%alloc(ndofn+1,nsubdivisions+1,rbuffer,'rbuffer','AverageToOneDimension')
   endif
  
   !Minimum and maximum subdivision contribution for every processor, we send only this range
   !of subdivisions
   gminmaxsub(1) = nsubdivisions+2
   gminmaxsub(2) = -1
  
   do ielem = 1,nelem
      call Mesh%ElementLoad(ielem,e)    
      call e%elmdel
      
      dvol = e%detjm
      
      !Now I need to add the corresponding part to each subelement contribution
      maxcoord = maxval(e%elcod(Avgdime,1:e%pnode))
      maxnode  = maxloc(e%elcod(Avgdime,1:e%pnode),1)
      maxpoin  = e%lnods(maxnode)
      
      mincoord = minval(e%elcod(Avgdime,1:e%pnode))
      minnode  = minloc(e%elcod(Avgdime,1:e%pnode),1)
      minpoin  = e%lnods(minnode)
      
      !I need to add to all the subdivisions
      minsub = ceiling((mincoord-1e-8-range(1,Avgdime))/sublength)+1
      maxsub = floor((maxcoord+1e-8-range(1,Avgdime))/sublength)+1   
      if (maxsub > nsubdivisions+1) maxsub = nsubdivisions+1
      insub = maxsub-minsub+1
      
      gminmaxsub(1) = min(gminmaxsub(1),minsub)
      gminmaxsub(2) = max(gminmaxsub(2),maxsub)
      
      !Count how many of lcal nodes and add only the proportional part
      local_nodes = 0
      do inode = 1,e%pnode
         if (e%lnods(inode) <= npoinLocal) then
            local_nodes = local_nodes+1
         endif
      enddo
      
      dvol = dvol/insub*(1.0_rp*local_nodes)/(1.0_rp*e%pnode)
      
      do isub = minsub,maxsub
         subcoord = range(1,Avgdime)+(isub-1)*sublength
         
         alpha = (subcoord-mincoord)/(maxcoord-mincoord)
         
         !Mass
         AuxSubArray(1,isub) = AuxSubArray(1,isub) + dvol
         !RHS
         AuxSubArray(2:1+ndofn,isub) = AuxSubArray(2:1+ndofn,isub)+dvol*(array(:,maxpoin)*alpha+array(:,minpoin)*(1-alpha))
      enddo
      
   enddo
   
   !We gather how many results root needs to receive
   call MPI_Gather( gminmaxsub, 2, MPI_INTEGER4, rGminmaxsub,2, MPI_INTEGER4, Mesh%MPIroot, Mesh%MPIcomm,ierr); 

   !Now we send the results root needs and root receives
   call MPI_ISEND(AuxSubArray(1:ndofn+1,gminmaxsub(1):gminmaxsub(2)),(ndofn+1)*(gminmaxsub(2)-gminmaxsub(1)+1), MPI_REAL8, Mesh%MPIroot, mtag1, Mesh%MPIcomm,irequest1, ierr) 
   !And receive
   if (Mesh%MPIrank == Mesh%MPIroot) then
      do irank = 0,Mesh%MPIsize-1
         call MPI_RECV(rbuffer,(ndofn+1)*(rGminmaxsub(2,irank+1)-rGminmaxsub(1,irank+1)+1),MPI_REAL8,irank,mtag1,Mesh%MPIcomm,status,ierr)
         AuxSubArray(:,rGminmaxsub(1,irank+1):rGminmaxsub(2,irank+1)) = AuxSubArray(:,rGminmaxsub(1,irank+1):rGminmaxsub(2,irank+1)) + &
                                                                  rbuffer(:,1:rGminmaxsub(2,irank+1)-rGminmaxsub(1,irank+1)+1)
      enddo
   endif
   !We all wait for the send to be complete
   call MPI_WAIT(irequest1, status, ierr)   
   
   if (Mesh%MPIrank == Mesh%MPIroot) then
      !We divide by the mass to have the mean value
      do isub = 1,nsubdivisions+1
         AuxSubArray(2:1+ndofn,isub) = AuxSubArray(2:1+ndofn,isub)/AuxSubArray(1,isub)
      enddo
      
      write(a%lun_outpu_graf,101) adjustl(trim(wopos)),istep,ttime, adjustl(trim(wopos)),istep,ttime
      do isub = 1,nsubdivisions+1
         write(a%lun_outpu_graf,*) (isub-1)*sublength, AuxSubArray(2:1+ndofn,isub)
      enddo
   endif
   
   !Memory deallocations
   call Memor%dealloc(ndofn+1,nsubdivisions+1,AuxSubArray,'AuxSubArray','AverageToOneDimension')
   call Mesh%ElementDeAlloc(e,Memor,'ForceClosedRule','vmass')
   if (Mesh%MPIrank == Mesh%MPIroot) then
      call Memor%dealloc(2,Mesh%MPIsize,rGminmaxsub,'rGminmaxsub','AverageToOneDimension')
      call Memor%dealloc(ndofn+1,nsubdivisions+1,rbuffer,'rbuffer','AverageToOneDimension')
   endif
   
   
   
   
   
   
!    call MPI_BCAST(AuxSubArray, (ndofn+1)*(nsubdivisions+1), MPI_REAL8, Mesh%MPIroot, a%MPIcomm, ierr)
!    !Now we put the Averaged value in each of the nodes
!    call Mesh%GetNpoin(npoin)
!    
!    do ipoin = 1,npoin
!       call Mesh%GetPointCoord(ipoin,coord)
!       subcoord = coord(AvgDime)
!       
!       isub = (subcoord-range(1,AvgDime))/(range(2,AvgDime)-range(1,AvgDime))*nsubdivisions+1
!       if (isub > nsubdivisions) isub = nsubdivisions
!       alpha = (subcoord-(range(1,AvgDime)+(isub-1)*sublength))/(sublength)
!       
!       AvgArray(1:ndofn,ipoin) = AuxSubArray(2:ndofn+1,isub+1)*alpha + AuxSubArray(2:ndofn+1,isub)*(1-alpha)
!    enddo
!    
!    !Finally we communicate the ghosts
!    call Mesh%ArrayCommunicator%GhostCommunicate(ndofn,AvgArray(1:ndofn,1:npoin))

   101 format ('# Graph: ''',a22, i6,f12.6,''' ',/,&
         '#'                                                    ,/,&
         '# X: ''Distance'' Y: ''',a22, i6,f12.6,''' '     ,/,&
         '#'                                                    ,/,&
         '# Units: {} {} ')
   
end subroutine

subroutine AverageToOneDimension1D(a,AvgDime,array,wopos,istep,ttime,Mesh,Memor)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_Postpr
   implicit none
   class(PostprFile) :: a
   integer(ip) :: AvgDime
   real(rp), target :: array(:)
   character(*), intent(in) :: wopos
   integer(ip) :: istep
   real(rp) :: ttime
   type(FemMesh)   :: Mesh
   type(MemoryMan) :: Memor
   
   interface 
      subroutine AverageToOneDimension(a,AvgDime,array,wopos,istep,ttime,Mesh,Memor)
         use typre
         use Mod_Memor
         use Mod_Mesh
         use Mod_Postpr
         implicit none
         class(PostprFile) :: a
         integer(ip) :: AvgDime
         real(rp) :: array(:,:)
         character(*), intent(in)  :: wopos
         integer(ip) :: istep
         real(rp) :: ttime
         type(FemMesh)   :: Mesh
         type(MemoryMan) :: Memor
      end subroutine
   
   end interface
   
   real(rp), pointer :: aux_array(:,:) => NULL()
   integer(ip) :: array_size
   
   array_size = size(array)
      
   aux_array(1:1,1:array_size) => array
   
   call AverageToOneDimension(a,AvgDime,aux_array,wopos,istep,ttime,Mesh,Memor)
   
end subroutine

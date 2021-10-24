module Mod_HangingNodes
   implicit none
   
contains

   subroutine ComputeMnodeForHangingNodesMeshes(a)
      use typre
      use Mod_Element
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ielem,nelem
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: pnode
      
      call a%GetNelem(nelem)
      
      a%mnode = maxval(a%nnode(1:a%nelty))
      
      call a%ElementAlloc(e,a%Memor,'DefaultRule','SetHangingMnode')
      
      elements : do ielem = 1,a%nelem
         !Load Element
         call a%ElementLoad(ielem,e)    
         
         call HangingComputePnodeAndFillIwaAndSeen(a,e,pnode)
         
         call HangingCleanIwaAndSeen(a,pnode)
         
         a%mnode = max(a%mnode,pnode)
      enddo elements
      
      call a%ElementDeAlloc(e,a%Memor,'DefaultRule','SetHangingMnode')
      
   end subroutine
   
   subroutine HangingComputePnodeAndFillIwaAndSeen(a,e,pnode)
      use typre
      use Mod_Element
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      class(FiniteElement) :: e
      integer(ip) :: pnode
      integer(ip) :: inode,ipoin,iphang,jhpos,jnode,jpoin
      
      pnode = 0
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         iphang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
         if (iphang == 0) then
            if (a%seenHanging(ipoin) == 0) then
               pnode = pnode + 1
               a%iwaHanging(pnode) = ipoin
               a%seenHanging(ipoin) = pnode
            endif
         else
            jhpos = a%pHangingList(ipoin)-1
            do jnode = 1,iphang
               jpoin = a%lHangingList(jhpos+jnode)
               if (a%seenHanging(jpoin) == 0) then
                  pnode = pnode + 1
                  a%iwaHanging(pnode) = jpoin
                  a%seenHanging(jpoin) = pnode
               endif
            enddo
         endif
      enddo
   end subroutine
   
   subroutine HangingCleanIwaAndSeen(a,pnode)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: pnode
      integer(ip) :: inode,ipoin
      
      !I need to clean seenHanging
      do inode = 1,pnode
         ipoin = a%iwaHanging(inode)
         a%seenHanging(ipoin) = 0
      enddo
   end subroutine
   
end module   

subroutine AssemblyToArrayHanging(a,e,ndofn,elrhs,array)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   class(FiniteElement) :: e
   integer(ip) :: ndofn
   real(rp) :: elrhs(ndofn,*),array(ndofn,*)
   
   integer(ip) :: inode,ipoin,phang,jspos,jpoin,jhang
   real(rp) :: jweight
   
   do inode = 1,e%pnode
      ipoin = e%lnods(inode)
      phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
      if (phang == 0) then
         array(:,ipoin) = array(:,ipoin) + elrhs(:,inode)
      else
         jspos = a%pHangingList(ipoin)-1
         do jhang = 1,phang
            jpoin = a%lHangingList(jspos+jhang)
            jweight = a%rHangingList(jspos+jhang)

            array(:,jpoin) = array(:,jpoin) + elrhs(:,inode)*jweight
         enddo
      endif
   enddo
   
end subroutine

subroutine PrepareHangingMatrices(a,e,elmat,elrhs)
   use typre
   use Mod_Element
   use Mod_Mesh
   use Mod_HangingNodes
   implicit none
   class(FemMesh), target :: a
   class(FiniteElement) :: e
   integer(ip) :: ndofn
   real(rp) ::  elmat(:,:,:,:),elrhs(:,:)
   integer(ip) :: HangCount
   integer(ip), pointer :: auxlnods(:) => NULL()
   integer(ip) :: inode,ipoin,iphang,ihpos
   integer(ip) :: jnode,jpoin,jphang,jhpos
   integer(ip) :: kpoin,knode,nodek,lpoin,lnode,nodel
   real(rp)    :: kweight,lweight
   integer(ip) :: pnode, lnods(1:e%mnode)
   
   real(rp) :: whelmat(size(elmat,1),size(elmat,2),size(elmat,3),size(elmat,4))
   real(rp) :: whelrhs(size(elrhs,1), size(elrhs,2))
   
   ndofn = size(elmat,1)
   
   !Check if there are no hanging nodes and assembly normally
   HangCount = 0
   do inode = 1,e%pnode
      ipoin = e%lnods(inode)
      iphang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
      if (iphang /= 0) HangCount = HangCount+1
   enddo

   !If there were no hanging nodes, assembly normally
   if (HangCount == 0) then
      !We can do a normal assembly
      return
   else
      !Modify the matrix and put it in whelmat, whelrhs
      !first we need to count the number of nodes I need to assembly to
      !but do not repeat nodes: I use seen and iwa

      !First of all we make sure that the element pointers and dimensions will be reset at the next element load
      e%ielty0 = -1

      !Compute Pnode and fill Iwa and Seen for later assembly
      call HangingComputePnodeAndFillIwaAndSeen(a,e,pnode)

      !Set whe data
      pnode = pnode
      lnods(1:pnode) = a%iwaHanging(1:pnode)

      !Matrices to zero
      whelrhs(:,1:pnode) = 0.0_rp
      whelmat(:,1:pnode,:,1:pnode) = 0.0_rp

      !begin the reassembly of the hanging nodes
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         !I only assembly localPoints rows
         !if (ipoin <= a%npoinLocal) then
         iphang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
         !no hanging nodes for node i
         if (iphang == 0) then
            nodek = a%seenHanging(ipoin)
            kweight = 1.0_rp

            !jnode loop for helmat
            call AssemblyLoop

            !node i has hanging nodes
         else
            ihpos = a%pHangingList(ipoin)-1
            do knode = 1,iphang
               kpoin = a%lHangingList(ihpos+knode)
               nodek = a%seenHanging(kpoin)
               kweight = a%rHangingList(ihpos+knode)

               !jnode loop for helmat
               call AssemblyLoop
            enddo
         endif

         !endif
      enddo

      call HangingCleanIwaAndSeen(a,pnode)

      !New nodes and number of nodes
      e%pnode = pnode

      !Swap whe and e
      !auxlnods => e%lnods
      e%lnods = lnods
      !whe%lnods => auxlnods

      elmat(1:ndofn,1:e%pnode,1:ndofn,1:e%pnode) = whelmat(1:ndofn,1:e%pnode,1:ndofn,1:e%pnode)
      elrhs(1:ndofn,1:e%pnode) = whelrhs(1:ndofn,1:e%pnode)

   endif
   
   e%ielty0 = -1

contains
      
   subroutine AssemblyLoop
      implicit none

      !First assembly rhs
      whelrhs(1:ndofn,nodek) = whelrhs(1:ndofn,nodek) + elrhs(1:ndofn,inode)*kweight

      !Now the AssemblyLoop
      do jnode = 1,e%pnode
         jpoin = e%lnods(jnode)
         jphang = a%pHangingList(jpoin+1)-a%pHangingList(jpoin)
         if (jphang == 0) then
            nodel = a%seenHanging(jpoin)
            lweight = 1.0_rp

            whelmat(1:ndofn,nodek,1:ndofn,nodel) = whelmat(1:ndofn,nodek,1:ndofn,nodel) + elmat(1:ndofn,inode,1:ndofn,jnode)*kweight*lweight
         else
            jhpos = a%pHangingList(jpoin)-1
            do lnode = 1,jphang
               lpoin = a%lHangingList(jhpos+lnode)
               nodel = a%seenHanging(lpoin)
               lweight = a%rHangingList(jhpos+lnode)

               whelmat(1:ndofn,nodek,1:ndofn,nodel) = whelmat(1:ndofn,nodek,1:ndofn,nodel) + elmat(1:ndofn,inode,1:ndofn,jnode)*kweight*lweight
            enddo
         endif
      enddo
   end subroutine   

end subroutine

subroutine BuildHangingGraph(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip), allocatable :: iwa(:),iwa2(:)
   type(i1p),allocatable :: i1wa(:)
   integer(ip), allocatable :: peaux(:)
   
   
   integer(ip) :: ielem,ispos,pnode,inode,ipoin,phang,ihpos,jpoin,jhang
   integer(ip) :: totalcoun, icoun,khang,kpoin,iaux2,elemi,iaux,jhpos,jnode
   integer(ip) :: myaux(4)
   
   !Create a%pelpo
   call a%Memor%alloc(a%npoin+1,a%pelpo,'pelpo','BuildHangingGraph')
   
   a%pelpo=0
   a%pelpo(1)=1
   do ielem=1,a%nelem
      !Only for last level elements
         call getlnods(ielem,pnode,ispos,a%pnods,a%lnods,a%nelty,a%nnode)
         do inode=1,pnode
            ipoin=a%lnods(ispos+inode)
            
            !We always add the point
            a%pelpo(ipoin+1)=a%pelpo(ipoin+1)+1
            
            !If it is a hanging node, we need to add the parent contributors
            phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
            if (phang /= 0) then
               ihpos = a%pHangingList(ipoin)-1
               do jhang = 1,phang
                  jpoin = a%lHangingList(ihpos+jhang)
                  a%pelpo(jpoin+1) = a%pelpo(jpoin+1)+1
               enddo   
            endif     
         end do
   end do
   
   !Compress a%pelpo
   do ipoin=1,a%npoin
      a%pelpo(ipoin+1)=a%pelpo(ipoin+1)+a%pelpo(ipoin)
   end do
   
   !Create a%lelpo
   call a%Memor%alloc(a%pelpo(a%npoin+1)-1,a%lelpo,'lelpo','BuildHangingGraph')
   call a%Memor%alloc(a%npoin,peaux,'peaux','BuildHangingGraph')

   peaux=-1
   do ielem=1,a%nelem
         call getlnods(ielem,pnode,ispos,a%pnods,a%lnods,a%nelty,a%nnode)
         do inode=1,pnode
            ipoin=a%lnods(ispos+inode)
            
            !We always add the point
            peaux(ipoin)=peaux(ipoin)+1
            a%lelpo(a%pelpo(ipoin)+peaux(ipoin))=ielem
            
            !If it is hanging, also add the recursive parent contribution
            phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
            if (phang /= 0) then
               ihpos = a%pHangingList(ipoin)-1
               do jhang = 1,phang
                  jpoin = a%lHangingList(ihpos+jhang)
                  peaux(jpoin)=peaux(jpoin)+1
                  a%lelpo(a%pelpo(jpoin)+peaux(jpoin))=ielem
               enddo   
            endif     
         end do
   end do
   call a%Memor%dealloc(a%npoin,peaux,'peaux','BuildHangingGraph')
   
   !Now we go for the graph
   call a%Memor%alloc(a%npoin,iwa,'iwa','BuildHangingGraph')
   call a%Memor%alloc(a%npoin,iwa2,'iwa2','BuildHangingGraph')
   call a%Memor%alloc(a%npoin,i1wa,'i1wa','BuildHangingGraph')
   
   !first to count and fill
   iwa = 0
   totalcoun = 0
   do ipoin = 1,a%npoin
      icoun = 0
      iwa(ipoin) = ipoin
      do elemi = a%pelpo(ipoin),a%pelpo(ipoin+1)-1
         ielem = a%lelpo(elemi)
         call getlnods(ielem,pnode,ispos,a%pnods,a%lnods,a%nelty,a%nnode)
         do jnode = 1,pnode
            jpoin = a%lnods(ispos+jnode)
            !Add all nodes, including hanging nodes 
            !even if they will not be assembled, hanging nodes need to be in the graph:
            !This is important in order to have a symmetric graph (required for calling PARMETIS in load rebalancing)
            if (iwa(jpoin) /= ipoin) then
               icoun = icoun+1
               totalcoun = totalcoun+1
               iwa(jpoin) = ipoin
               iwa2(icoun) = jpoin
            endif
            
            !If hanging, also add the parents
            phang = a%pHangingList(jpoin+1)-a%pHangingList(jpoin)
               
            !If it is a hanging node, I need to connect the point with all 
            !the final parents
            if (phang /= 0) then
               jhpos = a%pHangingList(jpoin)-1
               do khang = 1,phang
                  kpoin = a%lHangingList(jhpos+khang)
                     if (iwa(kpoin) /= ipoin) then
                     icoun = icoun+1
                     totalcoun = totalcoun+1
                     iwa(kpoin) = ipoin
                     iwa2(icoun) = kpoin
                  endif
               enddo   
            endif     
         enddo
      enddo
      allocate(i1wa(ipoin)%l(icoun))
      i1wa(ipoin)%l = iwa2(1:icoun)
   enddo
   
   !Compress the graph
   call a%Memor%alloc(a%npoin+1,a%ia,'ia','BuildHangingGraph')
   call a%Memor%alloc(totalcoun,a%ja,'ja','BuildHangingGraph')
   
   iaux2 = 1
   do ipoin = 1,a%npoin
      a%ia(ipoin) = iaux2
      if (associated(i1wa(ipoin)%l)) then
         iaux = size(i1wa(ipoin)%l)
         a%ja(iaux2:iaux2+iaux-1) = i1wa(ipoin)%l
         deallocate(i1wa(ipoin)%l)
         iaux2 = iaux2+iaux
      endif   
   enddo
   a%ia(a%npoin+1) = iaux2
   
   call a%Memor%dealloc(a%npoin,iwa,'iwa','BuildHangingGraph')
   call a%Memor%dealloc(a%npoin,iwa2,'iwa2','BuildHangingGraph')
   call a%Memor%dealloc(a%npoin,i1wa,'i1wa','BuildHangingGraph')
end subroutine

subroutine InterpolateHangingValues(a,ndofn,array)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   real(rp) :: array(ndofn,*)
   integer(ip) :: ipoin,phang,ihpos,jnode,jpoin
   real(rp) :: jweight
   
   do ipoin = 1,a%npoin
      phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
      if (phang /= 0) then
         array(:,ipoin) = 0.0_rp
         ihpos = a%pHangingList(ipoin)-1
         do jnode = 1,phang
            jpoin = a%lHangingList(ihpos+jnode)
            jweight = a%rHangingList(ihpos+jnode)
            array(:,ipoin) = array(:,ipoin) + jweight*array(:,jpoin)
         enddo
      endif
   enddo
end subroutine

subroutine InterpolateHangingValueSingle(a,ipoin,ndofn,array)
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   real(rp) :: array(ndofn,*)
   integer(ip) :: ipoin
   integer(ip) :: phang,ihpos,jnode,jpoin
   real(rp) :: jweight
   
   phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
   if (phang /= 0) then
      array(:,ipoin) = 0.0_rp
      ihpos = a%pHangingList(ipoin)-1
      do jnode = 1,phang
         jpoin = a%lHangingList(ihpos+jnode)
         jweight = a%rHangingList(ihpos+jnode)
         
         array(:,ipoin) = array(:,ipoin) + jweight*array(:,jpoin)
      enddo
   endif
   
end subroutine

subroutine CheckLogicalParents(a,ipoin,array,TrueOrFalse)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ipoin
   logical :: array(*)
   logical :: TrueOrFalse
   integer(ip) :: phang,jnode,jpoin,ihpos
   
   phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
   if (phang /= 0) then
      TrueOrFalse = .true.
      
      ihpos = a%pHangingList(ipoin)-1
      do jnode = 1,phang
         jpoin = a%lHangingList(ihpos+jnode)
         if (array(jpoin) .eqv. .false.) then
            TrueOrFalse = .false.
            return
         endif
      enddo
   else  
      TrueOrFalse = .false.
   endif
   
!    if (TrueOrFalse .eqv. .true.) then
!       write(*,*)
!    endif
  
end subroutine

subroutine AssemblyHangingNodesDiag(a,ndofn,LinearSystem,Memor)
   use typre
   use Mod_ParallelSystemInterface
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   class(ParallelSystemInterface) :: LinearSystem
   type(MemoryMan) :: Memor
   real(rp) :: elmat(ndofn,1,ndofn,1)
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: idofn,ipoin, phang
   real(rp)    :: diag(ndofn)
   
   !Manual allocation of the element e
   allocate(e)
   call Memor%palloc(1_ip,e%lnods,'e%lnods','AssemblyHangingNodesDiag')
   e%pnode = 1
   e%mnode = 1

   !Prepare elmat
   elmat = 0.0_rp
   do idofn = 1,ndofn
      elmat(idofn,1,idofn,1) = 1
   enddo

   !Assembly HangingNodes
   do ipoin = 1,a%npoin
      phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
      if (phang /= 0) then
         e%lnods(1) = ipoin
         call LinearSystem%AssemblyElmat(e,elmat)
      endif
   enddo   

   call Memor%pdealloc(1,e%lnods,'e%lnods','AssemblyHangingNodesDiag')
   deallocate(e)

end subroutine

subroutine HangingGetParents(a,ipoin,phang,lhang)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh), target :: a
   integer(ip) :: ipoin,phang
   integer(ip), pointer :: lhang(:)
   integer(ip) :: ispos
   
   phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
   ispos = a%pHangingList(ipoin)-1
   lhang => a%lhangingList(ispos+1:ispos+phang)
   
end subroutine

subroutine AssemblyHangingNodesDiagToZero(a,ndofn,LinearSystem,Memor)
   use typre
   use Mod_ParallelSystemInterface
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   class(ParallelSystemInterface) :: LinearSystem
   type(MemoryMan) :: Memor
   real(rp) :: elmat(ndofn,1,ndofn,1),elrhs(ndofn,1)
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: idofn,ipoin, phang
   real(rp)    :: diag(ndofn)
   
   !Manual allocation of the element e
   allocate(e)
   call Memor%palloc(1_ip,e%lnods,'e%lnods','AssemblyHangingNodesDiagToZero')
   e%pnode = 1
   e%mnode = 1

   !Prepare elmat
   elmat = 0.0_rp
   elrhs = 0.0_rp
   do idofn = 1,ndofn
      elmat(idofn,1,idofn,1) = -1
   enddo

   !Assembly HangingNodes
   do ipoin = 1,a%npoin
      phang = a%pHangingList(ipoin+1)-a%pHangingList(ipoin)
      if (phang /= 0) then
         e%lnods(1) = ipoin
         call LinearSystem%Assembly(e,elmat,elrhs)
      endif
   enddo   

   call Memor%pdealloc(1,e%lnods,'e%lnods','AssemblyHangingNodesDiagToZero')
   deallocate(e)

end subroutine

subroutine lev_GetRefinementCriteria(a,markel)
   use typre
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: markel(*)

   interface
      subroutine LayerRefinementCriteria(a,level,markel,Outlayers,GeneralRefinementLevels,InterfaceRefinementLevels,GlobalRefinementLevels,Offset)
         use typre
         use Mod_PhysicalProblem
         use Mod_Element
         use Mod_LinkedList
         implicit none
         class(PhysicalProblem) :: a
         real(rp) :: level(:)
         integer(ip) :: markel(*)

         integer(ip) :: OutLayers, GeneralRefinementLevels,InterfaceRefinementLevels,GlobalRefinementLevels
         real(rp) :: Offset
      end subroutine

   end interface



   call a%Timer%Refine%Tic

   call LayerRefinementCriteria(a,a%level(:,1),markel,a%OutLayers,a%GeneralRefinementLevels,a%InterfaceRefinementLevels,0,0.0_rp)

   call a%Timer%Refine%Toc

end subroutine


subroutine LayerRefinementCriteria(a,level,markel,Outlayers,InRefinementLevels,InterfaceRefinementLevels,GlobalRefinementLevels,Offset)
   use typre
   use Mod_PhysicalProblem
   use Mod_Element
   use Mod_LinkedList
   implicit none
   class(PhysicalProblem) :: a
   real(rp) :: level(:)
   integer(ip) :: markel(*)

   integer(ip) :: OutLayers, InRefinementLevels,InterfaceRefinementLevels,GlobalRefinementLevels
   real(rp)    :: Offset

   integer(ip) :: nelem,npoin
   integer(ip), allocatable :: RefinementLevel(:)
   class(FiniteElement), pointer :: e => NULL()

   integer(ip) :: ielem, inode,ipoin,maxtrue

   type(LinkedListHeader) :: HeaderPointsToDo
   type(LinkedListStorage) :: StoragePointsToDo

   integer(ip), allocatable :: auxlevel1(:), auxlevel2(:)
   integer(ip) :: nrefinementLayers,ilayer



   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin)

   call a%Memor%alloc(npoin,RefinementLevel,'Refinementlevel','lev_GetRefinementCriteria')
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_GetRefinementCriteria')


   allocate(auxlevel1(npoin))
   allocate(auxlevel2(npoin))
   auxlevel1 = 0
   auxlevel2 = 0

   !Initialize
   do ipoin = 1,npoin
      if (level(ipoin) >= Offset) then
         auxlevel1(ipoin) = 1
      else
         auxlevel1(ipoin) = 0
      endif
   enddo

    !call a%FilePostpr%postpr(auxlevel1,'auxlevel1pre0',a%istep,a%ctime,a%Mesh)

   !Mark the interface layer
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)

      if (maxval(level(e%lnods(1:e%pnode))) >= Offset .and. minval(level(e%lnods(1:e%pnode))) < Offset) then
         auxlevel1(e%lnods(1:e%pnode)) = 2
      endif
   enddo
   call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,auxlevel1)

   !call a%FilePostpr%postpr(auxlevel1,'auxlevel1pre',a%istep,a%ctime,a%Mesh)


   !Number of layers to refine
   nrefinementLayers = OutLayers
   do ilayer = 1,nrefinementLayers
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         maxtrue = maxval(auxlevel1(e%lnods(1:e%pnode)))
         if (maxtrue  == 2) then
            auxlevel2(e%lnods(1:e%pnode)) = maxtrue
         elseif (maxtrue == 1) then
            do inode = 1,e%pnode
               if (auxlevel2(e%lnods(inode)) <= 1) then
                  auxlevel2(e%lnods(inode)) = 1
               endif
            enddo
         endif
      enddo
      call move_alloc(auxlevel2,auxlevel1)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,auxlevel1)
      allocate(auxlevel2(npoin))
      auxlevel2 = 0
   enddo

    !call a%FilePostpr%postpr(auxlevel1,'auxlevel1',a%istep,a%ctime,a%Mesh)

   !We Refine two levels in the fluid 1 subdomain, we unrefine in the fluid 2 subdomain
   call a%Refiner%GetLevel(RefinementLevel)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      if (maxval(RefinementLevel(e%lnods(1:e%pnode))) < GlobalRefinementLevels) then
         markel(ielem) = 1
      elseif (maxval(auxlevel1(e%lnods(1:e%pnode))) == 1 .and. maxval(RefinementLevel(e%lnods(1:e%pnode))) < InRefinementLevels + GlobalRefinementLevels) then
         markel(ielem) = 1
      elseif (maxval(auxlevel1(e%lnods(1:e%pnode))) == 1 .and. maxval(RefinementLevel(e%lnods(1:e%pnode))) > InRefinementLevels + GlobalRefinementLevels) then
         markel(ielem) = -1
      elseif (maxval(auxlevel1(e%lnods(1:e%pnode))) == 2 .and. maxval(RefinementLevel(e%lnods(1:e%pnode))) < InRefinementLevels + InterfaceRefinementLevels + GlobalRefinementLevels) then
         markel(ielem) = 1
      elseif (maxval(auxlevel1(e%lnods(1:e%pnode))) == 0 .and. maxval(RefinementLevel(e%lnods(1:e%pnode))) > GlobalRefinementLevels) then
         markel(ielem) = -1
      endif
   enddo

   deallocate(auxlevel1)
   deallocate(auxlevel2)


   call a%Memor%dealloc(npoin,RefinementLevel,'Refinementlevel','lev_GetRefinementCriteria')
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','lev_GetRefinementCriteria')











end subroutine

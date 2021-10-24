subroutine lev_CutElementsAndListByLayers(a)
!DESCRIPTION
!   This subroutine get the element status: cut=1 and non cut=0
!-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   use Mod_CutMesh
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem), target :: a

   interface
      subroutine lev_LayerOut(a, NelemWithAtLeastOneauxPointDone,ElementDone,PointDone,ElementListByLayers)
         use typre
         use Mod_LevelSet
         use Mod_Element
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: NelemWithAtLeastOneauxPointDone
         logical :: ElementDone(*),PointDone(*)
         integer(ip) :: ElementListByLayers(*)
      end subroutine

   end interface

   class(FiniteElement), pointer :: e => NULL()

   integer(ip) :: icomp,ipoin,jpoin,npoin,ielem,nelem,nnode,inode,elemi,pelpo,helem,jelem
   integer(ip) :: elemh,hpoin,elemj

   !auxiliar variables
   integer(ip) :: auxlelemcut,auxnnode1,auxnnode2,kelem,existin,existout
   logical, allocatable     :: ElementDone(:),PointDone(:)
   !surface subroutine
   integer(ip) :: ndime,auxAdd,fnode,fpoin,gpoin,gnode
   real(rp)    :: dsurf,hclen,rdinv,cutcriteria
   real(rp), allocatable :: ellev(:,:)

   integer(ip) :: NelemWithAtLeastOneauxPointDone, ElementStatus
   real(rp), pointer :: level(:)=> Null()

   level => a%level(:,1)


   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)

   rdinv=1.0_rp/real(ndime)

   !Definition of arrays needed by the level set

   !allocate
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_CutElementsAndListByLayers')
   call a%Memor%alloc(e%mnode,1,ellev,'ellev','lev_CutElementsAndListByLayers')
   call a%Memor%alloc(nelem,ElementDone,'ElementDone','lev_CutElementsAndListByLayers')
   call a%Memor%alloc(npoin,PointDone,'PointDone','lev_CutElementsAndListByLayers')

   !Initialization
   ElementDone = .false.
   PointDone=.false.

   call a%CutMesh%ComputeIntersectionPoints(level)

   !Fluid definition
   !List of point cuts by the intreface and point status to define the sign of then
   call a%CutMesh%ComputeElementType(level)
   call a%CutMesh%ComputePoinType(level)
   !Only testing
   call a%CutMesh%ComputeSubelements(ndime)
   call a%CutMesh%ComputeEnrichBubble(ndime)


   kelem=0_ip

   !to check if the element has two fluids
   ielem=1
   do ielem=1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)

      call a%CutMesh%GetElementType(ielem,ElementStatus)

      if (ElementStatus == 0.0_rp) then
         !obtaining surface and intersection points
         dsurf=0.0_rp
         call a%CutMesh%GetSurfaceIntersection(ielem,e,dsurf)


         !Derivative (and detjm) at the center of gravity
         call e%elmdcg

         !Element length
         hclen=(e%weicg*e%detjm)**rdinv
         hclen=hclen/100000_rp

         if(dsurf>0.0_rp)then

            !List of cut elements initially not full
            kelem=kelem + 1

            a%ElementListByLayers(kelem) = ielem
            ElementDone(ielem) = .true.

            !Load Element
            call a%Mesh%ElementLoad(ielem,e)

            do gnode=1,e%pnode
               gpoin=e%lnods(gnode)
               PointDone(gpoin) = .true.
            end do
          end if

      end if
   end do

   NelemWithAtLeastOneauxPointDone = kelem

   call lev_LayerOut(a,NelemWithAtLeastOneauxPointDone,ElementDone,PointDone,a%ElementListByLayers)


   !deallocate
   call a%Memor%dealloc(e%mnode,1,ellev,'ellev','lev_CutElementsAndListByLayers')
   call a%Memor%dealloc(nelem,ElementDone,'ElementDone','lev_CutElementsAndListByLayers')

   call a%Memor%dealloc(npoin,PointDone,'PointDone','lev_CutElementsAndListByLayers')
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_CutElementsAndListByLayers')

end subroutine



subroutine lev_LayerOut(a, NelemWithAtLeastOneauxPointDone,ElementDone,PointDone,ElementListByLayers)
   use typre
   use Mod_LevelSet
   use Mod_Element
   use Mod_HangingNodes
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: NelemWithAtLeastOneauxPointDone
   logical :: ElementDone(*),PointDone(*)
   integer(ip) :: ElementListByLayers(*)

   class(FiniteElement), pointer :: e => NULL()
   logical, allocatable :: auxPointDone(:)
   real(rp), allocatable :: ellev(:,:)
   integer(ip), pointer     :: lelpo(:) => NULL()

   integer(ip) :: icomp,ipoin,jpoin,npoin,ielem,nelem,nnode,inode,elemi,pelpo,helem,jelem
   integer(ip) :: elemh,hpoin,elemj

   integer(ip) :: ndime,auxAdd,fnode,fpoin,gpoin,gnode,npoinInterface,kelem
   integer(ip), allocatable :: auxelementlist(:)
   integer(ip) :: kfl_Hanging

   interface
      subroutine CheckLogicalParents(a,ipoin,array,TrueOrFalse)
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh) :: a
      integer(ip) :: ipoin
      logical :: array(*)
      logical :: TrueOrFalse
      end subroutine

   end interface


   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)

   allocate(auxelementlist(nelem))
   auxelementlist= 0

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_CutElementsAndListByLayers')
   call a%Memor%alloc(e%mnode,1,ellev,'ellev','lev_CutElementsAndListByLayers')


   call a%Memor%alloc(npoin,auxPointDone,'auxPointDone','lev_CutElementsAndListByLayers')
   call a%Mesh%GetHanging(kfl_Hanging)

   !To complete the list of element layers
   !Do not do it if no element is cut
   if (NelemWithAtLeastOneauxPointDone > 0) then
      kelem = NelemWithAtLeastOneauxPointDone
      auxPointDone = .false.
      do elemh=1,nelem
         helem = ElementListByLayers(elemh)
         if (helem < 1 .or. helem > nelem) then
            return
            !We will complete it in the following passes

!             do ielem = 1,nelem
!                if (auxelementlist(ielem) == 0) then
!                   write(*,*) 'asdfasdf'
!                endif
!             enddo
!             call runend('Error in lev_layerOut')
         endif
         auxelementlist(helem) = 1
         !Load Element
         call a%Mesh%ElementLoad(helem,e)

         do inode = 1, e%pnode
            hpoin=e%lnods(inode)

            if(auxPointDone(hpoin) .eqv. .false.)then
               auxPointDone(hpoin)=.true.

               call a%Mesh%GetLelpo(hpoin,pelpo,lelpo)

               do elemj=1,pelpo
               jelem = lelpo(elemj)

                  auxAdd=0
                  call a%Mesh%ElementLoad(jelem,e)
                  do fnode =1,e%pnode
                     fpoin =e%lnods(fnode)

                     !Check for hanging nodes
                     if (kfl_Hanging == 1) then
                        if (PointDone(fpoin) .eqv. .false.) then
                           call CheckLogicalParents(a%Mesh,fpoin,PointDone,PointDone(fpoin))
                        endif
                     endif

                     if(PointDone(fpoin) .eqv. .true.)then
                        auxAdd = auxAdd +1
                     end if
                  end do

                  ! in linear elements we add elements when only one node is not calculated
                  if((ElementDone(jelem) .eqv. .false.) .and. (auxAdd>=ndime))then

                     ElementDone(jelem) = .true.
                     kelem=kelem+1
                     ElementListByLayers(kelem)=jelem

                     do gnode=1,e%pnode
                        gpoin=e%lnods(gnode)
                        PointDone(gpoin) = .true.
                     end do

                  end if
               end do
            end if
         end do
      end do
   endif

   call a%Memor%dealloc(npoin,auxPointDone,'auxPointDone','lev_CutElementsAndListByLayers')
   call a%Memor%dealloc(e%mnode,1,ellev,'ellev','lev_CutElementsAndListByLayers')

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_CutElementsAndListByLayers')

end subroutine



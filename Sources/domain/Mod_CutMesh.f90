module Mod_CutMesh
!-----------------------------------------------------
!
! This module contains the information of the cut elements for the levelset problem
!
!-----------------------------------------------------
   use typre
   use def_parame
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   private
   public CutMesh

   type CutElement
      integer(ip)           :: elemStatus = -20,ialone
      integer(ip)           :: ngauss_minus,ngauss_plus
      real(rp), allocatable :: weigp(:)
      integer(ip)           :: ninters,nsube
      real(rp), allocatable :: LocalIntersectionPoints(:,:)
      real(rp), allocatable :: xloc(:,:)
      real(rp), allocatable :: kEnrichBub(:)
   end type

   type CutMesh
      class(FemMesh), pointer :: Mesh => NULL()
      type(MemoryMan) , pointer :: Memor => NULL()

      integer(ip), allocatable :: poinStatus(:)
      integer(ip)              :: ncutelem,npointInterface
      type(CutElement), allocatable :: CutElementPoints(:)


contains
      procedure :: allocCutMesh
      procedure :: ComputeIntersectionPoints
      procedure :: GetLocalIntersectionPoints
      procedure :: GetRealIntersectionPoints
      procedure :: GetSurfaceIntersection
      procedure :: GetNinters
      procedure :: ComputeElementType
      procedure :: ComputePoinType
      procedure :: GetElementType
      procedure :: GetPointType
      procedure :: GetNCutElements
      procedure :: GetNgaussSide
      procedure :: GetInterfacePoints
      procedure :: ComputeSubelements
      procedure :: ComputeEnrichBubble
      procedure :: GetEnrichBuble
      procedure :: GetWeigpCut
      procedure :: GetXlocCut
      procedure :: GetIalone
      procedure :: deallocCutMesh
      procedure :: deallocCutElement
      procedure :: GetExternalNormal

   end type

contains

   subroutine allocCutMesh(a,Memor,Mesh)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh)     :: a
      class(FemMesh),target :: Mesh
      type(MemoryMan), target   :: Memor

      integer(ip)           :: npoin,ndime,nelem

      a%Mesh => Mesh
      a%Memor => Memor

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)

      call a%Memor%alloc(npoin,a%poinStatus,'poinStatus','allocCutMesh')

      !CutElement type
      allocate(a%CutElementPoints(nelem))
      call a%Memor%allocObj(0,'CutElementPoints','allocCutMesh',(nelem)*1_ip)

   end subroutine



   subroutine ComputeIntersectionPoints(a,array)
      use Mod_Memor
      implicit none
      class(CutMesh) :: a
      class(FiniteElement), pointer :: e => NULL()

      real(rp), intent(in)     :: array(*)
      integer(ip)              :: inode,ielem,jnode,idime,ipoin,ndime,auxinter,nelem,jpoin
      real(rp)                 :: lratio(4),tolerance,chale(1) !4 maximum number of intersection points in 3d
      real(rp), pointer        :: coord1(:) => NULL(), coord2(:) => NULL()
      integer(ip)              :: countmemor,auxdim,auxcount
      integer(ip)              :: side_minus1_blending,side_plus1_blending,side_minus1_level,side_plus1_level,icount
      real(rp)                 :: LocalIntersectionPoint(3), IntersectionPoint(3)

      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ComputeIntersectionPoints')


      !Number of intersection points 2 in 2d and 3 o 4 in 3d taking the size 4 always
      if(ndime==2)then
         auxdim=2
      elseif(ndime==3)then
         auxdim=4
      end if

      countmemor=0

      do ielem=1,nelem

         !Load Element
         call a%Mesh%ElementLoad(ielem,e)

         auxinter=0
         auxcount=0

         if(ndime==2)then
            if(e%pnode==3)then
               icount = count(array(e%lnods(1:e%pnode))>0.0_rp)
               if (icount > 0 .and. icount < e%pnode) then

                  do inode = 1,e%pnode
                     do jnode = 1+inode,e%pnode

                     ipoin=e%lnods(inode)
                     jpoin=e%lnods(jnode)

                        if((((array(ipoin) <= 0.0_rp) .and. (array(jpoin) > 0.0_rp))).or. &
                        (((array(ipoin) > 0.0_rp) .and. (array(jpoin) <= 0.0_rp))))then


                           if(auxcount==0)then
                              auxcount=auxcount+1
                              allocate(a%CutElementPoints(ielem)%LocalIntersectionPoints(ndime,auxdim))
                              countmemor=countmemor+1
                           end if

                           auxinter = auxinter + 1

                           !ratio to define the intersection points
                           lratio(auxinter) = array(ipoin)/(array(ipoin) - array(jpoin))
                           !We need to take it from elcod just in case we are using ALE
                           coord1 => e%elcod(:,inode)
                           coord2 => e%elcod(:,jnode)
                           !call a%Mesh%GetPointCoord(ipoin,coord1)
                           !call a%Mesh%GetPointCoord(jpoin,coord2)

                           IntersectionPoint(1:e%ndime) = coord1(:) + (coord2(:) - coord1(:))*lratio(auxinter)
                           call e%isoparinv(IntersectionPoint,a%CutElementPoints(ielem)%LocalIntersectionPoints(:,auxinter))
                        endif

                     end do
                  end do
               endif


            else
               call runend('lev_Surface : level set only for linear elements')
            end if

         elseif(ndime==3)then

            if(e%pnode==4)then

               auxcount=0
               do inode = 1,e%pnode
                  do jnode = 1+inode,e%pnode

                  ipoin=e%lnods(inode)
                  jpoin=e%lnods(jnode)

                     if((((array(ipoin) <= 0.0_rp) .and. (array(jpoin) > 0.0_rp))).or. &
                     (((array(ipoin) > 0.0_rp) .and. (array(jpoin) <= 0.0_rp))))then

                        auxinter = auxinter + 1

                        if(auxcount==0)then
                           auxcount=auxcount+1
                           allocate(a%CutElementPoints(ielem)%LocalIntersectionPoints(ndime,auxdim))
                           countmemor=countmemor+1
                        end if

                        !ratio to define the intersection points
                        lratio(auxinter) = array(ipoin)/(array(ipoin) - array(jpoin))
                        !We need to take it from elcod just in case we are using ALE
                        coord1 => e%elcod(:,inode)
                        coord2 => e%elcod(:,jnode)
                        !call a%Mesh%GetPointCoord(ipoin,coord1)
                        !call a%Mesh%GetPointCoord(jpoin,coord2)

                        IntersectionPoint(1:e%ndime) = coord1(:) + (coord2(:) - coord1(:))*lratio(auxinter)
                        call e%isoparinv(IntersectionPoint,a%CutElementPoints(ielem)%LocalIntersectionPoints(:,auxinter))

                     end if

                  end do
               end do
            else
               call runend('lev_Surface: level set only for linear tetrahedra in 3d ')
            end if
         end if

         a%CutElementPoints(ielem)%ninters=auxinter

      end do

      !CutElement Variables
      call a%Memor%allocObj(0,'ipoints','allocCutMesh',(countmemor*ndime*auxdim)*rp)

      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputeIntersectionPoints')

   end subroutine

   subroutine GetLocalIntersectionPoints(a,ielem,points)
      implicit none
      class(CutMesh), target :: a

      integer(ip), intent(in) :: ielem
      real(rp), pointer   :: points(:,:)

      points=>a%CutElementPoints(ielem)%LocalIntersectionPoints(:,:)

   end subroutine

   subroutine GetRealIntersectionPoints(a,ielem,e,points)
      implicit none
      class(CutMesh), target :: a

      integer(ip), intent(in) :: ielem
      class(FiniteElement) :: e
      real(rp)   :: points(:,:)

      integer(ip) :: iipoin

      do iipoin = 1,a%CutElementPoints(ielem)%ninters
         call e%isopar(e%ndime,a%CutElementPoints(ielem)%LocalIntersectionPoints(:,iipoin),e%elcod,points(:,iipoin))
      enddo

   end subroutine

   subroutine GetSurfaceIntersection(a,ielem,e,dsurf)
      implicit none
      class(CutMesh) :: a
      class(FiniteElement) :: e

      integer(ip), intent(in) :: ielem
      real(rp), intent(out)   :: dsurf
      real(rp)                :: t1(e%ndime),t2(e%ndime),cprod(e%ndime),rdsurf
      integer(ip)             :: idime

      integer(ip) :: iipoin
      real(rp) :: IntersectionPoints(e%ndime,4)

      !Compute the intersection points in real world coordinates
      call a%GetRealIntersectionPoints(ielem,e,IntersectionPoints)

      if(e%ndime==2)then
         if(a%CutElementPoints(ielem)%ninters == 2)then !linear triangle or Cuadrilateral
            dsurf=0.0_rp
            do idime = 1,e%ndime
               dsurf = dsurf + (IntersectionPoints(idime,1) - IntersectionPoints(idime,2))**2.0_rp
            end do
            dsurf=sqrt(dsurf)
         else
            write(*,*) 'Element number: ', ielem, ', NumberOfIntersectionPoints: ', a%CutElementPoints(ielem)%ninters
            call runend('lev_Surface: wrong number of intersection points')
         endif
      elseif(e%ndime==3)then

         if(a%CutElementPoints(ielem)%ninters == 3)then !linear tetrahedra

            dsurf=0.0_rp

            t1= IntersectionPoints(:,2)-IntersectionPoints(:,1)
            t2= IntersectionPoints(:,3)-IntersectionPoints(:,1)

            call vecpro(t1,t2,cprod,3)
            call vecnor(cprod,e%ndime,dsurf,2)
            dsurf = dsurf*0.5_rp

         elseif(a%CutElementPoints(ielem)%ninters==4)then

            rdsurf=0.0_rp
            dsurf=0.0_rp

            t1= IntersectionPoints(:,1)-IntersectionPoints(:,2)
            t2= IntersectionPoints(:,1)-IntersectionPoints(:,3)

            call vecpro(t1,t2,cprod,3)
            call vecnor(cprod,3,rdsurf,2)

            dsurf= dsurf + abs(rdsurf)*0.5_rp

            rdsurf=0.0_rp

            t1= IntersectionPoints(:,1)-IntersectionPoints(:,4)
            t2= IntersectionPoints(:,1)-IntersectionPoints(:,3)

            call vecpro(t1,t2,cprod,3)
            call vecnor(cprod,3,rdsurf,2)

            dsurf= dsurf + abs(rdsurf)*0.5_rp

         else
            call runend('lev_Surface : wrong number of intersection')
         end if

      endif

   end subroutine

   subroutine GetNinters(a,ielem,nipoin)
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ielem
      integer(ip), intent(out) :: nipoin

      nipoin=0_ip
      nipoin=a%CutElementPoints(ielem)%ninters

   end subroutine

   subroutine ComputeElementType(a,array)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh) :: a
      real(rp), intent(in)    :: array(*)


      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ielem,nelem,existin,existout,inode,ipoin,ndime
      integer(ip) :: auxngaus_minus,auxngaus_plus
      real(rp)    :: dsurf
      !CutElement
      integer(ip) :: countmemor,auxdim,auxdim1,idime
      real(rp)    ::chale(1)

      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ComputeElementType')

      !Initialization
      a%ncutelem=0
      countmemor=0

      auxdim =(ndime-1)*(ndime-1)+2
      auxdim1=(ndime+1)*auxdim

      !to check if the element has two fluids
      ielem=1
      do ielem=1,nelem

         !Load Element
         call a%Mesh%ElementLoad(ielem,e)

         existin=0_ip
         existout=0_ip

         auxngaus_plus=0
         auxngaus_minus=0

         a%CutElementPoints(ielem)%ngauss_plus =0
         a%CutElementPoints(ielem)%ngauss_minus=0

         do inode = 1,e%pnode
            !fluid definition
            ipoin=e%lnods(inode)
            if(array(ipoin)>0.0_rp)then
               existin = 1_ip
               auxngaus_plus=auxngaus_plus+1
            elseif(array(ipoin)<=0.0_rp)then
               existout= 1_ip
               auxngaus_minus=auxngaus_minus+1
            end if
         end do

         if(existin == 1_ip .and. existout == 0_ip)then
            a%CutElementPoints(ielem)%elemStatus = 1_ip
         elseif(existin == 0_ip .and. existout == 1_ip)then
            a%CutElementPoints(ielem)%elemStatus = -1_ip

         !Define the element type
         !Point and element status definition between fluid 1 and 2
         elseif(existin == 1_ip .and. existout == 1_ip)then

            a%CutElementPoints(ielem)%elemStatus=0_ip

            !allocate the CutElement type variables
            allocate(a%CutElementPoints(ielem)%weigp(auxdim1))
            countmemor = 1 + countmemor

            !number of Gauss points definition
            if(ndime==2)then !only for linear triangles
               if(auxngaus_plus==1)then
                  a%CutElementPoints(ielem)%ngauss_plus=3
               elseif(auxngaus_plus==2)then
                  a%CutElementPoints(ielem)%ngauss_plus=6
               end if

               if(auxngaus_minus==1)then
                  a%CutElementPoints(ielem)%ngauss_minus=3
               elseif(auxngaus_minus==2)then
                  a%CutElementPoints(ielem)%ngauss_minus=6
               endif
            elseif(ndime==3)then !only for linear tetrahedra
               if(a%CutElementPoints(ielem)%ninters==3)then
                  if(auxngaus_plus==1)then
                     a%CutElementPoints(ielem)%ngauss_plus=4
                  elseif(auxngaus_plus==2)then
                     a%CutElementPoints(ielem)%ngauss_plus=8
                  elseif(auxngaus_plus==3)then
                     a%CutElementPoints(ielem)%ngauss_plus=12
                  end if

                  if(auxngaus_minus==1)then
                     a%CutElementPoints(ielem)%ngauss_minus=4
                  elseif(auxngaus_minus==2)then
                     a%CutElementPoints(ielem)%ngauss_minus=8
                  elseif(auxngaus_minus==3)then
                     a%CutElementPoints(ielem)%ngauss_minus=12
                  endif
               elseif(a%CutElementPoints(ielem)%ninters==4)then
                  if(auxngaus_minus==2 .and. auxngaus_plus==2)then
                     a%CutElementPoints(ielem)%ngauss_plus=12
                     a%CutElementPoints(ielem)%ngauss_minus=12
                  end if
               end if
            end if

            !obtaining surface and intersection points
            dsurf=0.0_rp
            call GetSurfaceIntersection(a,ielem,e,dsurf)

            if(dsurf>0.0_rp)then
               !number of cut elements
               a%ncutelem =1_ip + a%ncutelem
            end if
         else
            call runend('This should not be happening')

         end if
      end do

      !CutElement Variables
      call a%Memor%allocObj(0,'weigp','ComputeElementType',(countmemor*auxdim1)*rp)

      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputeElementType')

!      write(*,*) 'Number of cut element Mod_CutMesh'
!      write(*,*) a%ncutelem
   end subroutine

   subroutine GetNCutElements(a,ncutelem)
      class(CutMesh) :: a
      integer(ip), intent(out) :: ncutelem

      ncutelem = a%ncutelem

   end subroutine

   subroutine GetNgaussSide(a,ielem,ngauss_plus,ngauss_minus)
      class(CutMesh) :: a

      integer(ip), intent(in) :: ielem
      integer(ip), intent(out)   :: ngauss_plus,ngauss_minus

      ngauss_minus=0
      ngauss_plus =0
      ngauss_minus=a%CutElementPoints(ielem)%ngauss_minus
      ngauss_plus =a%CutElementPoints(ielem)%ngauss_plus

   end subroutine


   subroutine ComputePoinType(a,array)
      use typre
      implicit none
      class(CutMesh) :: a
      class(FiniteElement), pointer :: e => NULL()

      real(rp), intent(in)    :: array(*)
      integer(ip) :: ipoin,jpoin,npoin,auxpoininterface
      integer(ip) :: ielem,nelem,inode
      integer(ip) :: idime
      real(rp)    :: chale(1)

      integer(ip), allocatable :: pointList(:)



      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ComputePoinType')
      call a%Memor%alloc(npoin,pointList,'pointList','ComputePointType')


      auxpoininterface=0
      jpoin=0_ip
      !Initialization of the pointList
      pointList=0

      !poinStatus in cut elements
      ielem=1
      do ielem=1,nelem

         !Load Element
         call a%Mesh%ElementLoad(ielem,e)


         if(a%CutElementPoints(ielem)%elemStatus==0_ip)then !Interface
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               if(pointList(ipoin)==0_ip)then
                  pointList(ipoin)=1_ip
                  if(array(ipoin)<=0.0_rp)then
                     a%poinStatus(ipoin)=-1_ip
                  elseif(array(ipoin)>0.0_rp)then
                     a%poinStatus(ipoin)=1_ip
                  end if
               end if
            end do
         end if
      end do

      !poinStatus in no-Cutelements
      ielem=1
      do ielem=1,nelem
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         if(a%CutElementPoints(ielem)%elemStatus==1_ip .or. a%CutElementPoints(ielem)%elemStatus==-1_ip)then !Interface
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               if(pointList(ipoin)==0_ip)then
                  pointList(ipoin)=1_ip
                  if(array(ipoin)<=0.0_rp)then
                     a%poinStatus(ipoin)=-2_ip
                  elseif(array(ipoin)>0.0_rp)then
                     a%poinStatus(ipoin)=2_ip
                  end if
               end if
            end do
         end if
      end do

     call a%Memor%dealloc(npoin,pointList,'pointList','ComputePointType')
     call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputePoinType')

   end subroutine

   subroutine GetInterfacePoints(a,npointInterface)
      use typre
      implicit none
      class(CutMesh) :: a
      integer(ip), intent(out) :: npointInterface

      npointInterface=a%npointInterface

   end subroutine

   subroutine GetElementType(a,ielem,eStatus)
      use typre
      implicit none
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ielem
      integer(ip), intent(out) :: eStatus

      eStatus =a%CutElementPoints(ielem)%elemStatus

   end subroutine

   subroutine GetPointType(a,ipoin,pStatus)
      use typre
      implicit none
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ipoin
      integer(ip), intent(out) :: pStatus

      pStatus=a%poinStatus(ipoin)

   end subroutine

   subroutine ComputeSubelements(a,ndime)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh)   :: a
      class(FiniteElement), pointer :: e => NULL()


      integer(ip), intent(in)    :: ndime
      integer(ip)              :: inode,auxdim
      integer(ip)              :: existin,existout,auxnode,numnod(3),ielem,nelem
      integer(ip)              :: ialone,ipoin,isube,igaus,auxgaus,knode,auxpnode
      integer(ip)              :: icoun,icount(3),auxsubelcod(2),iedge,nodeedge(2,2),itype,auxsubelement
      integer(ip)              :: auxdim1,auxdim2
      real(rp), allocatable    :: subelcod(:,:,:),gpcod(:),iweigp(:)
      integer(ip), allocatable :: edge(:)
      !to compute iweigp
      integer(ip)              :: sublinea
      real(rp)                 :: subdetjm,subcartd(ndime,ndime+1)
      real(rp)                 :: subxjacm(ndime,ndime)
      real(rp)                 :: subxjaci(ndime,ndime)
      real(rp)                 :: dvol,subdvol
      integer(ip)              :: auxcount,countmemor

      integer(ip) :: iipoin
      real(rp), allocatable :: IntersectionPoints(:,:)

      call a%Mesh%GetNelem(nelem)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ComputeSubelements')

      auxdim=ndime+1 !number of nodes per subelement (always is a triangle or tetrahedra even in multilineal elements)
      auxdim1=((ndime-1)*(ndime-1)+2) !number of subelements
      auxdim2=auxdim*((ndime-1)*(ndime-1)+2) !total number of gauss point per cut element

      call a%Memor%alloc(ndime,auxdim,auxdim1,subelcod,'subelcod','ComputeSubelements')
      call a%Memor%alloc(auxdim1,edge,'edge','ComputeSubelements')
      call a%Memor%alloc(ndime,gpcod,'gpcod','ComputeSubelements')
      call a%Memor%alloc(auxdim2,iweigp,'iweigp','ComputeSubelements')
      call a%Memor%alloc(e%ndime,4_ip,IntersectionPoints,'IntersectionPoints','ComputeSubelements')

      a%CutElementPoints%nsube=0
      countmemor=0

      do ielem=1,nelem

         !initializations
         iweigp(:) = 0

         if(a%CutElementPoints(ielem)%elemStatus == 0)then !cut lement

            call a%Mesh%ElementLoad(ielem,e)

            allocate(a%CutElementPoints(ielem)%xloc(ndime,auxdim2))
            countmemor=countmemor+1

            !Compute the intersection points in real world coordinates
            call a%GetRealIntersectionPoints(ielem,e,IntersectionPoints)

            if(ndime==2)then !only for triangles
               if(e%pnode==3)then

                  a%CutElementPoints(ielem)%nsube=3

                  edge=0
                  !face localization
                  if((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(2)))<0) edge(1)=1
                  if((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(3)))<0) edge(2)=1
                  if((a%poinStatus(e%lnods(2))*a%poinStatus(e%lnods(3)))<0) edge(3)=1

                  !subelements definition (First or last the alone node)
                  !anti-clockwise direction
                  auxsubelement=0
                  ialone=0
                  if(edge(1)==1 .and. edge(2)==1)then

                     ialone=1
                     ipoin=e%lnods(ialone)
                     if(a%poinStatus(ipoin)==-1)then
                        subelcod(:,1,1)= e%elcod(:,1)
                        subelcod(:,2,1)= IntersectionPoints(:,1)
                        subelcod(:,3,1)= IntersectionPoints(:,2)

                        subelcod(:,1,2)= e%elcod(:,2)
                        subelcod(:,2,2)= IntersectionPoints(:,2)
                        subelcod(:,3,2)= IntersectionPoints(:,1)

                        subelcod(:,1,3)= e%elcod(:,2)
                        subelcod(:,2,3)= e%elcod(:,3)
                        subelcod(:,3,3)= IntersectionPoints(:,2)
                     elseif(a%poinStatus(ipoin)== 1)then
                        subelcod(:,1,3)= e%elcod(:,1)
                        subelcod(:,2,3)= IntersectionPoints(:,1)
                        subelcod(:,3,3)= IntersectionPoints(:,2)

                        subelcod(:,1,1)= e%elcod(:,2)
                        subelcod(:,2,1)= IntersectionPoints(:,2)
                        subelcod(:,3,1)= IntersectionPoints(:,1)

                        subelcod(:,1,2)= e%elcod(:,2)
                        subelcod(:,2,2)= e%elcod(:,3)
                        subelcod(:,3,2)= IntersectionPoints(:,2)
                     end if

                  elseif(edge(1)==1 .and. edge(3)==1)then

                     ialone=2
                     ipoin=e%lnods(ialone)

                     if(a%poinStatus(ipoin)==-1)then
                        subelcod(:,1,1)= e%elcod(:,2)
                        subelcod(:,2,1)= IntersectionPoints(:,2)
                        subelcod(:,3,1)= IntersectionPoints(:,1)

                        subelcod(:,1,2)= e%elcod(:,1)
                        subelcod(:,2,2)= IntersectionPoints(:,1)
                        subelcod(:,3,2)= IntersectionPoints(:,2)

                        subelcod(:,1,3)= e%elcod(:,3)
                        subelcod(:,2,3)= e%elcod(:,1)
                        subelcod(:,3,3)= IntersectionPoints(:,2)
                     elseif(a%poinStatus(ipoin)==1)then
                        subelcod(:,1,3)= e%elcod(:,2)
                        subelcod(:,2,3)= IntersectionPoints(:,2)
                        subelcod(:,3,3)= IntersectionPoints(:,1)

                        subelcod(:,1,1)= e%elcod(:,1)
                        subelcod(:,2,1)= IntersectionPoints(:,1)
                        subelcod(:,3,1)= IntersectionPoints(:,2)

                        subelcod(:,1,2)= e%elcod(:,3)
                        subelcod(:,2,2)= e%elcod(:,1)
                        subelcod(:,3,2)= IntersectionPoints(:,2)
                     end if

                  elseif(edge(2)==1 .and. edge(3)==1)then

                     ialone=3
                     ipoin=e%lnods(ialone)

                     if(a%poinStatus(ipoin)==-1)then
                        subelcod(:,1,1)= e%elcod(:,3)
                        subelcod(:,2,1)= IntersectionPoints(:,1)
                        subelcod(:,3,1)= IntersectionPoints(:,2)

                        subelcod(:,1,2)= e%elcod(:,1)
                        subelcod(:,2,2)= IntersectionPoints(:,2)
                        subelcod(:,3,2)= IntersectionPoints(:,1)

                        subelcod(:,1,3)= e%elcod(:,1)
                        subelcod(:,2,3)= e%elcod(:,2)
                        subelcod(:,3,3)= IntersectionPoints(:,2)
                     elseif(a%poinStatus(ipoin)==1)then
                        subelcod(:,1,3)= e%elcod(:,3)
                        subelcod(:,2,3)= IntersectionPoints(:,1)
                        subelcod(:,3,3)= IntersectionPoints(:,2)

                        subelcod(:,1,1)= e%elcod(:,1)
                        subelcod(:,2,1)= IntersectionPoints(:,2)
                        subelcod(:,3,1)= IntersectionPoints(:,1)

                        subelcod(:,1,2)= e%elcod(:,1)
                        subelcod(:,2,2)= e%elcod(:,2)
                        subelcod(:,3,2)= IntersectionPoints(:,2)
                     end if

                  end if

                  a%CutElementPoints(ielem)%ialone=ialone

                  !Cartesian derivatives and Jacobian at center of gravity
                  call e%elmdcg

                  !Element length at center of gravity
                  call e%elmlen

                  auxgaus=3*e%pgaus
                  isube=1
                  do isube=1,3

                     !detjm computation
                     call e%SubElmdetjm(subelcod(:,:,isube))


                     !we only need the volume to ponderate iweigp

                     do igaus=1,e%pgaus
                        e%igaus=igaus

                        dvol = e%weigp(e%igaus)*abs(e%detjm)
                        subdvol = e%weigp(e%igaus)*abs(e%subdetjm)

                        call e%interpg(ndime,subelcod(:,:,isube),gpcod)
                        knode=e%pgaus*(isube-1)+e%igaus

                        call e%isoparinv(gpcod,a%CutElementPoints(ielem)%xloc(:,knode))

                        iweigp(knode)=e%weigp(e%igaus)*subdvol/dvol
                        !CutElement array
                        a%CutElementPoints(ielem)%weigp(knode)=iweigp(knode)

                     end do
                  end do

               elseif(e%pnode==4)then
                  call runend('ComputeSubelements only for linear triangles in 2d')
               end if

            elseif(ndime==3)then !only for tetrahedra

               if(e%pnode==4)then

                  edge=0
                  ialone=0
                  !face localization
                  if((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(2)))<0) edge(1)=1
                  if((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(3)))<0) edge(2)=1
                  if((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(4)))<0) edge(3)=1
                  if((a%poinStatus(e%lnods(2))*a%poinStatus(e%lnods(3)))<0) edge(4)=1
                  if((a%poinStatus(e%lnods(2))*a%poinStatus(e%lnods(4)))<0) edge(5)=1
                  if((a%poinStatus(e%lnods(3))*a%poinStatus(e%lnods(4)))<0) edge(6)=1


                  !in tetrahedra we can obtain 3 or 4 intersection points

                  if(a%CutElementPoints(ielem)%ninters==3)then

                     a%CutElementPoints(ielem)%nsube=4
                     auxsubelement=0
                     ialone=0

                     !find the alone element (first or last subelement)
                     !anti-clockwise direction between intersection points and positive right hand rule for the remaining node
                     if(edge(1)==1 .and. edge(2)==1 .and. edge(3)==1)then !node 1 alone
                        ialone=1
                        ipoin=e%lnods(ialone)
                        if(a%poinStatus(ipoin)==-1)auxsubelement=0
                        if(a%poinStatus(ipoin)== 1)auxsubelement=3
                        subelcod(:,1,1+auxsubelement)= IntersectionPoints(:,1)
                        subelcod(:,2,1+auxsubelement)= IntersectionPoints(:,3)
                        subelcod(:,3,1+auxsubelement)= IntersectionPoints(:,2)
                        subelcod(:,4,1+auxsubelement)= e%elcod(:,1)
                     elseif(edge(1)==1 .and. edge(4)==1 .and. edge(5)==1)then !node 2 alone
                        ialone=2
                        ipoin=e%lnods(ialone)
                        if(a%poinStatus(ipoin)==-1)auxsubelement=0
                        if(a%poinStatus(ipoin)== 1)auxsubelement=3
                        subelcod(:,1,1+auxsubelement)= IntersectionPoints(:,1)
                        subelcod(:,2,1+auxsubelement)= IntersectionPoints(:,2)
                        subelcod(:,3,1+auxsubelement)= IntersectionPoints(:,3)
                        subelcod(:,4,1+auxsubelement)= e%elcod(:,2)
                     elseif(edge(2)==1 .and. edge(4)==1 .and. edge(6)==1)then !node 3 alone
                        ialone=3
                        ipoin=e%lnods(ialone)
                        if(a%poinStatus(ipoin)==-1)auxsubelement=0
                        if(a%poinStatus(ipoin)== 1)auxsubelement=3
                        subelcod(:,1,1+auxsubelement)= IntersectionPoints(:,1)
                        subelcod(:,2,1+auxsubelement)= IntersectionPoints(:,3)
                        subelcod(:,3,1+auxsubelement)= IntersectionPoints(:,2)
                        subelcod(:,4,1+auxsubelement)= e%elcod(:,3)
                     elseif(edge(3)==1 .and. edge(5)==1 .and. edge(6)==1)then !node 4 alone
                        ialone=4
                        ipoin=e%lnods(ialone)
                        if(a%poinStatus(ipoin)==-1)auxsubelement=0
                        if(a%poinStatus(ipoin)== 1)auxsubelement=3
                        subelcod(:,1,1+auxsubelement)= IntersectionPoints(:,1)
                        subelcod(:,2,1+auxsubelement)= IntersectionPoints(:,2)
                        subelcod(:,3,1+auxsubelement)= IntersectionPoints(:,3)
                        subelcod(:,4,1+auxsubelement)= e%elcod(:,4)
                     end if

                     a%CutElementPoints(ielem)%ialone=ialone

                     !The other side with 3 subelements (the order depends on the elment cut always first the -1 fluid subelements)

                     icoun=0
                     do inode=1,e%pnode
                        if(inode /= ialone)then
                           icoun=icoun+1
                           ipoin=e%lnods(inode)
                           if(a%poinStatus(ipoin)==-1)then

                              if(icoun==1)then
                                 numnod(1)=inode
                                 subelcod(:,1,1) = e%elcod(:,inode)
                              elseif(icoun==2)then
                                 numnod(2)=inode
                                 subelcod(:,2,1) = e%elcod(:,inode)
                                 subelcod(:,1,2) = e%elcod(:,inode)
                                 subelcod(:,1,3) = e%elcod(:,inode)
                              elseif(icoun==3)then
                                 numnod(3)=inode
                                 subelcod(:,3,1) = e%elcod(:,inode)
                                 subelcod(:,2,3) = e%elcod(:,inode)
                              end if

                              !now using the intersection points

                              subelcod(:,4,1) = IntersectionPoints(:,1)
                              subelcod(:,2,2) = IntersectionPoints(:,1)
                              subelcod(:,3,2) = IntersectionPoints(:,2)
                              subelcod(:,4,2) = IntersectionPoints(:,3)
                              subelcod(:,3,3) = IntersectionPoints(:,1)
                              subelcod(:,4,3) = IntersectionPoints(:,3)


                           elseif(a%poinStatus(ipoin)==1)then

                              if(icoun==1)then
                                 numnod(1)=inode
                                 subelcod(:,1,2) = e%elcod(:,inode)
                              elseif(icoun==2)then
                                 numnod(2)=inode
                                 subelcod(:,2,2) = e%elcod(:,inode)
                                 subelcod(:,1,3) = e%elcod(:,inode)
                                 subelcod(:,1,4) = e%elcod(:,inode)
                              elseif(icoun==3)then
                                 numnod(3)=inode
                                 subelcod(:,3,2) = e%elcod(:,inode)
                                 subelcod(:,2,4) = e%elcod(:,inode)
                              end if
                              !now using the intersection points

                              subelcod(:,4,2) = IntersectionPoints(:,1)
                              subelcod(:,2,3) = IntersectionPoints(:,1)
                              subelcod(:,3,3) = IntersectionPoints(:,2)
                              subelcod(:,4,3) = IntersectionPoints(:,3)
                              subelcod(:,3,4) = IntersectionPoints(:,1)
                              subelcod(:,4,4) = IntersectionPoints(:,3)

                           end if
                        end if
                     end do

                     !Cartesian derivatives and Jacobian at center of gravity
                     call e%elmdcg

                     !Element length at center of gravity
                     call e%elmlen

                     auxgaus=4*e%pgaus
                     do isube=1,4

                        !detjm computation
                        call e%SubElmdetjm(subelcod(:,:,isube))

                        do igaus=1,e%pgaus
                           e%igaus=igaus

                           dvol = e%weigp(e%igaus)*abs(e%detjm)
                           subdvol = e%weigp(e%igaus)*abs(e%subdetjm)

                           call e%interpg(ndime,subelcod(:,:,isube),gpcod)
                           knode=e%pgaus*(isube-1)+e%igaus

                           call e%isoparinv(gpcod,a%CutElementPoints(ielem)%xloc(:,knode))

                           iweigp(knode)=e%weigp(e%igaus)*subdvol/dvol
                           !CutElement array
                           a%CutElementPoints(ielem)%weigp(knode)=iweigp(knode)

                        end do
                     end do

                  elseif(a%CutElementPoints(ielem)%ninters==4)then

                     a%CutElementPoints(ielem)%nsube=6

                     !auxiliar variables to do both sides of the cut element
                     auxsubelcod(1)=-1
                     auxsubelcod(2)=1
                     auxsubelement=0


                     do itype=1,2
                        !we have two subdomains with 3 element each one
                        !Always the -1 fluid first
                        if(auxsubelcod(itype)==-1) auxsubelement=0
                        if(auxsubelcod(itype)==1) auxsubelement=3

                        icoun=0
                        icount=0
                        do inode=1,e%pnode
                           ipoin=e%lnods(inode)
                           if(a%poinStatus(ipoin) == auxsubelcod(itype))then ! side definition
                              icoun=icoun+1
                              if(icoun==1)then
                                 numnod(1)=inode
                                 subelcod(:,1,1+auxsubelement)= e%elcod(:,inode)
                                 subelcod(:,1,3+auxsubelement)= e%elcod(:,inode)
                              elseif(icoun==2)then
                                 numnod(2)=inode
                                 subelcod(:,1,2+auxsubelement)= e%elcod(:,inode)
                                 subelcod(:,2,3+auxsubelement)= e%elcod(:,inode)
                              end if
                           end if
                        end do

                        icount(1)=1
                        icount(2)=1
                        icount(3)=2

                        !add the 2nd and the 3rd intersection-node in 1rst and 2nd subelements
                        do iedge=1,6
                           if(edge(iedge)/=0)then
                              select case(iedge)

                                 case(1)
                                    do inode=1,2
                                       if(numnod(inode)==1 .or. numnod(inode)==2)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:1)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 1)then
                                             nodeedge(icount(inode)-1,inode) = 2
                                          elseif(numnod(inode) == 2)then
                                             nodeedge(icount(inode)-1,inode) = 1
                                          end if
                                       end if
                                    end do

                                 case(2)
                                    do inode=1,2
                                       if(numnod(inode)==1 .or. numnod(inode)==3)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:2)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 1)then
                                             nodeedge(icount(inode)-1,inode) = 3
                                          elseif(numnod(inode) == 3)then
                                             nodeedge(icount(inode)-1,inode) = 1
                                          end if
                                       end if
                                    end do

                                 case(3)
                                    do inode=1,2
                                       if(numnod(inode)==1 .or. numnod(inode)==4)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:3)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 1)then
                                             nodeedge(icount(inode)-1,inode) = 4
                                          elseif(numnod(inode) == 4)then
                                             nodeedge(icount(inode)-1,inode) = 1
                                          end if
                                       end if
                                    end do

                                 case(4)
                                    do inode=1,2
                                       if(numnod(inode)==2 .or. numnod(inode)==3)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:4)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 2)then
                                             nodeedge(icount(inode)-1,inode) = 3
                                          elseif(numnod(inode) == 3)then
                                             nodeedge(icount(inode)-1,inode) = 2
                                          end if
                                       end if
                                    end do

                                 case(5)
                                    do inode=1,2
                                       if(numnod(inode)==2 .or. numnod(inode)==4)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:5)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 2)then
                                             nodeedge(icount(inode)-1,inode) = 4
                                          elseif(numnod(inode) == 4)then
                                             nodeedge(icount(inode)-1,inode) = 2
                                          end if
                                       end if
                                    end do

                                 case(6)
                                    do inode=1,2
                                       if(numnod(inode)==3 .or. numnod(inode)==4)then
                                          icount(inode) = icount(inode) + 1
                                          subelcod(:,icount(inode),inode+auxsubelement) = IntersectionPoints(:,sum(edge(1:6)))
                                          nodeedge(icount(inode)-1,inode) = iedge
                                          if(numnod(inode) == 3)then
                                             nodeedge(icount(inode)-1,inode) = 4
                                          elseif(numnod(inode) == 4)then
                                             nodeedge(icount(inode)-1,inode) = 3
                                          end if
                                       end if
                                    end do

                              end select
                           end if
                        end do


                        subelcod(:,4,2+auxsubelement) = subelcod(:,2,1+auxsubelement)
                        subelcod(:,3,3+auxsubelement) = subelcod(:,2,1+auxsubelement)
                        if(nodeedge(1,1)==nodeedge(1,2))then
                           subelcod(:,4,1+auxsubelement) = subelcod(:,3,2+auxsubelement)
                           subelcod(:,4,3+auxsubelement) = subelcod(:,3,2+auxsubelement)
                        elseif(nodeedge(1,1)==nodeedge(2,2))then
                           subelcod(:,4,1+auxsubelement) = subelcod(:,2,2+auxsubelement)
                           subelcod(:,4,3+auxsubelement) = subelcod(:,2,2+auxsubelement)
                        else
                           call runend('ComputeSubelements: wrong intersection in 3d')
                        end if
                     end do


                     a%CutElementPoints(ielem)%ialone=0

                     !Cartesian derivatives and Jacobian at center of gravity
                     call e%elmdcg

                     !Element length at center of gravity
                     call e%elmlen

                     auxgaus=6*e%pgaus
                     do isube=1,6

                        !detjm computation
                        call e%SubElmdetjm(subelcod(:,:,isube))

                        do igaus=1,e%pgaus
                           e%igaus=igaus

                           dvol = e%weigp(e%igaus)*abs(e%detjm)
                           subdvol = e%weigp(e%igaus)*abs(e%subdetjm)

                           call e%interpg(ndime,subelcod(:,:,isube),gpcod)
                           knode=e%pgaus*(isube-1)+e%igaus

                           call e%isoparinv(gpcod,a%CutElementPoints(ielem)%xloc(:,knode))

                           iweigp(knode)=e%weigp(e%igaus)*subdvol/dvol
                           !CutElement array
                           a%CutElementPoints(ielem)%weigp(knode)=iweigp(knode)

                        end do
                     end do

                  else
                     call runend('ComputeSubelements: wrong number of intersection points')
                  end if

               else
                  call runend('ComputeSubelements: only for linear tetrahedra in 3d')
               endif

            end if

         end if
      end do

      !CutElement Variables
      call a%Memor%allocObj(0,'xloc','ComputeSubelements',(countmemor*ndime*auxdim2)*rp)

      call a%Memor%dealloc(ndime,auxdim,auxdim1,subelcod,'subelcod','ComputeSubelements')
      call a%Memor%dealloc(auxdim1,edge,'edge','ComputeSubelements')
      call a%Memor%dealloc(ndime,gpcod,'gpcod','ComputeSubelements')
      call a%Memor%dealloc(auxdim2,iweigp,'iweigp','ComputeSubelements')
      call a%Memor%dealloc(e%ndime,4,IntersectionPoints,'IntersectionPoints','ComputeSubelements')

      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputeSubelements')

   end subroutine

   subroutine ComputeEnrichBubble(a,ndime)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh)   :: a
      class(FiniteElement), pointer :: e => NULL()


      integer(ip), intent(in)  :: ndime
      integer(ip)              :: inode,auxdim
      integer(ip)              :: ielem,nelem
      integer(ip)              :: ialone,ipoin,isube,igaus,auxgaus,knode,auxpnode
      !enrich bubble element
      integer(ip)              :: nipoin,kpoin,aux_12_34,aux_13_24,aux_14_23
      real(rp), allocatable    :: kxloc(:,:),kweigp(:)
      integer(ip)              :: auxcount,countmemor
      real(rp), allocatable    :: gpcod(:)
      integer(ip), allocatable :: auxlist(:)
      real(rp)                 :: tolerance,knorm
      real(rp)                 :: tolerance3d

      call a%Mesh%GetNelem(nelem)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ComputeEnrichBubble')

      call a%Memor%alloc(ndime,e%mnode,kxloc,'kxloc','ComputeEnrichBubble')
      call a%Memor%alloc(ndime,gpcod,'gpcod','ComputeEnrichBubble')
      call a%Memor%alloc(e%mnode,kweigp,'kweigp','ComputeEnrichBubble')
      call a%Memor%alloc(e%mnode,auxlist,'auxlist','ComputeEnrichBubble')

      countmemor=0

      tolerance=1e-16
      tolerance3d=1e-10

      do ielem=1,nelem

         if(a%CutElementPoints(ielem)%elemStatus == 0)then !cut lement

            call a%Mesh%ElementLoad(ielem,e)

            !Enriched bubble functions
            allocate(a%CutElementPoints(ielem)%kEnrichBub(e%mnode))
            countmemor=countmemor+1

            if(ndime==2)then !only for triangles
               if(e%pnode==3)then

                  ialone=a%CutElementPoints(ielem)%ialone
                  !Computation of enrich bubble functions
                  !Compute shape functions at the intersection points
                  nipoin=a%CutElementPoints(ielem)%ninters
                  do kpoin = 1,nipoin
                     !call e%isoparinv(a%CutElementPoints(ielem)%ipoints(:,kpoin),kxloc(:,kpoin))
                     kweigp(kpoin) = 1.0_rp/nipoin
                  end do

                  !The rutine give the needed shape functions associated
                  call e%SetParticularGaussPoints(a%Memor,nipoin,a%CutElementPoints(ielem)%LocalIntersectionPoints,kweigp)

                  auxlist=0
                  auxlist(1)=ialone
                  if(ialone==1)then
                     auxlist(2)=2
                     auxlist(3)=3
                  elseif(ialone==2)then
                     auxlist(2)=1
                     auxlist(3)=3
                  elseif(ialone==3)then
                     auxlist(2)=1
                     auxlist(3)=2
                  end if

                  a%CutElementPoints(ielem)%kEnrichBub(:)=0.0_rp

                  if(abs(e%shape(auxlist(1),1))<tolerance)then
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))=0.0_rp
                  else
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))=1.0_rp/e%shape(auxlist(1),1)
                  end if

                  if(abs(e%shape(auxlist(2),1))<tolerance)then
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(2))=0.0_rp
                  else
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(2))=1.0_rp/e%shape(auxlist(2),1)
                  end if

                  if(abs(e%shape(auxlist(3),2))<tolerance)then
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))=0.0_rp
                  else
                     a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))=a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))*(e%shape(auxlist(1),2)/e%shape(auxlist(3),2))
                  end if

               elseif(e%pnode==4)then
                  call runend('ComputeEnrichBubble only for linear triangles in 2d')
               end if

            elseif(ndime==3)then !only for tetrahedra

               if(e%pnode==4)then

                  !in tetrahedra we can obtain 3 or 4 intersection points

                  if(a%CutElementPoints(ielem)%ninters==3)then

                     ialone=a%CutElementPoints(ielem)%ialone

                     !Computation of enrich bubble functions
                     !Compute shape functions at the intersection points
                     nipoin=a%CutElementPoints(ielem)%ninters
                     do kpoin = 1,nipoin
                        !call e%isoparinv(a%CutElementPoints(ielem)%ipoints(:,kpoin),kxloc(:,kpoin))
                        kweigp(kpoin) = 1.0_rp/nipoin
                     end do

                     !The rutine give the needed shape functions associated
                     call e%SetParticularGaussPoints(a%Memor,nipoin,a%CutElementPoints(ielem)%LocalIntersectionPoints,kweigp)

                     auxlist=0
                     auxlist(1)=ialone
                     if(ialone==1)then
                        auxlist(2)=2
                        auxlist(3)=3
                        auxlist(4)=4
                     elseif(ialone==2)then
                        auxlist(2)=1
                        auxlist(3)=3
                        auxlist(4)=4
                     elseif(ialone==3)then
                        auxlist(2)=1
                        auxlist(3)=2
                        auxlist(4)=4
                     elseif(ialone==4)then
                        auxlist(2)=1
                        auxlist(3)=2
                        auxlist(4)=3
                     end if

                     if(abs(e%shape(auxlist(1),1))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))=0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))=1.0_rp/e%shape(auxlist(1),1)
                     end if

                     if(abs(e%shape(auxlist(2),1))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(2))=0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(2))=1.0_rp/e%shape(auxlist(2),1)
                     end if

                     if(abs(e%shape(auxlist(3),2))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))=0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))=a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))*(e%shape(auxlist(1),2)/e%shape(auxlist(3),2))
                     end if

                     if(abs(e%shape(auxlist(4),3))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(4))=0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(4))=a%CutElementPoints(ielem)%kEnrichBub(auxlist(1))*(e%shape(auxlist(1),3)/e%shape(auxlist(4),3))
                     end if

                  elseif(a%CutElementPoints(ielem)%ninters==4)then

                     !Definition of the configuration of the cut
                     !Case 1: node 1-2 in the same side and 3-4 in the other
                     !Case 2: node 1-3 in the same side and 2-4 in the other
                     !Case 3: node 1-4 in the same side and 2-3 in the other
                     aux_12_34=0
                     aux_13_24=0
                     aux_14_23=0

                     if(((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(2)))>0) .and. ((a%poinStatus(e%lnods(3))*a%poinStatus(e%lnods(4)))>0))then
                        aux_12_34=1
                     elseif(((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(3)))>0) .and. ((a%poinStatus(e%lnods(2))*a%poinStatus(e%lnods(4)))>0))then
                        aux_13_24=1
                     elseif(((a%poinStatus(e%lnods(1))*a%poinStatus(e%lnods(4)))>0) .and. ((a%poinStatus(e%lnods(2))*a%poinStatus(e%lnods(3)))>0))then
                        aux_14_23=1
                     end if

                     !Computation of enrich bubble functions
                     !Compute shape functions at the intersection points
                     nipoin=a%CutElementPoints(ielem)%ninters
                     do kpoin = 1,nipoin
                        !call e%isoparinv(a%CutElementPoints(ielem)%ipoints(:,kpoin),kxloc(:,kpoin))
                        kweigp(kpoin) = 1.0_rp/nipoin
                     end do

                     !The rutine give the needed shape functions associated
                     call e%SetParticularGaussPoints(a%Memor,nipoin,a%CutElementPoints(ielem)%LocalIntersectionPoints,kweigp)

                     auxlist=0
                     if(aux_12_34==1)then
                        auxlist(1)=1
                        auxlist(2)=3
                        auxlist(3)=2
                        auxlist(4)=4
                     elseif(aux_13_24==1)then
                        auxlist(1)=1
                        auxlist(2)=2
                        auxlist(3)=3
                        auxlist(4)=4
                     elseif(aux_14_23==1)then
                        auxlist(1)=1
                        auxlist(2)=1
                        auxlist(3)=4
                        auxlist(4)=3
                     end if

                     if(abs(e%shape(auxlist(1),1))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(1)) = 0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(1)) = 1.0_rp/e%shape(auxlist(1),1)
                     endif
                     if(abs(e%shape(auxlist(2),1))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(2)) = 0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(2)) = 1.0_rp/e%shape(auxlist(2),1)
                     endif
                     if(abs(e%shape(auxlist(3),3))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(3)) = 0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))= a%CutElementPoints(ielem)%kEnrichBub(auxlist(2))*(e%shape(auxlist(2),3)/e%shape(auxlist(3),3))
                     endif
                     if(abs(e%shape(auxlist(4),4))<tolerance3d)then
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(4)) = 0.0_rp
                     else
                        a%CutElementPoints(ielem)%kEnrichBub(auxlist(4))= a%CutElementPoints(ielem)%kEnrichBub(auxlist(3))*(e%shape(auxlist(3),4)/e%shape(auxlist(4),4))
                     endif


                  else
                     call runend('ComputeEnrichBubble: wrong number of intersection points')
                  end if

               else
                  call runend('ComputeEnrichBubble: only for linear tetrahedra in 3d')
               endif

            end if

         end if
      end do

      !CutElement Variables
      call a%Memor%allocObj(0,'kEnrichBub','ComputeEnrichBubble',(countmemor)*rp)

      call a%Memor%dealloc(ndime,gpcod,'gpcod','ComputeEnrichBubble')
      call a%Memor%dealloc(ndime,size(kxloc,2),kxloc,'kxloc','ComputeEnrichBubble')
      call a%Memor%dealloc(size(kweigp,1),kweigp,'kweigp','ComputeEnrichBubble')
      call a%Memor%dealloc(size(auxlist,1),auxlist,'auxlist','ComputeEnrichBubble')


      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputeEnrichBubble')

   end subroutine

   subroutine GetEnrichBuble(a,ielem,pnode,Kshape)
      use typre
      implicit none
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ielem,pnode
      real(rp), intent(out)    :: kshape(pnode)

      kshape(1:pnode)=a%CutElementPoints(ielem)%kEnrichBub(1:pnode)

   end subroutine

   subroutine GetWeigpCut(a,ielem,ndime,weigp)
      use typre
      implicit none
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ielem,ndime
      real(rp), intent(out)    :: weigp((ndime+1)*((ndime-1)*(ndime-1)+2))

      weigp(:)=a%CutElementPoints(ielem)%weigp(:)

   end subroutine

   subroutine GetIalone(a,ielem,ialone)
      use typre
      implicit none
      class(CutMesh) :: a

      integer(ip), intent(in)  :: ielem
      integer(ip), intent(out)    :: ialone

      ialone = a%CutElementPoints(ielem)%ialone

   end subroutine

   subroutine GetXlocCut(a,ielem,ndime,xloc)
      use typre
      implicit none
      class(CutMesh) :: a
      integer(ip), intent(in)  :: ielem,ndime
      real(rp), intent(out)    :: xloc(ndime,(ndime+1)*((ndime-1)*(ndime-1)+2))

      xloc(:,:)=a%CutElementPoints(ielem)%xloc(:,:)

   end subroutine

   subroutine deallocCutMesh(a)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh)   :: a

      integer(ip)       :: auxdim,auxdim1,npoin,ndime,nelem,auxsize

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)

      auxdim=(ndime-1)*(ndime-1)+2
      auxdim1=(ndime+1)*auxdim

      !status element: cut or not cut
      call a%Memor%dealloc(size(a%PoinStatus),a%poinStatus,'poinStatus','allocCutMesh')

      !CutElement type
      auxsize = size(a%CutElementPoints)
      call a%Memor%deallocObj(0,'CutElementPoints','allocCutMesh',auxsize*1_ip)
      deallocate(a%CutElementPoints)

   end subroutine

   subroutine deallocCutElement(a)
      use typre
      use Mod_Memor
      implicit none
      class(CutMesh)   :: a
      class(FiniteElement), pointer :: e => NULL()

      integer(ip)       :: auxdim,auxdim1,auxdim2,auxdim3,ndime,ielem
      integer(ip)       :: countmemor

      logical :: iskenrichbubble = .false.

      call a%Mesh%GetNdime(ndime)

      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','deallocCutElement')

      auxdim=(ndime+1) !number of nodes in linear elements 2d and 3d

      !number of intersection points
      if(ndime==2)then
         auxdim1=2
      elseif(ndime==3)then
         auxdim1=4
      endif
      auxdim2=(ndime-1)*(ndime-1)+2 !number of subelements
      auxdim3=auxdim*auxdim2 !total number of gaus points in a cut element (linear elements)



      countmemor=0
      do ielem=1,size(a%CutElementPoints,1)  !Should be nelem, but if adaptive refinement nelem has already changed, this is safe.

         if(a%CutElementPoints(ielem)%elemStatus==0)then

            countmemor=countmemor+1

            if (allocated(a%CutElementPoints(ielem)%weigp)) then

            else
               write(*,*) 'basdfasdsdf'

               write(*,*) 'elementStatus',ielem, a%CutElementPoints(ielem)%elemStatus

            endif

            deallocate(a%CutElementPoints(ielem)%weigp)
            deallocate(a%CutElementPoints(ielem)%LocalIntersectionPoints)
            deallocate(a%CutElementPoints(ielem)%xloc)
            if (allocated(a%CutElementPoints(ielem)%kEnrichBub)) then
               deallocate(a%CutElementPoints(ielem)%kEnrichBub)
               iskenrichbubble = .true.
            endif
         end if
      end do

      call a%Memor%deallocObj(0,'weigp','deallocCutElement',(countmemor*auxdim3)*rp)
      call a%Memor%deallocObj(0,'ipoints','deallocCutElement',(countmemor*ndime*auxdim1)*rp)
      call a%Memor%deallocObj(0,'xloc','deallocCutElement',(countmemor*ndime*auxdim3)*rp)
      if (iskenrichbubble) call a%Memor%deallocObj(0,'kEnrichBub','deallocCutElement',(countmemor)*rp)
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','deallocCutElement')

   end subroutine

   subroutine GetExternalNormal(a,ielem,e,UnitOrthogonalVector)
      implicit none
      class(CutMesh) :: a
      integer(ip) :: ielem
      class(FiniteElement), target :: e

      real(rp):: points(e%ndime,10)

      integer(ip) :: nipoin,ipoin,jpoin,jnode
      real(rp)    :: toleralce
      real(rp)    :: boundaryVector(3),hipo,tmp
      real(rp)    :: UnitBoundaryVector(3),UnitOrthogonalVector(3),BoundaryCenter(3),checkvector(3)
      !3d case
      real(rp)    :: boundaryVector2(3),hipo2,localx3(3),hipo3,projCheck
      integer(ip) :: ialone,idime,jdime,inode !node alone


      real(rp), pointer :: coord(:) => NULL()
      integer(ip) :: ipStatus, jpStatus

      toleralce=1e-6


      !Initializations
      boundaryVector=0.0_rp
      hipo=0.0_rp
      boundaryVector2=0.0_rp
      hipo2=0.0_rp
      hipo3=0.0_rp
      UnitBoundaryVector=0.0_rp
      UnitOrthogonalVector=0.0_rp
      localx3=0.0_rp
      BoundaryCenter=0.0_rp
      checkvector=0.0_rp
      projCheck=0.0_rp
      !using as a reference points(:,1)

      call a%GetREalIntersectionPoints(ielem,e,points)
      call a%GetNinters(ielem,nipoin)

      do idime=1,e%ndime
         boundaryVector(idime) = (points(idime,1) - points(idime,2))
         hipo = hipo + boundaryVector(idime)*boundaryVector(idime)
      end do

      hipo=sqrt(hipo)

      !known direction
      if(hipo>toleralce)then
         UnitBoundaryVector   = boundaryVector/hipo
      elseif(hipo<=toleralce)then
         UnitBoundaryVector = 0.0_rp
      end if


      if(e%ndime==2)then
         !in 2d the orthogonal vector is directly
         UnitOrthogonalVector(1)=-UnitBoundaryVector(2)
         UnitOrthogonalVector(2)=UnitBoundaryVector(1)
      elseif(e%ndime==3)then

         do idime=1,e%ndime
            boundaryVector2(idime) = boundaryVector2(idime) + (points(idime,1) - points(idime,3))
            hipo2 = hipo2 + boundaryVector2(idime)*boundaryVector2(idime)
         end do
         hipo2=sqrt(hipo2)

         if(hipo2>toleralce)then
            localx3=boundaryVector2/hipo2
         elseif(hipo2<=toleralce)then
            localx3=0.0_rp
         endif

         !We use the classical vectorial product rule
         UnitOrthogonalVector(1)= UnitBoundaryVector(2)*localx3(3)-UnitBoundaryVector(3)*localx3(2)
         UnitOrthogonalVector(2)= -(UnitBoundaryVector(1)*localx3(3)-UnitBoundaryVector(3)*localx3(1))
         UnitOrthogonalVector(3)= UnitBoundaryVector(1)*localx3(2)-UnitBoundaryVector(2)*localx3(1)

         do idime=1,e%ndime
            hipo3=hipo3+UnitOrthogonalVector(idime)*UnitOrthogonalVector(idime)
         end do
         hipo3=sqrt(hipo3)

         if(hipo3>toleralce)then
            UnitOrthogonalVector=UnitOrthogonalVector/hipo3
         elseif(hipo3<=toleralce)then
            UnitOrthogonalVector=0.0_rp
         end if
      end if

      !we need check if the normal is outer or inner. We need the outer

      do ipoin=1,nipoin
         do idime=1,e%ndime
            BoundaryCenter(idime)=BoundaryCenter(idime) + points(idime,ipoin)
         end do
      end do
      !Geometrical gravity center
      BoundaryCenter=BoundaryCenter/nipoin

      call a%GetIalone(ielem,ialone)
      if(ialone/=0)then
         ipoin=e%lnods(ialone)
         coord => e%elcod(:,ialone)

         do idime=1,e%ndime
            checkvector(idime)= checkvector(idime) + (coord(idime)-BoundaryCenter(idime))
         end do

         projCheck   = dot_product(checkvector,UnitOrthogonalVector)

         call a%GetPointType(ipoin,ipStatus)

         if(projCheck>0.0_rp .and. ipStatus>= 1) UnitOrthogonalVector=-UnitOrthogonalVector
         if(projCheck<0.0_rp .and. ipStatus < 1) UnitOrthogonalVector=-UnitOrthogonalVector
      elseif(ialone==0)then

         do inode=1,e%pnode
            jpoin=e%lnods(inode)
            call a%GetPointType(jpoin,jpStatus)
            if(jpStatus < 1) jnode = inode
         end do

         coord => e%elcod(:,jnode)


         do idime=1,e%ndime
            checkvector(idime)= checkvector(idime) + (coord(idime)-BoundaryCenter(idime))
         end do

         projCheck   = dot_product(checkvector,UnitOrthogonalVector)

         if(projCheck<0.0_rp) UnitOrthogonalVector=-UnitOrthogonalVector

      end if

   end subroutine

end module

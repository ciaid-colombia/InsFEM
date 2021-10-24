module Mod_Element
   use typre
   implicit none
   private
   public FiniteElement

   type FiniteElement

      integer(ip), pointer :: lboel(:) => NULL(), lnods(:) => NULL()
      integer(ip), allocatable :: lnodb(:)
      real(rp), allocatable :: elcod(:,:),bocod(:,:)

      !Force closed rule
      integer(ip) :: kfl_ClosedRule = 0

      !Is already computed linearderivatives
      logical :: kfl_IsAlreadyComputedLinearDerivatives = .false.

      !Dimensions
      integer(ip) :: &
         ndime,&           !Dimension of the problem? element?
         mnode,&           !Maximum Number of nodes of any element of the mesh
         mnodb,&           !Maximum Number of nodes in the boundary of any element of the mesh
         mface,&           !Maximum Number of faces of any element of the mesh
         ntens,&
         mlapl

      !Periodic Boundary Conditions
      logical :: kfl_perio = .false. !Element has periodic boundary conditions

      !Current element info
      integer(ip) :: &
         pnode,&        !Number of nodes of the current element
         pnodb,&        !Number of nodes in the boundary of the current elemetn
         pgaus,&        !# of gauss points current element
         pgaub,&        !# of gauss points current element in the boundary
         pface,&        !Number of faces of the current element
         iface,&        !Current face
         linea,&        !1:current element is 1st order (i.e. P1)
         ielty,&        !type of element
         ielty0,&       !Previous element type, for avoiding resetting the pointers
         iblty,&        !type of boundary element
         iblty0,&       !Previous type of boundary element
         igaus,&        !Current gauss point
         igaub,&        !Current gauss point, on the boundary
         llapl,&
         npol,&         !polynomial order of integration
         npol4,&        !npol^4
         npoinLocal     !Points number of the local mesh

      !Shape functions and derivatives
      real(rp), pointer  :: &
            shape(:,:) => NULL(),&               !Shape functions
            deriv(:,:,:) => NULL(),&
            heslo(:,:,:) => NULL(),&
            weigp(:) => NULL()                   !Weigth of the gauss points
      !Shape functions and derivatives center of gravity
      real(rp), pointer  :: &
            shacg(:) => NULL(),&
            dercg(:,:) => NULL()
      !Boundary shape functions and derivatives
      real(rp), pointer  :: &
            shapb(:,:) => NULL(),&
            derib(:,:,:) => NULL(),&
            weigb(:) => NULL()
      real(rp), pointer  :: &
            shaga(:,:) => NULL(),&
            shagb(:,:) => NULL()
      real(rp) :: &
            hnatu,&
            weicg

      !Element jacobians and derivatives at gauss points and center of gravity
      real(rp), pointer ::&
         cartd(:,:) => NULL(),&                 ! Cartesian derivatives
         cartb(:,:) => NULL(),&                 ! Cartesian derivatives in the boundary
         xjaci(:,:) => NULL(),&                 ! Inverse of jacobian
         xjacm(:,:) => NULL(),&                 ! Jacobian
         tragl(:,:) => NULL(),&                 ! Inverse of jacobian at c.o.g
         hleng(:) => NULL()                     ! Element length at c.o.g
      real(rp), pointer ::&
         hessi(:,:) => NULL()                   ! Hessian

      !Element values at Gauss points
      real(rp) ::&
         detjm                                 ! Jacobian determinant

      !Boundary jacobians and derivatives
      real(rp), pointer ::&
         baloc(:,:) => NULL()                  ! Local base of boundary element

      !Boundary values at Gauss points
      real(rp) ::&
         eucta                        ! Norm of baloc

      !Particular Gauss points
      integer(ip) :: mParticularGaussPoints = 0

      real(rp) :: subdetjm

      !Particular Shape functions and derivatives

      real(rp), allocatable  :: &
            Particularshape(:,:),&               !Shape functions
            Particularderiv(:,:,:),&
            Particularheslo(:,:,:),&
            Particularweigp(:),&                 !Weigth of the gauss points
            Particularcartd(:,:),&
            Particularhessi(:,:)

       integer(ip) :: &
            mPnode                      !Modify pnode used in enriched elements


contains

      procedure :: alloc  => ElementAlloc
      procedure :: dealloc=> ElementDealloc
      procedure :: elmdel => ElementElmdel
      procedure :: elmder => ElementElmder
      procedure :: elmdcg => ElementElmdcg
      procedure :: SubElmdetjm => SubElementDetjm
      procedure :: elmhes => ElementElmhes
      procedure :: elmderb=> ElementElmderb
      procedure :: gather => ElementGather
      procedure :: gatherb=> BoundaryGather
      procedure :: elmlen => ElementElmlen
      procedure :: bounor => ElementBounor
      procedure :: elenor => ElementElenor
      procedure :: interpc
      procedure :: interpg
      procedure :: interpb
      procedure :: hessian
      procedure :: gradient
      procedure :: gradientb
      procedure :: divergence
      procedure :: isopar
      procedure :: isoparinv
      procedure :: ForceClosedRule
      procedure :: DefaultRule
      procedure :: ComputeShafun
      procedure :: SetParticularGaussPoints
      procedure :: SetEnrichedElement
      procedure :: DeallocParticularGaussPoints
      procedure :: GaussToNodes
      procedure :: GaussToNodesB
      procedure :: GetWeightFactor
   end type

contains

   subroutine ElementAlloc(e,Memor,runam)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElement) :: e
      character(*) :: runam
      type(MemoryMan) :: Memor

      call Memor%palloc(e%mnode,e%lnods,'e%lnods',runam)
      call Memor%alloc(e%ndime,e%mnode,e%elcod,'e%elcod',runam)

      call Memor%alloc(e%mnodb,e%lnodb,'e%lnodb',runam)
      call Memor%alloc(e%ndime,e%mnodb,e%bocod,'e%bocod',runam)
      call Memor%palloc(e%ndime,e%ndime,e%baloc,'e%baloc',runam)

      call Memor%palloc(e%ndime,e%mnode,e%cartd,'e%cartd',runam)
      call Memor%palloc(e%ndime,e%mnode,e%cartb,'e%cartb',runam)
      call Memor%palloc(e%ndime,e%ndime,e%xjaci,'e%xjaci',runam)
      call Memor%palloc(e%ndime,e%ndime,e%xjacm,'e%xjacm',runam)
      call Memor%palloc(e%ndime,e%ndime,e%tragl,'e%tragl',runam)
      call Memor%palloc(e%ndime        ,e%hleng,'e%hleng',runam)

      if(e%mlapl>0) then
         call Memor%palloc(e%ntens,e%mnode*e%mlapl        ,e%hessi,'e%hessi',runam)
      end if

   end subroutine



   subroutine ElementDealloc(e,Memor,runam)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor
      character(*) :: runam

      call Memor%pdealloc(size(e%lnods,1),e%lnods,'e%lnods',runam)
      call Memor%dealloc(e%ndime,size(e%elcod,2),e%elcod,'e%elcod',runam)

      call Memor%dealloc(e%mnodb,e%lnodb,'e%lnodb',runam)
      call Memor%dealloc(e%ndime,e%mnodb,e%bocod,'e%bocod',runam)
      call Memor%pdealloc(e%ndime,e%ndime,e%baloc,'e%baloc',runam)

      call Memor%pdealloc(e%ndime,size(e%cartd,2),e%cartd,'e%cartd',runam)
      call Memor%pdealloc(e%ndime,size(e%cartb,2),e%cartb,'e%cartb',runam)
      call Memor%pdealloc(e%ndime,e%ndime,e%xjaci,'e%xjaci',runam)
      call Memor%pdealloc(e%ndime,e%ndime,e%xjacm,'e%xjacm',runam)
      call Memor%pdealloc(e%ndime,e%ndime,e%tragl,'e%tragl',runam)
      call Memor%pdealloc(e%ndime        ,e%hleng,'e%hleng',runam)

      if(e%mlapl>0) then
         call Memor%pdealloc(e%ntens,size(e%hessi,2),e%hessi,'e%hessi',runam)
      end if

      if(e%mParticularGaussPoints /= 0)then
         call e%DeallocParticularGaussPoints(Memor)
         e%mParticularGaussPoints = 0
      end if
   end subroutine

   subroutine ElementElmdel(e)
      implicit none
      class(FiniteElement) :: e

      e%kfl_IsAlreadyComputedLinearDerivatives = .true.
      if (e%linea == 1) call elmdel(e%pnode,e%ndime,e%elcod,e%cartd,e%detjm,e%linea,e%xjacm,e%xjaci)
   end subroutine

   subroutine ElementElmdcg(e)
      implicit none
      class(FiniteElement) :: e

      if (e%linea == 1) then
         if (e%kfl_IsAlreadyComputedLinearDerivatives .eqv. .false.) then
            call ElementElmdel(e)
         endif
      else
         call elmder(e%pnode,e%ndime,e%dercg,e%elcod,e%cartd,e%detjm,e%xjacm,e%xjaci)
      endif

   end subroutine

   subroutine SubElementDetjm(e,subelcod)
      implicit none
      class(FiniteElement) :: e
      real(rp)  :: subelcod(e%ndime,e%pnode)

      if (e%linea == 1) then
         call elemdetjm(e%pnode,e%ndime,subelcod,e%subdetjm)
      else
         call elmder(e%pnode,e%ndime,e%dercg,subelcod,e%cartd,e%subdetjm,e%xjacm,e%xjaci)
      endif

   end subroutine

   subroutine ElementElmder(e)
      implicit none
      class(FiniteElement) :: e

      if (e%linea == 1) then
         if (e%kfl_IsAlreadyComputedLinearDerivatives .eqv. .false.) then
            call ElementElmdel(e)
         endif
      else
         call elmder(e%pnode,e%ndime,e%deriv(1,1,e%igaus),e%elcod,e%cartd,e%detjm,e%xjacm,e%xjaci)
      endif

   end subroutine

   subroutine ElementElmhes(e)
      implicit none
      class(FiniteElement) :: e

      if (e%linea == 1) then
      else
         call elmhes(e%heslo(1,1,e%igaus),e%hessi,e%ndime,e%pnode,e%ntens,e%xjaci,e%deriv(1,1,e%igaus),e%elcod)
      endif

   end subroutine

   !Calculate cartesian derivates at boundary gauss points
   subroutine ElementElmderb(e)
      implicit none
      class(FiniteElement) :: e

      integer(ip) :: inode,idime,inodb,igaus
      real(rp)    :: aux1

      if (e%linea == 1) then
         if (e%kfl_IsAlreadyComputedLinearDerivatives .eqv. .false.) then
            call ElementElmdel(e)
         endif
         e%cartb = e%cartd
      else
          e%cartb=0.0_rp
          do igaus=1,e%pgaus
            e%igaus = igaus
            call e%elmder
            do inodb=1,e%pnodb
               aux1 = e%shapb(inodb,e%igaub)*e%shaga(igaus,e%lboel(inodb))
                do inode=1,e%pnode
                  do idime=1,e%ndime
                     e%cartb(idime,inode)=e%cartb(idime,inode)&
                        +aux1*e%cartd(idime,inode)
                  end do
               end do
            end do
         end do
      endif

   end subroutine

   subroutine ElementGather(e,ndime,elval,glval)
      implicit none
      class(FiniteElement)  :: e
      integer(ip)           :: ndime
      real(rp), intent(in)  :: glval(ndime,*)
      real(rp), intent(out) :: elval(ndime,e%pnode)
      integer(ip) :: inode,ipoin

      elval(:,:) = glval(:,e%lnods(1:e%pnode))

   end subroutine

   subroutine BoundaryGather(e,ndime,elval,glval)
      implicit none
      class(FiniteElement)  :: e
      integer(ip)           :: ndime
      real(rp), intent(in)  :: glval(ndime,*)
      real(rp), intent(out) :: elval(ndime,e%pnodb)
      integer(ip) :: inode,ipoin

      elval(:,:) = glval(:,e%lnodb(1:e%pnodb))

   end subroutine

   subroutine ElementBounor(e)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: cocog(e%ndime),produ
      integer(ip) :: idime

      !Evaluates the tangent vectors
      call mbmabt(e%baloc,e%bocod,e%derib(1,1,e%igaub),e%ndime,e%ndime-1,e%pnodb)

      !Computes the normal vector
      if(e%ndime==2) then
         e%baloc(1,2)= e%baloc(2,1)
         e%baloc(2,2)=-e%baloc(1,1)
      else
         call vecpro(e%baloc(1,1),e%baloc(1,2),e%baloc(1,3),3)
      end if

      !Computes the Euclidian norm of the normal vector
      call vecnor(e%baloc(1,e%ndime),e%ndime,e%eucta,2)

      !Coordinates center of gravity
      do idime = 1,e%ndime
         cocog(idime) = sum(e%elcod(idime,1:e%pnode))/real(e%pnode)
      enddo

      !Normalize baloc
      do idime=1,e%ndime
         call vecuni(e%ndime,e%baloc(1,idime),produ)
      end do

      !Ensure tangential vectors orthogonality
      if (e%ndime == 3) then
         produ = dot_product(e%baloc(1:e%ndime,1),e%baloc(1:e%ndime,2))
         e%baloc(:,2) = e%baloc(:,2) - produ*e%baloc(:,1)
         call vecuni(e%ndime,e%baloc(1,2),produ)
      endif

      !Calculate n.v
      produ=dot_product(cocog(1:e%ndime)-e%bocod(1:e%ndime,1),&
            e%baloc(1:e%ndime,e%ndime))

      !If n.v > 0, the normal is not the external one
      if(produ>0.0_rp) then
         e%baloc(:,1)     = -e%baloc(:,1)             !t1=-t1
         e%baloc(:,e%ndime) = -e%baloc(:,e%ndime)     !n =-n
      end if

   end subroutine

   subroutine ElementElenor(e)
      implicit none
      class(FiniteElement) :: e

      real(rp) :: cocog(e%ndime),produ
      integer(ip) :: idime

      !Evaluates the tangent vectors
      call mbmabt(e%baloc,e%bocod,e%derib(1,1,e%igaub),e%ndime,e%ndime-1,e%pnodb)

      !Computes the normal vector
      if(e%ndime==2) then
         e%baloc(1,2)= e%baloc(2,1)
         e%baloc(2,2)=-e%baloc(1,1)
      else
         call vecpro(e%baloc(1,1),e%baloc(1,2),e%baloc(1,3),3)
      end if

      !Computes the Euclidian norm of the tangent vector  (only a tangent (3D)?)
      call vecnor(e%baloc(1,e%ndime),e%ndime,e%eucta,2)

      !Coordinates center of gravity
      do idime = 1,e%ndime
         cocog(idime) = sum(e%elcod(idime,1:e%pnode))/real(e%pnode)
      enddo

      !Normalize baloc
      do idime=1,e%ndime
         call vecuni(e%ndime,e%baloc(1,idime),produ)
      end do

      !Ensure tangential vectors orthogonality
      if (e%ndime == 3) then
         produ = dot_product(e%baloc(1:e%ndime,1),e%baloc(1:e%ndime,2))
         e%baloc(:,2) = e%baloc(:,2) - produ*e%baloc(:,1)
         call vecuni(e%ndime,e%baloc(1,2),produ)
      endif

      !Calculate n.v
      produ=dot_product(cocog(1:e%ndime)-e%bocod(1:e%ndime,1),&
            e%baloc(1:e%ndime,e%ndime))

      !If n.v > 0, the normal is not the external one
      if(produ>0.0_rp) then
         e%baloc(:,1)     = -e%baloc(:,1)             !t1=-t1
         e%baloc(:,e%ndime) = -e%baloc(:,e%ndime)     !n =-n
      end if

   end subroutine

   subroutine ElementElmlen(e)
      implicit none
      class(FiniteElement) :: e

      real(rp)    :: enor0,h_tem
      integer(ip) :: idime,jdime

       if (e%linea == 0) then
         e%xjacm=matmul(e%elcod(:,1:e%pnode),transpose(e%dercg(:,1:e%pnode))) ! J^t
         call invmtx(e%xjacm,e%tragl,e%detjm,e%ndime) ! J^(-t)
      else
         e%tragl = e%xjaci
      endif

      !Element length
      do idime=1,e%ndime
         enor0=0.0_rp
         do jdime=1,e%ndime
            enor0=enor0+e%tragl(idime,jdime)**2
         end do
         enor0=sqrt(enor0)
         !call vecnor(tragl(1,idime),ndime,enor0,2)
         e%hleng(idime)=e%hnatu/enor0
      end do

      !Sort hleng: hleng(1)=max; hleng(ndime)=min
      do idime=1,e%ndime-1
         do jdime=idime+1,e%ndime
            if(e%hleng(jdime)>e%hleng(idime)) then  !swap
               h_tem       =e%hleng(jdime)
               e%hleng(jdime)=e%hleng(idime)
               e%hleng(idime)=h_tem
            end if
         end do
      end do

   end subroutine

   subroutine interpc(e,ndime,elvar,gpvar)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndime
      real(rp), intent(in)    :: elvar(ndime,e%pnode)
      real(rp), intent(out) :: gpvar(ndime)

      !integer(ip)             :: inode,idime

      gpvar = matmul(elvar,e%shacg(1:e%pnode))
      !do idime=1,e%ndime
      !   gpvar(idime)=0.0_rp
      !   do inode=1,e%pnode
      !      gpvar(idime) = gpvar(idime) + e%shacg(inode)*elvar(idime,inode)
      !   end do
      !end do
   end subroutine

   pure subroutine interpg(e,ndime,elvar,gpvar)
   !subroutine interpg(e,ndime,elvar,gpvar)
      implicit none
      class(FiniteElement), intent(in)  :: e
      !class(FiniteElement)              :: e
      integer(ip)         , intent(in)  :: ndime
      real(rp)            , intent(in)  :: elvar(ndime,e%pnode)
      real(rp)            , intent(out) :: gpvar(ndime)

      !integer(ip)             :: inode,idime

      gpvar = matmul(elvar,e%shape(1:e%pnode,e%igaus))
   end subroutine

   subroutine interpb(e,ndime,bovar,gpvar)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndime
      real(rp), intent(in)    :: bovar(ndime,e%pnodb)
      real(rp), intent(out) :: gpvar(ndime)

      !integer(ip)             :: inodb,idime

      gpvar = matmul(bovar,e%shapb(1:e%pnodb,e%igaub))
   end subroutine

   subroutine hessian(e,ndofn,elvar,hessvar)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndofn
      real(rp) :: elvar(ndofn,e%pnode), hessvar(ndofn,e%ntens)

      hessvar = matmul(elvar,transpose(e%hessi(:,1:e%pnode)))
   end subroutine

   subroutine gradient(e,ndofn,elvar,grvar)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndofn
      real(rp) :: elvar(ndofn,e%pnode), grvar(ndofn,e%ndime)

      grvar = matmul(elvar,transpose(e%cartd(:,1:e%pnode)))
   end subroutine

   subroutine gradientb(e,ndofn,elvar,grvar)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndofn
      real(rp) :: elvar(ndofn,e%pnode), grvar(ndofn,e%ndime)

      grvar = matmul(elvar,transpose(e%cartb(:,1:e%pnode)))
   end subroutine

   subroutine divergence(e,elvar,divvar)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: elvar(e%ndime,e%pnode), divvar

      integer(ip) :: idime,inode

      divvar = 0.0_rp
      do inode=1,e%pnode
         do idime=1,e%ndime
             divvar = divvar + e%cartd(idime,inode)*elvar(idime,inode)
          end do
      end do
   end subroutine

   subroutine isopar(e,ncomp,xloc,elval,xval)
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: ncomp
      real(rp) :: xloc(e%ndime), xval(ncomp),elval(ncomp,*)

      if (e%linea == 1) then
         call linear_isopar(e,ncomp,xloc,elval,xval)
      else
         call runend('FiniteElement: isopar only ready for linear elements')
      endif
   end subroutine

   subroutine linear_isopar(e,ncomp,xloc,elval,xval)
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: ncomp
      real(rp) :: xloc(e%ndime), xval(ncomp),elval(ncomp,*)


      real(rp) :: xjacm(ncomp,e%ndime)
      real(rp) :: b(ncomp),denom
      integer(ip) :: i

      if (e%ndime == 2) then
         xjacm(:,1) = elval(:,2)-elval(:,1);
         xjacm(:,2) = elval(:,3)-elval(:,1);
         b(:) = elval(:,1);
         xval(:) = xjacm(:,1)*xloc(1)+xjacm(:,2)*xloc(2)+b(:);

      else if (e%ndime == 3) then
         xjacm(:,1) = elval(:,2)-elval(:,1)
         xjacm(:,2) = elval(:,3)-elval(:,1)
         xjacm(:,3) = elval(:,4)-elval(:,1)

         b(:) = elval(:,1)
         xval(:) = xjacm(:,1)*xloc(1)+xjacm(:,2)*xloc(2)+xjacm(:,3)*xloc(3)+b(:);
      endif

   end subroutine

   subroutine isoparinv(e,xglob,xloc)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp) :: xloc(e%ndime), xglob(e%ndime)

      if (e%linea == 1) then
         call linear_isoparinv(e,xglob,xloc)
      else
         call nonlinear_isoparinv(e,xglob,xloc)
      endif


   end subroutine

   subroutine linear_isoparinv(e,xglob,xloc)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp) :: xloc(e%ndime), xglob(e%ndime)

      real(rp) :: xjacm(e%ndime,e%ndime)
      real(rp) :: b(e%ndime),denom,detjm, xjaci(e%ndime,e%ndime)
      integer(ip) :: i

      if (e%ndime == 2) then
         xjacm(:,1) = e%elcod(:,2)-e%elcod(:,1);
         xjacm(:,2) = e%elcod(:,3)-e%elcod(:,1);

         b(:) = xglob(:)-e%elcod(:,1);
         detjm = xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1);
         denom = 1.0_rp/detjm;
         xloc(1) = (b(1)*xjacm(2,2)-b(2)*xjacm(1,2))*denom;
         xloc(2) = (-b(1)*xjacm(2,1)+b(2)*xjacm(1,1))*denom;
      else if (e%ndime == 3) then
         xjacm(:,1) = e%elcod(:,2)-e%elcod(:,1)
         xjacm(:,2) = e%elcod(:,3)-e%elcod(:,1)
         xjacm(:,3) = e%elcod(:,4)-e%elcod(:,1)

         !inverse matrix
         xjaci(1,1)  = xjacm(2,2)*xjacm(3,3) - xjacm(3,2)*xjacm(2,3)
         xjaci(2,1)  =-xjacm(2,1)*xjacm(3,3) + xjacm(3,1)*xjacm(2,3)
         xjaci(3,1)  = xjacm(2,1)*xjacm(3,2) - xjacm(3,1)*xjacm(2,2)
         xjaci(2,2) =  xjacm(1,1)*xjacm(3,3) - xjacm(3,1)*xjacm(1,3)
         xjaci(3,2) = -xjacm(1,1)*xjacm(3,2) + xjacm(1,2)*xjacm(3,1)
         xjaci(3,3) =  xjacm(1,1)*xjacm(2,2) - xjacm(2,1)*xjacm(1,2)
         xjaci(1,2) = -xjacm(1,2)*xjacm(3,3) + xjacm(3,2)*xjacm(1,3)
         xjaci(1,3) =  xjacm(1,2)*xjacm(2,3) - xjacm(2,2)*xjacm(1,3)
         xjaci(2,3) = -xjacm(1,1)*xjacm(2,3) + xjacm(2,1)*xjacm(1,3)
         detjm = xjacm(1,1)*xjaci(1,1) + xjacm(1,2)*xjaci(2,1)+ xjacm(1,3)*xjaci(3,1)
         denom=1.0_rp/detjm
         xjaci = xjaci*denom

         b(:) = xglob(:)-e%elcod(:,1)
         forall (i = 1:e%ndime)
            xloc(i) = dot_product(xjaci(i,:),b(:))
         end forall
      endif

   end subroutine

   subroutine nonlinear_isoparinv(e,xglob,xloc)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp) :: xloc(e%ndime),xglob(e%ndime)

      real(rp) :: xjaci(e%ndime,e%ndime)
      real(rp) :: xjacm(e%ndime,e%ndime)
      real(rp) :: shape(e%pnode)
      real(rp) :: deriv(e%ndime,e%pnode)
      real(rp) :: heslo(e%ntens,e%pnode)
      real(rp) :: f_prev(e%ndime),delta(e%ndime),error,norm,detjm,cartd(e%ndime,e%pnode)
      integer(ip) :: idime,jdime,inode,maxit,i

      xloc = 0.0_rp
      error = 1.0_rp
      maxit = 6
      i = 1

      do while (i<=maxit .and. error>0.0000000001)
         call shafun(xloc,e%ndime,e%pnode,e%ntens,shape,deriv,heslo)
         call elmder(e%pnode,e%ndime,deriv,e%elcod,cartd,detjm,xjacm,xjaci)

         f_prev = 0.0_rp
         do idime=1,e%ndime
            do inode=1,e%pnode
               f_prev(idime) = f_prev(idime) - shape(inode)*e%elcod(idime,inode)
            end do
            f_prev(idime) = f_prev(idime) + xglob(idime)
         end do

         call mbvab1(delta,xjaci,f_prev,e%ndime,e%ndime,error,norm)
         xloc = xloc + delta
         i = i + 1
      end do
   end subroutine

   subroutine ForceClosedRule(e)
      use typre
      implicit none
      class(FiniteElement) :: e
      !If this is set, a close quadrature is enforced
      !This is necessary for integrating the RHS of a smoothing
      !also the LHS (vmass)

      e%kfl_ClosedRule = 1
   end subroutine

   subroutine DefaultRule(e)
      use typre
      implicit none
      class(FiniteElement) :: e
      !The default integration rule for the problem is used

      e%kfl_ClosedRule = 0
   end subroutine

   subroutine ComputeShafun(e,xloc,newshape)
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp) :: xloc(e%ndime),shape(e%pnode),deriv(e%ndime,e%pnode),heslo(e%ntens,e%pnode)
      real(rp), intent(out) :: newshape(e%pnode)

      call shafun(xloc,e%ndime,e%pnode,e%ntens,shape,deriv,heslo)
      newshape = shape

   end subroutine

   subroutine SetParticularGaussPoints(e,Memor,pgaus,posgp,weigp)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElement), target :: e
      type(MemoryMan) :: Memor

      integer(ip), intent(in) :: pgaus
      real(rp), intent(in)    :: posgp(:,:),weigp(:)
      integer(ip)             :: igaus

      if(pgaus > e%mParticularGaussPoints)then

         if(e%mParticularGaussPoints /= 0)then

            call Memor%realloc(e%mnode,pgaus,e%Particularshape,'e%Particularshape','SetParticularGaussPoints')
            call Memor%realloc(e%ndime,e%mnode,pgaus,e%Particularderiv,'e%Particularderiv','SetParticularGaussPoints')
            call Memor%realloc(e%ntens,e%mnode,pgaus,e%Particularheslo,'e%Particularheslo','SetParticularGaussPoints')
            call Memor%realloc(pgaus,e%Particularweigp,'e%Particularweigp','SetParticularGaussPoints')

         else

            call Memor%alloc(e%mnode,pgaus,e%Particularshape,'e%Particularshape','SetParticularGaussPoints')
            call Memor%alloc(e%ndime,e%mnode,pgaus,e%Particularderiv,'e%Particularderiv','SetParticularGaussPoints')
            call Memor%alloc(e%ntens,e%mnode,pgaus,e%Particularheslo,'e%Particularheslo','SetParticularGaussPoints')
            call Memor%alloc(pgaus,e%Particularweigp,'e%Particularweigp','SetParticularGaussPoints')
            call Memor%alloc(e%ndime,e%mnode,e%Particularcartd,'e%Particularcartd','SetParticularGaussPoints')

         end if

         e%mParticularGaussPoints=pgaus
         !to enforce the reallocation of the enriched elements when we change mParticularGaussPoints
         e%mPnode=0
      end if

      e%pgaus = pgaus
      do igaus = 1,e%pgaus
         e%igaus = igaus
         call shafun(posgp(1,e%igaus),e%ndime,e%pnode,e%ntens,e%Particularshape(1,e%igaus),e%Particularderiv(1,1,e%igaus),e%Particularheslo(1,1,e%igaus))
      end do

      e%Particularweigp(1:pgaus) = weigp(1:pgaus)


      e%shape => e%Particularshape
      e%deriv => e%Particularderiv
      e%heslo => e%Particularheslo
      e%weigp => e%Particularweigp

      e%ielty0 = -1


   end subroutine

   subroutine SetEnrichedElement(e,Memor,newpnode,kshape,ngauss_minus,ngauss_plus,nodeStatus)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElement), target :: e
      type(MemoryMan) :: Memor
      integer(ip), intent(in) :: newpnode,ngauss_minus,ngauss_plus,nodeStatus(e%pnode)
      real(rp), intent(in)    :: kshape(e%pnode)
      integer(ip) :: inode,gauss_points(2,e%pnode)

      if (newpnode > e%mnode) then
         if (e%mPnode== 0) then
            call Memor%realloc(newpnode,e%pgaus,e%Particularshape,'e%Particularshape','SetParticularGaussPoints')
            call Memor%realloc(e%ndime,newpnode,e%pgaus,e%Particularderiv,'e%Particularderiv','SetParticularGaussPoints')
            call Memor%realloc(e%ntens,newpnode,e%pgaus,e%Particularheslo,'e%Particularheslo','SetParticularGaussPoints')
            call Memor%realloc(e%pgaus,e%Particularweigp,'e%Particularweigp','SetParticularGaussPoints')
            call Memor%realloc(e%ndime,newpnode,e%Particularcartd,'e%Particularcartd','SetParticularGaussPoints')
         end if
         e%mPnode=newpnode
      end if

      do inode = 1,e%pnode
         if (nodeStatus(inode) < 0) then
            gauss_points(1,inode) = ngauss_minus+1
            gauss_points(2,inode) = ngauss_minus+ngauss_plus
         else
            gauss_points(1,inode) = 1
            gauss_points(2,inode) = ngauss_minus
         endif
      enddo

      e%Particularshape(e%pnode+1,:)=0.0_rp
      e%Particularderiv(1:e%ndime,e%pnode+1,:)=0.0_rp
      e%Particularheslo(1:e%ntens,e%pnode+1,:)=0.0_rp

      !Enriched shape function contruction
      do inode = 1,e%pnode
         e%Particularshape(e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) = e%Particularshape(e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) + &
               e%Particularshape(inode,gauss_points(1,inode):gauss_points(2,inode))*kshape(inode)

         e%Particularderiv(1:e%ndime,e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) = e%Particularderiv(1:e%ndime,e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) + &
               e%Particularderiv(1:e%ndime,inode,gauss_points(1,inode):gauss_points(2,inode))*kshape(inode)

         e%Particularheslo(1:e%ntens,e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) = e%Particularheslo(1:e%ntens,e%pnode+1,gauss_points(1,inode):gauss_points(2,inode)) + &
               e%Particularheslo(1:e%ntens,inode,gauss_points(1,inode):gauss_points(2,inode))*kshape(inode)
      enddo

      e%shape => e%Particularshape
      e%deriv => e%Particularderiv
      e%heslo => e%Particularheslo
      e%weigp => e%Particularweigp
      e%ielty0 = -1
   end subroutine


   subroutine DeallocParticularGaussPoints(e,Memor)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElement) :: e
      type(MemoryMan) :: Memor

      call Memor%dealloc(size(e%Particularshape,1),e%mParticularGaussPoints,e%Particularshape,'e%Particularshape','DeallocParticularGaussPoints')
      call Memor%dealloc(e%ndime,size(e%Particularderiv,2),e%mParticularGaussPoints,e%Particularderiv,'e%Particularderiv','DeallocParticularGaussPoints')
      call Memor%dealloc(e%ntens,size(e%Particularheslo,2),e%mParticularGaussPoints,e%Particularheslo,'e%Particularheslo','DeallocParticularGaussPoints')
      call Memor%dealloc(e%mParticularGaussPoints,e%Particularweigp,'e%Particularweigp','DeallocParticularGaussPoints')
      call Memor%dealloc(e%ndime,size(e%Particularcartd,2),e%Particularcartd,'e%Particularcartd','DeallocParticularGaussPoints')
   end subroutine

   subroutine GaussToNodes(e,GaussArray,NodalArray)
      class(FiniteElement) :: e
      real(rp) :: GaussArray(:,:), NodalArray(:,:)
      real(rp) :: mass(e%pnode,e%pnode)
      real(rp) :: transposeNodal(size(NodalArray,2),size(NodalArray,1))
      integer(ip) :: ipiv(e%pnode),info
      integer(ip) :: igaus,inode,jnode,ndofn

      ndofn = size(GaussArray,1)
      mass = 0.0_rp
      NodalArray = 0.0_rp
      !Build local mass matrix
      do igaus = 1,e%pgaus
         e%igaus = igaus
         do inode = 1,e%pnode
            do jnode = 1,e%pnode
               mass(inode,jnode) = mass(inode,jnode) + e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            enddo
            NodalArray(:,inode) = NodalArray(:,inode) + e%shape(inode,e%igaus)*GaussArray(:,e%igaus)
         enddo
      enddo
      transposeNodal = transpose(NodalArray)
      CALL DGESV( e%pnode, ndofn, mass, e%pnode, ipiv, transposeNodal, e%pnode, info )
      NodalArray = transpose(transposeNodal)
   end subroutine

   subroutine GaussToNodesB(e,GaussArray,NodalArray)
      class(FiniteElement) :: e
      real(rp) :: GaussArray(:,:), NodalArray(:,:)
      real(rp) :: mass(e%pnodb,e%pnodb)
      real(rp) :: transposeNodal(size(NodalArray,2),size(NodalArray,1))
      integer(ip) :: ipiv(e%pnodb),info
      integer(ip) :: igaub,inodb,jnodb,ndofn

      ndofn = size(GaussArray,1)
      mass = 0.0_rp
      NodalArray = 0.0_rp
      !Build local mass matrix
      do igaub = 1,e%pgaub
         e%igaub = igaub
         do inodb = 1,e%pnodb
            do jnodb = 1,e%pnodb
               mass(inodb,jnodb) = mass(inodb,jnodb) + e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)
            enddo
            NodalArray(:,inodb) = NodalArray(:,inodb) + e%shapb(inodb,e%igaub)*GaussArray(:,e%igaub)
         enddo
      enddo
      transposeNodal = transpose(NodalArray)
      CALL DGESV( e%pnodb, ndofn, mass, e%pnodb, ipiv, transposeNodal, e%pnodb, info )
      NodalArray = transpose(transposeNodal)
   end subroutine

   subroutine GetWeightFactor(e,weightfactor)
      use typre
      Implicit none
      class(FiniteElement) :: e
      real(rp), intent(out) :: weightfactor

      weightfactor = 1.0-real(count(e%lnods(1:e%pnode)>e%npoinLocal))/real(e%pnode)

   end subroutine


end module

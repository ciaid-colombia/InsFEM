module Mod_nsi_Refine
   use typre
   use Mod_Element
   use Mod_NavierStokes
   use Mod_AdaptiveInterface
   implicit none
   private
   public ElementNavierStokesTransferer

   type, extends(AdaptiveRefinerElementTransferer) :: ElementNavierStokesTransferer
      type(r1p), allocatable :: &
            prsgs(:)                            !Pressure subgrid scales
      type(r3p), allocatable :: &
            vesgs(:)                            !Velocity subgrid scales
      class(NavierStokesProblem), pointer :: NSP => NULL()
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ncsgs

contains
      procedure :: CopyElement => ENT_CopyElement
      procedure :: NodallyInterpolateElement => ENT_NodallyInterpolateElement
      procedure :: NodallyInterpolateElementFromOldBuffer => ENT_NodallyInterpolateElementFromOldBuffer
      

      procedure :: OldElementToBuffer => ENT_OldElementToBuffer
      procedure :: GetOldElementBufferSize => ENT_GetOldElementBufferSize

      procedure :: GetNewElementBufferSize => ENT_GetNewElementBufferSize
      procedure :: BufferToNewElement => ENT_BufferToNewElement
      procedure :: NewElementToBuffer => ENT_NewElementToBuffer

      procedure :: SetProblem
      procedure :: AllocateArrays

   end type


contains

   subroutine SetProblem(a,NSP)
      implicit none
      class(ElementNavierStokesTransferer) :: a
      class(NavierStokesProblem), target :: NSP

      a%NSP => NSP
      
      a%ncsgs=1                  !Number of components in time for accuracy
      if(NSP%kfl_tacsg>0) a%ncsgs=2

   end subroutine

   subroutine AllocateArrays(a)
      implicit none
      class(ElementNavierStokesTransferer) :: a

      integer(ip) :: ielem,nelem,ndime,pgaus,pnode
      integer(ip) ::  prsgs_coun,vesgs_coun
      
      call a%NSP%Mesh%GetNelem(nelem)
      call a%NSP%Mesh%GetNdime(ndime)
    
      prsgs_coun = 0
      vesgs_coun = 0

      !For explicit time integration scheme:
      !position 1: stage increment
      !position 2: stage contribution
      !position 3: previous time step 
      
      call a%NSP%Memor%alloc(nelem,a%prsgs,'prsgs','nsi_refine')
      call a%NSP%Memor%alloc(nelem,a%vesgs,'vesgs','nsi_refine')
      do ielem=1,nelem
         call a%NSP%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%prsgs(ielem)%a(pgaus))
         allocate(a%vesgs(ielem)%a(ndime,a%ncsgs,pgaus))
         a%prsgs(ielem)%a = 0.0_rp
         a%vesgs(ielem)%a = 0.0_rp
         prsgs_coun = prsgs_coun + pgaus
         vesgs_coun = vesgs_coun + ndime*a%ncsgs*pgaus
      end do
      call a%NSP%Memor%allocObj(0,'prsgs%a','nsi_refine',prsgs_coun*rp)
      call a%NSP%Memor%allocObj(0,'vesgs%a','nsi_refine',vesgs_coun*rp)

      !Allocate element
      call a%NSP%Mesh%ElementAlloc(a%e,a%NSP%Memor,'DefaultRule','nsi_pr_elmope')
    
   end subroutine

   subroutine ENT_CopyElement(a,oldielem,newielem)
      implicit none
      class(ElementNavierStokesTransferer), target :: a
      integer(ip) :: oldielem,newielem

      a%prsgs(newielem)%a = a%NSP%prsgs(oldielem)%a 
      a%vesgs(newielem)%a = a%NSP%vesgs(oldielem)%a 

   end subroutine

   subroutine ENT_NodallyInterpolateElement(a,Refiner,oldelem,newelem)
      implicit none
      class(ElementNavierStokesTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) :: oldelem, newelem

      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)
      real(rp), pointer :: OldNodalArray(:,:) => NULL()
      integer(ip) :: ndime,pgaus,pnode

      call a%NSP%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      call a%NSP%Mesh%GetNdime(ndime)
      ndofn = ndime*a%ncsgs+1
      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 

      call Refiner%InterpolateElement(ndofn,a%NSP%NodalArrayDataForRefinement(oldelem)%a,NewNodalArray)

      call ENT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)
   end subroutine

   subroutine ENT_NodalArrayToSubscales(a,newelem,NewNodalArray)
      implicit none
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: newelem
      real(rp)  :: NewNodalArray(:,:)
      integer(ip) :: pgaus,igaus,ndofn,ndime,icomp

      real(rp)  :: gpNewNodalArray(200)

      call a%NSP%Mesh%GetNdime(ndime)
      call a%NSP%Mesh%ElementLoad(newelem,a%e)
   
      ndofn  = size(NewNodalArray,1) 
      do igaus = 1,a%e%pgaus
         a%e%igaus = igaus
         call a%e%interpg(ndofn,NewNodalArray,gpNewNodalArray(1:ndofn))
         do icomp = 1, a%ncsgs
             a%vesgs(newelem)%a(:,icomp,a%e%igaus) = gpNewNodalArray(1+(icomp-1)*ndime:icomp*ndime)
         enddo
         a%prsgs(newelem)%a(a%e%igaus) = gpNewNodalArray(a%ncsgs*ndime+1)
      enddo

   end subroutine

   subroutine ENT_NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
      implicit none
      class(ElementNavierStokesTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) ::  newelem
      real(rp), target :: OldBuffer(*)
      
      real(rp), pointer :: auxOldNodalArray(:,:) => NULL()

      integer(ip) :: ndofn,ndime,pnode,pgaus
      real(rp), allocatable :: NewNodalArray(:,:)

      call a%NSP%Mesh%GetNdime(ndime)
      call a%NSP%Mesh%GetElemArraySize(newelem,pnode,pgaus)

      ndofn = ndime*a%ncsgs+1
      auxOldNodalArray(1:ndofn,1:pnode) => OldBuffer(1:ndofn*pnode)

      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 
 
      call Refiner%InterpolateElement(ndofn,auxOldNodalArray,NewNodalArray)

      call ENT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)

   end subroutine

   subroutine ENT_OldElementToBuffer(a,ielem,Buffer)
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)
      

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

     
      icount = 0
      auxpointer(1:size(a%NSP%vesgs(ielem)%a)) => a%NSP%vesgs(ielem)%a
      icount = icount + size(a%NSP%vesgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%NSP%prsgs(ielem)%a)) => a%NSP%prsgs(ielem)%a
      icount = icount + size(a%NSP%prsgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer
     
   end subroutine

   subroutine ENT_NewElementToBuffer(a,ielem,Buffer)
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%vesgs(ielem)%a)) => a%vesgs(ielem)%a
      icount = icount + size(a%vesgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%prsgs(ielem)%a)) => a%prsgs(ielem)%a
      icount = icount + size(a%prsgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer
     
   end subroutine

   subroutine ENT_BufferToNewElement(a,ielem,Buffer)
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%vesgs(ielem)%a)) => a%vesgs(ielem)%a
      icount = icount + size(a%vesgs(ielem)%a)
      auxpointer = Buffer(1:icount)

      oldicount = icount

      auxpointer(1:size(a%prsgs(ielem)%a)) => a%prsgs(ielem)%a
      icount = icount + size(a%prsgs(ielem)%a)
      auxpointer =  Buffer(oldicount+1:icount)
     
   end subroutine

   subroutine ENT_GetOldElementBufferSize(a,ielem,BufferSize)
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount

      icount = 0
      icount = icount + size(a%NSP%vesgs(ielem)%a)
      icount = icount + size(a%NSP%prsgs(ielem)%a)
     
      BufferSize = icount
      
   end subroutine

   subroutine ENT_GetNewElementBufferSize(a,ielem,BufferSize)
      class(ElementNavierStokesTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount
      icount = 0
      icount = icount + size(a%vesgs(ielem)%a)
      icount = icount + size(a%prsgs(ielem)%a)
     
      BufferSize = icount
   end subroutine

end module

subroutine nsi_Prerefine(a)
   use typre
   use Mod_NavierStokes
   use Mod_phpRefineArrays
   use Mod_nsi_Refine
   use Mod_Element
   implicit none
   class(NavierStokesProblem) :: a

   integer(ip) :: nelem,ielem,ndime,icomp,ncsgs
   class(FiniteElement), pointer :: e => NULL()
   real(rp), target :: rwa_GaussArray(500)
   real(rp), pointer :: GaussArray(:,:) => NULL()

   if(a%kfl_trasg/=0) then

      ncsgs=2

      !Things to be done prior to refinement
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_Refine')
   
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      !Allocate data array
      call a%Memor%alloc(nelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','nsi_refine')
   
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         allocate(a%NodalArrayDataForRefinement(ielem)%a(ncsgs*ndime+1,e%pnode))
         GaussArray(1:ncsgs*ndime+1,1:e%pgaus) => rwa_GaussArray(1:(ncsgs*ndime+1)*e%pgaus)
         do icomp = 1,ncsgs
             GaussArray(1+(icomp-1)*ndime:icomp*ndime,1:e%pgaus) = a%vesgs(ielem)%a(:,icomp,1:e%pgaus)
         enddo
         GaussArray(ncsgs*ndime+1,:) = a%prsgs(ielem)%a
         call e%GaussToNodes(GaussArray,a%NodalArrayDataForRefinement(ielem)%a)
      enddo
   
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsi_Refine')

    end if

end subroutine

subroutine nsi_refine(a,itask)
   use typre
   use Mod_NavierStokes
   use Mod_phpRefineArrays
   use Mod_nsi_Refine
   
   implicit none
   class(NavierStokesProblem) :: a
   type(ElementNavierStokesTransferer) ::  ENT
   character(6) :: itask

   integer(ip) :: oldnpoin,oldnboun,newnpoin,newnboun,oldnelem,newnelem,icomp,ndime,ielem,iboun,pgaus
   integer(ip) :: visc_coun,vesgs_coun,prsgs_coun,pnode,ncsgs
   integer(ip) :: imat=1, auxsize
   real(rp), allocatable    :: auxveloc(:,:,:),auxpress(:,:),auxTausmo(:,:),auxbtraction(:,:),auxuavg(:,:),auxutrap(:,:)
   real(rp), allocatable    :: auxveloc_cp(:,:),auxpress_cp(:)
   integer(ip), allocatable :: auxkfl_fixrs(:),auxkfl_bours(:),auxkfl_weakbc(:)
   integer(ip), pointer     :: BoundaryNewToOld(:) => NULL(),iauxBoundaryMatch1(:) => NULL()
   type(i1p), pointer       :: iauxBoundaryMatch2(:) => NULL()
   real(rp), allocatable    :: auxturbi(:)
   real(rp), allocatable    :: auxtimav(:,:,:)
   real(rp), allocatable    :: auxtimap(:,:)
   real(rp), allocatable    :: auxrmsv(:,:)
   real(rp), allocatable    :: auxrmsp(:),auxTimaTractions(:,:,:)
   real(rp), allocatable    :: auxrepro(:,:),auxgrprj(:,:,:)

   !This subroutine transforms the tempe arrays from one mesh to the other
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%veloc,2)
   oldnboun = size(a%kfl_bours)

   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNboun(newnboun)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(newnelem)

   !Veloc
   call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxveloc,'veloc','nsi_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
      endif
   enddo
   call move_alloc(auxveloc,a%veloc)
   call a%Memor%deallocObj(0,'veloc','nsi_Refine',rp*ndime*oldnpoin*a%ncomp)

   !Pressure
   call a%Memor%alloc(newnpoin,a%ncomp,auxpress,'press','nsi_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      endif
   enddo
   call move_alloc(auxpress,a%press)
   call a%Memor%deallocObj(0,'press','nsi_Refine',rp*oldnpoin*a%ncomp)

   if(a%kfl_docoupconv)  then
       !Veloc
       call a%Memor%alloc(ndime,newnpoin,auxveloc_cp,'veloc_cp','nsi_Refine')
       do icomp = 1,a%ncomp
           if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(ndime,a%veloc_cp(:,:),auxveloc_cp(:,:))
           elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(ndime,a%veloc_cp(:,:),auxveloc_cp(:,:))
           endif
       enddo
       call move_alloc(auxveloc_cp,a%veloc_cp)
       call a%Memor%deallocObj(0,'veloc_cp','nsi_Refine',rp*ndime*oldnpoin)

       !Pressure
       call a%Memor%alloc(newnpoin,auxpress_cp,'press_cp','nsi_Refine')
       if (itask == 'Refine') then
           call a%Refiner%UpdateVariable(1_ip,a%press_cp(:),auxpress_cp(:))
       elseif (itask == 'Rebala') then
           call a%Refiner%RebalanceVariable(1_ip,a%press_cp(:),auxpress_cp(:))
       endif
       call move_alloc(auxpress_cp,a%press_cp)
       call a%Memor%deallocObj(0,'press_cp','nsi_Refine',rp*oldnpoin)

   endif

   !Btraction
   call a%Memor%alloc(ndime,newnpoin,auxbtraction,'btraction','nsi_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%btraction(:,:),auxbtraction(:,:))
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%btraction(:,:),auxbtraction(:,:))
   endif
   call move_alloc(auxbtraction,a%btraction)
   call a%Memor%deallocObj(0,'btraction','nsi_Refine',rp*oldnpoin*ndime)

   !Residual projection
   if (a%kfl_repro >= 1) then
      auxsize = size(a%repro,1)
      call a%Memor%alloc(auxsize,newnpoin,auxrepro,'repro','nsi_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      endif
      call move_alloc(auxrepro,a%repro)
      call a%Memor%deallocObj(0,'repro','nsi_Refine',rp*oldnpoin*auxsize)
   endif

   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      auxsize = size(a%grprj,1)
      call a%Memor%alloc(auxsize,ndime,newnpoin,auxgrprj,'grprj','nsi_Refine')
      do icomp = 1,auxsize
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         endif   
      enddo
      call move_alloc(auxgrprj,a%grprj)
      call a%Memor%deallocObj(0,'grprj','nsi_Refine',rp*auxsize*ndime*oldnpoin)
   end if

   !Tau Smoothing
   if (a%kfl_tausm >= 1) then
      call a%Memor%alloc(size(a%Tausmo,1),newnpoin,auxTausmo,'Tausmo','nsi_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(size(a%Tausmo,1),a%Tausmo,auxTausmo(:,:))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(size(a%Tausmo,1),a%Tausmo,auxTausmo(:,:))
      endif
      call move_alloc(auxTausmo,a%Tausmo)
      call a%Memor%deallocObj(0,'Tausmo','nsi_Refine',rp*oldnpoin*size(a%Tausmo,1))
   endif

   !Subgrid scales
   if(a%kfl_trasg/=0) then
      call ENT%SetProblem(a) 
      call ENT%AllocateArrays 

      if (itask == 'Refine') then
         !Don't call us, we'll call YOU!
         call a%Refiner%ElementalRefine(ENT)
      elseif (itask == 'Rebala') then
         call a%Refiner%ElementalRebalance(ENT)
      endif

      oldnelem = size(a%prsgs,1)
      vesgs_coun = 0
      prsgs_coun = 0
      do ielem=1,oldnelem
         pgaus = size(a%vesgs(ielem)%a,3)
         deallocate(a%vesgs(ielem)%a)
         deallocate(a%prsgs(ielem)%a)
         vesgs_coun = vesgs_coun + ndime*ENT%ncsgs*pgaus
         prsgs_coun = prsgs_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'vesgs%a','nsi_refine',vesgs_coun*rp)
      call a%Memor%deallocObj(0,'prsgs%a','nsi_refine',prsgs_coun*rp)
      call a%Memor%dealloc(oldnelem,a%vesgs,'vesgs','nsi_refine')
      call a%Memor%dealloc(oldnelem,a%prsgs,'prsgs','nsi_refine')

      call move_alloc(ENT%vesgs,a%vesgs)
      call move_alloc(ENT%prsgs,a%prsgs)

      call a%Mesh%ElementDeAlloc(ENT%e,a%Memor,'DefaultRule','nsi_Refine')
      
      if (itask == 'Refine') then
         !DeAllocate data array
         call a%Memor%dealloc(oldnelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','nsi_refine')
      end if 

   endif
   !Postprocess residual
   if (a%npp_stepi(18) /= 0) then
      !call runend('Postprocess residual: Adaptive refinement not ready to transfer info on gauss points')
   endif

   !Viscosity law
   if(a%MatProp(imat)%lawvi/=0)then
      if ((a%kfl_repro >= 1).or.(a%npp_stepi(5)==1)) then
         call runend('Viscosity law: not ready for adaptivity')
      endif
   endif

   !Dissipation postprocess
   if (a%kfl_dispa /=0) then
      call a%Memor%realloc(newnpoin,a%Dissipation,'dissipation','nsi_memall')
   endif

   !Vorticity postprocess
   if (a%npp_stepi(10)>0) then
      call a%Memor%realloc(size(a%vorti,1),newnpoin,a%vorti,'vorti','nsi_memall')
   endif

   !Vorticity postprocess
   if (a%npp_stepi(11)>0) then
      !Deallocate
      visc_coun = 0
      do ielem=1,size(a%Divergence,1)
         pgaus = size(a%Divergence(ielem)%a,1)
         deallocate(a%Divergence(ielem)%a)
         visc_coun = visc_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'Divergence%a','nsi_turnof',visc_coun*rp)
      call a%Memor%dealloc(size(a%Divergence,1),a%Divergence,'Divergence','nsi_turnof')

      !Allocate again
      visc_coun = 0
      call a%Memor%alloc(newnelem,a%Divergence,'Divergence','nsi_memall')
      do ielem=1,newnelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%Divergence(ielem)%a(pgaus))
         visc_coun = visc_coun + pgaus
      end do
      call a%Memor%allocObj(0,'Divergence%a','nsi_memall',visc_coun*rp)
   endif

       !Statistics of fields
   if (a%npp_stepi(24) /= 0) then
      call a%Memor%alloc(ndime,newnpoin,5,auxtimav,'timav','nsi_Refine')
      call a%Memor%alloc(newnpoin,     3,auxtimap,'timap','nsi_Refine')
      call a%Memor%alloc(ndime,newnpoin,auxrmsv,'rmsv','nsi_Refine')
      call a%Memor%alloc(newnpoin,auxrmsp,'rmsp','nsi_Refine')
      call a%Memor%alloc(newnpoin,auxturbi,'turbi','nsi_Refine')
      call a%Memor%alloc(ndime,newnpoin,2,auxTimaTractions,'auxTimaTractions','nsi_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime,a%timav(:,:,1),auxtimav(:,:,1))
         call a%Refiner%UpdateVariable(ndime,a%timav(:,:,2),auxtimav(:,:,2))
         call a%Refiner%UpdateVariable(1_ip,a%timap(:,1),auxtimap(:,1))
         call a%Refiner%UpdateVariable(1_ip,a%timap(:,2),auxtimap(:,2))
         call a%Refiner%UpdateVariable(ndime,a%rmsv,auxrmsv)
         call a%Refiner%UpdateVariable(1_ip,a%rmsp,auxrmsp)
         call a%Refiner%UpdateVariable(1_ip,a%turbi,auxturbi)
         call a%Refiner%UpdateVariable(ndime,a%TimaTractions(:,:,1),auxTimaTractions(:,:,1))
         call a%Refiner%UpdateVariable(ndime,a%TimaTractions(:,:,2),auxTimaTractions(:,:,2))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime,a%timav(:,:,1),auxtimav(:,:,1))
         call a%Refiner%RebalanceVariable(ndime,a%timav(:,:,2),auxtimav(:,:,2))
         call a%Refiner%RebalanceVariable(1_ip,a%timap(:,1),auxtimap(:,1))
         call a%Refiner%RebalanceVariable(1_ip,a%timap(:,2),auxtimap(:,2))
         call a%Refiner%RebalanceVariable(ndime,a%rmsv,auxrmsv)
         call a%Refiner%RebalanceVariable(1_ip,a%rmsp,auxrmsp)
         call a%Refiner%RebalanceVariable(1_ip,a%turbi,auxturbi)
         call a%Refiner%RebalanceVariable(ndime,a%TimaTractions(:,:,1),auxTimaTractions(:,:,1))
         call a%Refiner%RebalanceVariable(ndime,a%TimaTractions(:,:,2),auxTimaTractions(:,:,2))
      endif
      call move_alloc(auxtimav,a%timav)
      call move_alloc(auxtimap,a%timap)
      call move_alloc(auxrmsv,a%rmsv)
      call move_alloc(auxrmsp,a%rmsp)
      call move_alloc(auxturbi,a%turbi)
      call move_alloc(auxTimaTractions,a%TimaTractions)
      call a%Memor%deallocObj(0,'timav','nsi_Refine',rp*ndime*oldnpoin*2)
      call a%Memor%deallocObj(0,'timap','nsi_Refine',rp*oldnpoin*2)
      call a%Memor%deallocObj(0,'rmsv','nsi_Refine',rp*ndime*oldnpoin)
      call a%Memor%deallocObj(0,'rmsp','nsi_Refine',rp*oldnpoin)
      call a%Memor%deallocObj(0,'turbi','nsi_Refine',rp*oldnpoin)
       call a%Memor%deallocObj(0,'timatractions','nsi_Refine',rp*ndime*oldnpoin*2)
   endif

   !Split_OSS
   if (a%kfl_repro == 2) then
      !Deallocate
      visc_coun = 0
      do ielem=1,size(a%residualGraP2,1)
         pgaus = size(a%residualGraP2(ielem)%a(1,:),1)
         deallocate(a%residualGraP2(ielem)%a)
         visc_coun = visc_coun + ndime*pgaus
      end do
      call a%Memor%deallocObj(0,'residualGraP2%a','nsi_refine',visc_coun*rp)
      call a%Memor%dealloc(size(a%residualGraP2,1),a%residualGraP2,'Divergence','nsi_refine')

      !Allocate again
      visc_coun = 0
      call a%Memor%alloc(newnelem,a%residualGraP2,'residualGraP2','nsi_refine')
      do ielem=1,newnelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%residualGraP2(ielem)%a(ndime,pgaus))
         visc_coun = visc_coun + ndime*pgaus
      end do
      call a%Memor%allocObj(0,'residualGraP2%a','nsi_refine',visc_coun*rp)
   endif

   !Gauss point dissipation
   if (a%npp_stepi(20) /= 0) then
      call runend('Gauss point dissipation: Not ready for adaptivity')
   endif

   !Arrays for conditions on boundaries

   !kfl_fixrs
   call a%Memor%alloc(newnpoin,auxkfl_fixrs,'kfl_fixrs','nsi_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%kfl_fixrs,auxkfl_fixrs,'minab')
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%kfl_fixrs,auxkfl_fixrs)
   endif
   call move_alloc(auxkfl_fixrs,a%kfl_fixrs)
   call a%Memor%deallocObj(0,'kfl_fixrs','nsi_Refine',ip*oldnpoin)

   !kfl_bours
   if (itask == 'Refine') then
      call a%Mesh%GetBoundaryNewToOld(BoundaryNewToOld)
      !kfl_bours
      call a%Memor%alloc(newnboun,auxkfl_bours,'kfl_bours','nsi_Refine')
      do iboun = 1,newnboun
         if (BoundaryNewToOld(iboun) > 0) then
            auxkfl_bours(iboun) = a%kfl_bours(BoundaryNewToOld(iboun))
         else
            auxkfl_bours(iboun) = 0
         endif
      enddo
      call move_alloc(auxkfl_bours,a%kfl_bours)
      call a%Memor%deallocObj(0,'kfl_bours','php_Refine',ip*oldnboun)

   elseif (itask == 'Rebala') then

      call a%Mesh%GetRebalanceInverseBoundaryMatches(iauxBoundaryMatch1,iauxBoundaryMatch2)
      call a%Memor%alloc(newnboun,auxkfl_bours,'kfl_bours','php_Refine')
      call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,1,a%kfl_bours,auxkfl_bours)
      call move_alloc(auxkfl_bours,a%kfl_bours)
      call a%Memor%deallocObj(0,'kfl_bours','php_Refine',oldnboun*ip)
   endif

   if (a%kfl_StabilizeFreeSurface /= 0) then
      call php_RefineArrays_b(a,itask,a%StabilizeFreeSurface_velocityGradient,'StabilizeFreeSurface_VelocityGradient')
      call php_RefineArrays_b(a,itask,a%StabilizeFreeSurface_PressureGradientAndForce,'StabilizeFreeSurface_PressureGradientAndForce')
   endif

end subroutine

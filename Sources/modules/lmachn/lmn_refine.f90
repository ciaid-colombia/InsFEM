module Mod_lmn_Refine
   use typre
   use Mod_Element
   use Mod_AdaptiveInterface
   use Mod_LowMach
   implicit none
   private
   public ElementLowMachTransferer

   type, extends(AdaptiveRefinerElementTransferer) :: ElementLowMachTransferer
      type(r1p), allocatable :: &
            prsgs(:)                            !Pressure subgrid scales
      type(r2p), allocatable :: &
            tesgs(:)                            !Temperature subgrid scales

      type(r3p), allocatable :: &
            vesgs(:)                            !Velocity subgrid scales
      class(LowMachProblem), pointer :: LMP => NULL()
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ncsgs

contains
      procedure :: CopyElement => ELT_CopyElement
      procedure :: NodallyInterpolateElement => ELT_NodallyInterpolateElement
      procedure :: NodallyInterpolateElementFromOldBuffer => ELT_NodallyInterpolateElementFromOldBuffer
      

      procedure :: OldElementToBuffer => ELT_OldElementToBuffer
      procedure :: GetOldElementBufferSize => ELT_GetOldElementBufferSize

      procedure :: GetNewElementBufferSize => ELT_GetNewElementBufferSize
      procedure :: BufferToNewElement => ELT_BufferToNewElement
      procedure :: NewElementToBuffer => ELT_NewElementToBuffer

      procedure :: SetProblem
      procedure :: AllocateArrays

   end type


contains

   subroutine SetProblem(a,LMP)
      implicit none
      class(ElementLowMachTransferer) :: a
      class(LowMachProblem), target :: LMP

      a%LMP => LMP
      
      a%ncsgs=1                  !Number of components in time for accuracy
      if(LMP%kfl_tacsg>0) a%ncsgs=LMP%ncomp-1

   end subroutine

   subroutine AllocateArrays(a)
      implicit none
      class(ElementLowMachTransferer) :: a

      integer(ip) :: ielem,nelem,ndime,pgaus,pnode
      integer(ip) ::  prsgs_coun,vesgs_coun,tesgs_coun
      
      call a%LMP%Mesh%GetNelem(nelem)
      call a%LMP%Mesh%GetNdime(ndime)
    
      prsgs_coun = 0
      vesgs_coun = 0
      tesgs_coun = 0

      !For explicit time integration scheme:
      !position 1: stage increment
      !position 2: stage contribution
      !position 3: previous time step 
      
      call a%LMP%Memor%alloc(nelem,a%prsgs,'prsgs','lmn_refine')
      call a%LMP%Memor%alloc(nelem,a%vesgs,'vesgs','lmn_refine')
      call a%LMP%Memor%alloc(nelem,a%tesgs,'tesgs','lmn_refine')
      do ielem=1,nelem
         call a%LMP%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%prsgs(ielem)%a(pgaus))
         allocate(a%vesgs(ielem)%a(ndime,a%ncsgs,pgaus))
         allocate(a%tesgs(ielem)%a(a%ncsgs,pgaus))
         a%prsgs(ielem)%a = 0.0_rp
         a%vesgs(ielem)%a = 0.0_rp
         a%tesgs(ielem)%a = 0.0_rp
         prsgs_coun = prsgs_coun + pgaus
         vesgs_coun = vesgs_coun + ndime*a%ncsgs*pgaus
         tesgs_coun = tesgs_coun + a%ncsgs*pgaus
      end do
      call a%LMP%Memor%allocObj(0,'prsgs%a','lmn_refine',prsgs_coun*rp)
      call a%LMP%Memor%allocObj(0,'vesgs%a','lmn_refine',vesgs_coun*rp)
      call a%LMP%Memor%allocObj(0,'tesgs%a','lmn_refine',tesgs_coun*rp)

      !Allocate element
      call a%LMP%Mesh%ElementAlloc(a%e,a%LMP%Memor,'DefaultRule','lmn_pr_elmope')
    
   end subroutine

   subroutine ELT_CopyElement(a,oldielem,newielem)
      implicit none
      class(ElementLowMachTransferer), target :: a
      integer(ip) :: oldielem,newielem

      a%prsgs(newielem)%a = a%LMP%prsgs(oldielem)%a 
      a%vesgs(newielem)%a = a%LMP%vesgs(oldielem)%a 
      a%tesgs(newielem)%a = a%LMP%tesgs(oldielem)%a 

   end subroutine

   subroutine ELT_NodallyInterpolateElement(a,Refiner,oldelem,newelem)
      implicit none
      class(ElementLowMachTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) :: oldelem, newelem

      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)
      real(rp), pointer :: OldNodalArray(:,:) => NULL()
      integer(ip) :: ndime,pgaus,pnode

      call a%LMP%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      call a%LMP%Mesh%GetNdime(ndime)
      ndofn = (ndime+1)*a%ncsgs+1
      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 

      call Refiner%InterpolateElement(ndofn,a%LMP%NodalArrayDataForRefinement(oldelem)%a,NewNodalArray)

      call ELT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)
   end subroutine

   subroutine ELT_NodalArrayToSubscales(a,newelem,NewNodalArray)
      implicit none
      class(ElementLowMachTransferer) :: a
      integer(ip) :: newelem
      real(rp)  :: NewNodalArray(:,:)
      integer(ip) :: pgaus,igaus,ndofn,ndime,icomp

      real(rp)  :: gpNewNodalArray(200)

      call a%LMP%Mesh%GetNdime(ndime)
      call a%LMP%Mesh%ElementLoad(newelem,a%e)
   
      ndofn  = size(NewNodalArray,1) 
      do igaus = 1,a%e%pgaus
         a%e%igaus = igaus
         call a%e%interpg(ndofn,NewNodalArray,gpNewNodalArray(1:ndofn))
         do icomp = 1, a%ncsgs
             a%vesgs(newelem)%a(:,icomp,a%e%igaus) = gpNewNodalArray(1+(icomp-1)*ndime:icomp*ndime)
         enddo
         a%tesgs(newelem)%a(:,a%e%igaus) = gpNewNodalArray(a%ncsgs*ndime+1:a%ncsgs*(ndime+1))
         a%prsgs(newelem)%a(a%e%igaus) = gpNewNodalArray(a%ncsgs*(ndime+1)+1)
      enddo

   end subroutine

   subroutine ELT_NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
      implicit none
      class(ElementLowMachTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) ::  newelem
      real(rp), target :: OldBuffer(*)
      
      real(rp), pointer :: auxOldNodalArray(:,:) => NULL()

      integer(ip) :: ndofn,ndime,pnode,pgaus
      real(rp), allocatable :: NewNodalArray(:,:)

      call a%LMP%Mesh%GetNdime(ndime)
      call a%LMP%Mesh%GetElemArraySize(newelem,pnode,pgaus)

      ndofn = (ndime+1)*a%ncsgs+1
      auxOldNodalArray(1:ndofn,1:pnode) => OldBuffer(1:ndofn*pnode)

      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 
 
      call Refiner%InterpolateElement(ndofn,auxOldNodalArray,NewNodalArray)

      call ELT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)

   end subroutine

   subroutine ELT_OldElementToBuffer(a,ielem,Buffer)
      class(ElementLowMachTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)
      

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

     
      icount = 0
      auxpointer(1:size(a%LMP%vesgs(ielem)%a)) => a%LMP%vesgs(ielem)%a
      icount = icount + size(a%LMP%vesgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%LMP%tesgs(ielem)%a)) => a%LMP%tesgs(ielem)%a
      icount = icount + size(a%LMP%tesgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%LMP%prsgs(ielem)%a)) => a%LMP%prsgs(ielem)%a
      icount = icount + size(a%LMP%prsgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer
     
   end subroutine

   subroutine ELT_NewElementToBuffer(a,ielem,Buffer)
      class(ElementLowMachTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%vesgs(ielem)%a)) => a%vesgs(ielem)%a
      icount = icount + size(a%vesgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%tesgs(ielem)%a)) => a%tesgs(ielem)%a
      icount = icount + size(a%tesgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%prsgs(ielem)%a)) => a%prsgs(ielem)%a
      icount = icount + size(a%prsgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer
     
   end subroutine

   subroutine ELT_BufferToNewElement(a,ielem,Buffer)
      class(ElementLowMachTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%vesgs(ielem)%a)) => a%vesgs(ielem)%a
      icount = icount + size(a%vesgs(ielem)%a)
      auxpointer = Buffer(1:icount)

      oldicount = icount

      auxpointer(1:size(a%tesgs(ielem)%a)) => a%tesgs(ielem)%a
      icount = icount + size(a%tesgs(ielem)%a)
      auxpointer = Buffer(oldicount+1:icount)

      oldicount = icount

      auxpointer(1:size(a%prsgs(ielem)%a)) => a%prsgs(ielem)%a
      icount = icount + size(a%prsgs(ielem)%a)
      auxpointer =  Buffer(oldicount+1:icount)
     
   end subroutine

   subroutine ELT_GetOldElementBufferSize(a,ielem,BufferSize)
      class(ElementLowMachTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount

      icount = 0
      icount = icount + size(a%LMP%vesgs(ielem)%a)
      icount = icount + size(a%LMP%tesgs(ielem)%a)
      icount = icount + size(a%LMP%prsgs(ielem)%a)
     
      BufferSize = icount
      
   end subroutine

   subroutine ELT_GetNewElementBufferSize(a,ielem,BufferSize)
      class(ElementLowMachTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount
      icount = 0
      icount = icount + size(a%vesgs(ielem)%a)
      icount = icount + size(a%tesgs(ielem)%a)
      icount = icount + size(a%prsgs(ielem)%a)
     
      BufferSize = icount
   end subroutine

end module

subroutine lmn_Prerefine(a)
   use typre
   use Mod_LowMach
   use Mod_phpRefineArrays
   use Mod_lmn_Refine
   use Mod_Element
   implicit none
   class(LowMachProblem) :: a

   integer(ip) :: nelem,ielem,ndime,icomp,ncsgs
   class(FiniteElement), pointer :: e => NULL()
   real(rp), target :: rwa_GaussArray(500)
   real(rp), pointer :: GaussArray(:,:) => NULL()

   if(a%kfl_trasg/=0) then

      ncsgs=a%ncomp-1

      !Things to be done prior to refinement
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_Refine')
   
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      !Allocate data array
      call a%Memor%alloc(nelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','lmn_refine')
   
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         allocate(a%NodalArrayDataForRefinement(ielem)%a(ncsgs*(ndime+1)+1,e%pnode))
         GaussArray(1:ncsgs*(ndime+1)+1,1:e%pgaus) => rwa_GaussArray(1:(ncsgs*(ndime+1)+1)*e%pgaus)
         do icomp = 1,ncsgs
             GaussArray(1+(icomp-1)*ndime:icomp*ndime,1:e%pgaus) = a%vesgs(ielem)%a(:,icomp,1:e%pgaus)
         enddo
         GaussArray(ndime*ncsgs+1:ncsgs*(ndime+1),:) = a%tesgs(ielem)%a
         GaussArray(ncsgs*(ndime+1)+1,:) = a%prsgs(ielem)%a
         call e%GaussToNodes(GaussArray,a%NodalArrayDataForRefinement(ielem)%a)
      enddo
   
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lmn_Refine')

    end if

end subroutine

subroutine lmn_Refine(a,itask)
   use typre
   use Mod_LowMach
   use Mod_phpRefineArrays
   use Mod_lmn_Refine

   implicit none
   class(LowMachProblem) :: a
   type(ElementLowMachTransferer) ::  ELT
   character(6) :: itask
   
   integer(ip)              :: oldnpoin,newnpoin,oldnelem,newnelem,icomp,ndime
   real(rp), allocatable    :: auxPointwiseSource(:),auxinitemp(:),auxdensf(:)
   real(rp), allocatable    :: auxrepro(:,:),auxgrprj(:,:,:)
   real(rp), allocatable    :: auxveloc(:,:,:),auxpress(:,:),auxtempe(:,:)
   integer(ip)              :: ielem,pgaus,vesgs_coun,tesgs_coun,prsgs_coun,pnode,auxsize
   
   !Old Dimensions
   oldnpoin = size(a%veloc,2)
   
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(newnelem)

   !Velocity
   call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxveloc,'veloc','lmn_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
      endif   
   enddo
   call move_alloc(auxveloc,a%veloc)
   call a%Memor%deallocObj(0,'veloc','lmn_Refine',rp*ndime*oldnpoin*a%ncomp)
   
   !Pressure
   call a%Memor%alloc(newnpoin,a%ncomp,auxpress,'press','lmn_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      endif   
   enddo
   call move_alloc(auxpress,a%press)
   call a%Memor%deallocObj(0,'press','lmn_Refine',rp*oldnpoin*a%ncomp)
   
   !Temperature
   call a%Memor%alloc(newnpoin,a%ncomp,auxtempe,'tempe','lmn_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
      endif   
   enddo
   call move_alloc(auxtempe,a%tempe)
   call a%Memor%deallocObj(0,'tempe','lmn_Refine',rp*oldnpoin*a%ncomp)
   
   !Initial Temperature
   call a%Memor%alloc(newnpoin,auxinitemp,'itemp','lmn_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%itemp,auxinitemp)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%itemp,auxinitemp)
      endif   
   call move_alloc(auxinitemp,a%itemp)
   call a%Memor%deallocObj(0,'itemp','lmn_Refine',rp*oldnpoin)
  
   !Initial density field 
   if(a%npp_stepi(4)>0) then     
   call a%Memor%alloc(newnpoin,auxdensf,'densf','lmn_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%densf,auxdensf)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%densf,auxdensf)
      endif   
   call move_alloc(auxdensf,a%densf)
   call a%Memor%deallocObj(0,'densf','lmn_Refine',rp*oldnpoin)
   end if
   
   !PointwiseSource
   if (a%kfl_sourc == 2) then
      call a%Memor%alloc(newnpoin,auxPointwiseSource,'auxPointwiseSource','lmn_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%PointwiseSource,auxPointwiseSource)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%PointwiseSource,auxPointwiseSource)
      endif
      call move_alloc(auxPointwiseSource,a%PointwiseSource)
      call a%Memor%deallocObj(0,'PointwiseSource','lmn_Refine',rp*a%ndofn*oldnpoin)
   endif

   !Residual Projection
   if (a%kfl_repro >= 1) then
      auxsize = size(a%repro,1)
      call a%Memor%alloc(auxsize,newnpoin,auxrepro,'repro','lmn_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      endif   
      call move_alloc(auxrepro,a%repro)
      call a%Memor%deallocObj(0,'repro','lmn_Refine',rp*oldnpoin*auxsize)
   endif
   
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      auxsize = size(a%grprj,1)
      call a%Memor%alloc(auxsize,ndime,newnpoin,auxgrprj,'grprj','lmn_Refine')
      do icomp = 1,auxsize
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         endif   
      enddo
      call move_alloc(auxgrprj,a%grprj)
      call a%Memor%deallocObj(0,'grprj','lmn_Refine',rp*auxsize*ndime*oldnpoin)
   end if

   !Subgrid scales
   if(a%kfl_trasg/=0) then
      call ELT%SetProblem(a) 
      call ELT%AllocateArrays 

      if (itask == 'Refine') then
         !Don't call us, we'll call YOU!
         call a%Refiner%ElementalRefine(ELT)
      elseif (itask == 'Rebala') then
         call a%Refiner%ElementalRebalance(ELT)
      endif

      oldnelem = size(a%prsgs,1)
      vesgs_coun = 0
      tesgs_coun = 0
      prsgs_coun = 0
      do ielem=1,oldnelem
         pgaus = size(a%tesgs(ielem)%a,2)
         deallocate(a%vesgs(ielem)%a)
         deallocate(a%tesgs(ielem)%a)
         deallocate(a%prsgs(ielem)%a)
         vesgs_coun = vesgs_coun + ndime*ELT%ncsgs*pgaus
         tesgs_coun = tesgs_coun + ELT%ncsgs*pgaus
         prsgs_coun = prsgs_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'vesgs%a','lmn_refine',vesgs_coun*rp)
      call a%Memor%deallocObj(0,'tesgs%a','lmn_refine',tesgs_coun*rp)
      call a%Memor%deallocObj(0,'prsgs%a','lmn_refine',prsgs_coun*rp)
      call a%Memor%dealloc(oldnelem,a%vesgs,'vesgs','lmn_refine')
      call a%Memor%dealloc(oldnelem,a%tesgs,'tesgs','lmn_refine')
      call a%Memor%dealloc(oldnelem,a%prsgs,'prsgs','lmn_refine')

      call move_alloc(ELT%vesgs,a%vesgs)
      call move_alloc(ELT%tesgs,a%tesgs)
      call move_alloc(ELT%prsgs,a%prsgs)

      call a%Mesh%ElementDeAlloc(ELT%e,a%Memor,'DefaultRule','lmn_Refine')
      
      if (itask == 'Refine') then
         !DeAllocate data array
         call a%Memor%dealloc(oldnelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','lmn_refine')
      end if 

   endif
   
end subroutine   

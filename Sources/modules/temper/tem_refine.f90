module Mod_tem_Refine
   use typre
   use Mod_Temperature
   use Mod_AdaptiveInterface
   use Mod_Element
   implicit none
   private
   public ElementTemperatureTransferer

   type, extends(AdaptiveRefinerElementTransferer) :: ElementTemperatureTransferer
      type(r2p), allocatable :: &
            tesgs(:)                            !Temperature subgrid scales
      class(TemperatureProblem), pointer :: TEM => NULL()
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ncsgs

contains
      procedure :: CopyElement => ETT_CopyElement
      procedure :: NodallyInterpolateElement => ETT_NodallyInterpolateElement
      procedure :: NodallyInterpolateElementFromOldBuffer => ETT_NodallyInterpolateElementFromOldBuffer
      

      procedure :: OldElementToBuffer => ETT_OldElementToBuffer
      procedure :: GetOldElementBufferSize => ETT_GetOldElementBufferSize

      procedure :: GetNewElementBufferSize => ETT_GetNewElementBufferSize
      procedure :: BufferToNewElement => ETT_BufferToNewElement
      procedure :: NewElementToBuffer => ETT_NewElementToBuffer

      procedure :: SetProblem
      procedure :: AllocateArrays

   end type


contains

   subroutine SetProblem(a,TEM)
      implicit none
      class(ElementTemperatureTransferer) :: a
      class(TemperatureProblem), target :: TEM

      a%TEM => TEM
      
      a%ncsgs=1                  !Number of components in time for accuracy
      if(TEM%kfl_tacsg>0) a%ncsgs=2

   end subroutine

   subroutine AllocateArrays(a)
      implicit none
      class(ElementTemperatureTransferer) :: a

      integer(ip) :: ielem,nelem,ndime,pgaus,pnode
      integer(ip) :: tesgs_coun
      
      call a%TEM%Mesh%GetNelem(nelem)
      call a%TEM%Mesh%GetNdime(ndime)
    
      tesgs_coun = 0

      !For explicit time integration scheme:
      !position 1: stage increment
      !position 2: stage contribution
      !position 3: previous time step 
      
      call a%TEM%Memor%alloc(nelem,a%tesgs,'tesgs','tem_refine')
      do ielem=1,nelem
         call a%TEM%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%tesgs(ielem)%a(a%ncsgs,pgaus))
         a%tesgs(ielem)%a = 0.0_rp
         tesgs_coun = tesgs_coun + a%ncsgs*pgaus
      end do
      call a%TEM%Memor%allocObj(0,'tesgs%a','tem_refine',tesgs_coun*rp)

      !Allocate element
      call a%TEM%Mesh%ElementAlloc(a%e,a%TEM%Memor,'DefaultRule','tem_pr_elmope')
    
   end subroutine

   subroutine ETT_CopyElement(a,oldielem,newielem)
      implicit none
      class(ElementTemperatureTransferer), target :: a
      integer(ip) :: oldielem,newielem

      a%tesgs(newielem)%a = a%TEM%tesgs(oldielem)%a 

   end subroutine

   subroutine ETT_NodallyInterpolateElement(a,Refiner,oldelem,newelem)
      implicit none
      class(ElementTemperatureTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) :: oldelem, newelem

      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)
      real(rp), pointer :: OldNodalArray(:,:) => NULL()
      integer(ip) :: ndime,pgaus,pnode

      call a%TEM%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      call a%TEM%Mesh%GetNdime(ndime)
      ndofn = a%ncsgs
      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 

      call Refiner%InterpolateElement(ndofn,a%TEM%NodalArrayDataForRefinement(oldelem)%a,NewNodalArray)

      call ETT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)
   end subroutine

   subroutine ETT_NodalArrayToSubscales(a,newelem,NewNodalArray)
      implicit none
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: newelem
      real(rp)  :: NewNodalArray(:,:)
      integer(ip) :: pgaus,igaus,ndofn,ndime,icomp

      real(rp)  :: gpNewNodalArray(200)

      call a%TEM%Mesh%GetNdime(ndime)
      call a%TEM%Mesh%ElementLoad(newelem,a%e)
   
      ndofn  = size(NewNodalArray,1) 
      do igaus = 1,a%e%pgaus
         a%e%igaus = igaus
         call a%e%interpg(ndofn,NewNodalArray,gpNewNodalArray(1:ndofn))
         a%tesgs(newelem)%a(:,a%e%igaus) = gpNewNodalArray(1:a%ncsgs)
      enddo

   end subroutine

   subroutine ETT_NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
      implicit none
      class(ElementTemperatureTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) ::  newelem
      real(rp), target :: OldBuffer(*)
      
      real(rp), pointer :: auxOldNodalArray(:,:) => NULL()

      integer(ip) :: ndofn,ndime,pnode,pgaus
      real(rp), allocatable :: NewNodalArray(:,:)

      call a%TEM%Mesh%GetNdime(ndime)
      call a%TEM%Mesh%GetElemArraySize(newelem,pnode,pgaus)

      ndofn = a%ncsgs
      auxOldNodalArray(1:ndofn,1:pnode) => OldBuffer(1:ndofn*pnode)

      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 
 
      call Refiner%InterpolateElement(ndofn,auxOldNodalArray,NewNodalArray)

      call ETT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)

   end subroutine

   subroutine ETT_OldElementToBuffer(a,ielem,Buffer)
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

     
      icount = 0

      auxpointer(1:size(a%TEM%tesgs(ielem)%a)) => a%TEM%tesgs(ielem)%a
      icount = icount + size(a%TEM%tesgs(ielem)%a)
      Buffer(1:icount) = auxpointer
     
   end subroutine

   subroutine ETT_NewElementToBuffer(a,ielem,Buffer)
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0

      auxpointer(1:size(a%tesgs(ielem)%a)) => a%tesgs(ielem)%a
      icount = icount + size(a%tesgs(ielem)%a)
      Buffer(1:icount) = auxpointer
     
   end subroutine

   subroutine ETT_BufferToNewElement(a,ielem,Buffer)
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0

      auxpointer(1:size(a%tesgs(ielem)%a)) => a%tesgs(ielem)%a
      icount = icount + size(a%tesgs(ielem)%a)
      auxpointer = Buffer(1:icount)

     
   end subroutine

   subroutine ETT_GetOldElementBufferSize(a,ielem,BufferSize)
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount

      icount = 0
      icount = icount + size(a%TEM%tesgs(ielem)%a)
     
      BufferSize = icount
      
   end subroutine

   subroutine ETT_GetNewElementBufferSize(a,ielem,BufferSize)
      class(ElementTemperatureTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount
      icount = 0
      icount = icount + size(a%tesgs(ielem)%a)
     
      BufferSize = icount
   end subroutine

end module

subroutine tem_Prerefine(a)
   use typre
   use Mod_Temperature
   use Mod_phpRefineArrays
   use Mod_tem_Refine
   use Mod_Element
   implicit none
   class(TemperatureProblem) :: a

   integer(ip) :: nelem,ielem,ndime,icomp,ncsgs
   class(FiniteElement), pointer :: e => NULL()
   real(rp), target :: rwa_GaussArray(500)
   real(rp), pointer :: GaussArray(:,:) => NULL()

   if(a%kfl_trasg/=0) then

      ncsgs=2

      !Things to be done prior to refinement
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','tem_Refine')
   
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      !Allocate data array
      call a%Memor%alloc(nelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','tem_refine')
   
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         allocate(a%NodalArrayDataForRefinement(ielem)%a(ncsgs,e%pnode))
         GaussArray(1:ncsgs,1:e%pgaus) => rwa_GaussArray(1:ncsgs*e%pgaus)
         GaussArray(1:ncsgs,:) = a%tesgs(ielem)%a
         call e%GaussToNodes(GaussArray,a%NodalArrayDataForRefinement(ielem)%a)
      enddo
   
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','tem_Refine')

    end if

end subroutine

subroutine tem_Refine(a,itask)
   use typre
   use Mod_Temperature
   use Mod_phpRefineArrays
   use Mod_tem_Refine
   implicit none
   class(TemperatureProblem) :: a
   type(ElementTemperatureTransferer) ::  ETT
   character(6) :: itask
   
   integer(ip) :: oldnpoin,newnpoin,newnboun,icomp,oldnelem,newnelem,ielem
   integer(ip) :: tesgs_coun,pgaus,auxsize
   real(rp), allocatable :: auxtempe(:,:),auxPointwiseSource(:),auxdissipation(:)
   real(rp), allocatable :: auxrepro(:),auxgrprj(:,:)
   
   !This subroutine transforms the tempe arrays from one mesh to the other 
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%tempe,1)

   
   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNboun(newnboun)
   
   !Tempe
   call a%Memor%alloc(newnpoin,a%ncomp,auxtempe,'tempe','tem_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
      endif   
   enddo
   call move_alloc(auxtempe,a%tempe)
   call a%Memor%deallocObj(0,'tempe','tem_Refine',rp*a%ndofn*a%ncomp*oldnpoin)
   
   !PointwiseSource
   if (a%kfl_sourc == 2) then
      call a%Memor%alloc(newnpoin,auxPointwiseSource,'auxPointwiseSource','tem_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%PointwiseSource,auxPointwiseSource)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%PointwiseSource,auxPointwiseSource)
      endif
      call move_alloc(auxPointwiseSource,a%PointwiseSource)
      call a%Memor%deallocObj(0,'PointwiseSource','tem_Refine',rp*a%ndofn*oldnpoin)
   endif
   
   !repro
   if (a%kfl_repro == 1) then
      call a%Memor%alloc(newnpoin,auxrepro,'auxrepro','tem_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(a%ndofn,a%repro,auxrepro)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(a%ndofn,a%repro,auxrepro)
      endif
      call move_alloc(auxrepro,a%repro)
      call a%Memor%deallocObj(0,'repro','tem_Refine',rp*a%ndofn*oldnpoin)
   endif

   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      auxsize = size(a%grprj,1)
      call a%Memor%alloc(auxsize,newnpoin,auxgrprj,'grprj','tem_Refine')
      do icomp = 1,auxsize
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(1_ip,a%grprj(icomp,:),auxgrprj(icomp,:))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(1_ip,a%grprj(icomp,:),auxgrprj(icomp,:))
         endif   
      enddo
      call move_alloc(auxgrprj,a%grprj)
      call a%Memor%deallocObj(0,'grprj','tem_Refine',rp*auxsize*oldnpoin)
   end if
   
   !dissipation
   if (a%kfl_dispa == 1) then
      call a%Memor%alloc(newnpoin,auxdissipation,'auxdissipation','tem_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(a%ndofn,a%dissipation,auxdissipation)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(a%ndofn,a%dissipation,auxdissipation)
      endif
      call move_alloc(auxdissipation,a%dissipation)
      call a%Memor%deallocObj(0,'dissipation','tem_Refine',rp*a%ndofn*oldnpoin)
   endif
   
   !Subgrid scales
   if(a%kfl_trasg/=0) then
      call ETT%SetProblem(a) 
      call ETT%AllocateArrays 

      if (itask == 'Refine') then
         !Don't call us, we'll call YOU!
         call a%Refiner%ElementalRefine(ETT)
      elseif (itask == 'Rebala') then
         call a%Refiner%ElementalRebalance(ETT)
      endif

      oldnelem = size(a%tesgs,1)
      tesgs_coun = 0
      do ielem=1,oldnelem
         pgaus = size(a%tesgs(ielem)%a,2)
         deallocate(a%tesgs(ielem)%a)
         tesgs_coun = tesgs_coun + ETT%ncsgs*pgaus
      end do
      call a%Memor%deallocObj(0,'tesgs%a','lmn_refine',tesgs_coun*rp)
      call a%Memor%dealloc(oldnelem,a%tesgs,'tesgs','lmn_refine')

      call move_alloc(ETT%tesgs,a%tesgs)

      call a%Mesh%ElementDeAlloc(ETT%e,a%Memor,'DefaultRule','lmn_Refine')
      
      if (itask == 'Refine') then
         !DeAllocate data array
         call a%Memor%dealloc(oldnelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','lmn_refine')
      end if 

   endif
   
end subroutine   

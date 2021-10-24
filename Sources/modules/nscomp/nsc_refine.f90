module Mod_nsc_Refine
   use typre
   use Mod_NSCompressible
   use Mod_AdaptiveInterface
   use Mod_Element
   implicit none
   private
   public ElementCompressibleTransferer

   type, extends(AdaptiveRefinerElementTransferer) :: ElementCompressibleTransferer
      type(r2p), allocatable :: &
            cosgs(:),&                          !Continuity subgrid scales
            ensgs(:)                            !Energy subgrid scales

      type(r3p), allocatable :: &
            mosgs(:)                            !Momentum subgrid scales
      class(NSCompressibleProblem), pointer :: NSC => NULL()
      class(FiniteElement), pointer :: e => NULL()
contains
      procedure :: CopyElement => ECT_CopyElement
      procedure :: NodallyInterpolateElement => ECT_NodallyInterpolateElement
      procedure :: NodallyInterpolateElementFromOldBuffer => ECT_NodallyInterpolateElementFromOldBuffer
      

      procedure :: OldElementToBuffer => ECT_OldElementToBuffer
      procedure :: GetOldElementBufferSize => ECT_GetOldElementBufferSize

      procedure :: GetNewElementBufferSize => ECT_GetNewElementBufferSize
      procedure :: BufferToNewElement => ECT_BufferToNewElement
      procedure :: NewElementToBuffer => ECT_NewElementToBuffer

      procedure :: SetNSC
      procedure :: AllocateArrays

   end type


contains

   subroutine SetNSC(a,NSC)
      implicit none
      class(ElementCompressibleTransferer) :: a
      class(NSCompressibleProblem), target :: NSC

      a%NSC => NSC
   end subroutine

   subroutine AllocateArrays(a)
      implicit none
      class(ElementCompressibleTransferer) :: a

      integer(ip) :: ielem,nelem,ndime,pgaus,pnode
      integer(ip) ::  cosgs_coun,mosgs_coun,ensgs_coun
      
      call a%NSC%Mesh%GetNelem(nelem)
      call a%NSC%Mesh%GetNdime(ndime)
    
      cosgs_coun = 0
      mosgs_coun = 0
      ensgs_coun = 0

      !For explicit time integration scheme:
      !position 1: stage increment
      !position 2: stage contribution
      !position 3: previous time step 
      
      call a%NSC%Memor%alloc(nelem,a%cosgs,'cosgs','nsc_refine')
      call a%NSC%Memor%alloc(nelem,a%mosgs,'mosgs','nsc_refine')
      call a%NSC%Memor%alloc(nelem,a%ensgs,'ensgs','nsc_refine')
      do ielem=1,nelem
         call a%NSC%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%cosgs(ielem)%a(ncsgs,pgaus))
         allocate(a%mosgs(ielem)%a(ndime,ncsgs,pgaus))
         allocate(a%ensgs(ielem)%a(ncsgs,pgaus))
         a%cosgs(ielem)%a = 0.0_rp
         a%mosgs(ielem)%a = 0.0_rp
         a%ensgs(ielem)%a = 0.0_rp
         cosgs_coun = cosgs_coun + ncsgs*pgaus
         mosgs_coun = mosgs_coun + ndime*ncsgs*pgaus
         ensgs_coun = ensgs_coun + ncsgs*pgaus
      end do
      call a%NSC%Memor%allocObj(0,'cosgs%a','nsc_refine',cosgs_coun*rp)
      call a%NSC%Memor%allocObj(0,'mosgs%a','nsc_refine',mosgs_coun*rp)
      call a%NSC%Memor%allocObj(0,'ensgs%a','nsc_refine',ensgs_coun*rp)

      !Allocate element
      call a%NSC%Mesh%ElementAlloc(a%e,a%NSC%Memor,'DefaultRule','nsc_pr_elmope')
    
   end subroutine

   subroutine ECT_CopyElement(a,oldielem,newielem)
      implicit none
      class(ElementCompressibleTransferer), target :: a
      integer(ip) :: oldielem,newielem

      a%cosgs(newielem)%a = a%NSC%cosgs(oldielem)%a 
      a%mosgs(newielem)%a = a%NSC%mosgs(oldielem)%a 
      a%ensgs(newielem)%a = a%NSC%ensgs(oldielem)%a 

   end subroutine

   subroutine ECT_NodallyInterpolateElement(a,Refiner,oldelem,newelem)
      implicit none
      class(ElementCompressibleTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) :: oldelem, newelem

      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)
      real(rp), pointer :: OldNodalArray(:,:) => NULL()
      integer(ip) :: ndime,pgaus,pnode

      call a%NSC%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      call a%NSC%Mesh%GetNdime(ndime)
      ndofn = (ndime + 2)*ncsgs
      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 

      call Refiner%InterpolateElement(ndofn,a%NSC%NodalArrayDataForRefinement(oldelem)%a,NewNodalArray)

      call ECT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)
   end subroutine

   subroutine ECT_NodalArrayToSubscales(a,newelem,NewNodalArray)
      implicit none
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: newelem
      real(rp)  :: NewNodalArray(:,:)
      integer(ip) :: pgaus,igaus,ndofn,ndime,icomp

      real(rp)  :: gpNewNodalArray(200)

      call a%NSC%Mesh%GetNdime(ndime)
      call a%NSC%Mesh%ElementLoad(newelem,a%e)
   
      ndofn  = size(NewNodalArray,1) 
      do igaus = 1,a%e%pgaus
         a%e%igaus = igaus
         call a%e%interpg(ndofn,NewNodalArray,gpNewNodalArray(1))
         a%cosgs(newelem)%a(:,a%e%igaus) = gpNewNodalArray(1:ncsgs)
         do icomp = 1, ncsgs
             a%mosgs(newelem)%a(:,icomp,a%e%igaus) = gpNewNodalArray(ncsgs+1+(icomp-1)*ndime:ncsgs+icomp*ndime)
         enddo
         a%ensgs(newelem)%a(:,a%e%igaus) = gpNewNodalArray(ncsgs*(ndime+1)+1:ncsgs*(ndime+2))
      enddo

   end subroutine

   subroutine ECT_NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
      implicit none
      class(ElementCompressibleTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) ::  newelem
      real(rp), target :: OldBuffer(*)
      
      real(rp), pointer :: auxOldNodalArray(:,:) => NULL()

      integer(ip) :: ndofn,ndime,pnode,pgaus
      real(rp), allocatable :: NewNodalArray(:,:)

      call a%NSC%Mesh%GetNdime(ndime)
      call a%NSC%Mesh%GetElemArraySize(newelem,pnode,pgaus)

      ndofn = (ndime+2)*ncsgs
      auxOldNodalArray(1:ndofn,1:pnode) => OldBuffer(1:ndofn*pnode)

      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp 
 
      call Refiner%InterpolateElement(ndofn,auxOldNodalArray,NewNodalArray)

      call ECT_NodalArrayToSubscales(a,newelem,NewNodalArray)

      deallocate(NewNodalArray)

   end subroutine

   subroutine ECT_OldElementToBuffer(a,ielem,Buffer)
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)
      

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

     
      icount = 0
      auxpointer(1:size(a%NSC%cosgs(ielem)%a)) => a%NSC%cosgs(ielem)%a
      icount = icount + size(a%NSC%cosgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount
     
      auxpointer(1:size(a%NSC%mosgs(ielem)%a)) => a%NSC%mosgs(ielem)%a
      icount = icount + size(a%NSC%mosgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%NSC%ensgs(ielem)%a)) => a%NSC%ensgs(ielem)%a
      icount = icount + size(a%NSC%ensgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer


   end subroutine

   subroutine ECT_NewElementToBuffer(a,ielem,Buffer)
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%cosgs(ielem)%a)) => a%cosgs(ielem)%a
      icount = icount + size(a%cosgs(ielem)%a)
      Buffer(1:icount) = auxpointer

      oldicount = icount
     
      auxpointer(1:size(a%mosgs(ielem)%a)) => a%mosgs(ielem)%a
      icount = icount + size(a%mosgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer

      oldicount = icount

      auxpointer(1:size(a%ensgs(ielem)%a)) => a%ensgs(ielem)%a
      icount = icount + size(a%ensgs(ielem)%a)
      Buffer(oldicount+1:icount) = auxpointer

   end subroutine

   subroutine ECT_BufferToNewElement(a,ielem,Buffer)
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: icount, oldicount
      real(rp), pointer :: auxpointer(:) => NULL()

      icount = 0
      auxpointer(1:size(a%cosgs(ielem)%a)) => a%cosgs(ielem)%a
      icount = icount + size(a%cosgs(ielem)%a)
      auxpointer =  Buffer(1:icount)

      oldicount = icount
     
      auxpointer(1:size(a%mosgs(ielem)%a)) => a%mosgs(ielem)%a
      icount = icount + size(a%mosgs(ielem)%a)
      auxpointer = Buffer(oldicount+1:icount)

      oldicount = icount

      auxpointer(1:size(a%ensgs(ielem)%a)) => a%ensgs(ielem)%a
      icount = icount + size(a%ensgs(ielem)%a)
      auxpointer = Buffer(oldicount+1:icount)

   end subroutine

   subroutine ECT_GetOldElementBufferSize(a,ielem,BufferSize)
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount

      icount = 0
      icount = icount + size(a%NSC%cosgs(ielem)%a)
      icount = icount + size(a%NSC%mosgs(ielem)%a)
      icount = icount + size(a%NSC%ensgs(ielem)%a)
     
      BufferSize = icount
      
   end subroutine

   subroutine ECT_GetNewElementBufferSize(a,ielem,BufferSize)
      class(ElementCompressibleTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      integer(ip) :: icount
      icount = 0
      icount = icount + size(a%cosgs(ielem)%a)
      icount = icount + size(a%mosgs(ielem)%a)
      icount = icount + size(a%ensgs(ielem)%a)
     
      BufferSize = icount
   end subroutine

end module

subroutine nsc_Prerefine(a)
   use typre
   use Mod_NSCompressible
   use Mod_phpRefineArrays
   use Mod_nsc_Refine
   use Mod_Element
   implicit none
   class(NSCompressibleProblem) :: a

   integer(ip) :: nelem,ielem,ndime,icomp
   class(FiniteElement), pointer :: e => NULL()
   real(rp), target :: rwa_GaussArray(500)
   real(rp), pointer :: GaussArray(:,:) => NULL()

   if(a%kfl_trasg/=0) then

      !Things to be done prior to refinement
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_Refine')
   
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      !Allocate data array
      call a%Memor%alloc(nelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','nsc_refine')
   
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         allocate(a%NodalArrayDataForRefinement(ielem)%a(ncsgs*(ndime+2),e%pnode))
         GaussArray(1:ncsgs*(ndime+2),1:e%pgaus) => rwa_GaussArray(1:ncsgs*(ndime+2)*e%pgaus)
         GaussArray(1:ncsgs,:) = a%cosgs(ielem)%a
         do icomp = 1,ncsgs
             GaussArray(ncsgs+1+(icomp-1)*ndime:ncsgs+(icomp)*ndime,1:e%pgaus) = a%mosgs(ielem)%a(:,icomp,1:e%pgaus)
         enddo
         GaussArray((ncsgs+(ndime*ncsgs))+1:size(GaussArray,1),:) = a%ensgs(ielem)%a
         call e%GaussToNodes(GaussArray,a%NodalArrayDataForRefinement(ielem)%a)
      enddo
   
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsc_Refine')

    end if

end subroutine

subroutine nsc_Refine(a,itask)
   use typre
   use Mod_NSCompressible
   use Mod_phpRefineArrays
   use Mod_nsc_Refine

   implicit none
   class(NSCompressibleProblem) :: a
   character(6) :: itask
   
   integer(ip) :: oldnpoin,oldnelem,newnpoin,icomp,ndime,ielem
   real(rp), allocatable :: auxrepro(:,:),auxgrprj(:,:,:)
   real(rp), allocatable :: auxtimad(:,:),auxtimav(:,:,:)
   real(rp), allocatable :: auxtimap(:,:),auxtimat(:,:)
   real(rp), allocatable :: auxrmsd(:),auxrmsv(:,:)
   real(rp), allocatable :: auxrmsp(:),auxrmst(:)
   real(rp), allocatable :: auxturbi(:)
   real(rp), allocatable :: auxforfl(:,:),auxhfxfl(:)
   real(rp), allocatable :: auxmomfl(:,:)
   real(rp), allocatable :: auxbtraction(:,:)
   integer(ip) :: auxsize,pgaus
   integer(ip) ::  cosgs_coun,mosgs_coun,ensgs_coun

   type(ElementCompressibleTransferer) ::  ECT
   
   !Old Dimensions
   oldnpoin = size(a%densf,1)
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNdime(ndime)

   !Residual Projection
   if (a%kfl_repro >= 1) then
      auxsize = size(a%repro,1)
      call a%Memor%alloc(auxsize,newnpoin,auxrepro,'repro','nsc_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(auxsize,a%repro(:,:),auxrepro(:,:))
      endif   
      call move_alloc(auxrepro,a%repro)
      call a%Memor%deallocObj(0,'repro','nsc_Refine',rp*oldnpoin*auxsize)
   endif
   
   !Gradient Orthogonal projection
   if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
      auxsize = size(a%grprj,1)
      call a%Memor%alloc(auxsize,ndime,newnpoin,auxgrprj,'grprj','nsc_Refine')
      do icomp = 1,auxsize
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(ndime,a%grprj(icomp,:,:),auxgrprj(icomp,:,:))
         endif   
      enddo
      call move_alloc(auxgrprj,a%grprj)
      call a%Memor%deallocObj(0,'grprj','nsc_Refine',rp*auxsize*ndime*oldnpoin)
   end if

   !Body forces and moments 
   if(a%kfl_outfm==1)then   
      if(a%npp_stepi(11) /= 0)then   
        call a%Memor%alloc(ndime,newnpoin,auxforfl,'forfl','nsc_Refine')
        call a%Memor%alloc(newnpoin,auxhfxfl,'hfxfl','nsc_Refine')
        auxsize = size(a%momfl,1)
        call a%Memor%alloc(auxsize,newnpoin,auxmomfl,'momfl','nsc_Refine')
        if (itask == 'Refine') then
           call a%Refiner%UpdateVariable(ndime,a%forfl,auxforfl)
           call a%Refiner%UpdateVariable(1_ip,a%hfxfl,auxhfxfl)
           call a%Refiner%UpdateVariable(auxsize,a%momfl,auxmomfl)
        elseif (itask == 'Rebala') then
           call a%Refiner%RebalanceVariable(ndime,a%forfl,auxforfl)
           call a%Refiner%RebalanceVariable(1_ip,a%hfxfl,auxhfxfl)
           call a%Refiner%RebalanceVariable(auxsize,a%momfl,auxmomfl)
        endif   
        call move_alloc(auxforfl,a%forfl)
        call move_alloc(auxhfxfl,a%hfxfl)
        call move_alloc(auxmomfl,a%momfl)
        call a%Memor%deallocObj(0,'forfl','nsc_Refine',rp*oldnpoin*ndime)
        call a%Memor%deallocObj(0,'hfxfl','nsc_Refine',rp*oldnpoin)
        call a%Memor%deallocObj(0,'momfl','nsc_Refine',rp*oldnpoin*auxsize)
   
      end if
   end if

   !Postprocess the added shock capturing diffusiviy at the gauss points
   if (a%npp_stepi(9) /= 0) then
      call runend('Gauss point dissipation: Not ready for adaptivity')
   endif

   !Postprocess the added subgrid diffusiviy at the gauss points
   if (a%npp_stepi(10) /= 0) then
      call runend('Gauss point dissipation: Not ready for adaptivity')
   endif

   !Statistics of fields
   if (a%npp_stepi(15) /= 0) then
      call a%Memor%alloc(newnpoin,    2,auxtimad,'timad','nsc_Refine')
      call a%Memor%alloc(ndime,newnpoin,2,auxtimav,'timav','nsc_Refine')
      call a%Memor%alloc(newnpoin,     2,auxtimap,'timap','nsc_Refine')
      call a%Memor%alloc(newnpoin,     2,auxtimat,'timat','nsc_Refine')
      call a%Memor%alloc(newnpoin,auxrmsd,'rmsd','nsc_Refine')
      call a%Memor%alloc(ndime,newnpoin,auxrmsv,'rmsv','nsc_Refine')
      call a%Memor%alloc(newnpoin,auxrmsp,'rmsp','nsc_Refine')
      call a%Memor%alloc(newnpoin,auxrmst,'rmst','nsc_Refine')
      call a%Memor%alloc(newnpoin,auxturbi,'turbi','nsc_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%timad(:,1),auxtimad(:,1))
         call a%Refiner%UpdateVariable(1_ip,a%timad(:,2),auxtimad(:,2))
         call a%Refiner%UpdateVariable(ndime,a%timav(:,:,1),auxtimav(:,:,1))
         call a%Refiner%UpdateVariable(ndime,a%timav(:,:,2),auxtimav(:,:,2))
         call a%Refiner%UpdateVariable(1_ip,a%timap(:,1),auxtimap(:,1))
         call a%Refiner%UpdateVariable(1_ip,a%timap(:,2),auxtimap(:,2))
         call a%Refiner%UpdateVariable(1_ip,a%timat(:,1),auxtimat(:,1))
         call a%Refiner%UpdateVariable(1_ip,a%timat(:,2),auxtimat(:,2))
         call a%Refiner%UpdateVariable(1_ip,a%rmsd,auxrmsd)
         call a%Refiner%UpdateVariable(ndime,a%rmsv,auxrmsv)
         call a%Refiner%UpdateVariable(1_ip,a%rmsp,auxrmsp)
         call a%Refiner%UpdateVariable(1_ip,a%rmst,auxrmst)
         call a%Refiner%UpdateVariable(1_ip,a%turbi,auxturbi)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%timad(:,1),auxtimad(:,1))
         call a%Refiner%RebalanceVariable(1_ip,a%timad(:,2),auxtimad(:,2))
         call a%Refiner%RebalanceVariable(ndime,a%timav(:,:,1),auxtimav(:,:,1))
         call a%Refiner%RebalanceVariable(ndime,a%timav(:,:,2),auxtimav(:,:,2))
         call a%Refiner%RebalanceVariable(1_ip,a%timap(:,1),auxtimap(:,1))
         call a%Refiner%RebalanceVariable(1_ip,a%timap(:,2),auxtimap(:,2))
         call a%Refiner%RebalanceVariable(1_ip,a%timat(:,1),auxtimat(:,1))
         call a%Refiner%RebalanceVariable(1_ip,a%timat(:,2),auxtimat(:,2))
         call a%Refiner%RebalanceVariable(1_ip,a%rmsd,auxrmsd)
         call a%Refiner%RebalanceVariable(ndime,a%rmsv,auxrmsv)
         call a%Refiner%RebalanceVariable(1_ip,a%rmsp,auxrmsp)
         call a%Refiner%RebalanceVariable(1_ip,a%rmst,auxrmst)
         call a%Refiner%RebalanceVariable(1_ip,a%turbi,auxturbi)
      endif   
      call move_alloc(auxtimad,a%timad)
      call move_alloc(auxtimav,a%timav)
      call move_alloc(auxtimap,a%timap)
      call move_alloc(auxtimat,a%timat)
      call move_alloc(auxrmsd,a%rmsd)
      call move_alloc(auxrmsv,a%rmsv)
      call move_alloc(auxrmsp,a%rmsp)
      call move_alloc(auxrmst,a%rmst)
      call move_alloc(auxturbi,a%turbi)
      call a%Memor%deallocObj(0,'timad','nsc_Refine',rp*oldnpoin*2)
      call a%Memor%deallocObj(0,'timav','nsc_Refine',rp*ndime*oldnpoin*2)
      call a%Memor%deallocObj(0,'timap','nsc_Refine',rp*oldnpoin*2)
      call a%Memor%deallocObj(0,'timat','nsc_Refine',rp*oldnpoin*2)
      call a%Memor%deallocObj(0,'rmsd','nsc_Refine',rp*oldnpoin)
      call a%Memor%deallocObj(0,'rmsv','nsc_Refine',rp*ndime*oldnpoin)
      call a%Memor%deallocObj(0,'rmsp','nsc_Refine',rp*oldnpoin)
      call a%Memor%deallocObj(0,'rmst','nsc_Refine',rp*oldnpoin)
      call a%Memor%deallocObj(0,'turbi','nsc_Refine',rp*oldnpoin)
   endif   
   
   !Btraction
   call a%Memor%alloc(ndime,newnpoin,auxbtraction,'btraction','nsc_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%btraction,auxbtraction)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%btraction,auxbtraction)
   endif   
   call move_alloc(auxbtraction,a%btraction)
   call a%Memor%deallocObj(0,'btraction','nsc_Refine',rp*oldnpoin*ndime)

    
   !Subgrid scales
   if(a%kfl_trasg/=0) then
      call ECT%SetNSC(a) 
      call ECT%AllocateArrays 

      if (itask == 'Refine') then
         !Don't call us, we'll call YOU!
         call a%Refiner%ElementalRefine(ECT)
      elseif (itask == 'Rebala') then
         call a%Refiner%ElementalRebalance(ECT)
      endif

     
      oldnelem = size(a%cosgs,1)
      cosgs_coun = 0
      mosgs_coun = 0
      ensgs_coun = 0
      do ielem=1,oldnelem
         pgaus = size(a%cosgs(ielem)%a,2)
         deallocate(a%cosgs(ielem)%a)
         deallocate(a%mosgs(ielem)%a)
         deallocate(a%ensgs(ielem)%a)
         cosgs_coun = cosgs_coun + ncsgs*pgaus
         mosgs_coun = mosgs_coun + ndime*ncsgs*pgaus
         ensgs_coun = ensgs_coun + ncsgs*pgaus
      end do
      call a%Memor%deallocObj(0,'cosgs%a','nsc_refine',cosgs_coun*rp)
      call a%Memor%deallocObj(0,'mosgs%a','nsc_refine',mosgs_coun*rp)
      call a%Memor%deallocObj(0,'ensgs%a','nsc_refine',ensgs_coun*rp)
      call a%Memor%dealloc(oldnelem,a%cosgs,'cosgs','nsc_refine')
      call a%Memor%dealloc(oldnelem,a%mosgs,'mosgs','nsc_refine')
      call a%Memor%dealloc(oldnelem,a%ensgs,'ensgs','nsc_refine')

      call move_alloc(ECT%cosgs,a%cosgs)
      call move_alloc(ECT%mosgs,a%mosgs)
      call move_alloc(ECT%ensgs,a%ensgs)

      call a%Mesh%ElementDeAlloc(ECT%e,a%Memor,'DefaultRule','nsc_Refine')
      
      if (itask == 'Refine') then
         !DeAllocate data array
         call a%Memor%dealloc(oldnelem,a%NodalArrayDataForRefinement,'NodalArrayDataForRefinement','nsc_refine')
      end if 

   endif


   call a%SpecificNSCompRefine(itask) 

end subroutine   

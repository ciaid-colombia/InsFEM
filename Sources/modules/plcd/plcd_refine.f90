module Mod_plcd_Refine
   use typre
   use Mod_PLCD
   use Mod_AdaptiveInterface
   use Mod_Element
   use Mod_plcd_Material
   implicit none
   private
   public ElementMaterialTransferer


   type, extends(AdaptiveRefinerElementTransferer) :: ElementMaterialTransferer

      type(MatElementArray), pointer :: auxEMD(:) => NULL()
      class(PLCDProblem), pointer :: PLCD
      class(FiniteElement), pointer :: e => NULL()
contains
      procedure :: CopyElement => EMT_CopyElement
      procedure :: NodallyInterpolateElement => EMT_NodallyInterpolateElement
      procedure :: NodallyInterpolateElementFromOldBuffer => EMT_NodallyInterpolateElementFromOldBuffer
      

      procedure :: OldElementToBuffer => EMT_OldElementToBuffer
      procedure :: GetOldElementBufferSize => EMT_GetOldElementBufferSize

      procedure :: GetNewElementBufferSize => EMT_GetNewElementBufferSize
      procedure :: BufferToNewElement => EMT_BufferToNewElement
      procedure :: NewElementToBuffer => EMT_NewElementToBuffer

      procedure :: SetPLCD

   end type



contains

   subroutine EMT_CopyElement(a,oldielem,newielem)
      implicit none
      class(ElementMaterialTransferer), target :: a
      integer(ip) :: oldielem,newielem

      class(PLCDMaterial), pointer :: Material
      class(ElementMaterialData), pointer :: OldEMD, NewEMD

      integer(ip) :: pnode,pgaus

      if (newielem == 0) then
         write(*,*) 'hi newielem'
      end if

      call a%PLCD%ElementMaterialsData(oldielem)%p%GetMaterialPointer(Material)
      call a%PLCD%Mesh%GetElemArraySize(newielem,pnode,pgaus)
      call Material%CreateElementData(pgaus,a%auxEMD(newielem)%p)
      if (a%PLCD%kfl_ElementRotators == 1) call runend('refine not ready for element rotators')
      OldEMD => a%PLCD%ElementMaterialsData(oldielem)%p
      NewEMD => a%auxEMD(newielem)%p

      if (a%PLCD%IsElementSizeRequiredByEMDs) then
         call a%PLCD%Mesh%ElementLoad(newielem,a%e)
         !Compute Element length
         call a%e%elmdcg
         call a%e%elmlen

         Call NewEMD%SetElementSize(a%e%hleng(1))
      endif

      call OldEMD%CopyTo(NewEMD)
   end subroutine

   subroutine EMT_NodallyInterpolateElement(a,Refiner,oldelem,newelem)
      implicit none
      class(ElementMaterialTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) :: oldelem, newelem

      class(PLCDMaterial), pointer :: Material
      class(ElementMaterialData), pointer :: OldEMD, NewEMD
      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)
      real(rp), pointer :: OldNodalArray(:,:)
      integer(ip) :: pgaus,pnode

      call a%PLCD%ElementMaterialsData(oldelem)%p%GetMaterialPointer(Material)
      call a%PLCD%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      if (.not. associated(a%auxEMD(newelem)%p)) then
         call Material%CreateElementData(pgaus,a%auxEMD(newelem)%p)
         if (a%PLCD%kfl_ElementRotators == 1) call runend('refine not ready for element rotators')
      endif
      OldEMD => a%PLCD%ElementMaterialsData(oldelem)%p
      NewEMD => a%auxEMD(newelem)%p

      call OldEMD%GetHistoryVariablesNdofn(ndofn)

      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp
      call OldEMD%GetHistoryNodalArrayPointer(OldNodalArray)

      call Refiner%InterpolateElement(ndofn,OldNodalArray,NewNodalArray)

      call a%PLCD%Mesh%ElementLoad(newelem,a%e)

      if (a%PLCD%IsElementSizeRequiredByEMDs) then
         !Compute Element length
         call a%e%elmdcg
         call a%e%elmlen

         Call NewEMD%SetElementSize(a%e%hleng(1))
      endif


      call NewEMD%AddToGaussHistoryFromNodes(a%e,NewNodalArray)
      deallocate(NewNodalArray)
   end subroutine
   
   subroutine EMT_NodallyInterpolateElementFromOldBuffer(a,Refiner,newelem,OldBuffer)
      implicit none
      class(ElementMaterialTransferer) :: a
      class(AdaptiveRefinerInterface) :: Refiner
      integer(ip) ::  newelem
      real(rp), target :: OldBuffer(*)
      
      real(rp), pointer :: auxOldNodalArray(:,:)

      class(PLCDMaterial), pointer :: Material
      class(ElementMaterialData), pointer :: OldEMD, NewEMD
      integer(ip) :: ndofn

      real(rp), allocatable :: NewNodalArray(:,:)

      integer(ip) :: pgaus,pnode

      
      NewEMD => a%auxEMD(newelem)%p
      call NewEMD%GetHistoryVariablesNdofn(ndofn)
      
      call a%PLCD%Mesh%GetElemArraySize(newelem,pnode,pgaus)
      allocate(NewNodalArray(ndofn,pnode))
      NewNodalArray = 0.0_rp
      
      auxOldNodalArray(1:ndofn,1:pnode) => OldBuffer(2:ndofn*pnode+1)
      
      call Refiner%InterpolateElement(ndofn,auxOldNodalArray,NewNodalArray)
      
      call NewEMD%AddToGaussHistoryFromNodes(a%e,NewNodalArray)
      deallocate(NewNodalArray)
   end subroutine

   subroutine SetPLCD(a,PLCD)
      implicit none
      class(ElementMaterialTransferer) :: a
      class(PLCDProblem), target :: PLCD

      a%PLCD => PLCD
   end subroutine

   subroutine EMT_OldElementToBuffer(a,ielem,Buffer)
      class(ElementMaterialTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: idnum
      class(PLCDMaterial), pointer :: Material => NULL()

      call a%PLCD%ElementMaterialsData(ielem)%p%GetMaterialPointer(Material)

      call Material%GetIdNum(idnum)
      Buffer(1) = idnum

      call a%PLCD%ElementMaterialsData(ielem)%p%ToBuffer(Buffer(2))
   end subroutine

   subroutine EMT_NewElementToBuffer(a,ielem,Buffer)
      class(ElementMaterialTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: idnum
      class(PLCDMaterial), pointer :: Material => NULL()

      call a%auxEMD(ielem)%p%GetMaterialPointer(Material)

      call Material%GetIdNum(idnum)
      Buffer(1) = idnum

      call a%auxEMD(ielem)%p%ToBuffer(Buffer(2))
   end subroutine

   subroutine EMT_BufferToNewElement(a,ielem,Buffer)
      class(ElementMaterialTransferer) :: a
      integer(ip) :: ielem
      real(rp) :: Buffer(*)

      integer(ip) :: idnum,pnode,pgaus
      class(PLCDMaterial), pointer :: Material => NULL()

      idnum = Buffer(1)

      call a%PLCD%Mesh%GetElemArraySize(ielem,pnode,pgaus)

      Material => a%PLCD%Materials(idnum)%p
      call Material%CreateElementData(pgaus,a%auxEMD(ielem)%p)
      if (a%PLCD%kfl_ElementRotators == 1) call runend('refine not ready for element rotators')

      if (a%PLCD%IsElementSizeRequiredByEMDs) then
         call a%PLCD%Mesh%ElementLoad(ielem,a%e)
         !Compute Element length
         call a%e%elmdcg
         call a%e%elmlen

         call a%auxEMD(ielem)%p%SetElementSize(a%e%hleng(1))
      endif

      call a%auxEMD(ielem)%p%InitializeFromBuffer(Buffer(2))
   end subroutine

   subroutine EMT_GetOldElementBufferSize(a,ielem,BufferSize)
      class(ElementMaterialTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      call a%PLCD%ElementMaterialsData(ielem)%p%GetBufferSize(BufferSize)
      BufferSize = BufferSize+1 !in order to store the material number (unavoidable)
   end subroutine

   subroutine EMT_GetNewElementBufferSize(a,ielem,BufferSize)
      class(ElementMaterialTransferer) :: a
      integer(ip) :: ielem
      integer(ip) :: BufferSize

      call a%auxEMD(ielem)%p%GetBufferSize(BufferSize)
      BufferSize = BufferSize+1 !in order to store the material number (unavoidable)
   end subroutine

end module


subroutine plcd_Prerefine(a)
   use typre
   use Mod_PLCD
   use Mod_phpRefineArrays
   use Mod_plcd_Refine
   use Mod_Element
   implicit none
   class(PLCDProblem) :: a

   integer(ip) :: nelem,ielem
   class(FiniteElement), pointer :: e => NULL()

   !Things to be done prior to refinement
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Refine')

   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      call a%ElementMaterialsData(ielem)%p%AllocateAndComputeHistoryNodalArray(e)
   enddo
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Refine')

end subroutine

subroutine plcd_refine(a,itask)
   use typre
   use Mod_PLCD
   use Mod_phpRefineArrays
   use Mod_plcd_Refine
   use Mod_Element
   use Mod_plcd_StrainGenerator
   use Mod_int2str
   use Mod_r1pElementAllocation
   use Mod_plcd_SIMP_TopologyOptimization
   use Mod_plcd_TD_TopologyOptimization
   implicit none
   class(PLCDProblem) :: a
   character(6) :: itask

   type(MatElementArray), pointer :: auxEMD(:)
   type(ElementMaterialTransferer) :: MaterialTransferer

   integer(ip) :: nelem,ielem
   class(FiniteElement), pointer :: e => NULL()

   integer(ip) :: oldielem,oldnelem,ndime,auxpost,auxpost2

   procedure(), pointer :: GetStrain => NULL()
   procedure(), pointer :: AddSphericalComponent => NULL()
   procedure(), pointer :: GetStressTensor => NULL()
   procedure(), pointer :: GetStrainTensor => NULL()
   procedure(), pointer :: GetStressVector => NULL()
   procedure(), pointer :: GetStrainVector => NULL()
   integer(ip) :: vsize,istage,newnpoin,oldnpoin,imaterial

   integer(ip), allocatable :: auxkfl_fixno(:,:)
   integer(ip), allocatable :: auxarray(:,:), auxarray2(:,:)

   integer(ip) :: npoin

   integer(ip) :: ihpos,ipoin,jpoin,jnode,phang
   real(rp) :: auxvalue,jweight

   call php_RefineArrays(a,itask,a%displacement,'Displacement')
   if (a%kfl_TransientProblem /= 0) then
      call php_RefineArrays(a,itask,a%Velocity,'Velocity')
      call php_RefineArrays(a,itask,a%Acceleration,'Acceleration')
   endif
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%InterpolateHangingValues(ndime,a%Displacement(:,:,1))


   do istage = 1,a%NumberOfStages
      !-------------------------------------------------------------------------
      !Nodal forces require a special treatment
      !We are going to define an auxiliary integer array
      call a%Memor%alloc(size(a%Stages(istage)%NodalForces,1),size(a%Stages(istage)%NodalForces,2),auxarray,'auxarray','php_Refine')
      where (a%Stages(istage)%NodalForces /= 0.0_rp) auxarray = 1_ip


      call php_RefineArrays_b(a,itask,a%Stages(istage)%NodalForces,'NodalForces')

      call a%Memor%alloc(size(a%Stages(istage)%NodalForces,1),size(a%Stages(istage)%NodalForces,2),auxarray2,'auxarray2','php_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(size(a%Stages(istage)%NodalForces,1),auxarray,auxarray2,'minim')
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(size(a%Stages(istage)%NodalForces,1),auxarray,auxarray2)
      endif
      where(auxarray2 == 0) a%Stages(istage)%NodalForces = 0.0_rp

      call a%Memor%dealloc(size(auxarray,1),size(auxarray,2),auxarray,'auxarray','php_Refine')
      call a%Memor%dealloc(size(auxarray2,1),size(auxarray2,2),auxarray2,'auxarray2','php_Refine')
      !----------------------------------------------------------------------------



      call php_RefineArrays_b(a,itask,a%Stages(istage)%bvess,'bvess')

      !kfl_fixno
      call a%Mesh%GetNpoin(newnpoin)
      oldnpoin = size(a%Stages(istage)%kfl_fixno,2)
      call a%Memor%alloc(a%ndofbc,newnpoin,auxkfl_fixno,'kfl_fixno','php_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(a%ndofbc,a%Stages(istage)%kfl_fixno,auxkfl_fixno,'minim')
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(a%ndofbc,a%Stages(istage)%kfl_fixno,auxkfl_fixno)
      endif
      call move_alloc(auxkfl_fixno,a%Stages(istage)%kfl_fixno)
      call a%Memor%deallocObj(0,'kfl_fixno','php_Refine',ip*a%ndofbc*oldnpoin)
   enddo

   !Element Material Data
   if (itask == 'Refine') then
      call TransferElementMaterialData
   elseif (itask == 'Rebala') then
      call RebalanceElementMaterialData
   endif

   !Stresses
   call GetVoigtsize(ndime,vsize)
   call AllocateGetStrain(ndime,GetStrain,a%kfl_LargeStrains)
   call AllocateAddSphericalComponent(ndime,AddSphericalComponent)
   call AllocateGetStressTensor(ndime,GetStressTensor)
   call AllocateGetStrainTensor(ndime,GetStrainTensor)
   call AllocateGetStressVector(ndime,GetStressVector)
   call AllocateGetStrainVector(ndime,GetStrainVector)
   call DeAllocateR2pElement(a%Mesh,vsize,a%Stress,a%Memor,'stress')
   call AllocateR2pElement(a%Mesh,vsize,a%Stress,a%Memor,'stress')
   !Strains
   if(a%PostprocessStrain) then
      call DeAllocateR2pElement(a%Mesh,vsize,a%Strain,a%Memor,'strain')
      call AllocateR2pElement(a%Mesh,vsize,a%Strain,a%Memor,'strain')
   endif



   !Internal and external forces and residual
   call a%Mesh%GetNpoin(npoin)
   !These will be recomputed, so we only reallocate them
   call a%Memor%realloc(a%ndofn,npoin,a%InternalForcesVector,'InternalForcesVector','php_Refine')
   call a%Memor%realloc(a%ndofn,npoin,a%ExternalForcesVector,'ExternalForcesVector','php_Refine')
   call a%Memor%realloc(a%ndofn,npoin,a%ResidualForcesVector,'ResidualForcesVector','php_Refine')
   a%InternalForcesVector = 0.0_rp
   a%ExternalForcesVector = 0.0_rp


   if (a%UseSmoothedDisplacementGradient) then
      call php_RefineArrays_b(a,itask,a%SmoothedDisplacementGradient,'SmoothedDisplacementGradient')
   endif

    !UseUPFormulation
   if (a%UseUPFormulation) then
      call php_RefineArrays(a,itask,a%Pressure,'pressure')
      call a%Mesh%InterpolateHangingValues(1,a%Pressure(:,1))

      call php_RefineArrays_b(a,itask,a%UPResidualProjection,'UPResidualProjection')

      if (a%UPStoreSubscales) then
         call a%Mesh%GetNdime(ndime)
         call DeallocateR2pElement(a%Mesh,ndime,a%UPSubscales,a%Memor,'UPSubscales')
         call AllocateR2pElement(a%Mesh,ndime,a%UPSubscales,a%Memor,'UPSubscales')
      endif

   endif

   !Topology Optimization
   if (a%kfl_TopologyOptimization == 1) then
      call AdaptiveSIMP(a)
   elseif (a%kfl_TopologyOptimization == 2) then
      call AdaptiveTD(a,itask)
   end if


   !Do the endite operations to compute residual etc
   !But do not postprocess results at each iteration
   auxpost = a%kfl_PostprocessMatDataAtEachIteration
   auxpost2 = a%kfl_PostprocessDisplacementAtEachIteration
   a%kfl_PostprocessMatDataAtEachIteration = 0
   a%kfl_PostprocessDisplacementAtEachIteration = 0
   call a%Mesh%GetNdime(ndime)
   a%unkno(1:ndime,:) = a%Displacement(:,:,1)
   if (a%UseUPFormulation) a%unkno(ndime+1,:) = a%Pressure(:,1)

   a%itera = 0_ip  !So that some things are not done when computing forces (SIMP)
   call a%SpecificEndite(10_ip)
   a%ResidualForcesVector = a%ExternalForcesVector-a%InternalForcesVector

   a%kfl_PostprocessMatDataAtEachIteration = auxpost
   a%kfl_PostprocessDisplacementAtEachIteration = auxpost2

contains
   subroutine TransferElementMaterialData
      call a%Mesh%GetNelem(nelem)
      allocate(auxEMD(nelem))
      call a%Memor%allocObj(0,'ElementMaterialsData','plcd_refine',nelem)
      MaterialTransferer%auxEMD => auxEMD

      call MaterialTransferer%SetPLCD(a)
      call a%Mesh%ElementAlloc(MaterialTransferer%e,a%Memor,'DefaultRule','plcd_Refine')

      !Don't call us, we'll call YOU!
      call a%Refiner%ElementalRefine(MaterialTransferer)

      call a%Mesh%GetNelem(nelem)
      do ielem = 1,nelem
         call auxEMD(ielem)%p%SpecificHistoryFinalRefinementUpdate
      enddo


      oldnelem = size(a%ElementMaterialsData,1)
      do oldielem = 1,oldnelem
         call a%ElementMaterialsData(oldielem)%p%DeallocateNodalArray
         call a%ElementMaterialsData(oldielem)%p%Finalize
         deallocate(a%ElementMaterialsData(oldielem)%p)
      enddo
      deallocate(a%ElementMaterialsData)
      call a%Memor%deallocObj(0,'ElementMaterialsData','plcd_refine',oldnelem)
      a%ElementMaterialsData => auxEMD

      call a%Mesh%ElementDeAlloc(MaterialTransferer%e,a%Memor,'DefaultRule','plcd_Refine')
   end subroutine




   subroutine RebalanceElementMaterialData
      call a%Mesh%GetNelem(nelem)
      allocate(auxEMD(nelem))
      call a%Memor%allocObj(0,'ElementMaterialsData','plcd_refine',nelem)
      MaterialTransferer%auxEMD => auxEMD

      call MaterialTransferer%SetPLCD(a)
      call a%Mesh%ElementAlloc(MaterialTransferer%e,a%Memor,'DefaultRule','plcd_Refine')

      call a%Refiner%ElementalRebalance(MaterialTransferer)
      oldnelem = size(a%ElementMaterialsData,1)
      do oldielem = 1,oldnelem
         call a%ElementMaterialsData(oldielem)%p%Finalize
         deallocate(a%ElementMaterialsData(oldielem)%p)
      enddo
      deallocate(a%ElementMaterialsData)
      call a%Memor%deallocObj(0,'ElementMaterialsData','plcd_refine',oldnelem)
      a%ElementMaterialsData => auxEMD

      call a%Mesh%ElementDeAlloc(MaterialTransferer%e,a%Memor,'DefaultRule','plcd_Refine')
   end subroutine







end subroutine

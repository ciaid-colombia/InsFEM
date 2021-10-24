module Mod_plcd_CompositeMaterial
   use typre
   use Mod_plcd
   use Mod_plcd_Material
   use Mod_Listen
   use Mod_Element
   use Mod_plcd_ReadMaterials
   use Mod_r1pElementAllocation
   implicit none
   private
   public CompositeMaterial
   
   type :: CompositeStress
        type(r2p), allocatable :: SimpleMaterialStress(:)
   end type
   
   type, extends(PLCDMaterial) :: CompositeMaterial
      integer(ip) :: NumberOfMaterials
      real(rp), allocatable :: VolumetricParticipations(:)
      type(MatArray), allocatable :: Materials(:)
      
      type(CompositeStress), allocatable :: CompositeStresses(:)
      
      
      character(5), allocatable :: AuxRead_MaterialTypeList(:)
      
   
contains
      procedure :: SpecificCreateElementData => CM_SpecificCreateElementData
      procedure :: SpecificReadData => CM_SpecificReadData
      procedure :: SpecificScatterData => CM_SpecificScatterData
      procedure :: SpecificSetup       => CM_SpecificSetup
      !Postprocess
      procedure :: SetupPostprocess    => CM_SetupPostprocess
      procedure :: DoPostprocess       => CM_DoPostprocess
      procedure :: IsElementSizeRequiredByEMDs => CM_IsElementSizeRequiredByEMDs
     
   end type
   
   

   
   type, extends(ElementMaterialData) :: CompositeMaterial_ED
      
      class(CompositeMaterial), pointer :: Material => NULL()
      
      type(MatElementArray), allocatable :: EMDs(:)
      real(rp), allocatable :: C(:,:,:), S(:,:)
      integer(ip) :: pgaus
         
contains   

      procedure :: ComputeHistoryAndConstitutiveTensor => CM_ComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress => CM_ComputeStress
      procedure :: GetConstitutiveTensorPointer => CM_GetConstitutiveTensorPointer
      procedure :: GetStressTensorPointer => CM_GetStressTensorPointer
      procedure :: SetElementSize => CM_SetElementSize
      procedure :: MoveHistoryVariablesToConverged => CM_MoveHistoryVariablesToConverged
      !Postprocess
      procedure :: ContributeToPostprocess => CM_ContributeToPostprocess
      procedure :: GetMaterialPointer
      procedure :: Finalize
      procedure :: CopyTo
      procedure :: SpecificHistoryAddFromGaussArray
      procedure :: SpecificHistoryGetToGaussArray
      procedure :: SpecificHistoryFinalRefinementUpdate
      procedure :: GetHistoryVariablesNdofn
      
      procedure :: GetBufferSize
      procedure :: ToBuffer
      procedure :: InitializeFromBuffer
      
      !For u-p elements
      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetSecantCNorm
      procedure :: SetPressureForComputeHistoryAndConstitutiveTensor
   end type
   
   

contains

   subroutine CM_SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(CompositeMaterial), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData
      
      type(CompositeMaterial_ED), pointer :: CM_ElementData
      
      integer(ip) :: igaus,imaterial
      
      allocate(CM_ElementData)
      CM_ElementData%Material => a
      
      allocate(CM_ElementData%EMDs(a%NumberOfMaterials))
      do imaterial = 1,a%NumberOfMaterials
         call a%Materials(imaterial)%p%SpecificCreateElementData(pgaus,CM_ElementData%EMDs(imaterial)%p)
      enddo
      
      allocate(CM_ElementData%C(a%CT%VoigtSize,a%CT%VoigtSize,pgaus))
      allocate(CM_ElementData%S(a%CT%VoigtSize,pgaus))
      CM_ElementData%pgaus = pgaus
      
      ElementData => CM_ElementData
   end subroutine

   
   subroutine CM_SpecificReadData(a,Listener)
      use Mod_Memor
      implicit none
      class(CompositeMaterial) :: a
      type(ListenFile) :: Listener
      
      type(MemoryMan) :: Memor
      
      call ReadMaterials_Root(Listener,a%NumberOfMaterials,a%Materials,a%AuxRead_MaterialTypeList,a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank,a%ndime,a%kfl_LargeStrains,Memor)
      
      if (Listener%words(1) == 'NUMBE' .and. Listener%words(3) == 'MATER') then
         allocate(a%VolumetricParticipations(a%NumberOfMaterials))
         a%VolumetricParticipations = 1.0_rp/a%NumberOfMaterials
      elseif (Listener%words(1) == 'VOLUM') then
         a%VolumetricParticipations = Listener%param(1:a%NumberOfMaterials)
      endif
      
      
   end subroutine
   
   subroutine CM_SpecificScatterData(a)
      use MPI
      use Mod_Memor
      class(CompositeMaterial) :: a
      
      integer(ip) :: ierr
      
      type(MemoryMan) :: Memor
      
      write(*,*) 'Memor should be passed to the material'
      write(*,*) 'right now Memor is not the correct one'
      read(*,*)
      call ReadMaterials_Scatter(a%NumberOfMaterials,a%Materials,a%AuxRead_MaterialTypeList,a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank,a%ndime,a%kfl_LargeStrains,Memor)
      
      if (a%MPIrank /= a%MPIroot) then
         allocate(a%VolumetricParticipations(a%NumberOfMaterials))
      endif
      call MPI_BCAST(a%VolumetricParticipations,a%NumberOfMaterials,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
      
   end subroutine
   
   subroutine CM_SpecificSetup(a)
      class(CompositeMaterial):: a
      
      integer(ip) :: imaterial
      
      do imaterial = 1,a%NumberOfMaterials
         call a%Materials(imaterial)%p%Setup
      enddo
   end subroutine
   
   
   subroutine CM_ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      class(CompositeMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      integer(ip) :: imaterial
      real(rp) :: strain(a%Material%CT%VoigtSize)
      
      !Parallel
      call a%Material%CT%GetStrain(gradDisp,strain)
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%ComputeHistoryAndConstitutiveTensor(igaus,gradDisp)
      enddo
   end subroutine
   
   subroutine CM_ComputeStress(a,igaus,gradDisp)
      class(CompositeMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      integer(ip) :: imaterial
      real(rp) :: strain(a%Material%CT%VoigtSize)
      
      !Parallel
      call a%Material%CT%GetStrain(gradDisp,strain)
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%ComputeStress(igaus,gradDisp)
      enddo
   end subroutine
   
   subroutine CM_GetStressTensorPointer(a,igaus,S)
      class(CompositeMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: S(:)
      
      integer(ip) :: imaterial
      
      a%S(:,igaus) = 0.0_rp
      do imaterial = 1,a%Material%NumberOfMaterials  
         call a%EMDs(imaterial)%p%GetStressTensorPointer(igaus,S)
         a%S(:,igaus) = a%S(:,igaus)+S*a%Material%VolumetricParticipations(imaterial)
      enddo
      S => a%S(1:a%Material%CT%VoigtSize,igaus)
   end subroutine
   
   subroutine CM_GetConstitutiveTensorPointer(a,igaus,C)
      class(CompositeMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)
      
      integer(ip) :: imaterial
      
      a%C(:,:,igaus) = 0.0_rp
      do imaterial = 1,a%Material%NumberOfMaterials  
         call a%EMDs(imaterial)%p%GetConstitutiveTensorPointer(igaus,C)
         a%C(:,:,igaus) = a%C(:,:,igaus)+C*a%Material%VolumetricParticipations(imaterial)
      enddo
      C => a%C(:,:,igaus)
      
   end subroutine
   
   subroutine   CM_SetElementSize(a,h)
      class(CompositeMaterial_ED) :: a
      real(rp) :: h
      
      integer(ip) :: imaterial
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%SetElementSize(h)
      enddo
      
   end subroutine
   
   subroutine CM_MoveHistoryVariablesToConverged(a)
      class(CompositeMaterial_ED) :: a
      
      integer(ip) :: imaterial
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%MoveHistoryVariablesToConverged
      enddo
      
   end subroutine
   
   subroutine CM_SetupPostprocess(a,Mesh,Memor)
      use Mod_Mesh
      use Mod_Memor
      use Mod_RpMeshAllocator
      implicit none
      class(CompositeMaterial) :: a
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      
      integer(ip) :: imaterial,nelem
      
      do imaterial = 1,a%NumberOfMaterials
         call a%Materials(imaterial)%p%SetupPostprocess(Mesh,Memor)
      enddo
      
      allocate(a%CompositeStresses(a%NumberOfMaterials))
      call Memor%allocObj(0,'CompositeStresses','CompoMat',a%NumberOfMaterials)
      call Mesh%GetNelem(nelem)
      do imaterial = 1,a%NumberOfMaterials
        call AllocateR2pElement(Mesh,a%CT%VoigtSize,a%CompositeStresses(imaterial)%SimpleMaterialStress,Memor,'SimpleStress')
        !call Memor%alloc(nelem,a%CompositeStresses(imaterial)%SimpleMaterialStress,'Simplemat','Compomat')
      enddo
      
      
      
   end subroutine
   
   subroutine CM_ContributeToPostprocess(a,ielem)
      class(CompositeMaterial_ED) :: a
      integer(ip) :: ielem
      
      integer(ip) :: igaus
      
      integer(ip) :: imaterial
      real(rp), pointer :: S(:)
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%ContributeToPostprocess(ielem)
         
         do igaus = 1,a%pgaus
            call a%EMDs(imaterial)%p%GetStressTensorPointer(igaus,S)
            a%Material%CompositeStresses(imaterial)%SimpleMaterialStress(ielem)%a(:,igaus) = S
         enddo
      enddo
      
   end subroutine
   
   subroutine CM_DoPostprocess(a,istep,ctime,FilePostpr,Mesh,Memor,iiter_char)
      use Mod_Mesh
      use Mod_Postpr
      use Mod_Memor
      use Mod_RpMeshAllocator  
      use Mod_int2str
      implicit none
      class(CompositeMaterial) :: a
      class(PostprFile) :: FilePostpr
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: istep
      real(rp) :: ctime
      character(5) :: iiter_char

      integer(ip) :: imaterial
      
      do imaterial = 1,a%NumberOfMaterials
         call a%Materials(imaterial)%p%DoPostprocess(istep,ctime,FilePostpr,Mesh,Memor,iiter_char)
         
         call FilePostpr%postgp(a%CompositeStresses(imaterial)%SimpleMaterialStress,'Stress//'//(adjustl(trim(a%Materials(imaterial)%p%MaterialName))),istep,ctime,Mesh,'voigtT')
         
      enddo
            
      do imaterial = 1,a%NumberOfMaterials
        call DeAllocateR2pElement(Mesh,a%CT%VoigtSize,a%CompositeStresses(imaterial)%SimpleMaterialStress,Memor,'SimpleStress')
      enddo
      deallocate(a%CompositeStresses)
      call Memor%deallocObj(0,'CompositeStresses','CompoMat',a%NumberOfMaterials)
      
   end subroutine
   
   subroutine GetMaterialPointer(a,Material)
      class(CompositeMaterial_ED), target :: a
      class(PLCDMaterial), pointer :: Material
      
      Material => a%Material
      
   end subroutine
   
   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(CompositeMaterial_ED) :: a
      integer(ip) :: ndofn
      
      integer(ip) :: ndofn1,imaterial
      
      ndofn = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%GetHistoryVariablesNdofn(ndofn1)
         ndofn = ndofn+ndofn1
      enddo
      
   end subroutine

   
   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(CompositeMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)
      
      integer(ip) :: imaterial,ndofn,ndofn0
      
      ndofn0 = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%GetHistoryVariablesNdofn(ndofn)
         call a%EMDs(imaterial)%p%SpecificHistoryGetToGaussArray(GaussArray(ndofn0+1:ndofn0+ndofn,:))
         
         
         ndofn0 = ndofn0+ndofn
      enddo
   end subroutine
   
   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(CompositeMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)
      
      integer(ip) :: imaterial,ndofn,ndofn0
      
      ndofn0 = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%GetHistoryVariablesNdofn(ndofn)
         call a%EMDs(imaterial)%p%SpecificHistoryAddFromGaussArray(GaussArray(ndofn0+1:ndofn0+ndofn,:))
         
         
         ndofn0 = ndofn0+ndofn
      enddo
      
     
      
   end subroutine
   
   subroutine CopyTo(a,EMD)
      use Mod_plcd_Material
      class(CompositeMaterial_ED) :: a
      class(ElementMaterialData), pointer :: EMD 
      
      integer(ip) :: imaterial
      
      select type(EMD)
      type is (CompositeMaterial_ED)
      
         do imaterial = 1,a%Material%NumberOfMaterials
            call a%EMDs(imaterial)%p%CopyTo(EMD%EMDs(imaterial)%p)
         enddo
         
      end select
      
   end subroutine
   
   subroutine Finalize(a)
      class(CompositeMaterial_ED) :: a
      
      integer(ip) :: imaterial
      
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%Finalize
      enddo
      deallocate(a%EMDs)
   
   end subroutine
   
   
   subroutine CM_IsElementSizeRequiredByEMDs(a,isRequired)
      implicit none
      class(CompositeMaterial), target :: a
      logical :: IsRequired
      
      integer(ip) :: imaterial
      logical :: auxisRequired
      
      do imaterial = 1,a%NumberOfMaterials
         call a%Materials(imaterial)%p%IsElementSizeRequiredByEMDs(auxisRequired)
         if (auxisRequired .eqv. .true.) isRequired = .true.
      enddo
   end subroutine  
   
    subroutine GetBufferSize(a,BufferSize)
      class(CompositeMaterial_ED) :: a
      integer(ip) :: BufferSize
      
      integer(ip) :: imaterial
      integer(ip) :: auxBufferSize
      
      BufferSize = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%GetBufferSize(auxBufferSize)
         BufferSize = BufferSize + auxBufferSize
      enddo

   end subroutine
   
   subroutine ToBuffer(a,Buffer)
      class(CompositeMaterial_ED) :: a
      real(rp) :: Buffer(*)
      
      integer(ip) :: imaterial
      integer(ip) :: BufferPos,BufferSize
      
      BufferPos = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%ToBuffer(Buffer(BufferPos))
         call a%EMDs(imaterial)%p%GetBufferSize(BufferSize)
         BufferPos = BufferPos + BufferSize
      enddo
   end subroutine
   
   subroutine InitializeFromBuffer(a,Buffer)
      class(CompositeMaterial_ED) :: a
      real(rp) :: Buffer(*)
      
      integer(ip) :: imaterial
      integer(ip) :: BufferPos,BufferSize
      
      BufferPos = 0
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%InitializeFromBuffer(Buffer(BufferPos))
         call a%EMDs(imaterial)%p%GetBufferSize(BufferSize)
         BufferPos = BufferPos + BufferSize
      enddo
      
   end subroutine
   
   
   !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(CompositeMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: invK
      
      call runend('Not implemented')
      !call a%Material%CT%GetInverseVolumetricDeformationModulus(invK)
   end subroutine
   
   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(CompositeMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: SecantCNorm

      call runend('Not implemented')
      !call a%Material%CT%GetCNorm(SecantCNorm)
      
   end subroutine
  
   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(CompositeMaterial_ED) :: a
      real(rp) :: Pressure
      
      integer(ip) :: imaterial
      
      !Set the pressure for compound material
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%SetPressureForComputeHistoryAndConstitutiveTensor(pressure)
      enddo
      
      
      
      
   end subroutine
   
   subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(CompositeMaterial_ED) :: a
      
      integer(ip) :: imaterial
      
      !Set the pressure for compound material
      do imaterial = 1,a%Material%NumberOfMaterials
         call a%EMDs(imaterial)%p%SpecificHistoryFinalRefinementUpdate
      enddo
      
      
   end subroutine
   
   
end module
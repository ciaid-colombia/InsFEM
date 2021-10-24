module Mod_plcd_Material
   use typre
   use Mod_Listen
   use Mod_MPIObject
   use Mod_plcd_ConstitutiveTensor
   use Mod_plcd_ConstitutiveTensorFactory
   use Mod_plcd_StrainGenerator
   use Mod_Element
   use Mod_plcd_Rotation
   use Mod_plcd_RotationFactory
   use Mod_plcd_Rotators
   use Mod_Memor
   implicit none
   private
   public PLCDMaterial, ElementMaterialData

   type, extends(MPIObject), abstract :: PLCDMaterial
      integer(ip) :: ndime
      integer(ip) :: idnum
      character(len=5) :: MaterialName

      real(rp) :: density = 1.0_rp
      integer(ip) :: TangentTensorType = 2 !0: Initial  1: secant  2: tangent
      character(5) :: CTType = 'NOTDE', CTType2 = 'NOTDE', RerSysType = 'GLOBA'
      real(rp), pointer :: C(:,:)
      class(ConstitutiveTensor), pointer :: CT => NULL()
      class(Rotation), pointer :: Rotator => NULL()
      
      Type(Rotator_Not) :: DonotRotate

      logical :: UseUPFormulation = .false.
      
      integer(ip) :: kfl_LargeStrains = 0 !0: Linear Elastic, 1: Updated Formulation NonLinear Elastic, 2: Total Formulation Nonlinear Elastic
    
contains
      procedure                                 :: SetNdime
      procedure                                 :: SetFlagLargeStrains
      procedure                                 :: ReadData
      procedure                                 :: ScatterData
      procedure(SpecificReadData), deferred     :: SpecificReadData
      procedure(SpecificScatterData), deferred  :: SpecificScatterData
      procedure                                 :: CreateElementData
      procedure(SpecificCreateElementData), deferred    :: SpecificCreateElementData
      procedure                                 :: Setup
      procedure(SpecificSetup), deferred        :: SpecificSetup

      procedure                                 :: SetupPostprocess
      procedure                                 :: DoPostprocess
      procedure(isElementSizeRequiredByEMDs), deferred :: IsElementSizeRequiredByEMDs
      procedure                                 :: SetIdnum
      procedure                                 :: SetMaterialName
      procedure                                 :: GetIdnum

      !for u-p elements
      procedure                                 :: SetUseUPFormulation
   end type

   
   type, abstract :: ElementMaterialData
      !For adaptive mesh refinement
      real(rp), allocatable :: HistoryNodalArray(:,:)
      class(Rotation), pointer :: Rotator => NULL()

contains
      procedure(ComputeStress), deferred        :: ComputeStress
      procedure(ComputeHistoryAndConstitutiveTensor), deferred :: ComputeHistoryAndConstitutiveTensor
      procedure(GetStressTensorPointer), deferred  :: GetStressTensorPointer
      procedure(GetConstitutiveTensorPointer), deferred :: GetConstitutiveTensorPointer
      procedure                                 :: MoveHistoryVariablesToConverged
      procedure(GetMaterialPointer), deferred   :: GetMaterialPointer
      procedure                                 :: AllocateRotator
      procedure                                 :: SpecificReadElementData

      procedure                                 :: ContributeToPostprocess

      !Adaptivity
      procedure                                 :: GetHistoryNodalArrayPointer
      procedure                                 :: AddToGaussHistoryFromNodes
      procedure                                 :: GetHistoryVariablesNdofn
      procedure                                 :: SpecificHistoryGetToGaussArray
      procedure                                 :: SpecificHistoryAddFromGaussArray
      procedure                                 :: SpecificHistoryFinalRefinementUpdate
      procedure                                 :: CopyTo
      procedure                                 :: AllocateAndComputeHistoryNodalArray
      procedure                                 :: DeallocateNodalArray
      procedure                                 :: Finalize
      procedure                                 :: GetBufferSize
      procedure                                 :: InitializeFromBuffer
      procedure                                 :: ToBuffer
      procedure                                 :: GetComparisonCriteriaRatio

      !Optional data which needs some kind of interface
      procedure                                 :: SetElementSize

      !For u-p elements
      procedure                                 :: GetInverseVolumetricDeformationModulus
      procedure                                 :: GetSecantCNorm
      procedure                                 :: SetPressureForComputeHistoryAndConstitutiveTensor
   end type


   abstract interface
      subroutine SpecificCreateElementData(a,pgaus,ElementData)
         import
         implicit none
         class(PLCDMaterial), target :: a
         integer(ip) :: pgaus
         class(ElementMaterialData), pointer :: ElementData
      end subroutine

      subroutine SpecificReadData(a,Listener)
         import
         implicit none
         class(PLCDMaterial) :: a
         type(ListenFile) :: Listener
      end subroutine

      subroutine SpecificScatterData(a)
         import
         implicit none
         class(PLCDMaterial) :: a

      end subroutine

      subroutine SpecificSetup(a)
         import
         implicit none
         class(PLCDMaterial) :: a

      end subroutine

      subroutine GetStressTensorPointer(a,igaus,S)
         import
         implicit none
         class(ElementMaterialData), target :: a
         integer(ip) :: igaus
         real(rp), pointer :: S(:)
      end subroutine

      subroutine GetConstitutiveTensorPointer(a,igaus,C)
         import
         implicit none
         class(ElementMaterialData), target:: a
         integer(ip) :: igaus
         real(rp), pointer :: C(:,:)
      end subroutine

      subroutine ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
         import
         implicit none
         class(ElementMaterialData), target :: a
         integer(ip) :: igaus
         real(rp) :: gradDisp(:,:)
      end subroutine
      
      subroutine ComputeStress(a,igaus,gradDisp)
         import
         implicit none
         class(ElementMaterialData), target :: a
         integer(ip) :: igaus
         real(rp) :: gradDisp(:,:)
      end subroutine

      subroutine GetMaterialPointer(a,Material)
         import
         implicit none
         class(ElementMaterialData), target :: a
         class(PLCDMaterial), pointer :: Material
      end subroutine

      subroutine IsElementSizeRequiredByEMDs(a,isRequired)
         import
         implicit none
         class(PLCDMaterial), target :: a
         logical :: isRequired
      end subroutine

   end interface

contains

   subroutine CreateElementData(a,pgaus,ElementData)
         implicit none
         class(PLCDMaterial), target :: a
         integer(ip) :: pgaus
         class(ElementMaterialData), pointer :: ElementData
         
         call a%SpecificCreateElementData(pgaus,ElementData)
         
         ElementData%Rotator => a%DonotRotate
         
   end subroutine
   
   subroutine SetNdime(a,ndime)
      implicit none
      class(PLCDMaterial) :: a
      integer(ip) :: ndime

      a%ndime = ndime
   end subroutine
   
   subroutine SetFlagLargeStrains(a,flag)
      implicit none
      class(PLCDMaterial) :: a
      integer(ip) :: flag
      
      a%kfl_LargeStrains = flag
   end subroutine

   subroutine ReadData(a,Listener)
      implicit none
      class(PLCDMaterial), target :: a
      type(ListenFile) :: Listener
      real(rp) :: auxangle(3)
      real(rp), pointer :: new(:)
      type(MemoryMan) :: Memor
      
      a%Rotator => a%DonotRotate

      do while(Listener%words(1)/='ENDMA')
         call Listener%listen('plcd_reaphy')
         if(Listener%words(1) == 'DENSI') then
            a%density = Listener%param(1)
         
         elseif(Listener%words(1) == 'TANGT') then
            if (Listener%words(2) == 'INITI') then
               a%TangentTensorType = 0
            elseif (Listener%words(2) == 'SECAN') then
               a%TangentTensorType = 1
            elseif (Listener%words(2) == 'TANGE') then
               a%TangentTensorType = 2
            endif
         elseif(Listener%words(1) == 'CONST') then
            a%CTType = Listener%words(2)
            a%CTType2 = Listener%words(3)
            call AllocateConstitutiveTensor(a%ndime,a%kfl_LargeStrains,a%CTType,a%CTType2,a%CT)
            call a%CT%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
            call a%CT%ReadData(Listener)
         elseif(Listener%words(1) == 'REFSY' .and. Listener%words(2) == 'LOCAL') then
            a%RerSysType = Listener%words(2)
            call AllocatePLCDRotator(a%RerSysType,a%Rotator)
            call a%Rotator%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
            call a%Rotator%SetNdime(a%ndime)
            auxangle = Listener%param(2:4)
            call a%Rotator%SetAngle(auxangle,Memor)
         else
            call a%SpecificReadData(Listener)
         
         endif

      enddo
      !This is just so that in the external loop there are no problems
      Listener%words(1) = '     '
      
   end subroutine

   subroutine ScatterData(a)
      use MPI
      implicit none
      class(PLCDMaterial), target :: a

      integer(ip) :: ierr
      call MPI_BCAST(a%idnum,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%density,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%TangentTensorType,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%MaterialName,len(a%MaterialName),MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      
      call MPI_BCAST(a%CTType,len(a%CTType),MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%CTType2,len(a%CTType2),MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      
      call MPI_BCAST(a%RerSysType,len(a%RerSysType),MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)

      if (a%MPIrank /= a%MPIroot) then
         call AllocateConstitutiveTensor(a%ndime,a%kfl_LargeStrains,a%CTType,a%CTType2,a%CT)
         call a%CT%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
         a%Rotator => a%DonotRotate
         
         if(a%RerSysType == 'LOCAL') then
            call AllocatePLCDRotator(a%RerSysType,a%Rotator)
            call a%Rotator%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
            call a%Rotator%SetNdime(a%ndime)
         endif
      endif
      
      
      call a%CT%ScatterData
      call a%CT%CheckData
     
      call a%Rotator%ScatterData
            
      call a%SpecificScatterData

   end subroutine

   subroutine Setup(a)
      implicit none
      class(PLCDMaterial) :: a
      type(MemoryMan) :: Memor
      
      call a%CT%Build
      call a%CT%GetPointerToC(a%C)
      
      call a%Rotator%Initialize(Memor)
      
      call a%SpecificSetup
   end subroutine

   subroutine SetIdnum(a,idnum)
      implicit none
      class(PLCDMaterial) :: a
      integer(ip) :: idnum

      a%idnum = idnum
   end subroutine
   
   subroutine SetMaterialName(a,MaterialName)
      implicit none
      class(PLCDMaterial) :: a
      character(5) :: MaterialName

      a%MaterialName = MaterialName
   end subroutine

   subroutine GetIdnum(a,idnum)
      implicit none
      class(PLCDMaterial) :: a
      integer(ip) :: idnum

      idnum = a%idnum
   end subroutine

   subroutine SetElementSize(a,h)
      implicit none
      class(ElementMaterialData) :: a
      real(rp) :: h

   end subroutine

   subroutine MoveHistoryVariablesToConverged(a)
      class(ElementMaterialData) :: a
   end subroutine

   subroutine SetupPostprocess(a,Mesh,Memor)
      use Mod_Mesh
      use Mod_Memor
      use Mod_Element
      implicit none
      class(PLCDMaterial) :: a
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor


   end subroutine

   subroutine ContributeToPostprocess(a,ielem)
      class(ElementMaterialData) :: a
      integer(ip) :: ielem

   end subroutine

   subroutine DoPostprocess(a,istep,ctime,FilePostpr,Mesh,Memor,iiter_char)
      use Mod_Mesh
      use Mod_Postpr
      use Mod_Memor
      use Mod_Element
      implicit none
      class(PLCDMaterial) :: a
      class(PostprFile) :: FilePostpr
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: istep
      real(rp) :: ctime
      character(5) :: iiter_char

   end subroutine
   
   subroutine AllocateRotator(a,ndime,angle,Memor)
      implicit none
      class(ElementMaterialData), target :: a
      integer(ip), intent(in) :: ndime
      real(rp), intent(in) :: angle(:)
      type(MemoryMan), intent(in) :: Memor
      
      logical :: isZeroAngle
      
      isZeroAngle = .false.
      
      if (norm2(angle) == 0.0_rp) isZeroAngle = .true.
      
      if (.not. isZeroAngle) then
      
         call AllocatePLCDRotator('LOCAL',a%Rotator)
         call a%Rotator%SetNdime(ndime)
         call a%Rotator%SetAngle(angle,Memor)
         call a%Rotator%Initialize(Memor)
      endif
   end subroutine
   
   subroutine SpecificReadElementData(a,param)
      class(ElementMaterialData) :: a
      real(rp) :: param(*)
      
   end subroutine

   subroutine AllocateAndComputeHistoryNodalArray(a,e)
      class(ElementMaterialData) :: a
      class(FiniteElement) :: e

      real(rp), allocatable :: GaussArray(:,:)

      integer(ip) :: ndofn


      call a%GetHistoryVariablesNdofn(ndofn)

      allocate(a%HistoryNodalArray(ndofn,e%pnode))
      allocate(GaussArray(ndofn,e%pgaus))

      call a%SpecificHistoryGetToGaussArray(GaussArray)
      call e%GaussToNodes(GaussArray,a%HistoryNodalArray)


      deallocate(GaussArray)
   end subroutine

   subroutine DeallocateNodalArray(a)
      class(ElementMaterialData) :: a

      deallocate(a%HistoryNodalArray)
   end subroutine

   subroutine GetHistoryNodalArrayPointer(a,NodalArray)
      class(ElementMaterialData), target :: a
      real(rp), pointer :: NodalArray(:,:)

      NodalArray => a%HistoryNodalArray
   end subroutine

   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(ElementMaterialData) :: a
      integer(ip) :: ndofn

      call runend('GetHistoryVariablesNdofn not defined')
   end subroutine

   subroutine AddToGaussHistoryFromNodes(a,e,NodalArray)
      class(ElementMaterialData) :: a
      class(FiniteElement) :: e
      real(rp) :: NodalArray(:,:)

      real(rp), allocatable :: GaussArray(:,:)
      integer(ip) :: igaus,ndofn


      call a%GetHistoryVariablesNdofn(ndofn)
      allocate(GaussArray(ndofn,e%pgaus))

      if (ndofn > 0) then
         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%interpg(ndofn,NodalArray,GaussArray(:,e%igaus))
         enddo
      endif

      call a%SpecificHistoryAddFromGaussArray(GaussArray)

      deallocate(GaussArray)



   end subroutine

   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(ElementMaterialData) :: a
      real(rp) :: GaussArray(:,:)

      call runend('SpecificHistoryGetFromGaussArray not defined')
   end subroutine

   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(ElementMaterialData) :: a
      real(rp) :: GaussArray(:,:)

      call runend('SpecificHistoryPutInGaussArray not defined')
   end subroutine

   subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(ElementMaterialData) :: a


      call runend('SpecificFinalRefinementUpdate not defined')
   end subroutine

   subroutine CopyTo(a,EMD)
      class(ElementMaterialData) :: a
      class(ElementMaterialData), pointer :: EMD

      call runend('CopyTo not implemented for this specific element material data')
   end subroutine

   subroutine Finalize(a)
      class(ElementMaterialData) :: a

      call runend('Finalize not implemented for this element material data')
   end subroutine

   subroutine GetBufferSize(a,BufferSize)
      class(ElementMaterialData) :: a
      integer(ip) :: BufferSize

      call runend('GetBufferSize not implemented for this element material data')
   end subroutine

   subroutine ToBuffer(a,Buffer)
      class(ElementMaterialData) :: a
      real(rp) :: Buffer(*)

      call runend('ToBuffer not implemented for this element material data')
   end subroutine

   subroutine InitializeFromBuffer(a,Buffer)
      class(ElementMaterialData) :: a
      real(rp) :: Buffer(*)

      call runend('InitializeFromBuffer not implemented for this element material data')
   end subroutine

   subroutine GetComparisonCriteriaRatio(a,ComparisonRatio)
      class(ElementMaterialData) :: a
      real(rp) :: ComparisonRatio

      call runend('GetComparisonCriteriaRatio not implemented for this element material data')
   end subroutine

   !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(ElementMaterialData) :: a
      integer(ip) :: igaus
      real(rp) :: invK


      call runend('GetInverseVolumetricDeformationModulus not implemented for this element material data')
   end subroutine

   subroutine SetUseUPFormulation(a,UseUPFormulation)
      class(PLCDMaterial) :: a
      logical :: UseUPFormulation

      a%UseUPFormulation = UseUPFormulation
      if (a%UseUPFormulation) call a%CT%SetOnlyDeviatoric(.true.)
   end subroutine

   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(ElementMaterialData) :: a
      integer(ip) :: igaus
      real(rp) :: SecantCNorm

      call runend('GetSecantCNorm not implemented for this element material data')

   end subroutine

   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(ElementMaterialData) :: a
      real(rp) :: Pressure

      call runend('SetPressure Not implemented for this element material data')

   end subroutine
end module

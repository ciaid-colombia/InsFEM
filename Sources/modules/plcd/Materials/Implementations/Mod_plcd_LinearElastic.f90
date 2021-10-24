module Mod_plcd_LinearElastic
   use typre
   use Mod_plcd_Material
   use Mod_Listen
   use Mod_Element
   implicit none
   private
   public LinearElastic, LinearElastic_ED

   type, extends(PLCDMaterial) :: LinearElastic

contains
      procedure :: SpecificCreateElementData => LE_SpecificCreateElementData
      procedure :: SpecificReadData => LE_SpecificReadData
      procedure :: SpecificScatterData => LE_SpecificScatterData
      procedure :: SpecificSetup       => NULLSUB
      procedure :: IsElementSizeRequiredByEMDs => LE_IsElementSizeRequiredByEMDs


   end type

   type, extends(ElementMaterialData) :: LinearElastic_ED

      class(LinearElastic), pointer :: Material => NULL()
      real(rp), allocatable :: S(:,:)

contains

      procedure :: ComputeHistoryAndConstitutiveTensor => LE_ComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress => LE_ComputeStress
      procedure :: GetConstitutiveTensorPointer => LE_GetConstitutiveTensorPointer
      procedure :: GetStressTensorPointer => LE_GetStressTensorPointer
      procedure :: GetMaterialPointer

      !Adaptivity
      procedure :: GetHistoryVariablesNdofn
      procedure :: Finalize
      procedure :: CopyTo
      procedure :: SpecificHistoryAddFromGaussArray
      procedure :: SpecificHistoryGetToGaussArray
      procedure :: SpecificHistoryFinalRefinementUpdate

      procedure :: GetBufferSize
      procedure :: ToBuffer
      procedure :: InitializeFromBuffer

      !For u-p elements
      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetSecantCNorm
      procedure :: SetPressureForComputeHistoryAndConstitutiveTensor


   end type

contains

   subroutine LE_SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(LinearElastic), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData

      type(LinearElastic_ED), pointer :: LinEL_ElementData

      allocate(LinEL_ElementData)
      LinEL_ElementData%Material => a

      ElementData => LinEL_ElementData
      
      allocate(LinEl_ElementData%S(a%CT%VoigtSize,pgaus))
      
   end subroutine


   subroutine LE_SpecificReadData(a,Listener)
      class(LinearElastic) :: a
      type(ListenFile) :: Listener
   end subroutine

   subroutine LE_SpecificScatterData(a)
      class(LinearElastic) :: a
   end subroutine

   subroutine NULLSUB(a)
      class(LinearElastic):: a
   end subroutine

    subroutine LE_IsElementSizeRequiredByEMDs(a,isRequired)
      implicit none
      class(LinearElastic), target :: a
      logical :: IsRequired

      isRequired = .false.
   end subroutine


   !-------------------------------------------------------------------
   !EMDS

   subroutine LE_ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      class(LinearElastic_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)

   end subroutine
  

   subroutine LE_ComputeStress(a,igaus,gradDisp)
      class(LinearElastic_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)

      real(rp) :: strain(a%Material%CT%VoigtSize)

      call a%Material%CT%GetStrain(gradDisp,strain)
      a%S(1:a%Material%CT%VoigtSize,igaus) = matmul(a%Material%C,strain(1:a%Material%CT%VoigtSize))
   end subroutine

   subroutine LE_GetStressTensorPointer(a,igaus,S)
      class(LinearElastic_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: S(:)

      S => a%S(1:a%Material%CT%VoigtSize,igaus)
   end subroutine

   subroutine LE_GetConstitutiveTensorPointer(a,igaus,C)
      class(LinearElastic_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)

      C => a%Material%C
   end subroutine

   subroutine GetMaterialPointer(a,Material)
      class(LinearElastic_ED), target :: a
      class(PLCDMaterial), pointer :: Material

      Material => a%Material

   end subroutine

   !Adaptivity
   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(LinearElastic_ED) :: a
      integer(ip) :: ndofn

      ndofn = 0
   end subroutine

   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(LinearElastic_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(LinearElastic_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine CopyTo(a,EMD)
      class(LinearElastic_ED) :: a
      class(ElementMaterialData), pointer :: EMD
   end subroutine

   subroutine Finalize(a)
      class(LinearElastic_ED) :: a
      
      deallocate(a%S)

   end subroutine

   subroutine GetBufferSize(a,BufferSize)
      class(LinearElastic_ED) :: a
      integer(ip) :: BufferSize

      BufferSize = 0
   end subroutine

   subroutine ToBuffer(a,Buffer)
      class(LinearElastic_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   subroutine InitializeFromBuffer(a,Buffer)
      class(LinearElastic_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(LinearElastic_ED) :: a
      integer(ip) :: igaus
      real(rp) :: invK

      call a%Material%CT%GetInverseVolumetricDeformationModulus(invK)
   end subroutine

   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(LinearElastic_ED) :: a
      real(rp) :: SecantCNorm
      integer(ip) :: igaus

      call a%Material%CT%GetCNorm(SecantCNorm)

   end subroutine

   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(LinearElastic_ED) :: a
      real(rp) :: Pressure

      !The pressure is not needed for elastic materials

   end subroutine

    subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(LinearElastic_ED) :: a



   end subroutine

end module

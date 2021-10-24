module Mod_plcd_SIMPMaterial
   use typre
   use Mod_plcd_Material
   use Mod_plcd_LinearElastic
   implicit none
   private
   public  SIMPMaterial, SIMPED

   type, extends(LinearElastic) :: SIMPMaterial
      real(rp) :: p

      !For Postprocess
      type(r1p), allocatable :: ChiPostprocess(:)

contains
      procedure :: SpecificCreateElementData => SIMPSpecificCreateElementData
      procedure :: SetPValue

      procedure :: SetupPostprocess    => SIMP_SetupPostprocess
      procedure :: DoPostprocess       => SIMP_DoPostprocess
   end type

   type :: SIMPGPD
      real(rp) :: chi_p = 1.0_rp
   end type


   type, extends(LinearElastic_ED) :: SIMPED
      real(rp), allocatable :: SIMPC(:,:)

      type(SIMPGPD), allocatable :: GP(:)
      logical :: JustCreated = .true.

      class(SIMPMaterial), pointer :: MaterialS


contains
      procedure :: ComputeHistoryAndConstitutiveTensor => SIMPComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress  => SIMPComputeStress
      procedure :: GetConstitutiveTensorPointer => SIMPGetConstitutiveTensorPointer

      procedure :: SetChiValue
      procedure :: GetChiValue

      !Adaptivity
      procedure :: CopyTo
      procedure :: SpecificHistoryAddFromGaussArray
      procedure :: SpecificHistoryFinalRefinementUpdate
      procedure :: SpecificHistoryGetToGaussArray
      procedure :: GetHistoryVariablesNdofn
      procedure :: GetComparisonCriteriaRatio

      procedure :: GetBufferSize
      procedure :: ToBuffer
      procedure :: InitializeFromBuffer

      !Postprocess
      !Postprocess
      procedure :: ContributeToPostprocess => SIMP_ContributeToPostprocess

   end type


contains

   subroutine SIMPSpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(SIMPMaterial), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData

      class(SIMPED), pointer :: SIMPElementData
      type(SIMPGPD), pointer :: GP(:)
      integer(ip) :: igaus
      
      REAL(8), PARAMETER :: D_QNAN = &
TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)

      

      allocate(SIMPElementData)
      
      allocate(SimpELementData%S(a%CT%VoigtSize,pgaus))
      
      
      SIMPElementData%Material => a
      SIMPElementData%MaterialS => a
      allocate(SIMPElementData%GP(pgaus))

      GP => SIMPElementData%GP

      allocate(SIMPElementData%SIMPC(a%CT%VoigtSize,a%CT%VoigtSize))
      SIMPElementData%SIMPC = SIMPElementData%Material%C

!       do igaus = 1,pgaus
!          GP(igaus)%chi_p = d_qnan
!       enddo

      ElementData => SIMPElementData
   end subroutine

   subroutine SetPValue(a,p)
      class(SIMPMaterial) :: a
      real(rp) :: p
      a%p = p
   end subroutine

   subroutine SIMPComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      class(SIMPED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)

     a%SIMPC = a%Material%C*a%GP(igaus)%chi_p
   end subroutine

    subroutine SIMPComputeStress(a,igaus,gradDisp)
      class(SIMPED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)

      real(rp) :: strain(6)

      call a%Material%CT%GetStrain(gradDisp,strain)
      a%S(1:a%Material%CT%VoigtSize,igaus) = matmul(a%SIMPC,strain(1:a%Material%CT%VoigtSize))
   end subroutine

   subroutine SIMPGetConstitutiveTensorPointer(a,igaus,C)
      class(SIMPED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)

      C => a%SIMPC
   end subroutine


   subroutine GetChiValue(a,igaus,chi)
      class(SIMPED), target :: a
      integer(ip) :: igaus
      real(rp) :: chi


      chi = a%GP(igaus)%chi_p**(1.0_rp/a%MaterialS%p)
   end subroutine

   subroutine SetChiValue(a,igaus,chi)
      class(SIMPED), target :: a
      integer(ip) :: igaus
      real(rp) :: chi

      a%GP(igaus)%chi_p = chi**a%MaterialS%p
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !For adaptivity
     subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(SIMPED) :: a
      integer(ip) :: ndofn

      !In this subroutine the dimension of the history array to be communicated
      !in the adaptive process needs to be set
      !This subroutine is to be overwritten by each material

      ndofn = 1
   end subroutine

   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(SIMPED) :: a
      real(rp) :: GaussArray(:,:)

      integer(ip) :: igaus

      do igaus = 1,size(a%GP,1)
         GaussArray(1,igaus) = a%GP(igaus)%chi_p
      enddo


   end subroutine

   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(SIMPED) :: a
      real(rp) :: GaussArray(:,:)

      integer(ip) :: igaus

      !If this is being added from multiple elements, we need to set it to zero the first time
      if (a%JustCreated) then
         a%JustCreated = .false.
         a%GP(:)%chi_p = 0.0_rp
      endif

      do igaus = 1,size(a%GP,1)
         !Reference Tau
         a%GP(igaus)%chi_p = a%GP(igaus)%chi_p + GaussArray(1,igaus)
      enddo
   end subroutine

   subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(SIMPED) :: a

      integer(ip) :: igaus

   end subroutine

   subroutine Finalize(a)
      class(SIMPED) :: a
      
      deallocate(a%S)
      
      deallocate(a%GP)
      deallocate(a%SIMPC)
   end subroutine

   subroutine CopyTo(a,EMD)
      use Mod_plcd_Material
      class(SIMPED) :: a
      class(ElementMaterialData), pointer :: EMD

      integer(ip) :: igaus

      select type(EMD)
      type is (SIMPED)
         do igaus = 1,size(a%GP,1)
            EMD%GP(igaus)%chi_p = a%GP(igaus)%chi_p
         enddo
      end select

   end subroutine

   subroutine GetBufferSize(a,BufferSize)
      class(SIMPED) :: a
      integer(ip) :: BufferSize

      BufferSize = size(a%GP,1)
   end subroutine

   subroutine ToBuffer(a,Buffer)
      class(SIMPED) :: a
      real(rp) :: Buffer(*)

      integer(ip) :: igaus

      do igaus = 1,size(a%GP,1)
         Buffer(igaus) = a%GP(igaus)%chi_p
      enddo
   end subroutine

   subroutine InitializeFromBuffer(a,Buffer)
      class(SIMPED) :: a
      real(rp) :: Buffer(*)

      integer(ip) :: igaus

      do igaus = 1,size(a%GP,1)
         a%GP(igaus)%chi_p = Buffer(igaus)
      enddo
   end subroutine

   subroutine GetComparisonCriteriaRatio(a,ComparisonRatio)
      class(SIMPED) :: a
      real(rp) :: ComparisonRatio

      ComparisonRatio = maxval(a%GP(:)%chi_p**(1.0_rp/a%MaterialS%p))
   end subroutine


   !Postprocess
   subroutine SIMP_SetupPostprocess(a,Mesh,Memor)
      use Mod_Mesh
      use Mod_Memor
      use Mod_RpMeshAllocator
      implicit none
      class(SIMPMaterial) :: a
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor

      !Allocate array for gauss point info
      call AllocR1p(Memor,Mesh,a%ChiPostprocess)
   end subroutine

   subroutine SIMP_ContributeToPostprocess(a,ielem)
      class(SIMPED) :: a
      integer(ip) :: ielem

      integer(ip) :: igaus

      do igaus = 1,size(a%GP,1)
         a%MaterialS%ChiPostprocess(ielem)%a(igaus) = a%GP(igaus)%chi_p
      enddo
   end subroutine

   subroutine SIMP_DoPostprocess(a,istep,ctime,FilePostpr,Mesh,Memor,iiter_char)
      use Mod_Mesh
      use Mod_Postpr
      use Mod_Memor
      use Mod_RpMeshAllocator
      use Mod_int2str
      implicit none
      class(SIMPMaterial) :: a
      class(PostprFile) :: FilePostpr
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: istep
      real(rp) :: ctime
      character(5) :: iiter_char

      call FilePostpr%postgp(a%ChiPostprocess,'ChiValueGaussPoints'//iiter_char,istep,ctime,Mesh)
      !Deallocate array for GaussPoint info
      call DeallocR1p(Memor,Mesh,a%ChiPostprocess)

   end subroutine

end module



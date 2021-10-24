module Mod_plcd_ConstitutiveTensor
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none
   private
   public ConstitutiveTensor


   type, extends(MPIObject), abstract :: ConstitutiveTensor
      logical :: OnlyDeviatoricPart = .false.

      procedure(), NOPASS, pointer :: GetStrain => NULL()
      procedure(), NOPASS, pointer :: AddSphericalComponent => NULL()
      procedure(), NOPASS, pointer :: GetStrainTensor => NULL()
      procedure(), NOPASS, pointer :: GetStressTensor => NULL()
      procedure(), NOPASS, pointer :: GetStressVector => NULL()
      procedure(), NOPASS, pointer :: GetStrainVector => NULL()
      
      integer(ip) :: VoigtSize

contains
      procedure :: Initialize
      procedure(ReadData), deferred :: ReadData
      procedure(ScatterData), deferred :: ScatterData
      procedure(CheckData), deferred :: CheckData
      procedure(Build), deferred :: Build
      procedure(GetPointerToC), deferred :: GetPointerToC
      procedure(GetPointerToCDeviatoric), deferred :: GetPointerToCDeviatoric
      procedure(GetPointerToCSpherical), deferred :: GetPointerToCSpherical
      procedure(GetInverseVolumetricDeformationModulus), deferred :: GetInverseVolumetricDeformationModulus
      procedure(GetCNorm), deferred :: GetCNorm
      procedure :: SetOnlyDeviatoric
      procedure(ComputeStressNorm), deferred :: ComputeStressNorm
      procedure(ComputeMeanStress), deferred :: ComputeMeanStress
      procedure :: GetYoungModulusForDirection
      procedure :: GetPoissonRatio
      procedure(GetYoungModulusEquiv), deferred :: GetYoungModulusEquiv
      procedure :: GetYoungModulus
      procedure :: SetPoissonRatio
      procedure :: SetYoungModulus
      procedure :: SetYoungModulusForDirection


   end type

   abstract interface

      subroutine ReadData(a,Listener)
         import
         implicit none
         class(ConstitutiveTensor) :: a
         type(ListenFile) :: Listener
      end subroutine

      subroutine ScatterData(a)
         import
         implicit none
         class(ConstitutiveTensor) :: a
      end subroutine

      subroutine CheckData(a)
         import
         implicit none
         class(ConstitutiveTensor) :: a
      end subroutine
      
      subroutine Build(a)
         import
         implicit none
         class(ConstitutiveTensor) :: a
      end subroutine

      subroutine GetInverseVolumetricDeformationModulus(a,invK)
         import
         implicit none
         class(ConstitutiveTensor) :: a
         real(rp) :: invK
      end subroutine

      subroutine GetCNorm(a,Cnorm)
         import
         implicit none
         class(ConstitutiveTensor) :: a
         real(rp) :: CNorm
      end subroutine

      subroutine GetPointerToC(a,C)
         import
         implicit none
         class(ConstitutiveTensor), target :: a
         real(rp), pointer :: C(:,:)
      end subroutine

      subroutine GetPointerToCDeviatoric(a,C)
         import
         implicit none
         class(ConstitutiveTensor), target :: a
         real(rp), pointer :: C(:,:)
      end subroutine

      subroutine GetPointerToCSpherical(a,C)
         import
         implicit none
         class(ConstitutiveTensor), target :: a
         real(rp), pointer :: C(:,:)
      end subroutine

      subroutine ComputeStressNorm(a,Stress,StressNorm)
         import
         implicit none
         class(ConstitutiveTensor), target :: a
         real(rp) :: Stress(:),StressNorm
      end subroutine

      subroutine ComputeMeanStress(a,Stress,MeanStress)
         import
         implicit none
         class(ConstitutiveTensor), target :: a
         real(rp) :: Stress(:),MeanStress
      end subroutine

      subroutine GetYoungModulusEquiv(a,E)
         import
         implicit none
         class(ConstitutiveTensor) :: a
         real(rp) :: E
      end subroutine
   
      
   end interface


contains
   subroutine Initialize(a,ndime,LargeStrainsflag)
      use Mod_plcd_StrainGenerator
      implicit none
      class(ConstitutiveTensor) :: a
      integer(ip) :: ndime
      integer(ip) :: LargeStrainsflag

      call AllocateGetStrain(ndime,a%GetStrain,LargeStrainsflag)
      call GetVoigtsize(ndime,a%VoigtSize)
      call AllocateAddSphericalComponent(ndime,a%AddSphericalComponent)
      call AllocateGetStressTensor(ndime,a%GetStressTensor)
      call AllocateGetStrainTensor(ndime,a%GetStrainTensor)
      call AllocateGetStressVector(ndime,a%GetStressVector)
      call AllocateGetStrainVector(ndime,a%GetStrainVector)

   end subroutine

   subroutine SetOnlyDeviatoric(a,OnlyDeviatoricPart)
      implicit none
      class(ConstitutiveTensor) :: a
      logical :: OnlyDeviatoricPart

      a%OnlyDeviatoricPart = OnlyDeviatoricPart
   end subroutine

   subroutine GetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      call runend('GetYoungModulusForDirection not implemented')
   end subroutine

   subroutine GetPoissonRatio(a,nu)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: nu

      call runend('GetPoissonRatio not implemented')
   end subroutine
   
      subroutine GetYoungModulus(a,E)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: E

      call runend('GetYoungModulus not implemented')
   end subroutine
   
   subroutine SetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      call runend('SetYoungModulusForDirection not implemented')
   end subroutine

   subroutine SetPoissonRatio(a,nu)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: nu

      call runend('SetPoissonRatio not implemented')
   end subroutine
   
      subroutine SetYoungModulus(a,E)
      implicit none
      class(ConstitutiveTensor) :: a
      real(rp) :: E

      call runend('SetYoungModulus not implemented')
   end subroutine

end module

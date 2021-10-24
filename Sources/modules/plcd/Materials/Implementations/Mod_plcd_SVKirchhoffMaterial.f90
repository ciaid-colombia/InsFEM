module Mod_plcd_SVKirchhoffMaterial
   use typre
   use Mod_plcd_Material
   use Mod_Listen
   use Mod_Element
   implicit none
   private
   public SVKirchhoffMaterial, SVKirchhoffMaterial_ED

   type, extends(PLCDMaterial) :: SVKirchhoffMaterial
   
      real(rp) :: lambda, mu

contains
      procedure :: SpecificCreateElementData => SVK_SpecificCreateElementData
      procedure :: SpecificReadData => SVK_SpecificReadData
      procedure :: SpecificScatterData => SVK_SpecificScatterData
      procedure :: SpecificSetup       => NULLSUB
      procedure :: IsElementSizeRequiredByEMDs => SVK_IsElementSizeRequiredByEMDs

   end type
   
   type :: SVKirchhoffMaterial_GPD
      real(rp), allocatable :: C(:,:)
      real(rp), allocatable :: S(:) 
      real(rp), allocatable :: F_mat(:,:), F_spat(:,:)
      real(rp) :: J 
      
   end type

   type, extends(ElementMaterialData) :: SVKirchhoffMaterial_ED

      class(SVKirchhoffMaterial), pointer :: Material => NULL()
      type(SVKirchhoffMaterial_GPD), allocatable :: GP(:)


contains

      procedure :: ComputeHistoryAndConstitutiveTensor => SVK_ComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress => SVK_ComputeStress
      procedure :: GetConstitutiveTensorPointer => SVK_GetConstitutiveTensorPointer
      procedure :: GetStressTensorPointer => SVK_GetStressTensorPointer
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

   subroutine SVK_SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(SVKirchhoffMaterial), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData

      type(SVKirchhoffMaterial_ED), pointer :: SVK_ElementData
      type(SVKirchhoffMaterial_GPD), pointer :: GP(:)
      integer(ip) :: igaus

      allocate(SVK_ElementData)
      SVK_ElementData%Material => a
      
      allocate(SVK_ElementData%GP(pgaus))
      GP => SVK_ElementData%GP
      
      do igaus = 1,pgaus
         allocate(GP(igaus)%C(a%CT%VoigtSize,a%CT%VoigtSize))
         allocate(GP(igaus)%S(a%CT%VoigtSize))
         allocate(GP(igaus)%F_mat(a%ndime,a%ndime))
         allocate(GP(igaus)%F_spat(a%ndime,a%ndime))
         GP(igaus)%J=0.0_rp
      enddo
      
      ElementData => SVK_ElementData
   end subroutine


   subroutine SVK_SpecificReadData(a,Listener)
      class(SVKirchhoffMaterial) :: a
      type(ListenFile) :: Listener
      
   end subroutine

   subroutine SVK_SpecificScatterData(a)
      use MPI
      class(SVKirchhoffMaterial) :: a      

      call ComputeLameParameters(a)
      
   end subroutine

   subroutine NULLSUB(a)
      class(SVKirchhoffMaterial):: a
   end subroutine

    subroutine SVK_IsElementSizeRequiredByEMDs(a,isRequired)
      implicit none
      class(SVKirchhoffMaterial), target :: a
      logical :: IsRequired

      isRequired = .false.
   end subroutine
   
   subroutine ComputeLameParameters(a)
      implicit none
      class(SVKirchhoffMaterial) :: a
      real(rp) :: E,nu
      
      call a%CT%GetYoungModulus(E)
      call a%CT%GetPoissonRatio(nu)
      
      a%mu     = E/(2.0_rp + 2.0_rp*nu)
      a%lambda = nu*E/((1.0_rp + nu)*(1.0_rp - 2.0_rp*nu)) 
      
   end subroutine

   !-------------------------------------------------------------------
   !EMDS

   subroutine SVK_ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      class(SVKirchhoffMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      call ComputeDeformationGradientandJ(a,igaus,gradDisp)
      call SVK_ComputeConstitutiveTensor(a,igaus)

   end subroutine
   
   subroutine SVK_ComputeConstitutiveTensor(a, igaus)
      class(SVKirchhoffMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp), pointer :: J
      
      call GetJPointer(a,igaus,J) 
      a%GP(igaus)%C = (1.0_rp/J)*a%Material%C

   end subroutine

    subroutine SVK_ComputeStress(a,igaus,gradDisp)
      class(SVKirchhoffMaterial_ED), target :: a
      integer(ip) :: igaus, idime
      real(rp) :: gradDisp(:,:)
      
      real(rp), pointer :: F_mat(:,:) => NULL()
      real(rp), pointer :: F_spat(:,:) => NULL()
      
      !real(rp) :: Identity(a%Material%ndime,a%Material%ndime)
      real(rp) :: PiolaKirchhoffStress2(a%Material%ndime,a%Material%ndime)
      !real(rp) :: LeftCauchyTensor(a%Material%ndime,a%Material%ndime)
      real(rp) :: GreenLagrangeStrainTensor(a%Material%ndime,a%Material%ndime)
      real(rp) :: CauchyStress(a%Material%ndime,a%Material%ndime)
      
      real(rp) :: GreenLagrangeStrainVector(a%Material%CT%VoigtSize)
      
      real(rp) :: TraceE
      real(rp), pointer :: lambda, mu, J
      integer(ip), pointer :: ndime

      PiolaKirchhoffStress2 = 0.0_rp
      !LeftCauchyTensor = 0.0_rp
      GreenLagrangeStrainTensor = 0.0_rp
      CauchyStress = 0.0_rp
      GreenLagrangeStrainVector = 0.0_rp
      
      call GetDeformationGradientPointer(a,igaus,F_mat,F_spat)
      call GetJPointer(a,igaus,J)
      call GetLameParametersPointer(a,lambda,mu)
      
      ndime => a%Material%ndime

      !we compute Green strain tensor
      !GreenLagrangeStrainTensor= 0.5_rp*(matmul(transpose(F_mat),F_mat) - Identity)
      GreenLagrangeStrainTensor= 0.5_rp*(matmul(transpose(F_mat),F_mat)) 
      forall(idime = 1:ndime) GreenLagrangeStrainTensor(idime,idime) =  GreenLagrangeStrainTensor(idime,idime) - 0.5_rp

      !calculate PK2
      TraceE = 0.0_rp
      
      do idime = 1, ndime
         TraceE = TraceE + GreenLagrangeStrainTensor(idime,idime)
      end do
      
      !PiolaKirchhoffStress2 = lambda*TraceE*Identity + mu*2.0_rp*(GreenLagrangeStrainTensor)
      PiolaKirchhoffStress2 = mu*2.0_rp*(GreenLagrangeStrainTensor)
      forall(idime = 1:ndime) PiolaKirchhoffStress2(idime,idime) = PiolaKirchhoffStress2(idime,idime) + lambda*TraceE
      
      
      !Push forward of piolaK2 = kirchoff stress = J*cauchy
      !we solve for cauchy
      
      CauchyStress = (1.0_rp/J)*(matmul(F_mat,matmul(PiolaKirchhoffStress2,transpose(F_mat))))

      !get into voigt notation
      call a%Material%CT%GetStressVector(CauchyStress,a%GP(igaus)%S(1:a%Material%CT%VoigtSize))
      
   end subroutine

   subroutine SVK_GetStressTensorPointer(a,igaus,S)
      class(SVKirchhoffMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: S(:)

      S => a%GP(igaus)%S(1:a%Material%CT%VoigtSize)
   end subroutine

   subroutine SVK_GetConstitutiveTensorPointer(a,igaus,C)
      class(SVKirchhoffMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)

      C => a%GP(igaus)%C
   end subroutine

   subroutine GetMaterialPointer(a,Material)
      class(SVKirchhoffMaterial_ED), target :: a
      class(PLCDMaterial), pointer :: Material

      Material => a%Material

   end subroutine
   
   subroutine ComputeDeformationGradientandJ(a,igaus,gradDisp)
      class(SVKirchhoffMaterial_ED), target :: a
      integer(ip) :: idime, igaus
      integer(ip), pointer :: ndime
      real(rp) :: gradDisp(:,:)
      type(SVKirchhoffMaterial_GPD), pointer :: GP => NULL()
      
      ndime => a%Material%ndime
      GP => a%GP(igaus)

      !we transform gradDisp into deformation tensor F
      !u = x_sp - X_mat; grad(u) = ii - F(x)^-1
      GP%F_spat= - gradDisp
      forall(idime = 1:ndime) GP%F_spat(idime,idime) = GP%F_spat(idime,idime) + 1.0_rp
      
      !Compute F_mat as the inverse of F_spat and J as determinant
      call invmtx(GP%F_spat,GP%F_mat,GP%J,ndime)

      GP%J= 1.0_rp/GP%J
      
      if(GP%J .lt. epsilon(0.0_rp) ) then
         call runend('PLCD Base Elmope: Negative element determinant')
      end if
      
   end subroutine
   
   subroutine GetJPointer(a,igaus,J)
      class(SVKirchhoffMaterial_ED), target:: a
      real(rp), pointer :: J
      integer(ip) :: igaus

      J => a%GP(igaus)%J
   end subroutine
   
   subroutine GetDeformationGradientPointer(a,igaus,F_mat,F_spat)
      class(SVKirchhoffMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: F_mat(:,:), F_spat(:,:)

      F_mat => a%GP(igaus)%F_mat
      F_spat => a%GP(igaus)%F_spat
   end subroutine
   
   subroutine GetLameParametersPointer(a,lambda,mu)
      class(SVKirchhoffMaterial_ED), target:: a
      real(rp), pointer :: lambda, mu

      lambda => a%Material%lambda
      mu => a%Material%mu
   end subroutine
   
   !Adaptivity
   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(SVKirchhoffMaterial_ED) :: a
      integer(ip) :: ndofn

      ndofn = 0
   end subroutine

   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine CopyTo(a,EMD)
      class(SVKirchhoffMaterial_ED) :: a
      class(ElementMaterialData), pointer :: EMD

   end subroutine

   subroutine Finalize(a)
      class(SVKirchhoffMaterial_ED) :: a
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         deallocate(a%GP(igaus)%S)
         deallocate(a%GP(igaus)%C)
         deallocate(a%GP(igaus)%F_mat)
         deallocate(a%GP(igaus)%F_spat)
      enddo
      deallocate(a%GP)

   end subroutine

   subroutine GetBufferSize(a,BufferSize)
      class(SVKirchhoffMaterial_ED) :: a
      integer(ip) :: BufferSize

      BufferSize = 0
   end subroutine

   subroutine ToBuffer(a,Buffer)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   subroutine InitializeFromBuffer(a,Buffer)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(SVKirchhoffMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: invK

      call a%Material%CT%GetInverseVolumetricDeformationModulus(invK)
   end subroutine

   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: SecantCNorm
      integer(ip) :: igaus

      call a%Material%CT%GetCNorm(SecantCNorm)

   end subroutine

   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(SVKirchhoffMaterial_ED) :: a
      real(rp) :: Pressure

      !The pressure is not needed for elastic materials

   end subroutine

    subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(SVKirchhoffMaterial_ED) :: a



   end subroutine

end module

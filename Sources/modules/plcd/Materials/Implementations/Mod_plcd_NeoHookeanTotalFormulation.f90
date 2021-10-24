module Mod_plcd_NeoHookeanTotalFormulation
   use typre
   use Mod_plcd_Material
   use Mod_Listen
   use Mod_Element
   use Mod_plcd_ConstitutiveTensor
   implicit none
   private
   public NeoHookeanTotalMaterial, NeoHookeanTotalMaterial_ED

   type, extends(PLCDMaterial) :: NeoHookeanTotalMaterial
   
      real(rp) :: lambda, mu
      integer(ip), allocatable :: VoigtTransform(:,:)

contains
      procedure :: SpecificCreateElementData => NHT_SpecificCreateElementData
      procedure :: SpecificReadData => NHT_SpecificReadData
      procedure :: SpecificScatterData => NHT_SpecificScatterData
      procedure :: SpecificSetup       => NHT_SpecificSetup
      procedure :: IsElementSizeRequiredByEMDs => NHT_IsElementSizeRequiredByEMDs

   end type
   
   type :: NeoHookeanTotalMaterial_GPD
      real(rp), allocatable :: C(:,:)
      real(rp), allocatable :: S(:) 
      real(rp), allocatable :: F_mat(:,:)
      real(rp) :: J 
      
   end type

   type, extends(ElementMaterialData) :: NeoHookeanTotalMaterial_ED

      class(NeoHookeanTotalMaterial), pointer :: Material => NULL()
      type(NeoHookeanTotalMaterial_GPD), allocatable :: GP(:)

contains

      procedure :: ComputeHistoryAndConstitutiveTensor => NHT_ComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress => NHT_ComputeStress
      procedure :: GetConstitutiveTensorPointer => NHT_GetConstitutiveTensorPointer
      procedure :: GetStressTensorPointer => NHT_GetStressTensorPointer
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

   subroutine NHT_SpecificSetup(a)
      implicit none
      class(NeoHookeanTotalMaterial) :: a

   end subroutine   

   subroutine NHT_SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(NeoHookeanTotalMaterial), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData

      type(NeoHookeanTotalMaterial_ED), pointer :: NHT_ElementData
      type(NeoHookeanTotalMaterial_GPD), pointer :: GP(:)
      integer(ip) :: igaus

      allocate(NHT_ElementData)
      NHT_ElementData%Material => a
      
      allocate(NHT_ElementData%GP(pgaus))
      GP => NHT_ElementData%GP
      
      do igaus = 1,pgaus
         allocate(GP(igaus)%C(a%CT%VoigtSize,a%CT%VoigtSize))
         allocate(GP(igaus)%S(a%CT%VoigtSize))
         allocate(GP(igaus)%F_mat(a%ndime,a%ndime))
         GP(igaus)%J=0.0_rp
      enddo

      ElementData => NHT_ElementData
      
   end subroutine


   subroutine NHT_SpecificReadData(a,Listener)
      class(NeoHookeanTotalMaterial) :: a
      type(ListenFile) :: Listener
      
   end subroutine

   subroutine NHT_SpecificScatterData(a)
      use MPI
      class(NeoHookeanTotalMaterial) :: a
      
      allocate(a%VoigtTransform(3,a%CT%VoigtSize))
      call ConstructVoigtTransformation(a,a%VoigtTransform)
      call ComputeLameParameters(a)
      
   end subroutine

   subroutine NULLSUB(a)
      class(NeoHookeanTotalMaterial):: a
   end subroutine

   subroutine NHT_IsElementSizeRequiredByEMDs(a,isRequired)
      implicit none
      class(NeoHookeanTotalMaterial), target :: a
      logical :: IsRequired

      isRequired = .false.
   end subroutine
   
   subroutine ComputeLameParameters(a)
      implicit none
      class(NeoHookeanTotalMaterial) :: a
      real(rp) :: E,nu
      
      call a%CT%GetYoungModulus(E)
      call a%CT%GetPoissonRatio(nu)
      
      a%mu     = E/(2.0_rp + 2.0_rp*nu)
      a%lambda = nu*E/((1.0_rp + nu)*(1.0_rp - 2.0_rp*nu)) 
      
   end subroutine


   !-------------------------------------------------------------------
   !EMDS

   subroutine NHT_ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      class(NeoHookeanTotalMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      call ComputeDeformationGradientandJ(a,igaus,gradDisp)
      call NHT_ComputeConstitutiveTensor(a,igaus)
   end subroutine

   subroutine NHT_ComputeConstitutiveTensor(a,igaus)
      class(NeoHookeanTotalMaterial_ED), target :: a
      integer(ip) :: igaus, ivoigt, idime, jdime, ldime, kdime, jvoigt
      real(rp), pointer :: J => NULL()
      real(rp), pointer :: F_mat(:,:) => NULL()
      real(rp) :: lambda, mu
      real(rp) :: C(a%Material%ndime,a%Material%ndime,a%Material%ndime,a%Material%ndime)
      integer(ip), pointer :: ndime
      real(rp) :: RightCauchyTensor(a%Material%ndime,a%Material%ndime)
      
      ndime => a%Material%ndime
      
      call GetJPointer(a,igaus,J)
      call GetDeformationGradientPointer(a,igaus,F_mat)
      
      lambda = a%Material%lambda
      mu = a%Material%mu - a%Material%lambda*log(J)
      
      RightCauchyTensor = matmul(transpose(F_mat),F_mat)
      
      call invert(RightCauchyTensor,ndime,ndime)
      
      forall(idime=1:ndime,jdime=1:ndime,kdime=1:ndime,ldime=1:ndime)
         C(idime,jdime,kdime,ldime)=lambda*RightCauchyTensor(idime,jdime)*RightCauchyTensor(kdime,ldime)+ mu*(RightCauchyTensor(idime,kdime)*RightCauchyTensor(jdime,ldime)+RightCauchyTensor(idime,ldime)*RightCauchyTensor(jdime,kdime))
      end forall
      
      do ivoigt=1,a%Material%CT%VoigtSize
         do jvoigt=1,a%Material%CT%VoigtSize
            a%GP(igaus)%C(ivoigt,jvoigt)=C(a%Material%VoigtTransform(2,ivoigt),a%Material%VoigtTransform(3,ivoigt),a%Material%VoigtTransform(2,jvoigt),a%Material%VoigtTransform(3,jvoigt))
         enddo
      enddo

   end subroutine
   
   subroutine ComputeDeformationGradientandJ(a,igaus,gradDisp)
      class(NeoHookeanTotalMaterial_ED), target :: a
      integer(ip) :: idime, igaus
      integer(ip), pointer :: ndime
      real(rp) :: gradDisp(:,:)
      
      type(NeoHookeanTotalMaterial_GPD), pointer :: GP => NULL()
      
      ndime => a%Material%ndime
         
      GP => a%GP(igaus)

      !we transform gradDisp into deformation tensor F
     
      GP%F_mat = gradDisp
      forall (idime = 1:ndime) GP%F_mat(idime,idime) = GP%F_mat(idime,idime) + 1.0_rp 
      
      call deter(GP%F_mat,GP%J,ndime)
      
      if(GP%J .lt. epsilon(0.0_rp) ) then
         call runend('PLCD Base Elmope: Negative element determinant')
      end if
      
   end subroutine
   
   subroutine NHT_ComputeStress(a,igaus,gradDisp)
      class(NeoHookeanTotalMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)

      real(rp) :: GreenLagrangeStrain(a%Material%CT%VoigtSize)

      call a%Material%CT%GetStrain(gradDisp,GreenLagrangeStrain)
      a%GP(igaus)%S(1:a%Material%CT%VoigtSize) = matmul(a%GP(igaus)%C,GreenLagrangeStrain(1:a%Material%CT%VoigtSize))

   end subroutine

   subroutine NHT_GetStressTensorPointer(a,igaus,S)
      class(NeoHookeanTotalMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: S(:)

      S => a%GP(igaus)%S(1:a%Material%CT%VoigtSize)
   end subroutine

   subroutine NHT_GetConstitutiveTensorPointer(a,igaus,C)
      class(NeoHookeanTotalMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)

      C => a%GP(igaus)%C
   end subroutine

   subroutine GetMaterialPointer(a,Material)
      class(NeoHookeanTotalMaterial_ED), target :: a
      class(PLCDMaterial), pointer :: Material

      Material => a%Material
   end subroutine
   
   subroutine GetJPointer(a,igaus,J)
      class(NeoHookeanTotalMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: J

      J => a%GP(igaus)%J
   end subroutine
   
   subroutine GetDeformationGradientPointer(a,igaus,F_mat)
      class(NeoHookeanTotalMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: F_mat(:,:)

      F_mat => a%GP(igaus)%F_mat
   end subroutine
   
   subroutine GetLameParametersPointer(a,lambda,mu)
      class(NeoHookeanTotalMaterial_ED), target:: a
      real(rp), pointer :: lambda, mu

      lambda => a%Material%lambda
      mu => a%Material%mu
   end subroutine
   

   !Adaptivity
   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(NeoHookeanTotalMaterial_ED) :: a
      integer(ip) :: ndofn

      ndofn = 0
   end subroutine

   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)


   end subroutine

   subroutine CopyTo(a,EMD)
      class(NeoHookeanTotalMaterial_ED) :: a
      class(ElementMaterialData), pointer :: EMD

   end subroutine

   subroutine Finalize(a)
      class(NeoHookeanTotalMaterial_ED) :: a
      integer(ip) :: igaus

      do igaus = 1,size(a%GP,1)
         deallocate(a%GP(igaus)%S)
         deallocate(a%GP(igaus)%C)
         deallocate(a%GP(igaus)%F_mat)
      enddo
      deallocate(a%GP)
      
   end subroutine

   subroutine GetBufferSize(a,BufferSize)
      class(NeoHookeanTotalMaterial_ED) :: a
      integer(ip) :: BufferSize

      BufferSize = 0
   end subroutine

   subroutine ToBuffer(a,Buffer)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   subroutine InitializeFromBuffer(a,Buffer)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: Buffer(*)
   end subroutine

   !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(NeoHookeanTotalMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: invK

      call a%Material%CT%GetInverseVolumetricDeformationModulus(invK)
   end subroutine

   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: SecantCNorm
      integer(ip) :: igaus

      call a%Material%CT%GetCNorm(SecantCNorm)

   end subroutine

   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(NeoHookeanTotalMaterial_ED) :: a
      real(rp) :: Pressure

      !The pressure is not needed for elastic materials

   end subroutine

   subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(NeoHookeanTotalMaterial_ED) :: a

   end subroutine
   
   subroutine ConstructVoigtTransformation(a,VoigtTransform)
      class(NeoHookeanTotalMaterial) :: a
      integer(ip) :: VoigtTransform(3,a%CT%VoigtSize), idime
      
      forall(idime=1:a%ndime)
         VoigtTransform(:,idime)=idime
      endforall
      
      VoigtTransform(:,a%ndime+1)= (/a%ndime+1,1,2/)
      if (a%ndime == 3) then
         VoigtTransform(:,a%ndime+2)= (/a%ndime+2,1,3/)
         VoigtTransform(:,a%ndime+3)= (/a%ndime+3,2,3/)
      endif
      
   end subroutine

end module
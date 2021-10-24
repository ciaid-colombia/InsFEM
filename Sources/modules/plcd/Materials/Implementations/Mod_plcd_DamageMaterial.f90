module Mod_plcd_DamageMaterial 
   use typre
   use Mod_plcd_Material
   use Mod_Listen
   use Mod_Element
   use Mod_plcd_ComparisonCriteria
   use Mod_plcd_ComparisonCriteriaFactory
   implicit none
   private
   public DamageMaterial

   
   type, extends(PLCDMaterial) :: DamageMaterial
      
      real(rp) :: DamageLimit = 0.9999
      real(rp) :: Tau0, FractureEnergy,  UniaxialStressLimit
      character(5) :: CompareCriteria = 'SIMOJ'
      logical :: DeviatoricOnlyDamage = .false.
      
      type(r1p), allocatable :: DPostprocess(:),TauPostprocess(:),CompPostprocess(:)
      
      
      logical :: isReadComparisonCriteria = .false.
      class(ComparisonCriteria), pointer :: ComparisonCriteria
      procedure(ComputeDamage), pointer :: ComputeDamage => NULL()
      procedure(ComputeCT), pointer :: ComputeCT => NULL()

   
contains
      procedure :: SpecificCreateElementData => MD_SpecificCreateElementData
      procedure :: SpecificReadData => MD_SpecificReadData
      procedure :: SpecificScatterData => MD_SpecificScatterData
      procedure :: SpecificSetup       => MD_SpecificSetup
      
      !Postprocess
      procedure :: SetupPostprocess    => MD_SetupPostprocess
      procedure :: DoPostprocess       => MD_DoPostprocess
      procedure :: IsElementSizeRequiredByEMDs => MD_IsElementSizeRequiredByEMDs
      
   end type
   
   type :: DamageMaterial_GPD
      logical :: IsDamaged = .false.
      real(rp) :: ReferenceTau, CurrentTau
      real(rp) :: ReferenceD = 0.0_rp, CurrentD = 0.0_rp
      real(rp), pointer :: C(:,:) => NULL()
      real(rp), pointer :: S(:) => NULL()
      real(rp) :: ComparisonRatio = 0.0_rp
   end type

   
   type, extends(ElementMaterialData) :: DamageMaterial_ED
      
      logical :: JustCreated = .false. !Required for adaptivity
      class(DamageMaterial), pointer :: Material => NULL()
      
      
      real(rp) :: ElementSize = 1.0_rp   
      type(DamageMaterial_GPD), allocatable :: GP(:)
      
      !For UP elements
      real(rp) :: GpPressure
         
contains   

      procedure :: ComputeHistoryAndConstitutiveTensor => MD_ComputeHistoryAndConstitutiveTensor
      procedure :: ComputeStress => MD_ComputeStress
      procedure :: GetConstitutiveTensorPointer => MD_GetConstitutiveTensorPointer
      procedure :: GetStressTensorPointer => MD_GetStressTensorPointer
      procedure :: SetElementSize => MD_SetElementSize
      procedure :: MoveHistoryVariablesToConverged => MD_MoveHistoryVariablesToConverged
      
      !Postprocess
      procedure :: ContributeToPostprocess => MD_ContributeToPostprocess
      procedure :: GetMaterialPointer
      procedure :: Finalize
      
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
      
      !For u-p elements
      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetSecantCNorm
      procedure :: SetPressureForComputeHistoryAndConstitutiveTensor
      
   end type

   interface
      Subroutine ComputeDamage(a,ElementSize,GP)
         import
         implicit none
         class(DamageMaterial) :: a
         real(rp), intent(in) :: ElementSize   
         class(DamageMaterial_GPD) :: GP      
      end subroutine
      
      Subroutine ComputeCT(a,ElementSize,ElasticPredictor,GP)
         import
         implicit none
         class(DamageMaterial) :: a
         real(rp), intent(in) :: ElementSize
         real(rp), Intent(in) :: ElasticPredictor(:) 
         class(DamageMaterial_GPD) :: GP
      end subroutine

   end interface

contains

   subroutine MD_SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(DamageMaterial), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData
      
      type(DamageMaterial_ED), pointer :: MD_ElementData
      type(DamageMaterial_GPD), pointer :: GP(:)
      integer(ip) :: igaus
      
      allocate(MD_ElementData)
      MD_ElementData%JustCreated = .true.
      MD_ElementData%Material => a
      allocate(MD_ElementData%GP(pgaus))
      
      GP => MD_ElementData%GP
      
      do igaus = 1,pgaus
         GP(igaus)%ReferenceTau = a%Tau0
         GP(igaus)%CurrentTau = a%Tau0 
         GP(igaus)%C => a%C
         allocate(GP(igaus)%S(a%CT%VoigtSize))
      enddo
      
      ElementData => MD_ElementData
   end subroutine
   
   subroutine MD_SpecificReadData(a,Listener)
      class(DamageMaterial) :: a
      type(ListenFile) :: Listener
      
      if (Listener%words(1) == 'UNIAX') then
         a%UniaxialStressLimit = Listener%param(1)
      elseif (Listener%words(1) == 'FRACT') then
         a%FractureEnergy = Listener%param(1)
      elseif (Listener%words(1) == 'MAXDA' ) then
         a%DamageLimit = Listener%param(1)
      elseif (Listener%words(1) == 'COMPA' ) then
         !ComparisonCriteria
         a%CompareCriteria = Listener%words(2)   
         call ComparisonCriteriaFactory(a%CompareCriteria,a%ComparisonCriteria)
         call a%ComparisonCriteria%ReadData(Listener)
         call AllocateDamageAndTypeCT(a)
         a%isReadComparisonCriteria = .true.
      elseif (Listener%words(1) == 'TYPEO') then
         if(Listener%exists('DEVIA')) then
            a%DeviatoricOnlyDamage = .true.
         elseif(Listener%exists('FULL ')) then
            a%DeviatoricOnlyDamage = .false.
         endif
      endif   

      
   end subroutine
   
   subroutine MD_SpecificScatterData(a)
      use MPI
      class(DamageMaterial) :: a
      
      integer(ip) :: ierr
      
      call MPI_BCAST(a%UniaxialStressLimit,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%FractureEnergy,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%DamageLimit,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%CompareCriteria,5,MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%DeviatoricOnlyDamage,1,MPI_LOGICAL,a%MPIroot,a%MPIcomm,ierr)
      
      !Create the comparison Criteria
      if (a%MPIrank /= a%MPIroot .or. (a%IsReadComparisonCriteria .eqv. .false.)) then
         call ComparisonCriteriaFactory(a%CompareCriteria,a%ComparisonCriteria)
         call AllocateDamageAndTypeCT(a)
      endif
      call a%ComparisonCriteria%ScatterData
      
      !call a%ComparisonCriteria%ComputeTauValueFromUniaxialStressTest(a%UniaxialStressLimit,a%CT,a%Tau0)

   end subroutine
   
   subroutine MD_SpecificSetup(a)
      class(DamageMaterial) :: a
      
      call a%ComparisonCriteria%ComputeTauValueFromUniaxialStressTest(a%UniaxialStressLimit,a%CT,a%Tau0)

   end subroutine
         
   subroutine AllocateDamageAndTypeCT(a)
      class(DamageMaterial) :: a
   
      select case(a%TangentTensorType)
         case(0)
            a%ComputeCT => ComputeCT_Elastic
            
         case(1)
            a%ComputeCT => ComputeCT_Secant
            
         case(2)
            select case(a%CompareCriteria)
               case('SIMOJ')
                  a%ComputeCT => ComputeCT_Secant
            
               case('VONMI')
                  a%ComputeCT => ComputeCT_Secant !ComputeCT_Tangent_VonMises_Exp
               
            end select
         
         case default
            a%ComputeCT => ComputeCT_Secant
         
      end select
      
      select case(a%CompareCriteria)
         case('SIMOJ')
            a%ComputeDamage => ComputeDamage_Dev_Exp !ComputeDamage_Exp !ComputeDamage_SimoJu_Exp
            
         case('VONMI')
            a%ComputeDamage => ComputeDamage_Dev_Exp !ComputeDamage_Exp !ComputeDamage_VonMises_Exp
            
         case default
!            a%ComputeDamage => ComputeDamage_SimoJu_Exp
               
      end select
      
   end subroutine   
   
   subroutine NULLSUB(a)
      class(DamageMaterial):: a
   end subroutine
     
   subroutine MD_ComputeHistoryAndConstitutiveTensor(a,igaus,gradDisp)
      use Mod_plcd_ComparisonCriteria
      
      class(DamageMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      real(rp) :: strain(a%Material%CT%VoigtSize)
      real(rp) :: ElasticPredictor(a%Material%CT%VoigtSize)
      type(DamageMaterial_GPD), pointer :: GP => NULL()
      type(DamageMaterial), pointer :: Mat => NULL()
      
      GP => a%GP(igaus)
      Mat => a%Material
      
      call Mat%CT%GetStrain(gradDisp,strain)
      
      !Elastic predictor
      ElasticPredictor = matmul(Mat%C,strain)
      
      if (Mat%UseUPFormulation) call Mat%CT%AddSphericalComponent(-a%gpPressure,ElasticPredictor)
      
      !Not Converged Tau
      call Mat%ComparisonCriteria%Evaluate(ElasticPredictor,strain,Mat%CT,GP%CurrentTau)
      GP%ComparisonRatio = GP%CurrentTau/Mat%Tau0
      
      !Update damage parameter if necessary
      if (GP%CurrentTau > GP%ReferenceTau) then

         if (GP%isDamaged .eqv. .false.) then
            allocate(GP%C(Mat%CT%VoigtSize,Mat%CT%VoigtSize))
            GP%isDamaged = .true.
         endif
                    
         call Mat%ComputeDamage(a%ElementSize,GP)
         call Mat%ComputeCT(a%ElementSize,ElasticPredictor,GP)
         
      else
         GP%CurrentTau = GP%ReferenceTau
         GP%CurrentD = GP%ReferenceD
      endif
      
  end subroutine
  
   subroutine MD_ComputeStress(a,igaus,gradDisp)
      class(DamageMaterial_ED), target :: a
      integer(ip) :: igaus
      real(rp) :: gradDisp(:,:)
      
      real(rp) :: strain(a%Material%CT%VoigtSize)
      type(DamageMaterial_GPD), pointer :: GP => NULL()
      real(rp) :: ElasticPredictor(a%Material%CT%VoigtSize)
      
      real(rp), pointer :: CDev(:,:),CSph(:,:)
      
      call a%Material%CT%GetStrain(gradDisp,strain)
      
      GP => a%GP(igaus)
      
      !Normal damage
      if (.not. a%Material%DeviatoricOnlyDamage) then
         !Elastic predictor
         ElasticPredictor = matmul(a%Material%C,strain)
      
         !Final Stress 
         GP%S = ElasticPredictor*(1-GP%CurrentD)
         
      !DeviatoricOnlyDamage   
      else   
         call a%Material%CT%GetPointerToCDeviatoric(CDev)
         call a%Material%CT%GetPointerToCSpherical(CSph)
      
         GP%S = 0.0_rp
         
         !Spherical part, not damaged
         GP%S = matmul(CSph,strain) + (1-GP%CurrentD)*matmul(CDev,strain)
      endif
      
   end subroutine
   
   subroutine ComputeCTensorTangent(a,ElementSize,ElasticPredictor,GP)
      class(DamageMaterial), target :: a
      real(rp), intent(in) :: ElementSize
      real(rp), Intent(in) :: ElasticPredictor(:) 
      class(DamageMaterial_GPD) :: GP
     
      real(rp), pointer :: CDev(:,:),CSph(:,:)
      
      !Normal damage
      if (.not. a%DeviatoricOnlyDamage) then
         !C Tensor

         GP%C = a%C*(1-GP%CurrentD)
      
      !Deviatoric only damage
      else   
         call a%CT%GetPointerToCDeviatoric(CDev)
         call a%CT%GetPointerToCSpherical(CSph)
         
         GP%C = CSph + Cdev*(1-GP%CurrentD)
      endif
      
   end subroutine   
   
   subroutine MD_GetStressTensorPointer(a,igaus,S)
      class(DamageMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: S(:)
      
      S => a%GP(igaus)%S(1:a%Material%CT%VoigtSize)
   end subroutine
   
   subroutine MD_GetConstitutiveTensorPointer(a,igaus,C)
      class(DamageMaterial_ED), target:: a
      integer(ip) :: igaus
      real(rp), pointer :: C(:,:)
      
      C => a%GP(igaus)%C
   end subroutine
   
   subroutine MD_SetElementSize(a,h)
      class(DamageMaterial_ED) :: a
      real(rp) :: h
      
      a%ElementSize = h
   end subroutine
   
   subroutine MD_MoveHistoryVariablesToConverged(a)
      class(DamageMaterial_ED) :: a
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         if (a%GP(igaus)%CurrentTau > a%GP(igaus)%ReferenceTau) then
            a%GP(igaus)%ReferenceTau = a%GP(igaus)%CurrentTau
            a%GP(igaus)%ReferenceD = a%GP(igaus)%CurrentD
         endif
      enddo
   end subroutine
   
   subroutine MD_SetupPostprocess(a,Mesh,Memor)
      use Mod_Mesh
      use Mod_Memor
      use Mod_RpMeshAllocator
      implicit none
      class(DamageMaterial) :: a
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      
      !Allocate array for gauss point info
      call AllocR1p(Memor,Mesh,a%DPostprocess)
      call AllocR1p(Memor,Mesh,a%TauPostprocess)
      call AllocR1p(Memor,Mesh,a%CompPostprocess)
      
   end subroutine
   
   subroutine MD_ContributeToPostprocess(a,ielem)
      class(DamageMaterial_ED) :: a
      integer(ip) :: ielem
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         a%Material%DPostprocess(ielem)%a(igaus) = a%GP(igaus)%CurrentD
         a%Material%CompPostprocess(ielem)%a(igaus) = a%GP(igaus)%ComparisonRatio
         a%Material%TauPostprocess(ielem)%a(igaus) = a%GP(igaus)%CurrentTau
      enddo
   end subroutine
   
   subroutine MD_DoPostprocess(a,istep,ctime,FilePostpr,Mesh,Memor,iiter_char)
      use Mod_Mesh
      use Mod_Postpr
      use Mod_Memor
      use Mod_RpMeshAllocator   
      use Mod_int2str
      implicit none
      class(DamageMaterial) :: a
      class(PostprFile) :: FilePostpr
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: istep
      real(rp) :: ctime
      character(5) :: iiter_char
      
      call FilePostpr%postgp(a%DPostprocess,'Damage'//iiter_char,istep,ctime,Mesh)
      call FilePostpr%postgp(a%TauPostprocess,'Tau'//iiter_char,istep,ctime,Mesh)
      call FilePostpr%postgp(a%CompPostprocess,'ComparisonRatio'//iiter_char,istep,ctime,Mesh)
      
      !Deallocate array for GaussPoint info
      call DeallocR1p(Memor,Mesh,a%DPostprocess)
      call DeallocR1p(Memor,Mesh,a%TauPostprocess)
      call DeallocR1p(Memor,Mesh,a%CompPostprocess)
   end subroutine
   
   subroutine GetMaterialPointer(a,Material)
      class(DamageMaterial_ED), target :: a
      class(PLCDMaterial), pointer :: Material
      
      Material => a%Material
      
   end subroutine
   
   subroutine GetHistoryVariablesNdofn(a,ndofn)
      class(DamageMaterial_ED) :: a
      integer(ip) :: ndofn
      
      !In this subroutine the dimension of the history array to be communicated 
      !in the adaptive process needs to be set
      !This subroutine is to be overwritten by each material
      
      ndofn = 1
   end subroutine
   
   subroutine SpecificHistoryGetToGaussArray(a,GaussArray)
      class(DamageMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         GaussArray(1,igaus) = a%GP(igaus)%ReferenceTau
      enddo
      
      
   end subroutine
   
   subroutine SpecificHistoryAddFromGaussArray(a,GaussArray)
      class(DamageMaterial_ED) :: a
      real(rp) :: GaussArray(:,:)
      
      integer(ip) :: igaus
      
      !If this is being added from multiple elements, we need to set it to zero the first time
      if (a%JustCreated) then
         a%JustCreated = .false.
         a%GP(:)%ReferenceTau = 0.0_rp
      endif
      
      do igaus = 1,size(a%GP,1)
         !Reference Tau
         a%GP(igaus)%ReferenceTau = a%GP(igaus)%ReferenceTau + GaussArray(1,igaus) 
      enddo
   end subroutine
   
   subroutine SpecificHistoryFinalRefinementUpdate(a)
      class(DamageMaterial_ED) :: a
      real(rp) :: ElasPred(a%Material%CT%VoigtSize)
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         !Reference Tau
         a%GP(igaus)%ReferenceTau = max(a%GP(igaus)%ReferenceTau,a%Material%Tau0)
         !Current Tau
         a%GP(igaus)%CurrentTau = a%GP(igaus)%ReferenceTau
         
         !Damage parameter and C tensor
         if (a%GP(igaus)%ReferenceTau > a%Material%Tau0) then
            call a%Material%ComputeDamage(a%ElementSize,a%GP(igaus))
            call ComputeCT_Secant(a%Material,a%ElementSize,ElasPred,a%GP(igaus))
            a%GP(igaus)%ReferenceD = a%GP(igaus)%CurrentD
         endif
      enddo
   end subroutine
   
   
   
   subroutine CopyTo(a,EMD)
      use Mod_plcd_Material
      class(DamageMaterial_ED) :: a
      class(ElementMaterialData), pointer :: EMD 
      
      integer(ip) :: igaus
      
      select type(EMD)
      type is (DamageMaterial_ED)
         do igaus = 1,size(a%GP,1)
            EMD%GP(igaus)%ReferenceTau = a%GP(igaus)%ReferenceTau
            EMD%GP(igaus)%ReferenceD   = a%GP(igaus)%ReferenceD 
            EMD%GP(igaus)%CurrentTau   = a%GP(igaus)%CurrentTau 
            EMD%GP(igaus)%CurrentD     = a%GP(igaus)%CurrentD 
            
            if (EMD%GP(igaus)%CurrentD /= 0.0_rp) then
               EMD%GP(igaus)%IsDamaged = .true.
               allocate(EMD%GP(igaus)%C(a%Material%CT%VoigtSize,a%Material%CT%VoigtSize))
               EMD%GP(igaus)%C = a%GP(igaus)%C
            endif
         enddo
      end select
      
   end subroutine
   
   subroutine Finalize(a)
      class(DamageMaterial_ED) :: a
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         deallocate(a%GP(igaus)%S)
         if (a%GP(igaus)%IsDamaged) then
            deallocate(a%GP(igaus)%C)
         endif
      enddo
      deallocate(a%GP)
   
   end subroutine
   
   subroutine MD_IsElementSizeRequiredByEMDs(a,isRequired)
      implicit none
      class(DamageMaterial), target :: a
      logical :: IsRequired
      
      isRequired = .true.
   end subroutine   
   
   subroutine GetBufferSize(a,BufferSize)
      class(DamageMaterial_ED) :: a
      integer(ip) :: BufferSize
      
      BufferSize = size(a%GP,1)
   end subroutine
   
   subroutine ToBuffer(a,Buffer)
      class(DamageMaterial_ED) :: a
      real(rp) :: Buffer(*)
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         Buffer(igaus) = a%GP(igaus)%ReferenceTau
      enddo
   end subroutine
   
   subroutine InitializeFromBuffer(a,Buffer)
      class(DamageMaterial_ED) :: a
      real(rp) :: Buffer(*)
      real(rp) :: ElasPred(a%Material%CT%VoigtSize)
      
      integer(ip) :: igaus
      
      do igaus = 1,size(a%GP,1)
         a%GP(igaus)%ReferenceTau = Buffer(igaus)
         a%GP(igaus)%CurrentTau = a%GP(igaus)%ReferenceTau
         
         !Damage parameter and C tensor
         if (a%GP(igaus)%ReferenceTau > a%Material%Tau0) then
            call a%Material%ComputeDamage(a%ElementSize,a%GP(igaus))
            call ComputeCT_Secant(a%Material,a%ElementSize,ElasPred,a%GP(igaus))
            a%GP(igaus)%ReferenceD = a%GP(igaus)%CurrentD
         endif
      enddo
   end subroutine
   
   subroutine GetComparisonCriteriaRatio(a,ComparisonRatio)
      class(DamageMaterial_ED) :: a
      real(rp) :: ComparisonRatio
      
      ComparisonRatio = maxval(a%GP(:)%ComparisonRatio)
      
   end subroutine
   
    !For u-p elements
   subroutine GetInverseVolumetricDeformationModulus(a,igaus,invK)
      class(DamageMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: invK
      
      call a%Material%CT%GetInverseVolumetricDeformationModulus(invK)
      
      if (.not. a%Material%DeviatoricOnlyDamage) invk = invk/(1.0_rp-a%GP(igaus)%CurrentD) !It must be divided, since it is the inverse!!!!
   end subroutine
   
   subroutine GetSecantCNorm(a,igaus,SecantCNorm)
      class(DamageMaterial_ED) :: a
      integer(ip) :: igaus
      real(rp) :: SecantCNorm
      real(rp), pointer :: CDev(:,:), CSph(:,:)
      
      real(rp) :: EigenValues(size(a%Material%C,1))

      !Normal Damage
      if (.not. a%Material%DeviatoricOnlyDamage) then
         call a%Material%CT%GetCNorm(SecantCNorm)
         SecantCNorm = SecantCNorm*(1.0_rp-a%GP(igaus)%CurrentD)
      
      !Deviatoric Only damage 
      else
         call a%Material%CT%GetPointerToCDeviatoric(CDev)
         call a%Material%CT%GetPointerToCSpherical(CSph)
         
         EigenValues = 0.0_rp
         call MatrixEigenValues(size(CDev,1),CDev*(1-a%GP(igaus)%CurrentD)+CSph,Eigenvalues)
         SecantCNorm = maxval(abs(Eigenvalues))
      
      endif
      
   end subroutine
  
   subroutine SetPressureForComputeHistoryAndConstitutiveTensor(a,Pressure)
      class(DamageMaterial_ED) :: a
      real(rp) :: Pressure
      
      a%gpPressure = Pressure
      
   end subroutine

   
! Constitutive tensor options
   
   subroutine ComputeCT_Elastic(a,ElementSize,ElasticPredictor,GP)
      class(DamageMaterial) :: a
      real(rp), intent(in) :: ElementSize
      real(rp), intent(in) :: ElasticPredictor(:)
      class(DamageMaterial_GPD) :: GP
    
      GP%C = a%C
      
   end subroutine
   
   subroutine ComputeCT_Secant(a,ElementSize,ElasticPredictor,GP)
      class(DamageMaterial) :: a
      real(rp), intent(in) :: ElementSize
      real(rp), intent(in) :: ElasticPredictor(:)
      class(DamageMaterial_GPD) :: GP
      
      real(rp), pointer :: CDev(:,:), CSph(:,:)
      
      !Normal damage
      if (.not. a%DeviatoricOnlyDamage) then

         GP%C = a%C*(1-GP%CurrentD)
      
      !Deviatoric only damage
      else   
         call a%CT%GetPointerToCDeviatoric(CDev)
         call a%CT%GetPointerToCSpherical(CSph)
         
         GP%C = CSph + Cdev*(1-GP%CurrentD)
      endif
      
   end subroutine
   
   subroutine ComputeCT_Tangent_VonMises_Exp(a,ElementSize,ElasticPredictor,GP)
      class(DamageMaterial) :: a
      real(rp), intent(in) :: ElementSize
      real(rp), intent(in) :: ElasticPredictor(:)   
      class(DamageMaterial_GPD) :: GP
      
      real(rp), pointer :: CDev(:,:), CSph(:,:)
      real(rp) :: Aparam, gf, raux1, facTenTang, E0
      real(rp) :: dtau_de(a%CT%VoigtSize)
      integer(ip) :: is, js
      
      !Normal damage
      if (.not. a%DeviatoricOnlyDamage) then
      
         gf = a%FractureEnergy/ElementSize
      
         call a%CT%GetYoungModulusEquiv(E0)
      
         raux1 = gf*E0/(a%Tau0**2)-0.5_rp
         Aparam = 1/raux1

         dtau_de = ElasticPredictor/GP%CurrentTau
         facTenTang = ((a%Tau0 + Aparam*GP%CurrentTau)/GP%CurrentTau**2)*exp(Aparam*(1-GP%CurrentTau/a%Tau0))
         GP%C = a%C*(1-GP%CurrentD)
         
         Do is = 1, a%CT%VoigtSize
            Do js = 1, a%CT%VoigtSize
               GP%C(js,is) = GP%C(js,is) - facTenTang*ElasticPredictor(js)*dtau_de(is)
            End do
         End do
      
      !Deviatoric only damage
      else   
         call a%CT%GetPointerToCDeviatoric(CDev)
         call a%CT%GetPointerToCSpherical(CSph)
         
         GP%C = CSph + Cdev*(1-GP%CurrentD)
      endif

   end subroutine   
   
! Damage types

   subroutine ComputeDamage_Exp(a,ElementSize,GP)
      class(DamageMaterial) :: a
      real(rp), intent(in) :: ElementSize
      class(DamageMaterial_GPD) :: GP
      real(rp) :: E0, Aparam, gf, raux1
      
      !Fracture energy
      gf = a%FractureEnergy/ElementSize
      
      call a%CT%GetYoungModulusEquiv(E0)
      
      raux1 = E0*gf/(a%Tau0**2) - 0.5_rp
      Aparam = 1/raux1
      
      GP%CurrentD = 1.0_rp - a%Tau0/GP%CurrentTau*exp(Aparam*(1.0_rp-GP%CurrentTau/a%Tau0))
      GP%CurrentD = min(GP%CurrentD,a%DamageLimit)
      if (GP%CurrentD < 0.0_rp) call runend('Negative damage, fracture energy too low for current mesh')
      
   end subroutine
   
   subroutine ComputeDamage_Dev_Exp(a,ElementSize,GP)
      class(DamageMaterial) :: a
      real(rp), intent(in) :: ElementSize
      class(DamageMaterial_GPD) :: GP
      real(rp) :: Aparam, gf, raux1
      real(rp) :: Cnorm
      
      !Fracture energy
      gf = a%FractureEnergy/ElementSize
      
      !Cnorm is usually 2G
      call a%CT%GetCNorm(Cnorm)
      raux1 = gf*1.5_rp*Cnorm/(a%Tau0**2)-0.5_rp
      Aparam = 1/raux1
      
      GP%CurrentD = 1.0_rp - a%Tau0/GP%CurrentTau*exp(Aparam*(1.0_rp-GP%CurrentTau/a%Tau0))
      GP%CurrentD = min(GP%CurrentD,a%DamageLimit)
      if (GP%CurrentD < 0.0_rp) call runend('Negative damage, fracture energy too low for current mesh')
      
   end subroutine 
   
end module

subroutine plcd_memall(a)
   !------------------------------------------ -----------------------------
   !****f* PLCD/plcd_memall
   ! NAME
   !    plcd_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_Mesh
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   use Mod_PLCD
   use Mod_plcd_StrainGenerator
   use Mod_Element
   use Mod_r1pElementAllocation
   use Mod_plcd_SIMP_TopologyOptimization
   use Mod_plcd_TD_TopologyOptimization
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none

   class(PLCDProblem), target :: a
   class(FiniteElement), pointer :: e => NULL()

   type(TimeIntegratorDt2) :: Integrator
   integer(ip) :: nsteps,ndime,npoin,pgaus,nelem,ncomp,ielem,pnode,istage,nboun

   procedure(), pointer :: GetStrain => NULL()
   procedure(), pointer :: AddSphericalComponent => NULL()
   procedure(), pointer :: GetStressTensor => NULL()
   procedure(), pointer :: GetStrainTensor => NULL()
   procedure(), pointer :: GetStressVector => NULL()
   procedure(), pointer :: GetStrainVector => NULL()
   integer(ip) :: vsize,imaterial
   logical :: isRequired

   if (a%UseUPFormulation) then
      !We include a degree of freedom for the pressure
      a%ndofn = a%ndofn+1

   endif


   !Materials allocation
   call a%Mesh%GetNelem(nelem)
   !call a%Memor%alloc(nelem,a%AuxRead_ElementMaterialList,'AuxRead_ElementMaterialList','plcd_memall')

   allocate(a%ElementMaterialsData(nelem))
   call a%Memor%allocObj(0,'ElementMaterialsData','plcd_memall',nelem)


   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)

   call Integrator%Init(a%kfl_tsche_2nd_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)

   a%ncomp = nsteps + 1
   ncomp = a%ncomp

   call a%Memor%alloc(ndime,npoin,ncomp,a%Displacement, 'Displacement','plcd_memall')
   if (a%kfl_TransientProblem /= 0) then
      call a%Memor%alloc(ndime,npoin,ncomp,a%Velocity,'Velocity','plcd_memall')
      call a%Memor%alloc(ndime,npoin,ncomp,a%Acceleration,'Acceleration','plcd_memall')
   endif
   
   if (a%kfl_FSI /= 0) then
      call a%Memor%alloc(ndime,npoin,a%btraction, 'SolidTraction','plcd_memall')
   endif
   
   call a%Memor%alloc(a%ndofn,npoin,a%InternalForcesVector, 'InternalForcesVector','plcd_memall')
   call a%Memor%alloc(a%ndofn,npoin,a%ExternalForcesVector, 'ExternalForcesVector','plcd_memall')
   call a%Memor%alloc(a%ndofn,npoin,a%ResidualForcesVector, 'ResidualForcesVector','plcd_memall')


   !Stages data allocation
   do istage = 1,a%NumberOfStages
      call a%Memor%alloc(ndime,npoin,a%Stages(istage)%NodalForces,'NodalForces','plcd_memall')
      call a%Memor%alloc(ndime,npoin,a%Stages(istage)%kfl_fixno,'kfl_fixno','plcd_memall')
      call a%Memor%alloc(ndime,npoin,a%Stages(istage)%bvess,'bvess','plcd_memall')
      call a%Memor%alloc(nboun,a%Stages(istage)%kfl_fixbo,'kfl_fixbo','plcd_memall')
      call a%Memor%alloc(nboun,a%Stages(istage)%bvnat,'bvnat','plcd_memall')
      call a%Memor%alloc(nboun,a%Stages(istage)%kfl_funbo,'kfl_funbo','plcd_memall')
      a%Stages(istage)%bvnat_coun = 0
      !Default is no boundary condition
      a%Stages(istage)%kfl_fixno = -1
      a%Stages(istage)%kfl_fixbo = -1
   enddo

   a%cs => a%Stages(a%CurrentStage)
   a%css => a%cs%Substages(a%cs%CurrentSubstage)

   !Stresses
   call GetVoigtsize(ndime,vsize)
   call AllocateGetStrain(ndime,GetStrain,a%kfl_LargeStrains)
   call AllocateAddSphericalComponent(ndime,AddSphericalComponent)
   call AllocateGetStressTensor(ndime,GetStressTensor)
   call AllocateGetStrainTensor(ndime,GetStrainTensor)
   call AllocateGetStressVector(ndime,GetStressVector)
   call AllocateGetStrainVector(ndime,GetStrainVector)
   
   call AllocateR2pElement(a%Mesh,vsize,a%Stress,a%Memor,'stress')
   !Strains
   if (a%PostprocessStrain) then
      call AllocateR2pElement(a%Mesh,vsize,a%Strain,a%Memor,'strain')
   endif

   !Mean Compute
   if (a%ComputeMeanflag) then
      call a%Memor%alloc(vsize,a%MeanStress,'MeanStress','plcd_memall')
      call a%Memor%alloc(vsize,a%GMeanStress,'GMeanStress','plcd_memall')
      call a%Memor%alloc(vsize,vsize,a%MeanConst,'MeanConst','plcd_memall')
   endif


   if (a%UseSmoothedDisplacementGradient) then
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(ndime,ndime,npoin,a%SmoothedDisplacementGradient,'SmoothedDisplacementGradient','plcd_memall')
   endif

   !UseUPFormulation
   if (a%UseUPFormulation) then
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(npoin,ncomp,a%pressure,'pressure','plcd_memall')
      call a%Memor%alloc(ndime,npoin,a%UPResidualProjection,'UPResidualProjection','plcd_memall')
      if (a%RefinerErrorEstimator == 'SUBSC') a%UPStoreSubscales = .true.

      if (a%UPStoreSubscales) then
         call a%Mesh%GetNdime(ndime)
         call AllocateR2pElement(a%Mesh,ndime+1,a%UPSubscales,a%Memor,'UPSubscales')
      endif
   endif
   
    
   !We need to initialize materials here, before we enter the boundary conditions reading
   !Initially compute the constitutive tensor for each material
!    do imaterial = 1,a%NumberOfMaterials
!       if (a%UseUPFormulation) call a%Materials(imaterial)%p%SetUseUPFormulation(.true.)
!       call a%Materials(imaterial)%p%Setup
! 
!       !Check if any material needs the element size so that it is computed
!       call a%Materials(imaterial)%p%IsElementSizeRequiredByEMDs(IsRequired)
!       if (IsRequired) a%IsElementSizeRequiredByEMDs = .true.
!    enddo
   
   

   
end subroutine plcd_memall

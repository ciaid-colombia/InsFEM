subroutine plcd_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
!   use Mod_iofile
   use Mod_PLCD
   use Mod_Element
   use Mod_plcd_StrainGenerator
   use Mod_r1pElementAllocation
   use Mod_plcd_SIMP_TopologyOptimization
   use Mod_plcd_TD_TopologyOptimization
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem) :: a

   integer(ip) :: nelem,ndime,npoin,istage,ielem,vsize
   class(FiniteElement), pointer :: e => NULL()
   procedure(), pointer :: GetStrain => NULL()
   procedure(), pointer :: AddSphericalComponent => NULL()
   procedure(), pointer :: GetStressTensor => NULL()
   procedure(), pointer :: GetStrainTensor => NULL()
   procedure(), pointer :: GetStressVector => NULL()
   procedure(), pointer :: GetStrainVector => NULL()
   
   integer(ip) :: nsubstages

   deallocate(a%Materials)
   call a%Memor%deallocObj(0,'Materials','plcd_turnof',1*a%NumberOfMaterials)

   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)




   call a%Memor%dealloc(size(a%Displacement,1),size(a%Displacement,2),size(a%Displacement,3),a%Displacement, 'Displacement','plcd_turnof')
   if (a%kfl_TransientProblem /= 0) then
      call a%Memor%dealloc(size(a%Velocity,1),size(a%Velocity,2),size(a%Velocity,3),a%Velocity,'Velocity','plcd_turnof')
      call a%Memor%dealloc(size(a%Acceleration,1),size(a%Acceleration,2),size(a%Acceleration,3),a%Acceleration,'Acceleration','plcd_turnof')
   endif
   if (a%kfl_FSI /= 0) then
      call a%Memor%dealloc(size(a%btraction,1),npoin,a%btraction, 'SolidTraction','plcd_turnof')
   endif
   call a%Memor%dealloc(a%ndofn,npoin,a%InternalForcesVector, 'InternalForcesVector','plcd_turnof')
   call a%Memor%dealloc(a%ndofn,npoin,a%ExternalForcesVector, 'ExternalForcesVector','plcd_turnof')
   call a%Memor%dealloc(a%ndofn,npoin,a%ResidualForcesVector, 'ResidualForcesVector','plcd_turnof')

   call a%DeallocStages

   !Stresses
   call GetVoigtsize(ndime,vsize)
   call AllocateGetStrain(ndime,GetStrain,a%kfl_LargeStrains)
   call AllocateAddSphericalComponent(ndime,AddSphericalComponent)
   call AllocateGetStressTensor(ndime,GetStressTensor)
   call AllocateGetStrainTensor(ndime,GetStrainTensor)
   call AllocateGetStressVector(ndime,GetStressVector)
   call AllocateGetStrainVector(ndime,GetStrainVector)
   
   call DeAllocateR2pElement(a%Mesh,vsize,a%Stress,a%Memor,'stress')
   !Strains
   if(a%PostprocessStrain) then
      call DeAllocateR2pElement(a%Mesh,vsize,a%Strain,a%Memor,'strain')
   endif

   !Mean Compute
   if (a%ComputeMeanflag) then
      call a%Memor%dealloc(vsize,a%MeanStress,'MeanStress','plcd_turnof')
      call a%Memor%dealloc(vsize,a%GMeanStress,'GMeanStress','plcd_turnof')
      call a%Memor%dealloc(vsize,vsize,a%MeanConst,'MeanConst','plcd_turnof')
   endif

   do ielem = 1,nelem
      call a%ElementMaterialsData(ielem)%p%Finalize
      deallocate(a%ElementMaterialsData(ielem)%p)
   enddo
   deallocate(a%ElementMaterialsData)
   call a%Memor%deallocObj(0,'ElementMaterialsData','plcd_turnof',nelem)

   if (a%UseSmoothedDisplacementGradient) then
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(ndime,ndime,npoin,a%SmoothedDisplacementGradient,'SmoothedDisplacementGradient','plcd_turnof')
   endif

    !UseUPFormulation
   if (a%UseUPFormulation) then
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(size(a%pressure,1),size(a%pressure,2),a%pressure,'pressure','plcd_turnof')
      call a%Memor%dealloc(ndime,npoin,a%UPResidualProjection,'UPResidualProjection','plcd_turnof')


      if (a%UPStoreSubscales) then
         call a%Mesh%GetNdime(ndime)
         call DeallocateR2pElement(a%Mesh,ndime,a%UPSubscales,a%Memor,'UPSubscales')
      endif

   endif

   !TopologyOptimization
   if (a%kfl_TopologyOptimization == 1) then
      call TurnofSIMP(a)
    elseif (a%kfl_TopologyOptimization==2) then
      call TurnofTD(a)
   end if
   
   !Stochastic Optimization
   if (a%kfl_StochasticTopologyOptimization == 1) then
      call plcd_TDSTO_turnof(a)
   endif


end subroutine plcd_turnof

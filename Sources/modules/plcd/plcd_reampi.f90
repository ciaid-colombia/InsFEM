subroutine plcd_reampi(a)
   use typre
   use MPI
   use Mod_BroadCastBuffer
   use Mod_plcd_MaterialFactory
   use Mod_PLCD
   use Mod_plcd_ReadMaterials
   implicit none
   class(PLCDProblem) :: a

   integer(ip) :: ierr
   integer(ip) :: idnum,nsubstages,istage
   integer(ip) :: ndime, imaterial
   logical     :: isRequired

   call a%Mesh%GetNdime(ndime)
   call ReadMaterials_Scatter(a%NumberOfMaterials,a%Materials,a%AuxRead_MaterialTypeList,a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank,ndime,a%kfl_LargeStrains,a%Memor)

   !Stages Data
   call MPI_BCAST(a%NumberOfStages,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
   if (a%MPIrank /= a%MPIroot) then
      allocate(a%Stages(a%NumberOfStages))
      call a%Memor%allocObj(0,'Stages','plcd_memall',1*a%NumberOfMaterials)
   endif
   do istage = 1,a%NumberOfStages
      call MPI_BCAST(a%Stages(istage)%NumberOfSubStages,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      nsubstages = a%Stages(istage)%NumberOfSubStages
      if (a%MPIrank /= a%MPIroot) then
         allocate(a%Stages(istage)%Substages(nsubstages))
         call a%Memor%allocObj(0,'SubStages','plcd_memall',1*nsubstages)
      endif
      !Will only work as long as no allocatables in substage are present
      call MPI_BCAST(a%Stages(istage)%Substages,STORAGE_SIZE(a%Stages(istage)%Substages)/8,MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
    enddo

    !Numerical Treatment data
    call MPI_BCAST(a%UseSmoothedDisplacementGradient,1,MPI_LOGICAL,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%UseUPFormulation,1,MPI_LOGICAL,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%ErrorEstimatorTypeOfSubscales,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%UPStoreSubscales,1,MPI_LOGICAL,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%staco,size(a%staco),MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%kfl_LineSearch,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%LineSearchRadiusOfConvergence,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%LineSearchMaxIterations,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)

    !Postprocess strain

    call MPI_BCAST(a%PostprocessStrain,1,MPI_LOGICAL,a%MPIroot,a%MPIcomm,ierr)

    !Postprocess material data at each iteration
    call MPI_BCAST(a%kfl_PostprocessMatDataAtEachIteration,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    
    !Postprocess Displacement at each iteration
    call MPI_BCAST(a%kfl_PostprocessDisplacementAtEachIteration,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr) 
    
    !Large Strain information
    call MPI_BCAST(a%kfl_LargeStrains,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)

    !Topology Optimization information
    call MPI_BCAST(a%kfl_TopologyOptimization,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    if (a%kfl_TopologyOptimization == 1) then
      call a%SIMPData%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%SIMPData%ScatterData
    elseif (a%kfl_TopologyOptimization == 2) then
      call MPI_BCAST(a%kfl_StochasticTopologyOptimization,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call a%TDData%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%TDData%ScatterData
      if (a%kfl_StochasticTopologyOptimization == 1) then
         call a%STOData%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
         call a%STOData%ScatterData
      endif
    endif
    
    !Transient Problem information
    call MPI_BCAST(a%kfl_TransientProblem,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%Beta,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%Gamma,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
    
    !Gravity Force information
    call MPI_BCAST(a%kfl_GravityForce,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%gravity, size(a%gravity), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
    
    !Rotating Frame of Reference information
    call MPI_BCAST(a%kfl_RotatingFrame,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%angvelocitynorm2, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
    call MPI_BCAST(a%angvelocity, size(a%angvelocity), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)

    !Adaptive Mesh refinement
    call MPI_BCAST(a%NLayersRefinement,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%GeneralRefinementLevels,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    call MPI_BCAST(a%InterfaceRefinementLevels,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
    
    !We need to initialize materials here, before we enter the boundary conditions reading
    
   !Initially compute the constitutive tensor for each material
   do imaterial = 1,a%NumberOfMaterials

      if (a%UseUPFormulation) call a%Materials(imaterial)%p%SetUseUPFormulation(.true.)
      call a%Materials(imaterial)%p%Setup

      !Check if any material needs the element size so that it is computed
      call a%Materials(imaterial)%p%IsElementSizeRequiredByEMDs(IsRequired)
      if (IsRequired) a%IsElementSizeRequiredByEMDs = .true.
   enddo


end subroutine

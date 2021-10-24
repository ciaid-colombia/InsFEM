subroutine plcd_turnon(a)
   ! DESCRIPTION
   !    This routine starts the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use def_parame
!   use Mod_iofile
   use Mod_PLCD
   use Mod_Element
   use Mod_plcd_Material
   use Mod_plcd_SIMP_TopologyOptimization
   use Mod_plcd_TD_TopologyOptimization
   use Mod_plcd_TD_StochasticTopologyOptimization
   use MPI
   implicit none
   class(PLCDProblem), target :: a

   real(rp) :: time
   integer(ip) :: is, iss,ierr

   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: imaterial,ielem,nelem

   logical :: IsRequired
   
   interface
      subroutine plcd_ComputeInitialAcceleration(a)
      import
      implicit none
      class(PLCDProblem), target :: a
      
      end subroutine
   end interface

   !Setup stages and substages
   time = 0.0_rp
   do is = 1,a%NumberOfStages
      do iss = 1,a%Stages(is)%NumberOfSubstages
         a%Stages(is)%Substages(iss)%IniTime = time
         time = time + a%Stages(is)%Substages(iss)%TimeInterval
         a%Stages(is)%Substages(iss)%EndTime = time
      enddo
   enddo

!This is now done in plcd_memall, because the info is needed in read boundary conditions (for reading element info for each mat)
!    !Initially compute the constitutive tensor for each material
!    do imaterial = 1,a%NumberOfMaterials
! 
!       if (a%UseUPFormulation) call a%Materials(imaterial)%p%SetUseUPFormulation(.true.)
!       call a%Materials(imaterial)%p%Setup
! 
!       !Check if any material needs the element size so that it is computed
!       call a%Materials(imaterial)%p%IsElementSizeRequiredByEMDs(IsRequired)
!       if (IsRequired) a%IsElementSizeRequiredByEMDs = .true.
!    enddo


   !Allocate the data structures for the materials at each Gauss Point
!    call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
!    call a%Mesh%GetNelem(nelem)
!    do ielem = 1,nelem
!       !Load Element
!       call a%Mesh%ElementLoad(ielem,e)
!       imaterial = a%AuxRead_ElementMaterialList(ielem)
!       if (imaterial <= 0 .or. imaterial > a%NumberOfMaterials) call runend('Wrong material number')
!       call a%Materials(imaterial)%p%SpecificCreateElementData(e%pgaus,a%ElementMaterialsData(ielem)%p)
!    enddo
!    call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','plcd_Elmope')

   ! call a%Memor%dealloc(nelem,a%AuxRead_ElementMaterialList,'AuxRead_ElementMaterialList','plcd_memall')

   !Topology Optimization
   if (a%kfl_TopologyOptimization == 1) then
      call MemallSIMP(a)
   elseif (a%kfl_TopologyOptimization==2) then
      call MemallTD(a)
   end if
   
   !Large Strains, Updated Formulation
   if (a%kfl_LargeStrains == 1) then
      !So we can move the mesh for Updated Formulation
      call a%Mesh%SetALE(1_ip)
      call a%Mesh%SetDisplacements(a%Displacement(:,:,:))
   endif    
   
   !Transient Problem, Compute Initial Accelerations
   if (a%kfl_TransientProblem == 1) then
      call plcd_ComputeInitialAcceleration(a)
   endif

end subroutine plcd_turnon


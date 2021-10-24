module Mod_plcd_EnditeElmope
   use Mod_plcd_BaseElmope
   use Mod_plcd_PassElementSizeToEMDs
   use Mod_plcd_ComputeSmoothDisplacementGradients
   use Mod_plcd_UPFormulation
   use Mod_plcd_PostprocessStrain
   use Mod_plcd_ComputeMean
   use Mod_plcd_ExternalForces
   use Mod_plcd_SIMP_TopologyOptimization
   use Mod_plcd_TD_TopologyOptimization
   use Mod_plcd_TD_Derivative
   use Mod_int2str
   use Mod_plcd_LargeStrainsOperations
   use Mod_plcd_TransientProblem
   use Mod_plcd_RotatingFrame
   implicit none

contains

   subroutine SetPointers

      call ResetProcedureComposition
      call SetPointersAndHooksToNULLSUB

      !Initializations

      call SetPointersPassElementSizeToEMDs%Initialize
      call SetPointersComputeSmoothDisplacementGradients%Initialize
      call SetPointersUPFormulationEndite%Initialize
      call SetPointersPostprocessStrain%Initialize
      call SetPointersComputeMean%Initialize
      call SetPointersExternalForces%Initialize
      call SetPointersSIMP_TopologyOptimization%Initialize
      call SetPointersTD_TopologyOptimization%Initialize
      call SetPointersLargeStrainsEndite%Initialize
      call SetPointersTransientProblemEndite%Initialize
      call SetPointersRotatingFrameEndite%Initialize

      !SetPointers
      call SetPointersPassElementSizeToEMDs%Set
      call SetPointersComputeSmoothDisplacementGradients%Set
      call SetPointersUPFormulationEndite%Set
      call SetPointersPostprocessStrain%Set
      call SetPointersComputeMean%Set
      call SetPointersExternalForces%Set
      call SetPointersSIMP_TopologyOptimization%Set
      call SetPointersTD_TopologyOptimization%Set
      call SetPointersLargeStrainsEndite%Set
      call SetPointersTransientProblemEndite%Set
      call SetPointersRotatingFrameEndite%Set

      !Finalizations
      call SetPointersPassElementSizeToEMDs%Finalize
      call SetPointersComputeSmoothDisplacementGradients%Finalize
      call SetPointersUPFormulationEndite%Finalize
      call SetPointersPostprocessStrain%Finalize
      call SetPointersComputeMean%Finalize
      call SetPointersExternalForces%Finalize
      call SetPointersSIMP_TopologyOptimization%Finalize
      call SetPointersTD_TopologyOptimization%Finalize
      call SetPointersLargeStrainsEndite%Finalize
      call SetPointersTransientProblemEndite%Finalize
      call SetPointersRotatingFrameEndite%Finalize

   end subroutine

end module



subroutine plcd_EnditeElmope(b)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrix
   use Mod_plcd_BMatrixFactory
   use Mod_plcd_EnditeElmope
   implicit none
   class(PLCDProblem), target :: b

!   real(rp), pointer :: stress(:) => NULL()

   integer(ip) :: idofn,inode,ipoin,npoin, ielem2

   a=>b

   !SetPointers
   call SetPointers

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_EnditeElmope')

   call a%Memor%Alloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_EnditeElmope')
   call a%Memor%Alloc(e%ndime,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')

   call a%Memor%Alloc(a%ndofn,e%mnode,elInternalForces,'elInternalForces','plcd_EnditeElmope')
   call a%Memor%Alloc(a%ndofn,e%mnode,elExternalForces,'elExternalForces','plcd_EnditeElmope')

   call a%Memor%alloc(a%ndofn,e%mnode,GaussElInternalForces,'GaussElInternalForces','plcd_EnditeElmope')
   call a%Memor%alloc(a%ndofn,e%mnode,GaussElExternalForces,'GaussElExternalForces','plcd_EnditeElmope')

   call CreateBMatrix(a,a%Memor,BMat)
   call BMat%Alloc(e,a%Memor)


   !Initializations
   !GradDispHistory might be different than gradDisp, for instance if smoothed gradients are used
   gradDispHistory => gradDisp

!       if (a%istep >= 64) then
!       write(*,*) 'asdfasdfasdfñlkj'
!    endif
!    
!    
! !       if (a%itera > 0) then
!          call a%FilePostpr%postpr(a%Displacement(:,:,1),'DispPre',a%istep,a%ctime,a%Mesh)
!          
!          call PostprocessTopologicalDerivative(a%TopologicalDerivative,a%FilePostpr,a%itera+1,a%istep,a%ctime)
!          
!          call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
!          call plcd_PostprocessMaterialData(a,int2str(a%itera+1))
!          call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
!       endif
   !Hook
   call ProcHook%Initializations
   
!    if (a%itera > 0) then
!    
!       call PostprocessTopologicalDerivative(a%TopologicalDerivative,a%FilePostpr,a%itera+2,a%istep,a%ctime)
!       call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
!       call plcd_PostprocessMaterialData(a,int2str(a%itera+2))
!       call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
!    endif


   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)

      ElementMatData => a%ElementMaterialsData(ielem)%p

      call e%gather(e%ndime,eldisp,a%Displacement(:,:,1))

      elInternalForces = 0.0_rp
      elExternalForces = 0.0_rp

      !Hook
      call ProcHook%PreGauss

      !Compute linear derivatives
      call e%elmdel

      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus

         call e%elmder
         dvol = e%detjm*e%weigp(e%igaus)

         !ComputeDisplacementGradients
         call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)
         
         !ProcPointer
         call ProcHook%InGauss

         !ComputeDisplacementGradients
         !call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)

         call BMat%Setup(e)
         call ElementMatData%ComputeHistoryAndConstitutiveTensor(e%igaus,gradDispHistory)
         call ElementMatData%ComputeStress(e%igaus,gradDisp)

         !Internal Forces
         GaussElInternalForces = 0.0_rp
         call ElementMatData%GetStressTensorPointer(e%igaus,stress)
         call BMat%Bt_Times_Vector(stress,GaussElInternalForces(1:e%ndime,1:e%pnode))
         !call BMat%Bt_Times_Vector(stress,GaussElInternalForces)
         GaussElInternalForces = GaussElInternalForces*dvol

         ielem2 = ielem

         GaussElExternalForces = 0.0_rp
         !Hook
         call ProcHook%InGaussElmats

         elInternalForces(:,1:e%pnode) = elInternalForces(:,1:e%pnode) + GaussElInternalForces(:,1:e%pnode)
         elExternalForces(:,1:e%pnode) = elExternalForces(:,1:e%pnode) + GaussElExternalForces(:,1:e%pnode)

         a%Stress(ielem)%a(:,igaus) = stress
      enddo gauss_points
      !Hook
      call ProcHook%PostGaussElmats

      call a%Mesh%AssemblyToArray(e,a%ndofn,elInternalForces,a%InternalForcesVector)
      call a%Mesh%AssemblyToArray(e,a%ndofn,elExternalForces,a%ExternalForcesVector)

   enddo elements

   !Contribution to the external forces of the Dirichlet nodes
   call a%Mesh%GetNpoin(npoin)
   do ipoin = 1,npoin
      do idofn = 1,a%ndofbc
         if (a%cs%kfl_fixno(idofn,ipoin) == 1) then
            a%ExternalForcesVector(idofn,ipoin) = a%InternalForcesVector(idofn,ipoin)
         endif
      enddo
   enddo

   call BMat%DeAlloc(e,a%Memor)
   call DestroyBMatrix(a,a%Memor,BMat)

   call a%Memor%dealloc(a%ndofn,e%mnode,GaussElInternalForces,'GaussElInternalForces','plcd_EnditeElmope')
   call a%Memor%dealloc(a%ndofn,e%mnode,GaussElExternalForces,'GaussElExternalForces','plcd_EnditeElmope')


   call a%Memor%deAlloc(a%ndofn,e%mnode,elInternalForces,'elInternalForces','plcd_EnditeElmope')
   call a%Memor%deAlloc(a%ndofn,e%mnode,elExternalForces,'elExternalForces','plcd_EnditeElmope')


   call a%Memor%deAlloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_EnditeElmope')
   call a%Memor%deAlloc(e%ndime,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')

   !Hook
   call ProcHook%Finalizations

   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')


end subroutine

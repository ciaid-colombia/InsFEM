module Mod_plcd_Elmope
   use Mod_plcd_BaseElmope
   use Mod_plcd_HangingNodes
   use Mod_plcd_UPFormulation
   use Mod_plcd_LargeStrainsOperations
   use Mod_plcd_TransientProblem
   use Mod_plcd_RotatingFrame

contains

   subroutine SetPointers
      call ResetProcedureComposition

      call SetPointersAndHooksToNULLSUB

      !Initialize setpointers
      call SetPointersHangingNodes%Initialize
      call SetPointersUPFormulation%Initialize
      call SetPointersLargeStrains%Initialize
      call SetPointersTransientProblem%Initialize
      call SetPointersRotatingFrame%Initialize

      !SetPointers
      call SetPointersHangingNodes%Set
      call SetPointersUPFormulation%Set
      call SetPointersLargeStrains%Set
      call SetPointersTransientProblem%Set
      call SetPointersRotatingFrame%Set

      !Finalize SetPointers
      call SetPointersHangingNodes%Finalize
      call SetPointersUPFormulation%Finalize
      call SetPointersLargeStrains%Finalize
      call SetPointersTransientProblem%Finalize
      call SetPointersRotatingFrame%Finalize

   end subroutine
end module



subroutine plcd_Elmope(b)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrixFactory
   use Mod_plcd_Stages
   use Mod_plcd_elmdir
   use Mod_php_AssemblyVectorToSystem
   use Mod_plcd_Elmope
   use Mod_Debugging
   implicit none
   class(PLCDProblem), target :: b

   integer(ip) :: idofn, inode,jnode,ipoin,idime,jdime, ielem2

   !deb_PostprocessMatrix = 1

   a => b

   !SetPointers
   call SetPointers

   !This cannot be here
   !NodalForces => a%cs%NodalForces

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','plcd_elmope')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','plcd_elmope')

   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,GaussElmat,'GaussElmat','plcd_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,3,eldisp,'eldisp','plcd_elmope')
   call a%Memor%alloc(e%ndime,e%mnode,elNodalForces,'elNodalForces','plcd_elmope')

   !Hook
   call ProcHook%Initializations

   call CreateBMatrix(a,a%Memor,BMat)
   call BMat%Alloc(e,a%Memor)

   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem
      ielem2 = ielem

      !Load Element
      call a%Mesh%ElementLoad(ielem,e)

      !Hook
      call ProcHook%PreGauss

      elmat = 0.0_rp
      elrhs = 0.0_rp

      !Compute linear derivatives
      call e%elmdel

      ElementMatData => a%ElementMaterialsData(ielem)%p

      !Gathers
      !Displacements at the previous iteration
      call e%gather(e%ndime,eldisp(:,:,1),a%Displacement(:,:,1))

      !Gausspoint loop
      GaussPointLoop : do igaus = 1,e%pgaus
         e%igaus = igaus

         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         GaussElmat(:,1:e%pnode,:,1:e%pnode) = 0.0_rp

         call ElementMatData%GetConstitutiveTensorPointer(e%igaus,C)
         
         !Hook
         call ProcHook%InGauss

         call BMat%Setup(e)
         call BMat%Bt_Times_Matrix_Times_B(C,GaussElmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode))
         !call BMat%Bt_Times_Matrix_Times_B(C,GaussElmat)
         GaussElmat = GaussElmat*dvol

         !Hook
         call ProcHook%InGaussElmats

         !Elmat
         elmat(:,1:e%pnode,:,1:e%pnode) = elmat(:,1:e%pnode,:,1:e%pnode) + GaussElmat(:,1:e%pnode,:,1:e%pnode)
      enddo GaussPointLoop

      !Hook
      call ProcHook%PostGaussElmats

      !Hook
      call ProcHook%PreDirichlet

      !Dirichlet Boundary Conditions
      call plcd_elmdir(a,e,elmat,elrhs)

      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   enddo

   !Hook
   call ProcHook%Finalizations

   call BMat%DeAlloc(e,a%Memor)
   call DestroyBMatrix(a,a%Memor,BMat)

   call a%Memor%dealloc(e%ndime,e%mnode,3,eldisp,'eldisp','plcd_elmope')
   call a%Memor%dealloc(e%ndime,e%mnode,elNodalForces,'elNodalForces','plcd_elmope')

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_elmope')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','nsm_elmope')
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,GaussElmat,'GaussElmat','plcd_elmope')
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')



end subroutine plcd_Elmope

subroutine plcd_GetRefinementCriteria(b,markel)
   use typre
   use Mod_PLCD
   use Mod_Element
   use Mod_ZZErrorEstimator
   use Mod_plcd_BaseElmope
   use Mod_plcd_StrainGenerator
   use MPI
   implicit none
   class(PLCDProblem), target :: b
   integer(ip) :: markel(*)

   interface
      subroutine plcd_SubscalesRefCriteria(b,error,TotalEstimatedError)
      import
      implicit none
      class(PLCDProblem), target :: b
      real(rp) :: error(:),TotalEstimatedError
      end subroutine

      subroutine LayerRefinementCriteria(a,level,markel,Outlayers,InRefinementLevels,InterfaceRefinementLevels,GlobalRefinementLevels,Offset)
         use typre
         use Mod_PhysicalProblem
         use Mod_Element
         use Mod_LinkedList
         implicit none
         class(PhysicalProblem) :: a
         real(rp) :: level(:)
         integer(ip) :: markel(*)

         integer(ip) :: OutLayers, InRefinementLevels,InterfaceRefinementLevels, GlobalRefinementLevels
         real(rp) :: Offset
      end subroutine

   end interface


   real(rp), allocatable, target :: error(:)


   real(rp), pointer :: auxNodalChi(:,:)
   real(rp), pointer :: level(:)
   real(rp) :: meanchi,offset
   real(rp) :: TotalEstimatedError

   integer(ip) :: npoin,ndime
   
   integer(ip) :: elclass, istage
   real(rp) :: xmax, xmin, ymax, ymin


   a=> b


   call a%Timer%Refine%Tic
   call a%Mesh%GetNelem(nelem)
   call a%Memor%alloc(nelem,error,'error','plcd_GetRefinementCriteria')
   call a%Mesh%GetNdime(ndime)

   !First one is for debugging with load rebalancing
   if (.false.) then
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','MeshInit')   
      istage = mod(a%istep,4)
      
      
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         xmax = maxval(e%elcod(1,1:e%pnode))
         xmin = minval(e%elcod(1,1:e%pnode))
         ymax = maxval(e%elcod(2,1:e%pnode))
         ymin = minval(e%elcod(2,1:e%pnode))
         
         if (xmax < 100 .and. ymax < 40) then
            elclass = 1
         elseif (xmax < 100 .and. ymax >= 40) then
            elclass = 2
         elseif (xmax >= 100 .and. ymax < 40) then
            elclass = 3
         else
            elclass = 4
         endif
         
         if (istage == 0) then
            if (elclass == 1) then
               markel(ielem) = 1
            elseif (elclass == 2) then
               markel(ielem) = -1
            elseif (elclass == 3) then
               markel(ielem) = -1
            elseif (elclass == 4) then
               markel(ielem) = 1
            endif
               
               
         elseif (istage == 1) then
            if (elclass == 1) then
               markel(ielem) = 1
            elseif (elclass == 2) then
               markel(ielem) = 1
            elseif (elclass == 3) then
               markel(ielem) = -1
            elseif (elclass == 4) then
               markel(ielem) = -1
            endif
         
         elseif (istage == 2) then
         
            if (elclass == 1) then
               markel(ielem) = -1
            elseif (elclass == 2) then
               markel(ielem) = 1
            elseif (elclass == 3) then
               markel(ielem) = 1
            elseif (elclass == 4) then
               markel(ielem) = -1
            endif
         
         
         elseif (istage == 3) then
            if (elclass == 1) then
               markel(ielem) = -1
            elseif (elclass == 2) then
               markel(ielem) = -1
            elseif (elclass == 3) then
               markel(ielem) = 1
            elseif (elclass == 4) then
               markel(ielem) = 1
            endif
         endif
      enddo

      call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','MeshInit')   
      
      
   elseif (a%RefinerErrorEstimator /= 'LAYER') then
      if (a%RefinerErrorEstimator == 'ZZ   ') then
         call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%displacement(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'GRADI') then

         call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%displacement(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'SUBSC') then
         call plcd_SubscalesRefCriteria(a,error,TotalEstimatedError)

      elseif (a%RefinerErrorEstimator == 'MATCO') then
         !Material Comparison Criteria
         error = 0.0_rp
         call a%Mesh%GetNelem(nelem)
         do ielem = 1,nelem
            call a%ElementMaterialsData(ielem)%p%GetComparisonCriteriaRatio(error(ielem))
         enddo
      elseif (a%RefinerErrorEstimator == 'CHIGR') then
         if (a%kfl_TopologyOptimization == 1) then
            auxNodalChi(1:1,1:size(a%SIMPData%NodalChi)) => a%SIMPData%NodalChi
            call GradientErrorEstimator(a%Mesh,a%Memor,1_ip,auxNodalChi,error)
         elseif (a%kfl_TopologyOptimization == 2) then
            auxNodalChi(1:1,1:size(a%TDData%NodalChi)) => a%TDData%NodalChi
            call GradientErrorEstimator(a%Mesh,a%Memor,1_ip,auxNodalChi,error)
         else
            call runend('PLCD: CHI Gradient refinement criteria only makes sense when running topology optimization ')
         end if
      endif

      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
   elseif (a%RefinerErrorEstimator == 'LAYER') then
      if (a%kfl_TopologyOptimization == 2) then
         !We use the error array memory space as a level array (error is not needed)
         call a%Mesh%GetNpoin(npoin)
         allocate(level(npoin))
         
         level = 0.0_rp
         call a%Mesh%GetNpoin(npoin)
         meanchi = (a%TDData%chimin + a%TDData%chimax)*0.5_rp
         where (a%TDData%NodalChi > meanchi) level = 1.0_rp
         where (a%TDData%NodalChi <= meanchi) level = -1.0_rp
         call LayerRefinementCriteria(a,level,markel,a%NLayersRefinement,0_ip,a%InterfaceRefinementLevels,a%GeneralRefinementLevels,0.0_rp)
         deallocate(level)
      elseif (a%kfl_topologyOptimization == 1) then
         Offset = 0.5_rp*(a%SIMPData%Chimax+a%SimpData%Chimin)
         call LayerRefinementCriteria(a,a%SIMPData%NodalChi,markel,a%NLayersRefinement,0_ip,a%InterfaceRefinementLevels,a%GeneralRefinementLevels,Offset)
      end if
   end if

   call a%Memor%dealloc(nelem,error,'error','plcd_GetRefinementCriteria')

   call a%Timer%Refine%Toc
end subroutine

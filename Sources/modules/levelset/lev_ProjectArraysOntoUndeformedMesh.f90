subroutine lev_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
   use typre
   use Mod_Mesh
   use Mod_Levelset
   use Mod_MeshInterpolator
   implicit none
   class(LevelSetProblem) :: a
   type(Interpolator) :: Interp
   integer(ip) :: itask

   integer(ip) :: ndime,icomp

   interface
      subroutine lev_CutElementsAndListByLayers(a)
      use Mod_LevelSet
      implicit none
      class(LevelSetProblem) :: a
      end subroutine
   end interface

   if (itask == 1) then

!       call a%FilePostpr%postpr(a%level(:,1),'levelpreproj',a%istep,a%ctime,a%Mesh)

      !If we force no ALE, then we do not need to project, since the levelset is Eulerian
      if (a%kfl_ForceEulerianAdvection == 0) then
         do icomp = 1,size(a%level,2)
            call Interp%Interpolate(1_ip,a%level(:,icomp),a%level(:,icomp))
         enddo
      endif
      if (a%kfl_SmoothGradient == 1) then
         call Interp%Interpolate(size(a%SmoothGradient,1),a%SmoothGradient(:,:),a%SmoothGradient(:,:))
      endif

!       call a%FilePostpr%postpr(a%level(:,1),'levelpostproj',a%istep,a%ctime,a%Mesh)




   elseif (itask == 2) then
      if (a%kfl_ForceEulerianAdvection == 0) then
         !Cut Elements deallocation
         call a%CutMesh%deallocCutElement
         !Recomputations of Cut Elements
         call lev_CutElementsAndListByLayers(a)
      endif

   endif

end subroutine


subroutine lev_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
   use typre
   use Mod_Mesh
   use Mod_Levelset
   use Mod_Advector
   implicit none
   class(LevelSetProblem), target :: a
   type(Advector) :: Advect
   integer(ip) :: itask

   integer(ip) :: ndime,icomp
   real(rp), pointer :: auxarray(:,:) => NULL()

   interface
      subroutine lev_CutElementsAndListByLayers(a)
      use Mod_LevelSet
      implicit none
      class(LevelSetProblem) :: a
      end subroutine
   end interface

   if (itask == 1) then

       call a%FilePostpr%postpr(a%level(:,1),'levelpreproj',a%istep,a%ctime,a%Mesh)

      !If we force no ALE, then we do not need to project, since the levelset is Eulerian
      if (a%kfl_ForceEulerianAdvection == 0) then
         do icomp = 1,size(a%level,2)
            auxarray(1:1,1:size(a%level,1)) => a%level(:,icomp)
            call Advect%Advect(1_ip,auxarray,auxarray)
         enddo
      endif
      if (a%kfl_SmoothGradient == 1) then
         call Advect%Advect(size(a%SmoothGradient,1),a%SmoothGradient(:,:),a%SmoothGradient(:,:))
      endif

       call a%FilePostpr%postpr(a%level(:,1),'levelpostproj',a%istep,a%ctime,a%Mesh)




   elseif (itask == 2) then
      if (a%kfl_ForceEulerianAdvection == 0) then
         !Cut Elements deallocation
         call a%CutMesh%deallocCutElement
         !Recomputations of Cut Elements
         call lev_CutElementsAndListByLayers(a)
      endif

   endif

end subroutine

subroutine lev_begste(a)
   use typre
   use Mod_LevelSet
   use def_parame
   use Mod_TimeIntegrator
   use Mod_CutMesh
   implicit none
   class(LevelSetProblem) :: a

   integer(ip) :: ipoin,npoin,icomp,ndime

   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy
   logical :: isALE

  interface
      subroutine lev_coupling_outerr(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_CutElementsAndListByLayers(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_ComputeCuts(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

   end interface



   if(a%kfl_timei==1) then

      !call a%FilePostpr%postpr(a%level(:,1),'lev_PreBegste',a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%level(:,2),'lev2',a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%level(:,3),'lev3',a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%level(:,4),'lev4',a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%level(:,5),'lev5',a%istep,a%ctime,a%Mesh)
      !write(*,*) 'sizse the level: ', size(a%level,2)
      !write(*,*) 'ncomp: ', a%ncomp

      if(a%changeTimeIntegrator==1)then
         a%kfl_tsche_1st_current='BDF1 '
         call Integrator%Init(a%kfl_tsche_1st_current)
         call Integrator%GetNumberOfTimeSteps(nsteps)
         call Integrator%GetAccuracy(Accuracy)
      elseif(a%changeTimeIntegrator==0)then
         a%kfl_tsche_1st_current=a%kfl_tsche_1st_datafile
         call Integrator%Init(a%kfl_tsche_1st_current)
         call Integrator%GetNumberOfTimeSteps(nsteps)
         call Integrator%GetAccuracy(Accuracy)
      end if



      if (Accuracy == 1) then
         a%level(:,2) = a%level(:,3)
      !Higher order schemes
      elseif (Accuracy <= 3 .and. a%ncomp >= Accuracy + 2) then
         !call GetTimeProjection(Accuracy,1,npoin,a%ncomp,a%level,a%level(:,2))
         a%level(:,2) = a%level(:,3)
       !Default
      else
         a%level(:,2) = a%level(:,3)
      endif

      !Set the first one
      a%level(:,1) = a%level(:,2)

      !call a%FilePostpr%postpr(a%level(:,1),'lev_PostBegste',a%istep,a%ctime,a%Mesh)

   end if


   !Compute the cuts
   call lev_ComputeCuts(a)



end subroutine

subroutine lev_ComputeCuts(a)
    use typre
   use Mod_LevelSet
   use def_parame
   use Mod_TimeIntegrator
   use Mod_CutMesh
   implicit none
   class(LevelSetProblem) :: a


   integer(ip) :: ipoin,icomp

   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy
   logical :: isALE

  interface
      subroutine lev_coupling_outerr(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_CutElementsAndListByLayers(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

   end interface


   !Check if the coupling between modules is done
   call lev_coupling_outerr(a)

   !If we have fixed-mesh ALE, then we don't need to recompute the cuts, since the level set is not advected
   call a%Mesh%GetALE(isALE)
   if (isALE .and. a%kfl_ForceEulerianAdvection == 0) then
      return
   endif

   !First we deallocate the previous one
   call a%CutMesh%deallocCutElement

   !now we recompute the cuts
   call lev_CutElementsAndListByLayers(a)





end subroutine

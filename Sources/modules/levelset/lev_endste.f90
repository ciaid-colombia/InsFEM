subroutine lev_endste(a,itask)
   use typre
   use Mod_Memor
   use Mod_LevelSet
   use def_parame
   use Mod_TimeIntegrator
   use Mod_CutMesh
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask
   type(TimeIntegratorDt1) :: Integrator
   character(5) :: TimeScheme
   integer(ip) :: ielem,nelem,icomp,nsteps

   interface

      subroutine lev_CutElementsAndListByLayers(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_ReinitLevel(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

      subroutine lev_EndsteElmope(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

   end interface

   !Pre-CrankNicolson
   if (itask == 1) then


   !Post-CrankNicolson
   elseif (itask == 2) then

      !Update a%velocity and a%pressure
      if (a%ncomp > 3) then
          do icomp = a%ncomp, 4,-1
            a%level(:,icomp) = a%level(:,icomp-1)
         enddo
      endif

      a%level(:,3) = a%level(:,1)  ! Vn = Vn+1

      if(a%istep/a%NstepsReinitLevel==a%nReinit)then
         call a%CutMesh%deallocCutElement
         call lev_CutElementsAndListByLayers(a)

         a%nReinit = a%nReinit +1
         call lev_ReinitLevel(a)
         call a%FilePostpr%postpr(a%level(:,3),'LEVELR',a%istep,a%ctime,a%Mesh)

         !change of the time integrator order to do following time step
         if(a%kfl_tsche_1st_current=='BDF2 ')then
            a%changeTimeIntegrator=1
         end if
      end if

      call lev_EndsteElmope(a)

      if(a%changeTimeIntegrator==1)then
         !Return to the original time integration order
         a%changeTimeIntegrator=0
      end if





   end if

end subroutine

subroutine tem_begste(a)
   use typre
   use Mod_Temperature
   use def_parame
   use Mod_TimeIntegrator      
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: npoin,icomp
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: nsteps, Accuracy
   
  interface
      subroutine tem_coupling_outerr(a)
         use Mod_Temperature
         implicit none
         class(TemperatureProblem) :: a
      end subroutine
   
   end interface
   
   !!Check if the coupling between modules is done
   call tem_coupling_outerr(a)
   
   call a%Mesh%GetNpoin(npoin)
   !Only for Burgers equation which is non-linear
   !Other non-linearities should be added too
   if(a%kfl_timei==1 .and. a%kfl_advec == 9) then
      call Integrator%Init(a%kfl_tsche_1st_current)
      call Integrator%GetNumberOfTimeSteps(nsteps)
      call Integrator%GetAccuracy(Accuracy)
      
      if (Accuracy == 1) then
         a%tempe(:,2) = a%tempe(:,3)
      !Higher order schemes
      elseif (Accuracy <= 3 .and. a%ncomp >= Accuracy + 2) then
         call GetTimeProjection(Accuracy,1_ip,npoin,a%ncomp,a%tempe,a%tempe(:,2))          
      !Default
      else
         a%tempe(:,2) = a%tempe(:,3)
      end if
      a%tempe(:,1) = a%tempe(:,2)
   end if

end subroutine
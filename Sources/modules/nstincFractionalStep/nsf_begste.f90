subroutine nsf_begste(a)
   use typre
   use Mod_NavierStokes
   use Mod_NSFractionalStep
   use Mod_NavierStokes
   implicit none
   class(NSFractionalStepProblem) :: a

   interface
      subroutine nsi_begste(a)
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
      
      subroutine nsm_EndElmope(NSProblem,task)
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem), target :: NSProblem
         character(6) :: task
      end subroutine
      
      subroutine nsf_HydrostaticPressureNew(a)
         use Mod_NSFractionalStep
         implicit none
         class(NSFractionalStepProblem) :: a
       end subroutine  
   end interface
   
   !First do what needs to be done usually in NavierStokes
   call nsi_begste(a)

   !Now, if this is the first time step, and there is a gravity field, 
   !we need to generate a hydrostatic pressure field and then compute the residual
   !Do not do this if this is a restart problem
   !This is not done in turnon because Nstokes problem might be delayed

   !if (a%istep == 1_ip .and. (a%grnor /= 0.0_rp .or. a%kfl_cotem /= 0) .and. a%kfl_restar == 0_ip) then
   if (a%istep == 1_ip  .and. a%kfl_rstar == 0_ip) then
      !Compute an hydrostatic Pressure Field
      
!       if(a%kfl_colev==1 .and. a%kfl_fsurf==1)then
!          write(*,*) 'No hydrostatic pressure initialization'
!       else
      call nsf_HydrostaticPressureNew(a)
!       end if
      if(a%kfl_colev==1)then
         call a%FilePostpr%postpr(a%press(:,1),'P_Hydro',a%istep,a%ctime,a%Mesh)       
      endif

      
      if (a%grnor /= 0.0_rp .or. a%kfl_cotem /= 0) then
         !Compute the residual projection with the newly created hydrostatic pressure field
         call a%EndElmope('Endite')
      endif
   endif
   
end subroutine

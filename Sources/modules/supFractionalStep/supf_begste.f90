subroutine supf_begste(a)
   use typre
   use Mod_ThreeField
   use Mod_SUPFractionalStep
   implicit none
   class(SUPFractionalStepProblem) :: a

   interface
      subroutine sup_begste(a)
         use Mod_ThreeField
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine
      
!       subroutine nsm_EndElmope(NSProblem,task)
!          use Mod_NavierStokes
!          implicit none
!          class(NavierStokesProblem), target :: NSProblem
!          character(6) :: task
!       end subroutine
!       
!       subroutine nsf_HydrostaticPressureNew(a)
!          use Mod_NSFractionalStep
!          implicit none
!          class(NSFractionalStepProblem) :: a
!        end subroutine  
   end interface
   
   !First do what needs to be done usually in NavierStokes
   call sup_begste(a)

   !Now, if this is the first time step, and there is a gravity field, 
   !we need to generate a hydrostatic pressure field and then compute the residual
   !Do not do this if this is a restart problem
   !This is not done in turnon because Nstokes problem might be delayed

   !if (a%istep == 1_ip .and. (a%grnor /= 0.0_rp .or. a%kfl_cotem /= 0) .and. a%kfl_restar == 0_ip) then
!    if (a%istep == 1_ip  .and. a%kfl_restar == 0_ip) then
!       !Compute an hydrostatic Pressure Field
!       call nsf_HydrostaticPressureNew(a)
!    
!       if (a%grnor /= 0.0_rp .or. a%kfl_cotem /= 0) then
!          !Compute the residual projection with the newly created hydrostatic pressure field
!          call a%EndElmope('Endite')
!       endif
!    endif
   
end subroutine
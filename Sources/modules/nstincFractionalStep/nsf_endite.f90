subroutine nsf_endite(a,itask)
   !-----------------------------------------------------------------------
   ! NAME 
   !    nsi_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   
   use def_parame
   use Mod_NSFractionalStep
   implicit none
   class(NSFractionalStepProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ndime
   
   interface
      subroutine nsm_EndElmope(NSProblem,task)
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem), target :: NSProblem
         character(6) :: task
      end subroutine
   end interface

   call a%Timer%Endite%Tic

   select case(itask)

   case(1)

       !Update the subscale and/or compute residual projection
       !call a%EnditeElmope
       call a%EndElmope('Endite')

   case(2)

       !Compute convergence residual of the external iteration
       call a%Cvgunk(two)

       !Update u(n,i-1,*) <-- u(n,i,*)
       a%veloc(:,:,2) = a%veloc(:,:,1)
       ! Update p(n,i-1,*) <-- p(n,i,*)
       a%press(:,2) = a%press(:,1)

   case(4)

       call a%Cvgunk(4)

       !Update coupling values
       a%veloc_cp(:,:) = a%veloc(:,:,1)
       a%press_cp(:) = a%press(:,1)

   end select 
   
   call a%Timer%Endite%Toc

end subroutine nsf_endite

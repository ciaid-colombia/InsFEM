subroutine sldale_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* alemov/sldale_endite
   ! NAME 
   !    sld_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem),target :: a
   type(TimeIntegratorDt1) :: Integrator
   integer(ip) :: itask,ierr,npoin,nsteps
   integer(ip) :: ndime,idime
   logical :: AreFoldedALE 
   real(rp), allocatable :: valRHS(:,:)
   real(rp)     :: LHSdtinv,w
   real(rp), pointer     :: v_i2(:,:) => NULL(),v_i1(:,:) => NULL(),v_i(:,:) => NULL()
   
   select case(itask)

   case(1)

       a%Displacement(:,:,1) = a%unkno(:,:)

   case(2)

   case(4)

       call a%Mesh%SetALE(1_ip)

       call php_SetTimeIntegrator(a,Integrator,LHSdtinv,nsteps)
       call a%Mesh%GetNdime(ndime)
       call a%Mesh%GetNpoin(npoin)

       call a%Memor%alloc(ndime,npoin,valRHS,'valRHS','sldale_solite')

       call Integrator%GetRHS(ndime*npoin,a%Displacement(:,:,2:nsteps),valRHS)
       
       a%Velocity(:,:,1) = LHSdtinv*a%Displacement(:,:,1) - valRHS*a%dtinv

       call a%Memor%dealloc(ndime,npoin,valRHS,'valRHS','sldale_solite')


       !Remesh/Project only when folded elements
       if (a%kfl_RemeshingCriteria == 0) then
           !Check if the computed displacement leads to element folding
           call a%Mesh%SetDisplacements(a%Displacement)
           call a%Mesh%ComputeCheckFoldedALE
           call a%Mesh%GetCheckFoldedALE(AreFoldedALE)
           call a%Mesh%SetDisplacements(a%Displacement)

           if (AreFoldedALE .eqv. .true.) then
               a%DoRemesh = 1
           else
               a%DoRemesh = 0
           endif

       endif

       if (a%DoRemesh == 0) then

           !Dealloc and recompute ExtnorLpoty and Vmass
           call a%Mesh%DeallocExnorLpoty
           call a%Mesh%DeallocVmass

           !Recompute
           call a%Mesh%ComputeVmass
           call a%Mesh%ExtnorLpoty

       else
           if (a%MPIrank == a%MPIroot) then
               write(*,*) '****: sldale_solite: mesh overlap, run will stop;'
               write(*,*) '    : attempting to print final mesh displacement for debugging'
           endif
           call a%FilePostpr%postpr(a%Displacement(:,:,1),'ALEdispl_final',a%istep,a%ctime,a%Mesh)
           call a%FilePostpr%postpr(a%Velocity(:,:,1),'ALEveloc_final',a%istep,a%ctime,a%Mesh)
           call flush
           call MPI_Barrier(a%MPIcomm, ierr)

       endif

   end select  

end subroutine sldale_endite

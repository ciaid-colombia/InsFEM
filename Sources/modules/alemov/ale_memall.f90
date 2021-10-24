subroutine ale_memall(a)
   !-----------------------------------------------------------------------
   !> This routine allocates the arrays of displacements and velocities, as 
   !! well as the arrays of boundary conditions. The amount of allocated 
   !! components depends on the number of steps of the integration scheme. 
   !-----------------------------------------------------------------------
   use typre 
   use Mod_TimeIntegrator
   use Mod_Alemov
   use Mod_sldAlemov
   implicit none

   type(TimeIntegratorDt1) :: Integrator
   class(AlemovProblem) :: a
   integer(ip) :: npoin,nsteps,oldnpoin

   call a%Mesh%GetNpoin(npoin)
   
   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = nsteps
   
   call a%Memor%alloc(a%ndofbc,npoin,a%ncomp,a%Displacement,'Displacement','ale_memall')
   call a%Memor%alloc(a%ndofbc,npoin,a%ncomp,a%Velocity,'Velocity','ale_memall')
   call a%Memor%alloc(a%ndofbc,npoin,a%kfl_fixno0,'kfl_fixno0','ale_memall')

   select type (a)
   type is (sldAlemovProblem)
       call a%Memor%alloc(a%ndofbc,npoin,2,a%bdisp,'bdisp','ale_memall')
       call a%Memor%alloc(a%ndofbc,npoin,3,a%bres,'bres','ale_memall')
   end select

end subroutine ale_memall

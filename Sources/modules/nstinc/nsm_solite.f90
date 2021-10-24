subroutine nsm_solite(a,itask)
   !Specific operations for solite
   use typre
   use def_parame
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: itask
   
   integer(ip),save :: kfl_advec_old,kfl_timei_old,kfl_cotem_old
   real(rp),save    :: dtinv_old
   
   interface
      subroutine nsm_rotunk(a,itask)
         use typre
         import NavierStokesProblem
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip)                :: itask
      end subroutine
   end interface
   
   
   if (itask == 1) then
      
      if (a%kfl_TurbulentBodyForces == 1) then
         call a%TBF%ResetRandomSeed(a%istep)
      endif
      
      !PETSc usually departs from unkno = 0
      !write(*,*) 'nsi_begite: updunk not working'
      !if(kfl_algor==1) then      ! Monolithic
      !   call nsi_updunk(itask_updunk_m_vp_ig)
      !endif
   
      !If initial solution is Stokes: save original values
      if(a%kfl_inist==1) then
         kfl_advec_old = a%kfl_advec
         kfl_timei_old = a%kfl_timei
         kfl_cotem_old = a%kfl_cotem
         dtinv_old     = a%dtinv
         a%kfl_advec = 0
         a%kfl_timei = 0
         a%kfl_cotem = 0
         a%dtinv     = 0.0_rp
         if(a%kfl_rstar==1) write(*,*) 'Warning: This is a restart run with initial stokes ON for NSTINC, case might not converge'
      end if
   
   elseif (itask == 2) then
      !Tranform from local to global systems 
      call nsm_rotunk(a,two)
      
      !If initial solution is Stokes: recover original values
      if(a%kfl_inist==1) then
         a%kfl_inist = 0
         a%kfl_advec = kfl_advec_old 
         a%kfl_timei = kfl_timei_old 
         a%kfl_cotem = kfl_cotem_old 
         a%dtinv     = dtinv_old
         if(a%kfl_rstar==1) write(*,*) 'Warning: This is a restart run with initial stokes ON for NSTINC, case might not converge'
      end if
  endif

end subroutine nsm_solite

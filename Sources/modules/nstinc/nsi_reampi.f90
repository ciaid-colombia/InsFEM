subroutine nsi_reampi(a)
   use MPI
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: ierr,imat
   
   !Communicate nsi_reaphy
   CALL MPI_BCAST(a%doRobin, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_advec, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_visco, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotem, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)   
   CALL MPI_BCAST(a%kfl_cotur, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%kfl_SwitchOff, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%kfl_colev, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)   
   CALL MPI_BCAST(a%kfl_fsurf, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%kfl_fsurfLapla, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_fsurfDirichlet, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_fsurfDirichletMethod, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%kfl_SurfaceTension, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)   
   CALL MPI_BCAST(a%kfl_EnrichElem, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)  
   
   CALL MPI_BCAST(a%alfa_robin, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%fcons, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%fvins, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   !Frame of ReferenceAcceleration
   CALL MPI_BCAST(a%kfl_FORAcceleration, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_FORFunty, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)  
   CALL MPI_BCAST(a%FORAcceleration, size(a%FORAcceleration), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%FORParam, size(a%FORParam), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_FORAxesRotation, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%FORAxesAngularVeloc, size(a%FORAxesAngularVeloc), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_CoriolisForce, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%CoriolisW, 3, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)  
   
   
   !Physical properties   
   CALL MPI_BCAST(a%nmat, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)     
   do imat=1,a%nmat
   
      CALL MPI_BCAST(a%MatProp(imat)%densi, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%MatProp(imat)%visco, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)   
      CALL MPI_BCAST(a%MatProp(imat)%lawvi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)  
      
      if(a%MatProp(imat)%lawvi/=0)then
         CALL MPI_BCAST(a%MatProp(imat)%LawViParam, size(a%MatProp(imat)%LawViParam), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      end if      
   
   end do 
   
   CALL MPI_BCAST(a%grnor, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%gravi, size(a%gravi), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%boube, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%boutr, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%bougr, size(a%bougr), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%turbu, size(a%turbu), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   !Communicate nsi_reanut
   CALL MPI_BCAST(a%kfl_repro, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_repro_SkipFE, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_adapsgs, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_shock, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_wtemp, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_wlapl, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_trasg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_nolsg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_nolsgNewtonRaphson, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_bousg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   CALL MPI_BCAST(a%mtrit, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_tacsg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_stabm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_tausm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_hdifumin, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_penal, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_postBtract, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_StabilizeFreeSurface, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   CALL MPI_BCAST(a%staco, size(a%staco), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%bstco, size(a%bstco), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%shock, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%tosgs, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%relsg, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%penal, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%subrelax, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%StabilizeFreeSurface_Param, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
  
   !Communicate nsi_reaous
   CALL MPI_BCAST(a%kfl_dispa, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%Avg1DIdime, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   !Forces and moments
   CALL MPI_BCAST(a%kfl_outfm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_computeTractions, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if(a%kfl_outfm==1)then
      CALL MPI_BCAST(a%adimf, size(a%adimf), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%adimm, size(a%adimm), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%origm, size(a%origm), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr) 
   end if
   
   !For fixing CN
   if (a%kfl_tsche_1st_datafile == 'CNOBS') a%kfl_tsche_1st_datafile = 'CN   '
   
   CALL MPI_BCAST(a%StatisticsStartTime, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)

   !Elastic Dirichlet boundary conditions
   CALL MPI_BCAST(a%EB_Stiffness, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EB_ShearStiffness, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EB_Density, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EB_Damping, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EB_NDelaySteps, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Turbulent Inlet boundary conditions
   call MPI_BCAST(a%kfl_TurbulentInletBoundaryConditions, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_TurbulentInletBoundaryConditions == 1) then
      call a%TIBC%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%TIBC%SetMemor(a%Memor)
      call a%TIBC%ReadDataMPI
   endif
   
   call MPI_BCAST(a%kfl_ExitBodyForces, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_ExitBodyForces == 1) then
      call a%EBF%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%EBF%ReadDataMPI
      call a%EBF%SetDensity(a%MatProp(1)%densi)
   end if
   
   call MPI_BCAST(a%kfl_TurbulentBodyForces, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_TurbulentBodyForces == 1) then
      call a%TBF%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%TBF%ReadDataMPI
      call a%TBF%SetDensity(a%MatProp(1)%densi)
   endif
   
   CALL MPI_BCAST(a%kfl_ExtraInitialViscosity, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EIV_Visco, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%EIV_nsteps, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_PorousMedia, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_PorousMedia == 1) then
      call a%PME%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%PME%ReadDataMPI
   endif
   
   CALL MPI_BCAST(a%kfl_Plasma, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%kfl_Plasma == 1) then
      call a%PAC%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%PAC%ReadDataMPI
   endif
   
   call MPI_BCAST(a%kfl_RVE, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   !Domain decomposition
   CALL MPI_BCAST(a%kfl_doAitken, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%relax, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%relax_max, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)

   call a%NstincSpecificReampi

end subroutine

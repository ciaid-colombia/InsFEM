subroutine nsf_solite(a)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_NSFractionalStep
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   use def_parame
   use Mod_int2str
   use Mod_nsi_BouwalStats
   implicit none
   class(NSFractionalStepProblem) :: a
   
   integer(ip) :: ndime,kfl_perio, kfl_HangingNodes,iexternal
   
   interface 
      subroutine nsf_elmope_1st(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine

      subroutine nsf_elmope_1stnew(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
     
      subroutine nsf_elmope_2ndnew(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
      
      subroutine nsf_bouope_2nd(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
      
      subroutine nsf_elmope_3rd(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
      
      subroutine nsf_bouope_3rd(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
      
      subroutine nsf_bouope_1st(a)
         use Mod_NSFractionalStep
         use typre
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine
      
      subroutine nsm_rotunk(a,itask)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip)                :: itask
      end subroutine
      
      subroutine nsi_cvgunk(a,itask)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip), intent(in)    :: itask
      end subroutine
      
   end interface
   
   !This subroutine solves an iteration of the Navier-Stokes equations,
   !FRACTIONAL STEP FORM
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 1st step : Compute intermediate velocity   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ExternalLoop: do iexternal = 1,a%FractionalExternalIterations  
      
      a%kfl_goite = 1
      a%itera = 0
      
      do while(a%kfl_goite==1)
         
         if (a%kfl_StabilizeFreeSurface == 1) then
            call a%EndElmope('Endite')
         endif
         
         if (a%kfl_TurbulentBodyForces == 1) then
            call a%TBF%ResetRandomSeed(a%istep)
         endif
         
         !If wall law update statistiscs
         call nsi_BouwalStats(a,0)

         
         !Update inner iteration counter and write headings in the solver file.
         a%itera = a%itera + 1
         if (a%MPIrank == a%MPIroot) then
            if(a%itera==1) write(a%lun_solve,100) a%istep
            write(a%lun_solve,101) a%itera
         endif
         
         call a%Timer%BuildLinearSystem%Tic
         !Compute Linear System
         call a%LinearSystem%ToZero
         call nsf_elmope_1stnew(a)
         call nsf_bouope_1st(a)
         !Periodic Boundary Conditions
         call a%Mesh%GetPerio(kfl_perio)
         call a%Mesh%GetNdime(ndime)
         if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(ndime,a%LinearSystem,a%Memor)
         !Hanging nodes
         call a%Mesh%GetHanging(kfl_HangingNodes)
         if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(ndime,a%LinearSystem,a%Memor)
         
         
         call a%Timer%BuildLinearSystem%Toc
         
         call a%Timer%SolveLinearSystem%Tic
         
         call a%LinearSystem%Solve(a%unkno)
         
         !Ghost communicates
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%unkno)
         
         !HangingNodes
         call a%Mesh%GetHanging(kfl_HangingNodes)
         if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(ndime,a%unkno)
         
         !Periodic boundary conditions
         call a%Mesh%GetPerio(kfl_perio)
         if (kfl_perio == 1) call a%Mesh%MasterToSlave(ndime,a%unkno)
         
         call nsm_rotunk(a,two)
         
         call a%Timer%SolveLinearSystem%Toc
   
         !Compute convergence residual of the internal a%iteration 
         call nsi_cvgunk(a,one)

         !Update unknowns
         a%veloc(:,:,1) = a%unkno(1:ndime,:)*a%subrelax + a%veloc(:,:,1)*(1.0_rp-a%subrelax)
         
         !Update the subscale and/or compute residual projection
         !At this point, this is done at the end of the step
         !call a%EndElmope('Endite')
         if (a%kfl_StabilizeFreeSurface == 1) then
            call a%EndElmope('Endite')
         endif
         
         call a%EndElmope('Endste')
         
         
         !call a%FilePostpr%postpr(a%veloc(:,:,1),'VELOCINTERMEDIATE',a%istep,a%ctime,a%Mesh)
      end do
      
      if(a%kfl_tsche_1st_current == 'CN   ') then
         a%veloc(:,:,2) = a%veloc(:,:,1) ! int. u_n+1/2
         a%veloc(:,:,1) = 2.0_rp*a%veloc(:,:,1) - a%veloc(:,:,3)
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 2nd step : Compute pressure at n+1         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call a%Timer%BuildLinearSystem%Tic
      
      if (a%MPIrank == a%MPIroot) write(a%lun_solve,102)
      
      !Initialize values
      call a%LinearSystemP%ToZero
      
      call nsf_elmope_2ndnew(a)
      call nsf_bouope_2nd(a)         !for openflow, but it doesn't work
      !Periodic Boundary Conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(1_ip,a%LinearSystemP,a%Memor)
      !Hanging nodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(1_ip,a%LinearSystemP,a%Memor)
      
      call a%Timer%BuildLinearSystem%Toc
      
      call a%Timer%SolveLinearSystem%Tic
      call a%LinearSystemP%Solve(a%unknoP)
      
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%unknoP)
      
      !HangingNodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(1,a%unknoP)
      
      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(1,a%unknoP)
      
      !Update unknowns
      a%press(:,2) = a%press(:,1)
      a%press(:,1) = a%unknoP(1,:)
      call a%Timer%SolveLinearSystem%Toc

      !call a%FilePostpr%postpr(a%press(:,1),'PRESSINT'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%press(:,2),'PRESSPRE'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%press(:,1)-a%press(:,2),'pressDIF'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 3rd step : Compute end-of-step velocity    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call a%Timer%BuildLinearSystem%Tic
      
      if (a%MPIrank == a%MPIroot) write(a%lun_solve,103) 
      
      call a%LinearSystem%ToZero
      
      call nsf_elmope_3rd(a)
      call nsf_bouope_3rd(a)        !For open flow, but it doesn't work
      !Periodic Boundary Conditions
      call a%Mesh%GetPerio(kfl_perio)
      call a%Mesh%GetNdime(ndime)
      if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(ndime,a%LinearSystem,a%Memor)
      !Hanging nodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(ndime,a%LinearSystem,a%Memor)
      call a%Timer%BuildLinearSystem%Toc
   
      call a%Timer%SolveLinearSystem%Tic
      call a%LinearSystem%Solve(a%unkno)
      
      !GhostCommunicates
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%unkno)
      
      !HangingNodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(ndime,a%unkno)
      
      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(ndime,a%unkno)
      call nsm_rotunk(a,two)
      
      call a%Timer%SolveLinearSystem%Toc

      !Update unknowns
      
      a%veloc(:,:,1) = a%unkno(:,:)
      
      !call a%FilePostpr%postpr(a%veloc(:,:,1),'VELOCINT'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%veloc(:,:,2),'VELOCPRE'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      !call a%FilePostpr%postpr(a%veloc(:,:,1)-a%veloc(:,:,2),'VELOCDIF'//trim(adjustl(int2str(iexternal))),a%istep,a%ctime,a%Mesh)
      
      !Check residuals
      !call nsi_cvgunk(a,two)
      !write(*,*) 'Residuals: ', a%resid, a%resip
      
      !a%veloc(:,:,2) = a%veloc(:,:,1)
      
   enddo ExternalLoop      
   
   !Formats. 
   100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
   101 format('     INTERMEDIATE VELOCITY SOLVE, ITERATION NUMBER: ',i5)
   102 format(/,'   PRESSURE CORRECTION STEP')
   103 format(/,'   END-OF-STEP VELOCITY')

end subroutine nsf_solite

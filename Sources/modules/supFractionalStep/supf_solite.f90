subroutine supf_solite(a)
!    use Mod_PhysicalProblem
!    use Mod_NavierStokes
!    use Mod_ThreeField
   use Mod_SUPFractionalStep  
   use def_parame
   implicit none
   class(SUPFractionalStepProblem) :: a
   
   integer(ip) :: ndime,ntens
   
   interface 
   
      !--------------------------------------------------
      !Elmopes
      !Momentum
      subroutine supf_elmope_1st(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine

   
      subroutine supf_elmope_4rd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine
      
      subroutine supf_elmope_4rdY(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine       
      
      !Constituive
      subroutine supf_elmope_2nd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine
      
      subroutine supf_elmope_5th(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine
      
      !Continuity
      subroutine supf_elmope_3rd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine
      
      subroutine supf_elmope_hydro(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine      
      
      !-------------------------------------------
      !Bouopes
      !Momentum
      subroutine supf_bouope_1st(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine      

      subroutine supf_bouope_1stY(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine       
      
      subroutine supf_bouope_4rd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine

      !Constituive
      subroutine supf_bouope_2nd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine      
      
      subroutine supf_bouope_5th(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine
      
      !continuity
      subroutine supf_bouope_3rd(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine      
      
      
      
      !Endites
      subroutine supf_EndElmope_u(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a     
      end subroutine      
      
      subroutine supf_EndElmope_s(a)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a        
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
      
      subroutine sup_cvgunk(a,itask)
         use typre
         use Mod_ThreeField
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip), intent(in)    :: itask
      end subroutine
      
   
   
      subroutine supf_cvgunkU(a,itask)
         use typre
         import SUPFractionalStepProblem   
         implicit none
         class(SUPFractionalStepProblem) :: a
         integer(ip), intent(in)     :: itask
      end subroutine  

      subroutine supf_cvgunkS(a,itask)
         use typre
         import SUPFractionalStepProblem   
         implicit none
         class(SUPFractionalStepProblem) :: a
         integer(ip), intent(in)     :: itask
      end subroutine    
      
      subroutine supf_cvgunkUY(a,itask)
         use typre
         import SUPFractionalStepProblem 
         implicit none
         class(SUPFractionalStepProblem) :: a
         integer(ip), intent(in)     :: itask
      end subroutine   
      
      
   end interface
    
   !This subroutine solves an iteration of the Viscoelastic Navier-Stokes problem,
   !FRACTIONAL STEP FORM
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 0 step : Stokes problem to initialize the pressure !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(a%kfl_inist==1 .and. a%istep==1)then
      call a%Timer%BuildLinearSystem%Tic
      !Compute Linear System
      call a%LinearSystemST%ToZero
      call supf_elmope_hydro(a)
  
      call a%Mesh%GetNdime(ndime)
      call a%Timer%BuildLinearSystem%Toc
        
      !Solve the system
      call a%Timer%SolveLinearSystem%Tic      
      call a%LinearSystemST%Solve(a%unknoST)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime+1,a%unknoST)      
      call nsm_rotunk(a,two)      
      call a%Timer%SolveLinearSystem%Toc
   
         
      !Update unknowns
      a%veloc(:,:,1) = a%unknoST(1:ndime,:) 
      !Intermediate velocity
      a%veloc(:,:,2) = a%veloc(:,:,1)
      a%veloc(:,:,3) = a%veloc(:,:,2)
      
      !Update unknowns
      a%press(:,1) = a%unknoST(ndime+1,:) 
      !Intermediate velocity
      a%press(:,2) = a%press(:,1)
      a%press(:,3) = a%press(:,2)
         
      call a%Timer%Endite%Tic

      call a%Timer%Endite%Toc
      call a%FilePostpr%postpr(a%press(:,1),'P_Hydro',a%istep,a%ctime,a%Mesh) 
      call a%FilePostpr%postpr(a%veloc(:,:,1),'V_Hydro',a%istep,a%ctime,a%Mesh)
      
   end if
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 1st step : Compute intermediate velocity   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   do while(a%kfl_goite==1)
      
         !Update inner iteration counter and write headings in the solver file.
         a%itera = a%itera + 1
         if (a%MPIrank == a%MPIroot) then
            if(a%itera==1) write(a%lun_solve,100) a%istep
            write(a%lun_solve,101) a%itera
         endif
         
         call a%Timer%BuildLinearSystem%Tic
         !Compute Linear System
         call a%LinearSystem%ToZero
         call supf_elmope_1st(a)
         
         call supf_bouope_1st(a)     
         call a%Mesh%GetNdime(ndime)
         call a%Timer%BuildLinearSystem%Toc
         
         !Solve the system
         call a%Timer%SolveLinearSystem%Tic      
         call a%LinearSystem%Solve(a%unkno)
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%unkno)      
         call nsm_rotunk(a,two)      
         call a%Timer%SolveLinearSystem%Toc
   
         !Compute convergence residual of the internal a%iteration in the velocity case 
         call supf_cvgunkU(a,one)
         
         !Update unknowns
         a%veloc(:,:,1) = a%unkno(:,:)*a%subrelax + a%veloc(:,:,1)*(1.0_rp-a%subrelax) 
         !Intermediate velocity
         a%veloc(:,:,2) = a%veloc(:,:,1)
         
         !Update the subscale and/or compute residual projection
         !At this point, this is done at the end of the step
         
         call a%Timer%Endite%Tic
         
         call supf_EndElmope_u(a)        
         
         call a%Timer%Endite%Toc 

   end do 
   
   if(a%kfl_colev==1)then
      call a%FilePostpr%postpr(a%veloc(:,:,1),'V_Inter',a%istep,a%ctime,a%Mesh)
   end if
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !2nd step : Compute intermediate stress   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   do while(a%kfl_goites==1)
      
      !Update inner iteration counter and write headings in the solver file.
      a%iteras = a%iteras + 1
      if (a%MPIrank == a%MPIroot) then
         if(a%iteras==1) write(a%lun_solve,100) a%istep
         write(a%lun_solve,102) a%iteras
      endif
      
      call a%Timer%BuildLinearSystem%Tic
      !Compute Linear System
      call a%LinearSystemS%ToZero
      call supf_elmope_2nd(a)
      call supf_bouope_2nd(a)
      call a%Mesh%GetNdime(ndime)
      call a%Timer%BuildLinearSystem%Toc
      
      !Solve the system
      call a%Timer%SolveLinearSystem%Tic      
      call a%LinearSystemS%Solve(a%unknoS)
      call a%Mesh%GetNdime(ndime)
      ntens=(ndime-1)*(ndime-1)+2
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ntens,a%unknoS)      
      call nsm_rotunk(a,two)      
      call a%Timer%SolveLinearSystem%Toc
 
      !Compute convergence residual of the internal a%iteration in the velocity case 
      call supf_cvgunkS(a,one)
       
      !Update unknowns
      a%sigma(:,:,1) = a%unknoS(:,:)*a%subrelax + a%sigma(:,:,1)*(1.0_rp-a%subrelax)
      !Intermediate Stress
      a%sigma(:,:,2) = a%sigma(:,:,1)
      
      !Update the subscale and/or compute residual projection
      !At this point, this is done at the end of the step
      
      call a%Timer%Endite%Tic 
               
      call supf_EndElmope_s(a)     
      
      call a%Timer%Endite%Toc      
      
   end do   
   
   if(a%kfl_colev==1)then
      call a%FilePostpr%postpr(a%sigma(:,:,1),'S_Inter',a%istep,a%ctime,a%Mesh)
   end if   
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3rd step : Compute pressure at n+1         !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   call a%Timer%BuildLinearSystem%Tic
   
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_solve,103) 
   endif   
  
   !Initialize values
   call a%LinearSystemC%ToZero
   call supf_elmope_3rd(a)
   call a%Timer%BuildLinearSystem%Toc
   
   call a%Timer%SolveLinearSystem%Tic
   call a%LinearSystemC%Solve(a%unknoC)
   call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%unknoC)   
   !Update unknowns
   a%press(:,1) = a%unknoC(1,:)
   call a%Timer%SolveLinearSystem%Toc
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 4rd step : Compute end-of-step velocity    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(a%kfl_tsche_1st_datafile /= 'BDF3 ')then   
      call a%Timer%BuildLinearSystem%Tic
      
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_solve,104) 
      endif
      
      call a%LinearSystem%ToZero
      
      call supf_elmope_4rd(a)
      call a%Mesh%GetNdime(ndime)
      call a%Timer%BuildLinearSystem%Toc
   
      call a%Timer%SolveLinearSystem%Tic
      call a%LinearSystem%Solve(a%unkno)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%unkno)
      call nsm_rotunk(a,two)
      
      call a%Timer%SolveLinearSystem%Toc

      !Update unknowns
      a%veloc(:,:,1) = a%unkno(:,:)
   elseif(a%kfl_tsche_1st_datafile == 'BDF3 ')then 
      do while(a%kfl_goiteY==1)
         
         !Update inner iteration counter and write headings in the solver file.
         a%iteraY = a%iteraY + 1
         if (a%MPIrank == a%MPIroot) then   
            if(a%iteraY==1) write(a%lun_solve,100) a%istep
            write(a%lun_solve,106) a%iteraY
         endif
         
         call a%Timer%BuildLinearSystem%Tic
         !Compute Linear System
         call a%LinearSystem%ToZero
         call supf_elmope_4rdY(a)
         call supf_bouope_1stY(a)    
         call a%Mesh%GetNdime(ndime)
         call a%Timer%BuildLinearSystem%Toc
         
         !Solve the system
         call a%Timer%SolveLinearSystem%Tic      
         call a%LinearSystem%Solve(a%unkno)
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%unkno)      
         call nsm_rotunk(a,two)      
         call a%Timer%SolveLinearSystem%Toc
   
         !Compute convergence residual of the internal a%iteration in the velocity case 
         call supf_cvgunkUY(a,one)
         
         !Update unknowns
         !a%veloc(:,:,1) = a%unkno(:,:)*a%subrelax + a%veloc(:,:,1)*(1.0_rp-a%subrelax) 
         a%veloc(:,:,1) = a%unkno(:,:)         
         
         !Update the subscale and/or compute residual projection
         !At this point, this is done at the end of the step
         
         call a%Timer%Endite%Tic  
         
         call supf_EndElmope_u(a)
         
         call a%Timer%Endite%Toc         
      
      end do    
   
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 5th step : Compute end-of-step stress    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call a%Timer%BuildLinearSystem%Tic
   
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_solve,105) 
   endif
   
   call a%LinearSystemS%ToZero
   
   call supf_elmope_5th(a)
   call a%Mesh%GetNdime(ndime)
   call a%Timer%BuildLinearSystem%Toc
   ntens=(ndime-1)*(ndime-1)+2
   call a%Timer%SolveLinearSystem%Tic
   call a%LinearSystemS%Solve(a%unknoS)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%ArrayCommunicator%GhostCommunicate(ntens,a%unknoS)
   call nsm_rotunk(a,two)
   
   call a%Timer%SolveLinearSystem%Toc

   !Update unknowns
   a%sigma(:,:,1) = a%unknoS(:,:)   
   
   !Formats. 
   100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
   101 format('     INTERMEDIATE VELOCITY SOLVE, ITERATION NUMBER: ',i5)
   102 format('     INTERMEDIATE STRESS   SOLVE, ITERATION NUMBER: ',i5)   
   103 format(/,'   PRESSURE CORRECTION STEP')
   104 format(/,'   END-OF-STEP VELOCITY')
   105 format(/,'   END-OF-STEP STRESS')
   106 format('     FINAL VELOCITY SOLVE, ITERATION NUMBER: ',i5)   
end subroutine supf_solite

subroutine nsi_output(a,itask)
   !-----------------------------------------------------------------------
   !****f* Nstinc/nsi_output
   ! NAME 
   !    nsi_output
   ! DESCRIPTION
   !    End of a NSTINC time step 
   !    itask = 0  When timemarching is true. There is output or post-process
   !               of results if required.
   !    itask = 1  When timemarching is false. Output and/or post-process of
   !               results is forced if they have not been written previously.
   ! USES
   !    output
   !    a%FilePostpr%postpr
   ! USED BY
   !    nsi_endste (itask=1)
   !    nsi_turnof (itask=2)
   !***
   !-----------------------------------------------------------------------
   use typre
   use Mod_Postpr
   use Mod_Stream
   use Mod_NavierStokes
   use Mod_nsi_Statistics
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip), save     :: dopost(49)
   real(rp), pointer     :: coord(:) => NULL()
   integer(ip)           :: itask,ndime,ipoin,npoin,nelem
   real(rp)              :: auxcoord(3),auxveloc(3)
   real(rp), allocatable :: CoriolisVeloc(:,:)
   real(rp), allocatable :: error(:)
   real(rp)              :: TotalEstimatedError

   select case(itask)
   !End of a time step.
   case(0)
      a%pos_alrea=0
      !Tracking of points.
      if(a%nptra > 0) then
         call a%PointTracking
      end if
      if(a%kfl_exacs /= 0 .and. a%kfl_timei /= 0 ) call a%SpecificExaerr
   !End of the run.
   case(1)
      if(a%kfl_exacs/=0 .and. a%kfl_timei == 0) then
         call a%SpecificExaerr         
      end if     
   end select
  
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)

   !Do the actual postprocess
   !Velocity
   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%veloc(:,:,1),'Velocity',a%istep,a%ctime,a%Mesh)
      
      if (a%kfl_CoriolisForce == 1) then
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         call a%Memor%alloc(ndime,npoin,CoriolisVeloc,'Coriolis_Velocity','nsi_output')
         
         auxcoord(1:3) = 0.0_rp
         do ipoin = 1,npoin
            call a%Mesh%GetPointCoord(ipoin,coord)
            auxcoord(1:ndime) = coord(1:ndime)
            call vecpro(a%CoriolisW,auxcoord,auxveloc,3)
            CoriolisVeloc(:,ipoin) = a%veloc(:,ipoin,1) + auxveloc(1:ndime)
         enddo
         call a%FilePostpr%postpr(CoriolisVeloc,'Coriolis_Velocity',a%istep,a%ctime,a%Mesh)
         call a%Memor%dealloc(ndime,npoin,CoriolisVeloc,'Coriolis_Velocity','nsi_output')
      endif
   end if
   
   !Pressure
   if (dopost(2) == 1) then
      call a%FilePostpr%postpr(a%press(:,1),'Pressure',a%istep,a%ctime,a%Mesh)
   end if
   
   !Streamlines
   if (dopost(3) == 1) then
      call a%Mesh%GetNdime(ndime)
      call stream(a%veloc,ndime,a%istep,a%ctime,a%Mesh,a%Memor,a%FilePostpr,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPISize)
   end if
   
   !non Newtonian case
   if (dopost(5) == 1) then
     call a%FilePostpr%postgp(a%viscarray,'NN-Viscosity',a%istep,a%ctime,a%Mesh)
   end if
   
   !TAusmoothing Newtonian case
   if (dopost(9) == 1 .and. a%kfl_Tausm >= 1) then
     call a%FilePostpr%postpr(a%Tausmo(1,:),'Timom',a%istep,a%ctime,a%Mesh)
     call a%FilePostpr%postpr(a%Tausmo(2,:),'Tidiv',a%istep,a%ctime,a%Mesh)
   end if
   
   !Divergence
   if (dopost(11) == 1) then
      call a%FilePostpr%postgp(a%divergence,'Divergence',a%istep,a%ctime,a%Mesh)
   end if
   
   !Vorticity
   if (dopost(10) == 1) then
      call a%FilePostpr%postpr(a%vorti(:,:),'Vorticity',a%istep,a%ctime,a%Mesh)
   end if
   
   !Subscales
   if (dopost(14) == 1 .and. a%kfl_tacsg == 1) then    
      call a%FilePostpr%postgp(a%vesgs,'Velocity_SGS',a%istep,a%ctime,a%Mesh)
      if (a%kfl_repro >= 2 .and. a%kfl_tacsg == 1) then
         call a%FilePostpr%postgp(a%vesgs2,'Velocity_SGS2',a%istep,a%ctime,a%Mesh)
      end if
      !if (a%kfl_bousg == 1) then
      !   call a%FilePostpr%postgp(a%bvesgs,'Boundary_SGS',a%istep,a%ctime,a%Mesh,'BOUNDA')
      !end if
      if (a%RefinerErrorEstimator == 'SUBSC') then
         call a%Mesh%GetNelem(nelem)
         call a%Memor%alloc(nelem,error,'error','nsi_output')
         call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
         call a%FilePostpr%postgp(error,'ErrorSGS',a%istep,a%ctime,a%Mesh)
         call a%Memor%dealloc(nelem,error,'error','nsi_output')
      end if
   end if

   if (dopost(26) == 1) then
      call a%FilePostpr%postgp(a%prsgs,'Pressure_SGS',a%istep,a%ctime,a%Mesh)
   end if
   
   !Dissipation
   if (dopost(15) == 1) then
      call a%FilePostpr%postpr(a%Dissipation,'DISSI_NS',a%istep,a%ctime,a%Mesh)
   end if
   
   !Averaged velocity along one dimension
   if (dopost(16) == 1) then
      call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%veloc(:,:,1),'Averaged_Velocity',a%istep,a%ctime,a%Mesh,a%Memor)
   endif
   
   !Averaged dissipation along one dimension
   if (dopost(17) == 1) then
      call a%FilePostpr%postAvg1D(a%Avg1DIdime,a%dissipation,'Averaged_Dissipation_NSI',a%istep,a%ctime,a%Mesh,a%Memor)
   endif
   
   !Residual
   if (dopost(18) == 1) then
      call a%FilePostpr%postgp(a%residualU,'Residual_U',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(a%residualP,'Residual_P',a%istep,a%ctime,a%Mesh)
      if (a%kfl_repro == 2) call a%FilePostpr%postgp(a%residualGraP2,'Residual_GraP2',a%istep,a%ctime,a%Mesh)
   endif
   
   !Repro
   if (dopost(19) == 1) then
      call a%Mesh%GetNdime(ndime)
      call a%FilePostpr%postpr(a%repro(1:ndime,:),'Proj_residual_U',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%repro(ndime+1,:),'Proj_residual_P',a%istep,a%ctime,a%Mesh)
      if (a%kfl_repro == 2) call a%FilePostpr%postpr(a%repro(ndime+1:2*ndime+1,:),'Proj_residual_GraP2',a%istep,a%ctime,a%Mesh) 
   endif
            
   !Gauss Point Dissipations
    if (dopost(20) == 1) then
       call a%FilePostpr%postgp(a%GPDissipation,'GPDissipation',a%istep,a%ctime,a%Mesh)
   endif

    if (dopost(21) == 1) then
        call a%FilePostpr%postpr(a%qfac(:),'QFACT',a%istep,a%ctime,a%Mesh)
   endif

   if (dopost(22) == 1) then
       call a%FilePostpr%postpr(a%btraction,'BTraction',a%istep,a%ctime,a%Mesh)
   end if
   
   if (dopost(24) == 1) then
      call a%FilePostpr%postpr(a%timav(:,:,1),'Time_Average_Velocity',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%timap(:,1),'Time_Average_Pressure',a%istep,a%ctime,a%Mesh)
      call nsi_FinalStatistics(a)
      call a%FilePostpr%postpr(a%rmsv,'RMS_Velocity',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%rmsp,'RMS_Pressure',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%turbi,'Turbulence_Intensity',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%TimaTractions(:,:,1),'Time_Average_Tractions',a%istep,a%ctime,a%Mesh)
   endif
   
   if (dopost(25) == 1) then
       call a%FilePostpr%postpr(a%ExternalForcesArray,'ExternalForces',a%istep,a%ctime,a%Mesh)
   endif
   
   if(dopost(27) == 1) then
      call a%FilePostpr%postpr(a%AnalyticalVelocity(:,:,1),'AnalyticalVelocity',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%AnalyticalPressure(:,1),'AnalyticalPressure',a%istep,a%ctime,a%Mesh)
   endif

   call a%NstincSpecificOutput(itask)

end subroutine nsi_output

subroutine nsi_outerr(a)
!------------------------------------------------------------------------
! NAME 
!    nsi_outerr
! DESCRIPTION
!    This routine checks if there are errors and warnings
!------------------------------------------------------------------------
   use typre
   use Mod_int2str
   use Mod_TimeIntegrator
   use Mod_NavierStokes
   implicit none
   integer(ip) :: ierro=0,iwarn=0

   class(NavierStokesProblem) :: a
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: Accuracy
   
   if (a%kfl_trasg /= 0 .and. 'ForceClosedRule' .eq. a%EndLoopQuadrature) then
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outpu,*)
         write(a%lun_outpu,*) 'NstincProblem: If tracking of subgrid scales is used, quadratures at the end of step and end of iteration must be default'
         write(a%lun_outpu,*) 'Setting EndLoopQuadrature to Default'
         write(a%lun_outpu,*)
      endif
      a%EndLoopQuadrature = 'DefaultRule'
   endif
  
   if(a%kfl_trasg/=0) then
      !Set the Number of components necessary for the arrays
      call Integrator%Init(a%kfl_tsche_1st_datafile)
      call Integrator%GetAccuracy(Accuracy)
      if(Accuracy==1.and.a%kfl_tacsg==2) then
         if (a%MPIrank == a%MPIroot) then
            write(a%lun_outph,101) &
               'TIME ACCURACY OF Usg CANNOT BE GREATER THAN THAT OF Uh.'
         endif
         a%kfl_tacsg=1
      end if
   end if
   
   if((a%kfl_incnd/=0) .and. (a%kfl_local/=0) ) then
         ierro=ierro+1
         if (a%MPIrank == a%MPIroot) then
            write(a%lun_outph,100)&
               'INITIal conditions not thought to be used with local systems'
         endif
   end if
  
   if(a%kfl_cotur/=0.and.a%kfl_advec==0) then
      ierro=ierro+1
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outph,101)&
               'TURBULENCE MODELING REQUIRES A CONVECTIVE TERM'
         write(a%lun_outph,101)&
               'TURBULENCE MODELING REQUIRES A CONVECTIVE TERM'
      endif
   end if
   
   if (a%npp_stepi(14) /= 0 .and. a%kfl_trasg == 0) then
      ierro = ierro+1
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outph,100)&
               'SUBSCALE OUTPUT BUT NO SUBSCALES TRACKING, NAVIER-STOKES'
         write(a%lun_outph,100)&
               'SUBSCALE OUTPUT BUT NOT SUBSCALES TRACKING, NAVER-STOKES'
      endif
   endif
   
   if (a%npp_stepi(19) /= 0 .and. a%kfl_repro == 0) then
      a%npp_stepi(19) = 0
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outph,101)&
               'Residual projection Output, but no repro is computed, Navier-Stokes'
         write(a%lun_outph,101)&
               'Residual projection Output, but no repro is computed, Navier-Stokes'
      endif
   endif
   
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if
   
   if (a%kfl_fsurf == 1 .and. a%nmat /= 2) then
      call runend('Free surface but the number of materials is not 2')
   endif
   
   if (a%kfl_fsurf == 1 .and. a%kfl_confi ==  1) then
      call runend('Free surface: pressure should not be fixed')
   endif

!Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine nsi_outerr

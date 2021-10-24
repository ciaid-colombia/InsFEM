subroutine nsi_memall(a)
   !-----------------------------------------------------------------------
   !****f* Nstinc/nsi_memall
   ! NAME 
   !    nsi_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the
   !    NS equations
   ! USES
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   type(TimeIntegratorDt1)    :: Integrator
   integer(ip) :: pnode,pnodb,nsteps,ncomp,ncsgs,vesgs_coun,prsgs_coun,visc_coun,resid_coun,bvesgs_coun
   integer(ip) :: ndime,ielem,nelem,ipoin,npoin,npoinLocal,iboun,nboun,ibody,nbody,pgaus,pgaub
   integer(ip) :: imat=1
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetnpoinLocal(npoinLocal)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNbody(nbody)
   call a%Mesh%GetNboun(nboun)

   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = 1 + nsteps

   call a%Memor%alloc(ndime,npoin,a%ncomp,a%veloc,'veloc','nsi_memall')
   call a%Memor%alloc(npoin,a%ncomp      ,a%press,'press','nsi_memall')

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      bvesgs_coun = 0
      vesgs_coun = 0
      prsgs_coun = 0
      ncsgs=1
      if(a%kfl_tacsg>0) ncsgs=2
      call a%Memor%alloc(nelem,a%vesgs,'vesgs','nsi_memall')
      call a%Memor%alloc(nelem,a%prsgs,'prsgs','nsi_memall')
      if(a%kfl_repro == 2 .or. a%kfl_repro == 3) call a%Memor%alloc(nelem,a%vesgs2,'vesgs2','nsi_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%vesgs(ielem)%a(ndime,ncsgs,pgaus))
         allocate(a%prsgs(ielem)%a(pgaus))
         if (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
            allocate(a%vesgs2(ielem)%a(ndime,ncsgs,pgaus))
            a%vesgs2(ielem)%a = 0.0_rp
         end if
         !initialization
         a%vesgs(ielem)%a = 0.0_rp
         a%prsgs(ielem)%a = 0.0_rp
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         prsgs_coun = prsgs_coun + pgaus
      end do
      call a%Memor%allocObj(0,'vesgs%a','nsi_memall',vesgs_coun*rp)
      if (a%kfl_repro ==2 .or. a%kfl_repro == 3) call a%Memor%allocObj(0,'vesgs2%a','nsi_memall',vesgs_coun*rp)
      if (a%kfl_bousg == 1) then
         call a%Mesh%GetBounArraySize(iboun,pnodb,pgaub)
         call a%Memor%alloc(nboun,a%bvesgs,'bvesgs','nsi_memall')
         do iboun=1,nboun
            allocate(a%bvesgs(iboun)%a(ndime,pgaub))
            a%bvesgs(iboun)%a = 0.0_rp
            bvesgs_coun = bvesgs_coun + ndime*pgaus
         end do
         call a%Memor%allocObj(0,'bvesgs%a','nsi_memall',bvesgs_coun*rp)
      end if
   end if
   
   !Residual Projection
   a%ResidualSize = ndime+1
   if (a%kfl_repro >=1) then
      a%ResidualSize = ndime+1
      if (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
         a%ResidualSize = a%ResidualSize + ndime
      endif
      call a%Memor%alloc(a%ResidualSize,npoin,a%repro,'repro','nsi_memall')
   endif
   
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1 .or. a%kfl_error) then
      call a%Memor%alloc(ndime,ndime,npoin,a%grprj,'grprj','nsi_memall')
   end if
   
   !Tau Smoothing
   if (a%kfl_Tausm >= 1) then
      call a%Memor%alloc(2,npoin,a%Tausmo,'Tausmo','nsi_memall')
   endif
   
   !Postprocess the residual at the gauss points
   if (a%npp_stepi(18) /= 0) then
      resid_coun = 0
      call a%Memor%alloc(nelem,a%residualU,'residual','nsi_memall')
      call a%Memor%alloc(nelem,a%residualP,'residual','nsi_memall')
      
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%residualU(ielem)%a(ndime,pgaus))
         allocate(a%residualP(ielem)%a(pgaus))
         a%residualU(ielem)%a = 0.0_rp
         a%residualP(ielem)%a = 0.0_rp
         resid_coun = resid_coun + (ndime+1)*pgaus
      end do
      
      if (a%kfl_repro == 2) then 
         call a%Memor%alloc(nelem,a%residualGraP2,'residual','nsi_memall')
         do ielem = 1,nelem
            allocate(a%residualGraP2(ielem)%a(ndime,pgaus))
            resid_coun = resid_coun + (ndime)*pgaus
         end do
      endif
      
      call a%Memor%allocObj(0,'residual%a','nsi_memall',resid_coun*rp)
   endif
   
   !Physical Properties defined in type MatProperties
   !we need define the number of materials for the moment nmat=1   
  
   if(a%MatProp(imat)%lawvi/=0)then 
      if ((a%kfl_repro >= 1).or.(a%npp_stepi(5)==1)) then   
         visc_coun = 0
         call a%Memor%alloc(nelem,a%viscarray,'viscarray','nsi_memall')   
         do ielem=1,nelem
            call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
            allocate(a%viscarray(ielem)%a(pgaus))
            visc_coun = visc_coun + pgaus        
         end do
         call a%Memor%allocObj(0,'viscarray%a','nsi_memall',visc_coun*rp)
      end if
   end if 
   
   !Divergence
    if (a%npp_stepi(11)>=1) then   
      visc_coun = 0
      call a%Memor%alloc(nelem,a%divergence,'divergence','nsi_memall')   
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         allocate(a%Divergence(ielem)%a(pgaus))
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%allocObj(0,'Divergence%a','nsi_memall',visc_coun*rp)
   end if 
   
   if (a%kfl_dispa /=0) then
      call a%Memor%alloc(npoin,a%Dissipation,'dissipation','nsi_memall')
   endif

   if (a%npp_stepi(10)>0) then
      call a%Memor%alloc(3,npoin,a%vorti,'vorti','nsi_memall')
   endif

   if (a%npp_stepi(21)>0) then
      call a%Memor%alloc(npoin,a%qfac,'qfac','nsi_memall')
   endif
   
   if(a%kfl_outfm>=1)then   
      call a%Memor%alloc(nbody,a%lun_force,'lun_force','nsi_memall')
      call a%Memor%alloc(ndime,nbody,a%force,'force','nsi_memall')
      call a%Memor%alloc(3,nbody,a%momen,'momen','nsi_memall')
   end if
   
   !Gauss point dissipation
   if (a%npp_stepi(20) /= 0) then
      call a%Memor%alloc(nelem,a%GPDissipation,'GPDissipation','nsi_memall')   
      visc_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         allocate(a%GPDissipation(ielem)%a(pgaus))
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%allocObj(0,'GPDissipation%a','nsi_memall',visc_coun*rp)
   endif
   
   if (a%kfl_StabilizeFreeSurface == 1) then
      call a%Memor%alloc(ndime, ndime,npoin,a%StabilizeFreeSurface_VelocityGradient,'StabilizeFreeSurface_VelocityGradient','nsi_memall')
      call a%Memor%alloc(ndime,npoin,a%StabilizeFreeSurface_PressureGradientAndForce,'StabilizeFreeSurface_PressureGradientAndForce','nsi_memall')
   endif

    !Statistics of fields
   if (a%npp_stepi(24) /= 0) then
     call a%Memor%alloc(ndime,npoin,5,a%timav,'timav','nsi_memall')
     call a%Memor%alloc(npoin,     3,a%timap,'timap','nsi_memall')
     call a%Memor%alloc(ndime,npoin,a%rmsv,'rmsv','nsi_memall')
     call a%Memor%alloc(npoin,a%rmsp,'rmsp','nsi_memall')
     call a%Memor%alloc(npoin,a%turbi,'turbi','nsi_memall')
     call a%Memor%alloc(ndime,npoin,2,a%TimaTractions,'TimaTractions','nsi_memall')
   endif
   
   if (a%npp_stepi(25) /= 0) then  
      call a%Memor%alloc(ndime,npoin,a%ExternalForcesArray,'ExternalForces','nsi_memall')
   endif
   
   !Exact solution
   if(a%kfl_exacs/=0 .and. a%kfl_timei /=0 .and. a%npp_stepi(27) /=0) then
      call a%Memor%alloc(ndime,npoin,1,a%AnalyticalVelocity,'AnalyticalVelocity','nsi_memall')
      call a%Memor%alloc(npoin,1,a%AnalyticalPressure,'AnalyticalPressure','nsi_memall')
   endif

   !Boundary coupling
   call a%Memor%alloc(ndime,npoin,a%btraction,'btraction','nsi_memall')
   if (a%kfl_bousg == 1) then
      call a%Memor%alloc(npoin,a%btau,'btau','nsi_memall')
   end if
   if (a%kfl_docoupconv) then
      call a%Memor%alloc(ndime,npoin,a%veloc_cp,'veloc_cp','nsi_memall')
      call a%Memor%alloc(npoin,a%press_cp,'press_cp','nsi_memall')
      call a%Memor%alloc(a%ndofbc,npoin,2,a%velocDD,'velocDD','nsi_memall')
      call a%Memor%alloc(a%ndofbc,npoin,2,a%tractDD,'tractDD','nsi_memall')
      call a%Memor%alloc(a%ndofbc,npoin,3,a%bres,'bres','nsi_memall')
   end if
      
end subroutine nsi_memall

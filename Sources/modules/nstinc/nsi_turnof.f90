subroutine nsi_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the incompressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   interface
      subroutine nsi_FinalizeTurbulentInletBoundaryConditions(a)
         import
         implicit none
         class(NavierStokesProblem) :: a
      end subroutine
   end interface
   
   integer(ip) :: pnode,pnodb,ncomp,ncsgs,vesgs_coun,prsgs_coun,visc_coun,resid_coun,bvesgs_coun
   integer(ip) :: ndime,ielem,nelem,ipoin,npoin,npoinLocal,iboun,nboun,ibody,nbody,pgaus,pgaub
   integer(ip) :: imat=1
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetnpoinLocal(npoinLocal)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNbody(nbody)
   ncomp = a%ncomp
   
   if(allocated(a%PointwiseSource)) call a%Memor%dealloc(size(a%PointwiseSource,1),size(a%PointwiseSource,2),a%PointwiseSource,'PointwiseSource','nsi_turnof')

   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF NAVIER-STOKES RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF NAVIER-STOKES RUN * * *'
   endif
   
   call a%Memor%dealloc(ndime,npoin,a%ncomp,a%veloc,'veloc','nsi_turnof')
   call a%Memor%dealloc(npoin,a%ncomp      ,a%press,'press','nsi_turnof')
   
   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      ncsgs=1
      bvesgs_coun = 0
      vesgs_coun = 0
      prsgs_coun = 0
      if(a%kfl_tacsg>0) ncsgs=2
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%vesgs(ielem)%a)      
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         deallocate(a%prsgs(ielem)%a)      
         prsgs_coun = prsgs_coun + pgaus
         if(a%kfl_repro == 2 .or. a%kfl_repro == 3) deallocate(a%vesgs2(ielem)%a)
      end do
      call a%Memor%deallocObj(0,'vesgs%a','nsi_turnof',vesgs_coun*rp)
      call a%Memor%dealloc(nelem,a%vesgs,'vesgs','nsi_turnof')
      call a%Memor%deallocObj(0,'prsgs%a','nsi_turnof',prsgs_coun*rp)
      call a%Memor%dealloc(nelem,a%prsgs,'prsgs','nsi_turnof')
      if (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
         call a%Memor%deallocObj(0,'vesgs2%a','nsi_turnof',vesgs_coun*rp)
         call a%Memor%dealloc(nelem,a%vesgs2,'vesgs2','nsi_turnof')  
      end if
      if (a%kfl_bousg == 1) then
         call a%Mesh%GetBounArraySize(iboun,pnodb,pgaub)
         do iboun=1,nboun
            deallocate(a%bvesgs(iboun)%a)
            bvesgs_coun = bvesgs_coun + ndime*pgaus
         end do
         call a%Memor%deallocObj(0,'bvesgs%a','nsi_turnof',bvesgs_coun*rp)
         call a%Memor%dealloc(nboun,a%bvesgs,'bvesgs','nsi_turnof')  
      end if
   end if

   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%dealloc(size(a%repro,1),npoin,a%repro,'repro','nsi_turnof')
   endif
   
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1 .or. a%kfl_error) then
      call a%Memor%dealloc(ndime,ndime,npoin,a%grprj,'grprj','nsi_turnof')
   end if
   
   !Tau Smoothing
   if (a%kfl_Tausm >= 1) then
      call a%Memor%dealloc(2,npoin,a%Tausmo,'Tausmo','nsi_memall')
   endif
   
   !Postprocess the residual at the gauss points
   if (a%npp_stepi(18) /= 0) then
      resid_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%residualU(ielem)%a)
         deallocate(a%residualP(ielem)%a)
         resid_coun = resid_coun + (ndime+1)*pgaus
      end do
      
      if (a%kfl_repro == 2) then 
         do ielem = 1,nelem
            deallocate(a%residualGraP2(ielem)%a)
            resid_coun = resid_coun + (ndime)*pgaus
         ENDDO
         call a%Memor%dealloc(nelem,a%residualGraP2,'residual','nsi_memall')
         
      endif
      
      call a%Memor%deallocObj(0,'residual%a','nsi_turnof',resid_coun*rp)
      call a%Memor%dealloc(nelem,a%residualU,'residual','nsi_turnof')
      call a%Memor%dealloc(nelem,a%residualP,'residual','nsi_turnof')
   endif
   
   !Viscosity at gauss points
   if(a%MatProp(imat)%lawvi/=0)then 
      if ((a%kfl_repro >= 1).or.(a%npp_stepi(5)==1)) then      
      visc_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%viscarray(ielem)%a)
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'viscarray%a','nsi_turnof',visc_coun*rp)      
      call a%Memor%dealloc(nelem,a%viscarray,'viscarray','nsi_turnof')  
      end if      
   end if  
   
   !Divergence
   if (a%npp_stepi(11)>=1) then      
      visc_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%Divergence(ielem)%a)
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'Divergence%a','nsi_turnof',visc_coun*rp)      
      call a%Memor%dealloc(nelem,a%Divergence,'Divergence','nsi_turnof')  
   end if  
   
   if (a%kfl_dispa /=0) then
      call a%Memor%dealloc(npoin,a%Dissipation,'dissipation','nsi_turnof')
   endif

   if (a%npp_stepi(10)>0) then
      call a%Memor%dealloc(3,npoin,a%vorti,'vorti','nsi_memall')
   endif

   if (a%npp_stepi(21)>0) then
      call a%Memor%dealloc(npoin,a%qfac,'qfac','nsi_memall')
   endif   
   
   if(a%kfl_outfm==1)then   
      call a%Memor%dealloc(ndime,nbody,a%force,'force','nsi_turnof')
      call a%Memor%dealloc(3,nbody,a%momen,'momen','nsi_turnof')
      do ibody=1,nbody
         close(a%lun_force(ibody))
      end do
      call a%Memor%dealloc(nbody,a%lun_force,'lun_force','nsi_turnof')
   end if
  
   !Reference systems
   call a%Memor%dealloc(npoin,a%kfl_fixrs,'kfl_fixrs','nsi_turnof')
   call a%Memor%dealloc(nboun,a%kfl_bours,'kfl_bours','nsi_turnof')
   
   !Gauss point dissipation
   if (a%npp_stepi(20) /= 0) then
      visc_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%GPDissipation(ielem)%a)
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'GPDissipation%a','nsi_memall',visc_coun*rp)
      call a%Memor%dealloc(nelem,a%GPDissipation,'GPDissipation','nsi_memall')   
      
   endif
   
   !Free surface stabilization
   if (a%kfl_StabilizeFreeSurface == 1) then
      call a%Memor%dealloc(ndime, ndime,npoin,a%StabilizeFreeSurface_VelocityGradient,'StabilizeFreeSurface_VelocityGradient','nsi_memall')
      call a%Memor%dealloc(ndime,npoin,a%StabilizeFreeSurface_PressureGradientAndForce,'StabilizeFreeSurface_PressureGradientAndForce','nsi_memall')
   endif
   
   if (a%nptra/=0) then
      if (a%MPIrank == a%MPIroot) then
         call iofile(two,a%lun_trap,a%fil_trap,'NSI TRACKING OF POINTS')
      endif
   endif
   
   !Statistics of fields
   if (a%npp_stepi(24) /= 0) then
      call a%Memor%dealloc(ndime,npoin,5,a%timav,'timav','nsi_turnof')
      call a%Memor%dealloc(npoin,     3,a%timap,'timap','nsi_turnof')
      call a%Memor%dealloc(ndime,npoin,a%rmsv,'rmsv','nsi_turnof')
      call a%Memor%dealloc(npoin,a%rmsp,'rmsp','nsi_turnof')
      call a%Memor%dealloc(npoin,a%turbi,'turbi','nsi_turnof')
      call a%Memor%dealloc(ndime,npoin,2,a%TimaTractions,'TimaTractions','nsi_memall')
   endif
   
   if (a%kfl_ElasticBoundary == 1) then
      call a%Memor%dealloc(ndime,npoin,a%EB_Displacement,'EB_Displacement','nsi_memall')
   endif
   
   if (a%kfl_TurbulentInletBoundaryConditions == 1) call nsi_FinalizeTurbulentInletBoundaryConditions(a)

   
   if (a%npp_stepi(25) /= 0) then  
      call a%Memor%dealloc(ndime,npoin,a%ExternalForcesArray,'ExternalForces','nsi_memall')
   endif
   
   !Exact solution
   if(a%kfl_exacs/=0 .and. a%kfl_timei /=0 .and. a%npp_stepi(27) /=0) then
      call a%Memor%dealloc(ndime,npoin,1,a%AnalyticalVelocity,'AnayticalVelocity','nsi_turnof')
      call a%Memor%dealloc(npoin,1,a%AnalyticalPressure,'AnalyticalPressure','nsi_turnof')
   endif

   !Boundary coupling
   call a%Memor%dealloc(ndime,npoin,a%btraction,'btraction','nsi_turnof')
   if (a%kfl_bousg == 1) then
      call a%Memor%dealloc(npoin,a%btau,'btau','nsi_turnof')
   end if
   if (a%kfl_docoupconv) then
      call a%Memor%dealloc(ndime,npoin,a%veloc_cp,'veloc_cp','nsi_turnof')
      call a%Memor%dealloc(npoin,a%press_cp,'press_cp','nsi_turnof')
      call a%Memor%dealloc(a%ndofbc,npoin,2,a%velocDD,'velocDD','nsi_turnof')
      call a%Memor%dealloc(a%ndofbc,npoin,2,a%tractDD,'tractDD','nsi_turnof')
      call a%Memor%dealloc(a%ndofbc,npoin,3,a%bres,'bres','nsi_turnof')
   end if
   
   20 format(21(1x,e14.7))     

end subroutine nsi_turnof


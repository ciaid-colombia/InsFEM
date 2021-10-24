subroutine sup_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the incompressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem,ncomp,ncsgs,pgaus,iaux,jaux,ielem,pnode,nboun
   integer(ip) :: vesgs_coun,visc_coun,lamb_coun, newsigma_coun, ndim_aux, sisgs_coun, tau_coun
   integer(ip) :: auxtens,kfl_nonlinear,nbody,ibody
   
   !todo multy materials
   integer(ip) :: imat=1
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNbody(nbody)   
   ncomp = a%ncomp
   auxtens=(ndime-1)*(ndime-1)+2
   
   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF THREE-FIELD-NAVIER-STOKES RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF THREE-FIELD-NAVIER-STOKES RUN * * *'
   endif
   
   call a%Memor%dealloc(ndime,npoin,a%ncomp,a%veloc,'veloc','sup_turnof')
   call a%Memor%dealloc(npoin,a%ncomp      ,a%press,'press','sup_turnof')
   call a%Memor%dealloc(auxtens,npoin,a%ncomp,a%sigma,'sigma','sup_turnof')
   
   if(a%kfl_docoupconv) call a%Memor%dealloc(ndime,npoin,a%veloc_cp,'veloc_cp','sup_turnof')
   if(a%kfl_docoupconv) call a%Memor%dealloc(npoin,a%press_cp,'press_cp','sup_turnof')
   if(a%kfl_docoupconv) call a%Memor%dealloc(auxtens,npoin,a%sigma_cp,'sigma_cp','sup_turnof')
   
   call a%Memor%dealloc(ndime,npoin,a%btraction,'btraction','nsi_memall')
   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      ncsgs=1
      vesgs_coun = 0
      sisgs_coun = 0
      ndim_aux = 0 !counter to use dynamic split oss
      if(a%kfl_tacsg>0) ncsgs=2
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%vesgs(ielem)%a) 
         deallocate(a%sisgs(ielem)%a)
         if(a%kfl_repro>=2) then
            deallocate(a%vesgs2(ielem)%a)   
            deallocate(a%vesgs3(ielem)%a) 
            if (a%MatProp(imat)%lawvi<0 .and. a%kfl_repro==4) then
               deallocate(a%sisgs2(ielem)%a)   
               deallocate(a%sisgs3(ielem)%a) 
            end if
         end if   
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         sisgs_coun = sisgs_coun + auxtens*ncsgs*pgaus
      end do
      
      call a%Memor%deallocObj(0,'vesgs%a','sup_turnof',vesgs_coun*rp)
      call a%Memor%dealloc(nelem,a%vesgs,'vesgs','sup_turnof')
      
      
      call a%Memor%deallocObj(0,'sisgs%a','sup_turnof',sisgs_coun*rp)
      call a%Memor%dealloc(nelem,a%sisgs,'sisgs','sup_turnof')
      
      
      if(a%kfl_repro >= 2) then
        call a%Memor%deallocObj(0,'vesgs2%a','sup_turnof',vesgs_coun*rp)
        call a%Memor%dealloc(nelem,a%vesgs2,'vesgs2','sup_turnof')      
        call a%Memor%deallocObj(0,'vesgs3%a','sup_turnof',vesgs_coun*rp)
        call a%Memor%dealloc(nelem,a%vesgs3,'vesgs3','sup_turnof')  
        
         if (a%MatProp(imat)%lawvi<0 .and. a%kfl_repro==4) then
            call a%Memor%deallocObj(0,'sisgs2%a','sup_turnof',vesgs_coun*rp)
            call a%Memor%dealloc(nelem,a%sisgs2,'sisgs2','sup_turnof')      
            call a%Memor%deallocObj(0,'sisgs3%a','sup_turnof',vesgs_coun*rp)
            call a%Memor%dealloc(nelem,a%sisgs3,'sisgs3','sup_turnof')          
         end if
      end if
   end if
   
   
   if (a%kfl_repro >= 1) then
      iaux = auxtens+ndime+1
      call a%Memor%dealloc(iaux,npoin,a%repro,'repro','sup_turnof')
   endif   
   
   !Split Oss
   if(a%kfl_repro ==2 .or. a%kfl_repro ==3 .or. a%kfl_repro ==4)then
      call a%Memor%dealloc(ndime,npoin,a%reproSDiv,'reproSDiv','sup_turnof')
      call a%Memor%dealloc(ndime,npoin,a%reproUGradU,'reproUGradU','sup_turnof')
      call a%Memor%dealloc(ndime,npoin,a%reproGradP,'reproGradP','sup_turnof')
      call a%Memor%dealloc(1,npoin,a%reproDivU,'reproDivU','sup_turnof')
      if (a%kfl_repro ==4) then
         call a%Memor%dealloc(auxtens,npoin,a%reproGradU,'reproGradU','sup_turnof')
      end if
   end if   
         
   if (a%MatProp(imat)%lawvi < 0) then !Viscoelastic case
      call a%Memor%dealloc(auxtens,npoin,a%reproUGradS,'reproUGradS','sup_turnof')
      call a%Memor%dealloc(auxtens,npoin,a%reproSGradU,'reproSGradU','sup_turnof')   
      
      if (a%LogFormulation==1 .and. a%kfl_repro ==4) then
         call a%Memor%dealloc(auxtens,npoin,a%reproExpS,'reproExpS','sup_turnof')
      end if
      !non_linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)  
      if(kfl_nonlinear ==1)then
         call a%Memor%dealloc(ndime,npoin,a%reproLapla,'reproLapla','sup_turnof')
      end if
  
   endif 

   if(a%MatProp(imat)%lawvi/=0)then 
      if ((a%kfl_repro == 1).or.(a%npp_stepi(5)>=1).or.(a%kfl_shock == 1)) then         
      visc_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%viscarray(ielem)%a)
         visc_coun = visc_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'viscarray%a','sup_turnof',visc_coun*rp)      
      call a%Memor%dealloc(nelem,a%viscarray,'viscarray','sup_turnof')
      endif
   end if  
   
   !Discontinuity capturing
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 1)then
      iaux = auxtens*ndime
      call a%Memor%dealloc(iaux,npoin,a%reproGrad,'reproGrad','sup_turnof')    
   endif
   !DEVSS
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 3)then
      jaux = ndime*ndime
      call a%Memor%dealloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_turnof')      
   endif
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 2)then
      jaux = ndime*ndime
      call a%Memor%dealloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_turnof')      
   endif

   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 4)then
      jaux = ndime*ndime
      call a%Memor%dealloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_turnof')
      iaux = auxtens*ndime
      call a%Memor%dealloc(iaux,npoin,a%reproGrad,'reproGrad','sup_turnof')          
   endif
   
   if(a%kfl_exacs/=0.and.a%kfl_timei/=0)then
      call a%Memor%dealloc(ndime,npoin,1,a%exveloc,'exveloc','sup_turnof')
      call a%Memor%dealloc(npoin,1      ,a%express,'express','sup_turnof')
      call a%Memor%dealloc(auxtens,npoin,1,a%exsigma,'exsigma','sup_turnof')
   end if   
   
   call a%Memor%dealloc(npoin,a%kfl_fixrs,'kfl_fixrs','sup_turnof')
   call a%Memor%dealloc(nboun,a%kfl_bours,'kfl_bours','sup_turnof')
   
   if(a%kfl_outfm>=1)then   
      call a%Memor%dealloc(ndime,nbody,a%force,'force','nsi_turnof')
      call a%Memor%dealloc(3,nbody,a%momen,'momen','nsi_turnof')
      do ibody=1,nbody
         close(a%lun_force(ibody))
      end do
      call a%Memor%dealloc(nbody,a%lun_force,'lun_force','nsi_turnof')
   end if
  

   !Coupling with temperature model
   if (a%kfl_cotem_WLF==1 .or. a%kfl_cotem_Arrhenius==1) then
      lamb_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%lambdarray(ielem)%a)
         lamb_coun = lamb_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'lambdarray%a','sup_turnof',lamb_coun*rp)    
      call a%Memor%dealloc(nelem,a%lambdarray,'lambdarray','sup_turnof')
   end if 
   
   !Logarithmic model
   if (a%LogFormulation/=0) then
      call a%Memor%dealloc(auxtens,npoin,a%sigmaold,'sigmaold','sup_turnof')
      if(a%npp_stepi(28)>0) call a%Memor%dealloc(auxtens,npoin,a%psiReal,'psiReal','sup_turnof')
   
      lamb_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%alpha3array(ielem)%a)
         lamb_coun = lamb_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'alpha3array%a','sup_turnof',lamb_coun*rp)    
      call a%Memor%dealloc(nelem,a%alpha3array,'alpha3array','sup_turnof')
      
   end if
    
   !Tau Smoothing
   if (a%kfl_Tausm >= 1) then
      call a%Memor%dealloc(3,npoin,a%Tausmo,'Tausmo','sup_turnof')
   endif
      

   call a%Memor%dealloc(npoin,a%taudet,'taudet','sup_turnof')
   
   !Vorticity
   if (a%npp_stepi(10)>0) then
      call a%Memor%dealloc(3,npoin,a%vorti,'vorti','nsi_memall')
   endif
   
   !Taus
   if (a%npp_stepi(26)>=1) then
      tau_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         deallocate(a%tau_mom(ielem)%a)
         deallocate(a%tau_sig(ielem)%a)
         tau_coun = tau_coun + pgaus        
      end do
      call a%Memor%deallocObj(0,'tau_mom%a','sup_turnof',tau_coun*rp)    
      call a%Memor%dealloc(nelem,a%tau_mom,'tau_mom','sup_turnof')
      call a%Memor%deallocObj(0,'tau_sig%a','sup_turnof',tau_coun*rp)    
      call a%Memor%dealloc(nelem,a%tau_sig,'tau_sig','sup_turnof')
   end if   
  

end subroutine sup_turnof


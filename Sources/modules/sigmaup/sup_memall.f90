subroutine sup_memall(a)
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
   use Mod_Listen
   use Mod_Memor
   use Mod_MPIObject
   use Mod_Mesh
   use Mod_TimeIntegrator   
   use Mod_Postpr
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   use Mod_ThreeField
   implicit none
  
   class(ThreeFieldNSProblem), target :: a
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps   
   
   integer(ip) :: ndime,npoin,nelem,ncomp,ncsgs,pgaus,iaux,jaux,ielem,pnode,vesgs_coun,visc_coun,ntens,nbody, lamb_coun, newsigma_coun, tau1_coun
   integer(ip) :: kfl_nonlinear, sisgs_coun
   type(MemoryMan), pointer :: Memor
   integer(ip) :: imat=1
   
   Memor => a%Memor

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNbody(nbody)

   ntens=(ndime-1)*(ndime-1)+2    

   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = 1 + nsteps   
   
   call Memor%alloc(ndime,npoin,a%ncomp,a%veloc,'veloc','sup_memall')
   call Memor%alloc(npoin,a%ncomp      ,a%press,'press','sup_memall')
   call Memor%alloc(ntens,npoin,a%ncomp,a%sigma,'sigma','sup_memall')

   if (a%kfl_docoupconv) then
       call a%Memor%alloc(ndime,npoin,a%veloc_cp,'veloc_cp','sup_memall')
       call a%Memor%alloc(npoin,a%press_cp,'press_cp','sup_memall')
       call a%Memor%alloc(ntens,npoin,a%sigma_cp,'sigma_cp','sup_memall')
       a%veloc_cp = 0.0_rp
       a%press_cp = 0.0_rp
       a%sigma_cp = 0.0_rp
   endif

   call a%Memor%alloc(ndime,npoin,a%btraction,'btraction','sup_memall')
   
   !Residual Projection
   a%ResidualSize = ntens+ndime+1
   if (a%kfl_repro >=1) then
      call Memor%alloc(a%ResidualSize,npoin,a%repro,'repro','sup_memall')
   endif
   
   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      vesgs_coun = 0
      sisgs_coun = 0
      ncsgs=1
      if(a%kfl_tacsg>0) ncsgs=2
      call a%Memor%alloc(nelem,a%vesgs,'vesgs','sup_memall')
      call a%Memor%alloc(nelem,a%sisgs,'sisgs','sup_memall')
      
      if(a%kfl_repro >= 2) then
         call a%Memor%alloc(nelem,a%vesgs2,'vesgs2','sup_memall')
         call a%Memor%alloc(nelem,a%vesgs3,'vesgs3','sup_memall')
         
         if (a%kfl_repro==4 .and. a%MatProp(imat)%lawvi<0) then
            call a%Memor%alloc(nelem,a%sisgs2,'sisgs2','sup_memall')
            call a%Memor%alloc(nelem,a%sisgs3,'sisgs3','sup_memall')
         end if
         
      end if   
      
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%vesgs(ielem)%a(ndime,ncsgs,pgaus))
         allocate(a%sisgs(ielem)%a(ntens,ncsgs,pgaus))
         
         if(a%kfl_repro >= 2)then
            allocate(a%vesgs2(ielem)%a(ndime,ncsgs,pgaus))
            allocate(a%vesgs3(ielem)%a(ndime,ncsgs,pgaus))
            a%vesgs2(ielem)%a = 0.0_rp
            a%vesgs3(ielem)%a = 0.0_rp
            
            if (a%kfl_repro==4 .and. a%MatProp(imat)%lawvi<0) then
               allocate(a%sisgs2(ielem)%a(ntens,ncsgs,pgaus))
               allocate(a%sisgs3(ielem)%a(ntens,ncsgs,pgaus))
               a%sisgs2(ielem)%a = 0.0_rp
               a%sisgs3(ielem)%a = 0.0_rp
            end if
            
         end if
         !initialization
         a%vesgs(ielem)%a = 0.0_rp
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         
        
         a%sisgs(ielem)%a = 0.0_rp
         sisgs_coun = sisgs_coun + ntens*ncsgs*pgaus   
         
         
      end do
      call a%Memor%allocObj(0,'vesgs%a','sup_memall',vesgs_coun*rp)

      call a%Memor%allocObj(0,'sisgs%a','sup_memall',sisgs_coun*rp)
      
      if(a%kfl_repro >= 2) then
         call a%Memor%allocObj(0,'vesgs2%a','sup_memall',vesgs_coun*rp)
         call a%Memor%allocObj(0,'vesgs3%a','sup_memall',vesgs_coun*rp)
         if (a%MatProp(imat)%lawvi<0 .and. a%kfl_repro==4) then
            call a%Memor%allocObj(0,'sisgs2%a','sup_memall',sisgs_coun*rp)
            call a%Memor%allocObj(0,'sisgs3%a','sup_memall',sisgs_coun*rp)
         end if
      end if 
      
   end if
   
   !Residual Projection 
   if(a%kfl_repro ==2 .or. a%kfl_repro ==3 .or. a%kfl_repro ==4)then
      call Memor%alloc(ndime,npoin,a%reproSDiv,'reproSDiv','sup_memall')
      call Memor%alloc(ndime,npoin,a%reproUGradU,'reproUGradU','sup_memall')
      call Memor%alloc(ndime,npoin,a%reproGradP,'reproGradP','sup_memall')
      call Memor%alloc(1,npoin,a%reproDivU,'reproDivU','sup_memall')
      if (a%kfl_repro==4) then
         
         call Memor%alloc(ntens,npoin,a%reproGradU,'reproGradU','sup_memall')
      end if
   end if   
   
   if (a%MatProp(imat)%lawvi<0) then
      call Memor%alloc(ntens,npoin,a%reproUGradS,'reproUGradS','sup_memall')
      call Memor%alloc(ntens,npoin,a%reproSGradU,'reproSGradU','sup_memall')
      if (a%LogFormulation==1 .and. (a%kfl_repro ==4)) then
         call Memor%alloc(ntens,npoin,a%reproExpS,'reproExpS','sup_memall')
      end if
      !non_linear
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)      
      if(kfl_nonlinear == 1)then
         call Memor%alloc(ndime,npoin,a%reproLapla,'reproLapla','sup_memall')
      end if  
   endif   
   
   if(a%MatProp(imat)%lawvi/=0)then 
      if((a%kfl_repro == 1).or.(a%npp_stepi(5)>=1).or.(a%kfl_shock == 1))then
         visc_coun = 0
         call Memor%alloc(nelem,a%viscarray,'viscarray','sup_memall')   
         do ielem=1,nelem
            call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
            allocate(a%viscarray(ielem)%a(pgaus))
            visc_coun = visc_coun + pgaus        
         end do
         call Memor%allocObj(0,'viscarray%a','sup_memall',visc_coun*rp)     
      end if
   end if 
   
   !Discontinuity Capturing
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 1)then
      iaux = ntens*ndime
      call Memor%alloc(iaux,npoin,a%reproGrad,'reproGrad','sup_memall')  
   end if  
   !DEVSS
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 3)then
      jaux = ndime*ndime
      call Memor%alloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_memall')      
   end if  
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 2)then
      jaux = ndime*ndime
      call Memor%alloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_memall')      
   end if
   
   if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 4)then
      jaux = ndime*ndime
      call Memor%alloc(jaux,npoin,a%reproSGrad,'reproSGrad','sup_memall') 
      iaux = ntens*ndime
      call Memor%alloc(iaux,npoin,a%reproGrad,'reproGrad','sup_memall')        
   end if  
   
   if(a%kfl_exacs/=0.and.a%kfl_timei/=0)then
      call Memor%alloc(ndime,npoin,1,a%exveloc,'exveloc','sup_memall')
      call Memor%alloc(npoin,1      ,a%express,'express','sup_memall')
      call Memor%alloc(ntens,npoin,1,a%exsigma,'exsigma','sup_memall')
   end if
   
   
   if(a%kfl_outfm>=1)then   
      call a%Memor%alloc(nbody,a%lun_force,'lun_force','nsi_memall')
      call a%Memor%alloc(ndime,nbody,a%force,'force','nsi_memall')
      call a%Memor%alloc(3,nbody,a%momen,'momen','nsi_memall')
   end if
   
   !Temperature coupling   
   if (a%kfl_cotem_WLF==1 .or. a%kfl_cotem_Arrhenius==1) then
      lamb_coun = 0
      call Memor%alloc(nelem,a%lambdarray,'lambdarray','sup_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         allocate(a%lambdarray(ielem)%a(pgaus))
         lamb_coun = lamb_coun + pgaus        
      end do
      call Memor%allocObj(0,'lambdarray%a','sup_memall',lamb_coun*rp)     
   end if 
   
   !Logarithmic problem
   if (a%LogFormulation/=0) then
      call Memor%alloc(ntens,npoin,a%sigmaold,'sigmaold','sup_memall')
      if(a%npp_stepi(28)>0) call Memor%alloc(ntens,npoin,a%psiReal,'psiReal','sup_memall')
      lamb_coun = 0
      call Memor%alloc(nelem,a%alpha3array,'alpha3array','sup_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         allocate(a%alpha3array(ielem)%a(pgaus))
         lamb_coun = lamb_coun + pgaus        
      end do
      call Memor%allocObj(0,'alpha3array%a','sup_memall',lamb_coun*rp)  
   end if
   
   !Tau Smoothing
   if (a%kfl_Tausm >= 1) then
      call a%Memor%alloc(3,npoin,a%Tausmo,'Tausmo','sup_memall')
   endif

   call Memor%alloc(npoin,a%taudet,'taudet','sup_memall')
   
   !Vorticity
    if (a%npp_stepi(10)>0) then
      call a%Memor%alloc(3,npoin,a%vorti,'vorti','sup_memall')
   endif
   
   !Taus
   if (a%npp_stepi(26)>=1) then
      tau1_coun = 0
      call Memor%alloc(nelem,a%tau_mom,'tau_mom','sup_memall')
      call Memor%alloc(nelem,a%tau_sig,'tau_sig','sup_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
         allocate(a%tau_mom(ielem)%a(pgaus))
         allocate(a%tau_sig(ielem)%a(pgaus))
         tau1_coun = tau1_coun + pgaus        
      end do
      call Memor%allocObj(0,'tau_mom%a','sup_memall',tau1_coun*rp)    
      call Memor%allocObj(0,'tau_sig%a','sup_memall',tau1_coun*rp) 
   end if   

end subroutine sup_memall
  



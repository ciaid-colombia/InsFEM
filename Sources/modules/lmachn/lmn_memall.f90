subroutine lmn_memall(a)
   !-----------------------------------------------------------------------
   !****f* LMACHN/lmn_memall
   ! NAME 
   !    lmn_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_TimeIntegrator
   use Mod_LowMach
   implicit none
   
   class(LowMachProblem) :: a
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps
   integer(ip) :: ndime,npoin,nelem,ncsgs,nbody,pgaus,iaux,ielem,pnode,vesgs_coun,tesgs_coun,prsgs_coun,resid_coun
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNbody(nbody)
   
   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = 1 + nsteps

   call a%Memor%alloc(ndime,npoin,a%ncomp,a%veloc,'veloc','lmn_memall')
   call a%Memor%alloc(npoin,a%ncomp      ,a%press,'press','lmn_memall')
   call a%Memor%alloc(npoin,a%ncomp      ,a%tempe,'tempe','lmn_memall')
   call a%Memor%alloc(npoin              ,a%itemp,'itemp','lmn_memall')
   call a%Memor%alloc(      a%ncomp      ,a%pther,'pther','lmn_memall')
   if(a%npp_stepi(4)>0) then     
      call a%Memor%alloc(npoin            ,a%densf,'densf','lmn_memall')
   end if

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      vesgs_coun = 0
      tesgs_coun = 0
      prsgs_coun = 0
      ncsgs=1                  !Number of components in time for accuracy
      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1
      call a%Memor%alloc(nelem,a%vesgs,'vesgs','lmn_memall')
      call a%Memor%alloc(nelem,a%tesgs,'tesgs','lmn_memall')
      call a%Memor%alloc(nelem,a%prsgs,'prsgs','lmn_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%vesgs(ielem)%a(ndime,ncsgs,pgaus))
         a%vesgs(ielem)%a = 0.0_rp
         allocate(a%tesgs(ielem)%a(ncsgs,pgaus))
         a%tesgs(ielem)%a = 0.0_rp
         allocate(a%prsgs(ielem)%a(pgaus))
         a%prsgs(ielem)%a = 0.0_rp
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         tesgs_coun = tesgs_coun + ncsgs*pgaus
         prsgs_coun = prsgs_coun + pgaus
      end do
      call a%Memor%allocObj(0,'vesgs%a','lmn_memall',vesgs_coun*rp)
      call a%Memor%allocObj(0,'tesgs%a','lmn_memall',tesgs_coun*rp)
      call a%Memor%allocObj(0,'prsgs%a','lmn_memall',prsgs_coun*rp)
   end if

   !Residual Projection
   if (a%kfl_repro >=1) then
      iaux = ndime+2
      call a%Memor%alloc(iaux,npoin,a%repro,'repro','lmn_memall')
   endif
  
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      call a%Memor%alloc(ndime+1,ndime,npoin,a%grprj,'grprj','lmn_memall')
   end if

   !Postprocess the residual at the gauss points
   if (a%npp_stepi(11) /= 0) then
      resid_coun = 0
      call a%Memor%alloc(nelem,a%residualU,'residual','lmn_memall')
      call a%Memor%alloc(nelem,a%residualP,'residual','lmn_memall')
      call a%Memor%alloc(nelem,a%residualT,'residual','lmn_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%residualU(ielem)%a(ndime,pgaus))
         allocate(a%residualP(ielem)%a(pgaus))
         allocate(a%residualT(ielem)%a(pgaus))
         a%residualU(ielem)%a = 0.0_rp
         a%residualP(ielem)%a = 0.0_rp
         a%residualT(ielem)%a = 0.0_rp
         resid_coun = resid_coun + (ndime+1)*pgaus
      end do
      call a%Memor%allocObj(0,'residual%a','lmn_memall',resid_coun*rp)
   endif
   
   if(a%kfl_outfm==1)then   
      call a%Memor%alloc(nbody,a%lun_force,'lun_force','lmn_memall')
      call a%Memor%alloc(ndime,nbody,a%force,'force','lmn_memall')
      call a%Memor%alloc(3,nbody,a%momen,'momen','lmn_memall')
      call a%Memor%alloc(nbody,a%heatf,'heatf','lmn_memall')
   end if
end subroutine lmn_memall

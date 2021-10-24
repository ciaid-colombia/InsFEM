subroutine tem_memall(a)
   use typre
  
   use Mod_Memor
   use Mod_MPIObject
   use Mod_PhysicalProblem 
   use Mod_Mesh
   use Mod_Temperature
   use def_parame
   use Mod_TimeIntegrator
   use Mod_r1pElementAllocation
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem,tesgs_coun,ncsgs,pnode,pgaus,ielem, sig_coun

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps,nbody
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNbody(nbody)
   
   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = 1 + nsteps
   
   
   call a%Memor%alloc(npoin,a%ncomp,a%tempe,'tempe','tem_memall')
   
   if(a%kfl_trasg/=0) then
      call a%Memor%alloc(nelem,a%tesgs,'tesgs','tem_memall')
      tesgs_coun = 0
      ncsgs=2
      
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%tesgs(ielem)%a(ncsgs,pgaus))
         a%tesgs(ielem)%a = 0.0_rp
         tesgs_coun = tesgs_coun + ncsgs*pgaus
      end do
      call a%Memor%allocObj(0,'tesgs%a','tem_memall',tesgs_coun*rp)
   end if

   if(a%kfl_repro==1 .or. a%kfl_shock == 1 .or. a%kfl_shock == 2) then
      call a%Memor%alloc(npoin,a%repro,'repro','tem_memall')
   endif
  
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      call a%Memor%alloc(ndime,npoin,a%grprj,'grprj','tem_memall')
   end if
   
   if (a%kfl_shock == 3 .or. a%kfl_shock == 4) then
      call a%Memor%alloc(ndime,npoin,a%gradient,'gradient','tem_memall')
      call AllocateR2pElement(a%Mesh,ndime,a%GradientGaussPoints,a%Memor,'gradient')
   endif
   
   if (a%npp_stepi(8) /= 0) then
      call AllocateR1pElement(a%Mesh,a%ShockCapturingViscosity,a%Memor)
      
   endif
   
   if (a%kfl_dispa /= 0) then
      call a%Memor%alloc(npoin,a%dissipation,'dissipation','tem_memall')
   endif
   
   if(a%kfl_outfm==1)then   
      call a%Memor%alloc(nbody,a%lun_force,'lun_force','tem_memall')
      call a%Memor%alloc(nbody,a%heatf,'heatf','tem_memall')
   end if

   if (a%kfl_CouplingThreeField==1) then
      call a%Memor%alloc(ndime,ndime,npoin,a%SmoothedVelocityGradient,'SmoothedVelocityGradient','tem_memall')

      if (a%npp_stepi(9) /= 0) then   
         sig_coun = 0
         call a%Memor%alloc(nelem,a%sigmatermarray,'sigmatermarray','tem_memall')
         do ielem=1,nelem
            call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
            allocate(a%sigmatermarray(ielem)%a(pgaus))
            sig_coun = sig_coun + pgaus        
         end do
         call a%Memor%allocObj(0,'sigmatermarray%a','tem_memall',sig_coun*rp)   
         
         call a%Memor%alloc(npoin,a%ViscousDissipation,'ViscousDissipation','tem_memall')
      end if  
   end if
   
   if (a%NumberOfMaterials > 1) then
      call a%Memor%alloc(nelem,a%ElementMaterials,'ElementMaterials','tem_memall')
   endif

   
end subroutine

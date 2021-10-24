subroutine sup_bouope_Velocity(b)
!This subroutine computes the boundary contribution, incompressible Navier-Stokes equations
!-----------------------------------------------------------------------
   use typre
   use Mod_supm_BoundaryConditionsHooksAndSubroutines
   use Mod_NavierStokesBoundary
   implicit none
   class(ThreeFieldNSProblem), target  :: b
   real(rp)    :: density_wall,viscosity_wall
   real(rp)                :: dynpr
   real(rp)                :: dsurf0
   real(rp)                :: coorg,yplus
   integer(ip)             :: igaub,inodb,jnode,jnodb

   a=>b

   !to do multy materials
   imat=1
   
   if(  a%kfl_fixbo(iboun)==2  .or. & !Pressure
      & a%kfl_fixbo(iboun)==3  .or. & !Wall law
      & a%kfl_fixbo(iboun)==5  .or. & !Dynamic pressure
      & a%kfl_fixbo(iboun)==6  .or. & !Open flow
      & a%kfl_fixbo(iboun)==8  .or. & !Atmospheric stress
      & a%kfl_fixbo(iboun)==9  .or. & !Direct traction prescription
      & a%kfl_fixbo(iboun)==10 .or. & !Pressure prescribed in a weak form      
      & a%kfl_fixbo(iboun)==15 .or. & !NRBC \Gamma_O
      & a%kfl_fixbo(iboun)==16 .or. & !NRBC \Gamma_L
      & a%doRobin .eqv. .true.) then  !Robin boundary condition for FSI


      call ProcHook_PreGauss

      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)

      call e%elmdel
      
      !Gather operations
      call e%gatherb(e%ndime,bovel,a%veloc(:,:,1))
    
      call ProcHook_Gathers
      
      dsurf0 = 0.0_rp

      !Gauss-Point Loop
      GaussPointLoop: do igaub=1,e%pgaub
         e%igaub = igaub
   
         call ProcHook_ElmatsToZero

         !Derivatives at the boundary
         call e%elmderb
   
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         call e%interpb(e%ndime,bovel,gpvel(:,1))
         
         call ProcHook_InGauss

         !Velocity norm
         call vecnor(gpvel(:,1),e%ndime,gpvno,2)
         
         !Pressure: sig.n=-p n
         if(a%kfl_fixbo(iboun)==2) then
            call nsm_bouope_pressure

         !Wall law: sig.n=-rho*(U*^2)*u/|u|
         else if(a%kfl_fixbo(iboun)==3) then
            !Physical Parameters
            call a%GetPhysicalParameters(imat,density_wall,viscosity_wall)

            call nsm_bouope_WallLaw

            call nsm_bouwal(e,vel_wall_vec,viscosity_wall,density_wall,delta,wmatr,tract,yplus)
            
            a%ypmin=min(a%ypmin,yplus)
            a%ypmax=max(a%ypmax,yplus)
            a%ypmean=a%ypmean+yplus
            a%nmean_walllaw = a%nmean_walllaw + 1
            
            if(abs(a%grnor) /= 0.0_rp) call nsm_bopres(e,wmatr)
 
         !Dynamic pressure: sig.n=-(-1/2*rho*u^2) n
         else if(a%kfl_fixbo(iboun)==5) then
           dynpr=-0.5_rp*acden*gpvno*gpvno 
           tract(1:e%ndime)=-dynpr*e%baloc(1:e%ndime,e%ndime)
           
         !Open boundary: assemble 2*mu*Sym(grad(u).n + n·v p
         else if(a%kfl_fixbo(iboun)==6) then
            call nsm_bouopb(e,wmatr,acvis,a%fvins)

         !Atmospheric stress:  sig.n = - (p_0 + p_atm) n where p_0 is given
         else if(a%kfl_fixbo(iboun)==8) then
            tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*e%baloc(1:e%ndime,e%ndime)
            coorg = dot_product(gpcod(1:e%ndime),a%gravi(1:e%ndime))
            tract(1:e%ndime) = tract(1:e%ndime) - acden*a%grnor*coorg*e%baloc(1:e%ndime,e%ndime)

         !Prescribed traction: sig.n = tract_0 (tract_0 given)
         else if(a%kfl_fixbo(iboun) == 9) then
            tract(1:e%ndime) = a%bvnat(iboun)%a(1:e%ndime)
            
         !Prescribed Pressure in weak: 
         else if(a%kfl_fixbo(iboun) == 10) then
            call nsm_bouopb_p(e,wmatr,acvis,a%fvins)          
         endif
         
      call ProcHook_PostLoop
      end do GaussPointLoop

      call ProcHook_Finalizations

   endif
       
end subroutine sup_bouope_Velocity

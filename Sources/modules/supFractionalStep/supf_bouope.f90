subroutine supf_bouope(a,iboun,e,elmat,elrhs,ndofn,bcstart)
!This subroutine computes the boundary contribution, incompressible Navier-Stokes equations
!-----------------------------------------------------------------------
   use typre
   use Mod_SupExacso 
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_NavierStokesBoundary
   use Mod_nsm_Viscosity
   implicit none
   class(SUPFractionalStepProblem) :: a
   integer(ip)                :: iboun,bcstart,ndofn
   class(FiniteElement)        :: e
   type(SupExacso)  :: exacso    
   
   real(rp),intent(inout)     :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
   real(rp),intent(inout)     :: elrhs(ndofn,e%mnode)
!    real(rp),intent(inout)     :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)              ! Element matrices
!    real(rp),intent(inout)     :: elrhs(a%ndofn,e%mnode)

   !to do multy materials
   integer(ip) :: imat=1
   
   real(rp), external :: funcre

   real(rp)    :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)
   real(rp)    :: wrhsi(e%ndime+1,e%mnode)
   
   
   real(rp)    :: bovel(e%ndime,e%pnodb),elvel(e%ndime,e%pnode)
   real(rp)    :: gpcod(e%ndime),gpvel(e%ndime),gpvno,grvel(e%ndime,e%ndime)
   
   real(rp)    :: exveg(e%ndime,e%ndime),exvel(e%ndime)
   
   real(rp)    :: acden,acvis,vista,dynpr
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: tract(e%ndime)
   real(rp)    :: updbcn,coorg
   real(rp)    :: bosig((e%ndime-1)*(e%ndime-1)+2,e%pnodb),gpsig((e%ndime-1)*(e%ndime-1)+2)
   
   
   integer(ip) :: igaub,inodb,inode,idime,auxtens


   if(    a%kfl_fixbo(iboun)==2.or.&   !Pressure
      & a%kfl_fixbo(iboun)==3.or.&     !Wall law
      & a%kfl_fixbo(iboun)==5.or.&     !Dynamic pressure
      & a%kfl_fixbo(iboun)==6.or.&     !Open flow
      & a%kfl_fixbo(iboun)==8.or.&     !Atmospheric stress
      & a%kfl_fixbo(iboun)==9 .or.&    !Direct traction prescription
      & a%kfl_fixbo(iboun)==10 .or. &       !Pressure prescribed in a weak form
      & a%kfl_fixbo(iboun)==15 &
      ) then  

      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)
      
      call e%elmdel
      
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      !Gather operations
      call e%gatherb(e%ndime,bovel,a%veloc(:,:,1))
      call e%gatherb(auxtens,bosig,a%sigma(:,:,3))
      
      !Smagorinsky
      if (a%kfl_cotur == -1 .or. a%MatProp(imat)%lawvi /= 0) then
         call e%gather(e%ndime,elvel,a%veloc(:,:,1))
      endif
      
      dsurf0 = 0.0_rp
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         wmatr=0.0_rp
         wrhsi=0.0_rp
         tract=0.0_rp
   
         !Calculate exterior Normal
         call e%bounor
         
         !Derivatives at the boundary
         call e%elmderb         
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         call e%interpb(e%ndime,bovel,gpvel)
         
         call e%interpb(auxtens,bosig,gpsig)
         
         call e%gradientb(e%ndime,elvel,grvel)
         
         !Non-Newtonian viscosity
         if(a%MatProp(imat)%lawvi /= 0)then
            call nsi_vislaw(e,grvel,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis)                 
         end if       
         
         !Smagorinsky
         if(a%kfl_cotur == -1) then
            call nsm_smago(e,grvel,acden,a%turbu(1),vista)
            acvis = acvis + vista
         endif
         
         !Velocity norm
         call vecnor(gpvel,e%ndime,gpvno,2)

         !Pressure: sig.n=-p n
         if(a%kfl_fixbo(iboun)==2) then
            if(a%kfl_funbo(iboun)/=0) then
               updbcn=funcre(a%funpa(a%kfl_funbo(iboun))%a,&
               a%kfl_funty(a%kfl_funbo(iboun),2),&
               a%kfl_funty(a%kfl_funbo(iboun),1),a%bctime)
            else
               updbcn=1.0_rp
            end if
            tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*updbcn*e%baloc(1:e%ndime,e%ndime)
            
            tract(1) = tract(1) + e%baloc(1,e%ndime)*gpsig(1) + e%baloc(2,e%ndime)*gpsig(3)
            tract(2) = tract(2) + e%baloc(1,e%ndime)*gpsig(3) + e%baloc(2,e%ndime)*gpsig(2) 
            
         !Wall law: sig.n=-rho*(U*^2)*u/|u|
         else if(a%kfl_fixbo(iboun)==15) then
         
         !Exact solution
         call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)       
         call exacso%sup_GetVelocity(e%ndime,exvel,exveg) 

         
         tract(1) = exvel(1)
         tract(2) = exvel(2)
 
 
         !Dynamic pressure: sig.n=-(-1/2*rho*u^2) n
         else if(a%kfl_fixbo(iboun)==5) then
           dynpr=-0.5_rp*acden*gpvno*gpvno 
           tract(1:e%ndime)=-dynpr*e%baloc(1:e%ndime,e%ndime)
           
           
         !Open boundary: assemble 2*mu*Sym(grad(u).n + n·v p
         else if(a%kfl_fixbo(iboun)==6) then
            call nsm_bouopb(e,wrhsi,acvis,a%fvins)

 
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


         end if
         
          elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)= dsurf*wmatr + &
               elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)            
        
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,e%ndime
                  elrhs(bcstart+idime,inode)=elrhs(bcstart+idime,inode)+dsurf*(wrhsi(idime,inode)+tract(idime)*e%shapb(inodb,e%igaub))
            end do   
              
        end do
      end do
   endif
       


end subroutine supf_bouope

subroutine sup_bouope(b)
!This subroutine computes the boundary contribution, incompressible Navier-Stokes problem in three field form 
!-----------------------------------------------------------------------
   use typre
   use Mod_SupExacso 
   use Mod_NavierStokes
   use Mod_ThreeField
   use Mod_supm_BaseElmope
   use Mod_supm_BoundaryConditionsHooksAndSubroutines
   implicit none

   class(ThreeFieldNSProblem), target :: b  

   real(rp)    :: dsurf0
   integer(ip) :: igaub,inodb,auxiterv 
   real(rp)    :: dudy,dudz 
   
   interface
      subroutine sup_bouopb(e,wmatr,lawvi,Lawviparam,fvins)
         use typre
         use Mod_Element
         implicit none
         class(FiniteElement)        :: e
         real(rp),    intent(inout) :: wmatr(e%ndime+(e%ndime-1)*(e%ndime-1)+2+1,e%mnode,e%ndime+(e%ndime-1)*(e%ndime-1)+2+1,e%mnode)
         integer(ip), intent(in)    :: lawvi
         real(rp),    intent(in)    :: LawViParam(*),fvins       
      end subroutine
      
      subroutine sup_bouopb_p(e,wmatr,lawvi,Lawviparam,fvins)
         use typre
         use Mod_Element
         implicit none
         class(FiniteElement)        :: e
         real(rp),    intent(inout) :: wmatr(e%ndime+(e%ndime-1)*(e%ndime-1)+2+1,e%mnode,e%ndime+(e%ndime-1)*(e%ndime-1)+2+1,e%mnode)
         integer(ip), intent(in)    :: lawvi
         real(rp),    intent(in)    :: LawViParam(*),fvins      
      end subroutine
      
   end interface
   
   a=>b
   
   auxtens=(e%ndime-1)*(e%ndime-1)+2
   auxL=0.0_rp
   imat=1
    
   if( a%kfl_fixbo(iboun)==11.or. &    !Open Boundary using Stress field
      & a%kfl_fixbo(iboun)==12 .or. &  !Pressure prescribed in a weak form           
      & a%kfl_fixbo(iboun)==13 .or. &
      & a%kfl_fixbo(iboun)==15 .or. &      
      & a%LogFormulation/=0 )then                !Log conformation reformulation

      call ProcHook_PreGauss
      
      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)      
               
      beta=a%MatProp(imat)%LawViParam(2)
      call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)   
 
      call e%elmdel      
      
      call ProcHook_Gathers
      
      dsurf0 = 0.0_rp
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         wmatr=0.0_rp
         tract=0.0_rp
         tractexact = 0.0_rp
   
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf   
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         
         !Derivatives at the boundary
         call e%elmderb 

         !Open boundary: assemble 2*mu*Sym(grad(u).n + n·v p
         if (a%kfl_fixbo(iboun)==11)then
            call sup_bouopb(e,wmatr,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawviParam,a%fvins)
         else if(a%kfl_fixbo(iboun)==12)then
         !Prescribed Pressure in a weak form 
            call sup_bouopb_p(e,wmatr,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawviParam,a%fvins)
         !Wall law: sig.n=-rho*(U*^2)*u/|u|
         else if(a%kfl_fixbo(iboun)==15) then
            !Exact solution
            call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)       
            call exacso%sup_GetStress(e%ndime,a%LogFormulation,exsig,exsigr)  

            tractexact(1) = exsig(1)
            tractexact(2) = exsig(2)               
            tractexact(3) = exsig(3)  
          
         else if(a%kfl_fixbo(iboun) == 13)then
         
            !incremental scheme for Weissenberg number
            call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)           

            !n*sigma in the outflow boundary (Neumann Boundary condition) 
            !t1=n1*s11+n2*s12+n3*s13
            !t2=n1*s12+n2*s22+n3*s23
            !t3=n1*s13+n2*s23+n3*s33
            
            !Contracción 4:1
            !tract(1)=2.0_rp*lambda*(1.0_rp-beta)*acvis*(3.0_rp*gpcod(2))**2.0_rp
            !tract(2)=(1.0_rp-beta)*acvis*(-3.0_rp*gpcod(2))            
            
            !canal recto
            tract(1)=2.0_rp*lambda*(1.0_rp-beta)*acvis*(8.0_rp*gpcod(2))**2.0_rp
            tract(2)=(1.0_rp-beta)*acvis*(-8.0_rp*gpcod(2))
            
            !Cilindro
            !tract(1)=2.0_rp*lambda*(1.0_rp-beta)*acvis*(0.75_rp*gpcod(2))**2.0_rp
            !tract(2)=-3.0_rp*(1.0_rp-beta)*acvis*(-0.25_rp*gpcod(2))        
            
            !3d case 
            if(e%ndime==3)then             
               !4:1 contraction
            
               dudy=(3.0_rp)*gpcod(2)*((gpcod(3))**2 - 1.0_rp)
               dudz=(3.0_rp)*gpcod(3)*((gpcod(2))**2 - 1.0_rp)

               tract(1) = 2.0_rp*lambda*(1.0_rp-beta)*acvis*(dudy**2 + dudz**2)
               tract(2) = (1.0_rp-beta)*acvis*(dudy)
               tract(3) = (1.0_rp-beta)*acvis*(dudz)                        
            end if
         end if 
         
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,e%ndime
                  elrhs(bcstart+idime,inode)=elrhs(bcstart+idime,inode)+dsurf*((tract(idime))*e%shapb(inodb,e%igaub))
            end do   
        end do

         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,auxtens
                  elrhs(idime,inode)=elrhs(idime,inode)+dsurf*(tractexact(idime)*e%shapb(inodb,e%igaub))
            end do   
         end do        

      end do
      
      call ProcHook_Finalizations
   endif   


end subroutine sup_bouope

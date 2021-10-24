subroutine lmn_bouope(a,iboun,e,elmat,elrhs)
   use typre
   use Mod_Element
   use Mod_LowMach
   implicit none
   class(LowMachProblem)  :: a
   integer(ip)            :: iboun
   class(FiniteElement)    :: e
   real(rp),intent(inout) :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)              ! Element matrices
   real(rp),intent(inout) :: elrhs(a%ndofn,e%mnode)
   real(rp), external     :: funcre
   real(rp)    :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode)
   real(rp)    :: bovel(e%ndime,e%pnodb),botem(1,e%pnodb)
   real(rp)    :: gpcod(e%ndime),gpvel(e%ndime),gpvno
   real(rp)    :: acden,acvis,dynpr,gptem(1),gpden(1)
   real(rp)    :: dsurf,dsurf0,tract(e%ndime),updbcn,coorg
   integer(ip) :: igaub,inodb,inode,idime

   interface
      subroutine lmn_bouopb(a,e,wmatr,acvis,fvins)
      use typre
      use Mod_LowMach
      use Mod_Element
      implicit none
      class(LowMachProblem)      :: a
      class(FiniteElement)        :: e
      real(rp),    intent(inout) :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode)
      real(rp),    intent(in)    :: acvis,fvins
      end subroutine

      subroutine lmn_boucon(a,e,wmatr,acden)
         use typre
         use Mod_LowMach
         use Mod_Element
         implicit none
         class(LowMachProblem)      :: a
         class(FiniteElement)        :: e
         real(rp),    intent(inout) :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode)
         real(rp),    intent(in)    :: acden
      end subroutine
   end interface

      !Physical Parameters
      call a%GetPhysicalParameters(acvis=acvis)
      call e%elmdel
      
      !Gather operations
      call e%gatherb(e%ndime,bovel,a%veloc(:,:,1))
      call e%gatherb(1,botem,a%tempe(:,1))
      
      dsurf0 = 0.0_rp
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         wmatr=0.0_rp
         tract=0.0_rp
         !Interpolate density values 
         call e%interpb(1,botem,gptem)

         if (a%kfl_nolsg == 0) then
            gpden = a%pther(1)/a%sgasc/gptem
         else
            gpden = a%pther(1)/a%sgasc/(gptem + a%tesgs(e%lboel(e%pnodb+1))%a(1,e%igaub))
         end if
         acden = gpden(1)

         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         call e%interpb(e%ndime,bovel,gpvel)
         
         !Velocity norm
         call vecnor(gpvel,e%ndime,gpvno,2)
         
         !Derivatives at the boundary
         call e%elmderb
         
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
            
            
         !Wall law: sig.n=-rho*(U*^2)*u/|u|
         else if(a%kfl_fixbo(iboun)==3) then
            !nwall=nwall+1
            call lmn_bouwal(e,gpvel,acvis,acden,a%delta,wmatr,tract)

            if(abs(a%grnor) /= 0.0_rp) call lmn_bopres(e,wmatr)
 
 
         !Dynamic pressure: sig.n=-(-1/2*rho*u^2) n
         else if(a%kfl_fixbo(iboun)==5) then
           dynpr=-0.5_rp*acden*gpvno*gpvno 
           tract(1:e%ndime)=-dynpr*e%baloc(1:e%ndime,e%ndime)
           
           
         !Open boundary: assemble 2*mu*Sym(grad(u).n + n·v p
         else if(a%kfl_fixbo(iboun)==6) then
            call lmn_bouopb(a,e,wmatr,acvis,a%fvins)

 
         !Atmospheric stress:  sig.n = - (p_0 + p_atm) n where p_0 is given
         else if(a%kfl_fixbo(iboun)==8) then

            tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*e%baloc(1:e%ndime,e%ndime)
            coorg = dot_product(gpcod(1:e%ndime),a%gravi(1:e%ndime))
            tract(1:e%ndime) = tract(1:e%ndime) - acden*a%grnor*coorg*e%baloc(1:e%ndime,e%ndime)
            
         !Prescribed traction: sig.n = tract_0 (tract_0 given)
         else if(a%kfl_fixbo(iboun) == 9) then
            tract(1:e%ndime) = a%bvnat(iboun)%a(1:e%ndime)
         end if

         !Boundary term in the continuity eq  
         call lmn_boucon(a,e,wmatr,acden)

         elmat=elmat+dsurf*wmatr
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,e%ndime
                  elrhs(idime,inode)=elrhs(idime,inode)+dsurf*tract(idime)*e%shapb(inodb,e%igaub)
            end do
         end do
 
      end do 

end subroutine lmn_bouope

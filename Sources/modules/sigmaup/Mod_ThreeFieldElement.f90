module Mod_ThreeFieldElement
   use typre
   use Mod_Element
   use Mod_ThreeFieldElement_tools
   use Mod_ThreeField
   use Mod_supm_StressGenerator
   
contains

   ! Elemental components in the Three Field case
   ! Voigt notation 2d=>  s11 s22 s12  3d => s11 s22 s33 s23 s13 s12
   
  !*****************************************************************************************************************************
  !Three Field case no-elastic 
  !*****************************************************************************************************************************
   
   subroutine supm_elmbst(e,timom,tisig,dvolu,acvis,auxtens,auxoss,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: acvis,timom,tisig,dvolu
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      real(rp)                   :: mat(auxtens,e%mnode,auxtens,e%mnode)
      real(rp)                   :: aux,aux2,aux1

      aux1=1.0_rp/(2.0_rp*acvis)
      aux= (aux1*(1.0_rp-tisig*aux1*auxoss)*dvolu) !P_h(sigma_h)=0 !
      aux2=timom*dvolu

      call elmbst(e,e%ndime,auxtens,aux,aux2,mat)

      elmat = elmat + mat

   end subroutine supm_elmbst
   
   

   subroutine supm_elmbut(e,timom,tisig,dtinv,gpadv,dvolu,acden,acvis,auxtens,auxoss,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! -(gra_sym(u),T) +tau3(gra_sym(u),1/2mu*T) -tau1(Div(T),rho*a*grad(u)+rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: tisig,timom,dvolu,acvis,acden,dtinv,gpadv(e%pnode)  
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, itens, jdime
      real(rp)                   :: aux,aux1,tmp,aux2, tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)
      real(rp)                   :: mat(auxtens,e%mnode,e%ndime,e%mnode)
      
      
      tmp  = acden*dtinv
      aux2 = -timom*dvolu*auxoss
      aux = -(1.0_rp - (tisig/(2.0_rp*acvis))*auxoss)*dvolu
      
      call elmbut(e,e%ndime,auxtens,tmp,aux,aux2,acden,gpadv,mat)  

      elmat = elmat + mat
   
   end subroutine supm_elmbut
   
   
   subroutine supm_elmbpt(e,timom,dvolu,auxtens,elmat)
    !---------------------------------------------------------------------------
    !
    ! This routine computes the third lhs term for ASGS in constitutive equation
    !    -tau1*(gra(p),div(t))
    !
    !---------------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode, itens, jtens
      real(rp)                   :: aux, vector(auxtens), tensor(e%ndime,e%ndime)
      real(rp)                   :: mat(auxtens,e%mnode,1,e%mnode)

      aux=-timom*dvolu    

      call elmbpt(e,e%ndime,auxtens,aux,mat)

      elmat = elmat + mat

  end subroutine supm_elmbpt 
   
   
   subroutine supm_elmbsv(e,timom,tisig,dvolu,acvis,acden,gpadv,auxtens,auxoss,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in momentum equation
    !    (gra_sym(v),S) -tau3(gra_sym(v),1/2mu*S) -tau1(rho*a*grad(v),div(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: tisig,timom,dvolu,acvis,acden,gpadv(e%pnode) 
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode,jtens,idime
      real(rp)                   :: aux,aux1, tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)
      real(rp)                   :: mat(e%ndime,e%mnode,auxtens,e%mnode)

      aux  = (1.0_rp - (tisig/(2.0_rp*acvis))*auxoss)*dvolu 

      call elmbsv(e,e%ndime,auxtens,aux,timom,dvolu,acden,gpadv,mat)

      elmat = elmat + mat
 
   end subroutine supm_elmbsv
   
   subroutine supm_elmbuv2(e,tisig,dvolu,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in momentum equation
    !    tau3*(gra_sym(v),gra_sym(u))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: tisig,dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux
      real(rp)                   :: mat(e%ndime,e%mnode,e%ndime,e%mnode)

      aux = (tisig*dvolu) !aporte Newton
      
      call elmbuv2(e,e%ndime,aux,mat) 

      elmat = elmat + mat
 
   end subroutine supm_elmbuv2
   

   subroutine supm_elmbsq(e,timom,dvolu,auxtens,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs term for ASGS in momentum equation
    !    tau1*(gra(q),div(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
 
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elmat(1,e%mnode,auxtens,e%mnode)
      real(rp)                   :: aux
      real(rp)                   :: mat(1,e%mnode,auxtens,e%mnode)

      aux=-timom*dvolu   

      call elmbsq(e,e%ndime,auxtens,aux,mat)

      elmat = elmat + mat

  end subroutine supm_elmbsq 
   

   subroutine supm_elmrhc(e,timom,tisig,acvis,dvolu,elext,auxtens,elextS,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: elext(e%ndime),elextS(auxtens)
      real(rp),    intent(in)    :: dvolu,timom,tisig,acvis
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      real(rp)                   :: aux
      real(rp)                   :: rhs(auxtens,e%mnode)

      aux  = -timom*dvolu

      call elmrhc(e,e%ndime,auxtens,aux,tisig,acvis,dvolu,elext,elextS,rhs)

      elrhs = elrhs + rhs

   end subroutine supm_elmrhc
     
   subroutine supm_elmrhu(e,tidiv,tisig,auxtens,dvolu,elextC,elextS,elrhs)  
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: elextC(1),elextS(auxtens)
      real(rp),    intent(in)    :: dvolu,tidiv,tisig
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux,aux1
      real(rp)                   :: rhs(e%ndime,e%mnode)

      aux = tidiv*dvolu 
      aux1= tisig*dvolu

      call elmrhu(e,e%ndime,auxtens,aux,aux1,elextC,elextS,rhs)  

      elrhs = elrhs + rhs

   end subroutine supm_elmrhu

   subroutine supm_elmrhp(e,dvolu,elextC,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elextC(1)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,aux1
      real(rp)                   :: aux2

    
      do inode=1,e%pnode
         elrhs(1,inode) = e%shape(inode,e%igaus)*dvolu*elextC(1) &
               + elrhs(1,inode) 
      end do

   end subroutine supm_elmrhp

   !*****************************************************************************************************************************
   !Three Field viscoelastic 
   !*****************************************************************************************************************************
   
   subroutine sup_TimeIntegrationToElext(e,auxtens,Integrator,auxVE,dtinv,gpsig,elextS)
      use typre
      use Mod_Element
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      integer(ip), intent(in)    :: auxtens
      real(rp)                   :: gpsig(auxtens,*),elextS(auxtens)
      real(rp)                   :: auxVE,dtinv      
      real(rp)                   :: gprhs(auxtens)
      
      !Time integration
      call Integrator%GetRHS(auxtens,gpsig(:,2),gprhs)
      elextS = elextS + auxVE*dtinv*gprhs(1:auxtens)
      
      
   end subroutine sup_TimeIntegrationToElext
   
   !************************************************************************************************************************************
   !Galerkin terms 

   subroutine supm_elmbstGal(e,acvis,auxVE,auxG,auxPTT,dtinv,dvolu,auxtens,gpadv,grvel,gpsig,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmst matrix
    ! 1/2mu(S,T)+auxVE*(a*grad(S) - S*grad(a) - (grad(a))'*S,T) + ...
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,dtinv
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime),gpsig(auxtens)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens, jtens, itens2
      real(rp)                   :: aux1, aux
      real(rp)                   :: traza, tensor(e%ndime,e%ndime), tens_vel(auxtens,auxtens)   
      real(rp)                   :: vector(auxtens), value, tens_sig(auxtens,auxtens),value2
      
      aux=1.0_rp/(2.0_rp*acvis)
 
       do jnode=1,e%pnode   
          do inode=1,e%pnode
            aux1=e%shape(inode,e%igaus)*dvolu
               do itens=1, e%ndime
                  tensor(itens,:)= -2.0_rp*auxVE*grvel(itens,:)*e%shape(jnode,e%igaus)
               end do
               call GetStress_tt(e%ndime,tensor,tens_vel)
               vector(:)= auxVE*(auxG/acvis)*gpsig(:)*e%shape(jnode,e%igaus)
               call GetAdvStress_tt(e%ndime,auxPTT,vector,tens_sig)
        
               value=e%shape(jnode,e%igaus)*aux + auxVE*gpadv(jnode) + auxVE*dtinv*e%shape(jnode,e%igaus)
               call AddSphericalComponent_tt(e%ndime,value,tens_vel)
  
               elmat(1:auxtens,inode,:,jnode) = aux1*(tens_vel(1:auxtens,:)+tens_sig(1:auxtens,:)) + elmat(1:auxtens,inode,:,jnode)
         end do
      end do   

   end subroutine supm_elmbstGal   
   
   subroutine supm_elmbutGal(e,auxVE,dvolu,auxtens,beta,grvel,gpsig,grsig,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmut matrix 
    ! - ((1 - beta)* grad_sym(u), T) + auxVE*(u*grad(S) - S*grad(u) - (grad(u))'*S... ,T)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,beta,grvel(e%ndime,e%ndime),gpsig(auxtens),grsig(auxtens,e%ndime),auxVE  
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,itens
      real(rp)                   :: aux
      real(rp)                   :: grad_sig(auxtens,e%ndime), vector(e%ndime,e%ndime), grad_vel(auxtens,e%ndime)
      real(rp)                   :: tensor(auxtens,e%ndime)
      do jnode=1,e%pnode
         do inode=1,e%pnode
            aux=e%shape(inode,e%igaus)*dvolu  
        
            do itens=1, e%ndime
               tensor(:,itens)= -2.0_rp*auxVE*gpsig(:)*e%cartd(itens,jnode)
               vector(itens,itens) = -(1.0_rp - beta)*e%cartd(itens,jnode)
            end do
            grad_sig(1:auxtens,:)= auxVE*e%shape(jnode,e%igaus)*grsig(1:auxtens,:)
            call GetStress_ut(e%ndime,vector,grad_vel)
            call AddSigmaGradVel_ut(e%ndime,tensor,grad_sig)

            elmat(1:auxtens,inode,:,jnode) = aux*(grad_vel(1:auxtens,:)+grad_sig(1:auxtens,:)) + elmat(1:auxtens,inode,:,jnode)            
                                                
         end do  
      end do

   end subroutine supm_elmbutGal
   
   
   subroutine supm_elmrhcGal(e,auxVE,auxG,auxPTT,acvis,dvolu,auxtens,elextS,gpvel,grvel,gpsig,grsig,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmrhc matrix.
    !  
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: elextS(auxtens),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,gpsig(auxtens),grsig(auxtens,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, itens, jtens
      real(rp)                   :: aux1, value(auxtens)
      real(rp)                   :: tensor(e%ndime,auxtens),vect_GrSig(auxtens), tensGrVel(e%ndime,e%ndime)
      real(rp)                   :: tensor2(auxtens, auxtens), vect_GrVel(auxtens), vect_Sig(auxtens)
      
      
      !product GrSig*Vel 
      do itens=1, e%ndime
         tensor(itens,:) = auxVE*gpvel(itens)*grsig(:,itens)
         tensGrVel(itens,:)=-2.0_rp*auxVE*grvel(itens,:)
      end do    
      
      call GetGrSig_hc(e%ndime,tensor,vect_GrSig) 
      
      !product GrVel*Sig
      
      call GetGrVel_hc(e%ndime,tensGrVel,gpsig,vect_GrVel) 
      
      !product Sig*Sig
      do jtens=1, auxtens
         tensor2(jtens,:) = auxVE*(auxG/acvis)*gpsig(jtens)*gpsig(:)
      end do  
      call GetStress_hc(e%ndime,auxPTT,tensor2,vect_sig)
      
      value(:)=elextS(:)
      call AddValues_hc(e%ndime,value,vect_sig)
        
      do inode=1,e%pnode
         aux1=e%shape(inode,e%igaus)*dvolu
         elrhs(:,inode) = aux1*(vect_GrSig(:)+ vect_GrVel(:) + vect_sig(:)) + elrhs(:,inode) 
      end do

   end subroutine supm_elmrhcGal

   subroutine supm_elmbsvGal(e,dvolu,auxtens,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmbsv matrix.
    !  (S,gra_sym(v))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, idime, jtens
      real(rp)                   :: aux, tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)

      aux  = dvolu 
      do inode=1,e%pnode       
         do jnode=1,e%pnode
         
            do idime=1, e%ndime
               tensor(idime,idime)=e%cartd(idime,inode)*aux*e%shape(jnode,e%igaus)
            end do
            call GetStress_ut(e%ndime,tensor,tens_sol)

            do idime=1, e%ndime
               elmat(idime,inode,:,jnode) = tens_sol(:,idime) + elmat(idime,inode,:,jnode)
            end do
            
         end do  
      end do
 
   end subroutine supm_elmbsvGal
   
   
   subroutine supm_elmbuvGal(e,dvolu,acvis,beta,acden,gpadv,grvel,dtinv,linearization,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmbuv matrix.
    !    (2*beta*mu*gra_sym(u),gra_sym(v)) - (rho*a*gra(u), v)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,acden,gpadv(e%pnode),grvel(e%ndime,e%ndime),dtinv,acvis,beta
      integer(ip), intent(in)    :: linearization
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, idime, jdime
      real(rp)                   :: aux,aux1
      real(rp)                   :: tensor(e%ndime,e%ndime), value(e%ndime), cont

      aux  = dvolu*acden
      aux1 = dvolu*acvis*beta

      do inode=1,e%pnode       
         do jnode=1,e%pnode
            tensor=0.0_rp
            cont=0
            do idime=1,e%ndime
               tensor(idime,:)=aux*e%shape(inode,e%igaus)*(e%shape(jnode,e%igaus)*grvel(idime,:))*linearization 
               cont=cont+aux1*(e%cartd(idime,inode)*e%cartd(idime,jnode))
            end do
                 
            value(:)=aux*e%shape(inode,e%igaus)*(dtinv*e%shape(jnode,e%igaus)+ gpadv(jnode)) + cont   
            call AddSphericalComponent_uv(e%ndime,value,tensor)
            
            do idime=1,e%ndime
               elmat(idime,inode,:,jnode) =  tensor(idime,:)+ elmat(idime,inode,:,jnode)
            end do

         end do  
      end do
 
   end subroutine supm_elmbuvGal
   
   
   subroutine supm_elmbpvGal(e,dvolu,LCR,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmbpv matrix.
    ! -(p,div v) 
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,LCR
      real(rp)                   :: aux 
          
      do jnode=1,e%pnode
      
         aux=-dvolu*e%shape(jnode,e%igaus)    
         do inode=1,e%pnode
            elmat(:,inode,1,jnode) = aux*e%cartd(:,inode) + elmat(:,inode,1,jnode)  
         end do
         
      end do

   end subroutine supm_elmbpvGal 


   subroutine supm_elmrhuGal(e,acden,dvolu,elext,gpvel,grvel,linearization,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms in momentum equation
    !    
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime)
      integer(ip), intent(in)    :: linearization
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode, idime, Newton
      real(rp)                   :: aux, cont(e%ndime)

      cont=0.0_rp
      do inode=1,e%pnode
         aux = e%shape(inode,e%igaus)*dvolu

         if (linearization==1) then
            cont=0.0_rp
            do idime=1,e%ndime
               cont(:)=acden*(gpvel(idime)*grvel(:,idime))+cont(:)
            end do
         end if   
         
         elrhs(:,inode) = aux*(elext(:) + cont(:)) + elrhs(:,inode) 
      end do

   end subroutine supm_elmrhuGal 
   
   subroutine supm_elmbuqGal(e,dvolu,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the term in continuity equation
    ! (div(u),q)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux
      
      do inode=1,e%pnode  
         aux= e%shape(inode,e%igaus)*dvolu      
         do jnode=1,e%pnode                       
            elmat(1,inode,:,jnode) = aux*e%cartd(:,jnode) + elmat(1,inode,:,jnode)        
         end do
      end do

   end subroutine supm_elmbuqGal
   
   subroutine supm_elmbpqpena(e,acvis,epsi,dvolu,elmat)
    !-----------------------------------------------------------------------
    !
    !  (p,q)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: dvolu,acvis,epsi
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,k,aux1
      real(rp)                   :: aux,epsi_original

      k=6
      
      epsi_original = (10.0_rp**(-k))/acvis
      aux  = epsi_original*dvolu    
      
      do jnode=1,e%pnode  
         do inode=1,e%pnode    
         elmat(1,inode,1,jnode) = e%shape(inode,e%igaus)*aux*e%shape(jnode,e%igaus) + elmat(1,inode,1,jnode)
         end do  
      end do

   end subroutine supm_elmbpqpena    
   
   subroutine supm_elmrhpGal(e,dvolu,elextC,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms in continuity equation
    !  
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elextC(1)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode
    
      do inode=1,e%pnode
         elrhs(1,inode) = e%shape(inode,e%igaus)*dvolu*elextC(1) + elrhs(1,inode) 
      end do

   end subroutine supm_elmrhpGal     
   
   
 
   !**********************************************************************************************************************************
   !Estabilization terms of Momentum and Continuity equation
   subroutine supm_elmbstEst1(e,beta,timom,dvolu,auxtens,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes:
    ! tau1*(div(T),div(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu,beta
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens, jtens
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,e%ndime), tens_sol(auxtens,auxtens)
      
      aux= timom*dvolu*(1.0_rp-beta)      
      
      do jnode=1,e%pnode    
         do inode=1,e%pnode
         
            do itens=1,e%ndime
               tensor(itens,:)=e%cartd(itens,inode)*e%cartd(:,jnode)*aux 
            end do
            
            call GetStress_tt(e%ndime,tensor,tens_sol)
            
            do itens=1,auxtens
               elmat(itens,inode,:,jnode)= tens_sol(itens,:)+ elmat(itens,inode,:,jnode) 
            end do

         end do
      end do   

   end subroutine supm_elmbstEst1 
   
   
   subroutine supm_elmbutEst1(e,beta,acden,timom,dtinv,dvolu,auxtens,gpadv,grvel,auxoss,linearization,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! -tau1(Div(T),rho*a*grad(u) + rho*U*grad(a) +rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss,linearization
      real(rp),    intent(in)    :: timom,dvolu,acden,dtinv,gpadv(e%pnode),beta,grvel(e%ndime,e%ndime)  
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,itens

      real(rp)                   :: aux1,tmp,aux,aux2  
      real(rp)                   :: tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)
      real(rp)                   :: vector(e%ndime), tens_grvel(auxtens,e%ndime)

      tmp  = acden*dtinv
      aux= -timom*dvolu*(1.0_rp-beta)*auxoss       
      aux2= aux*acden
      
      do jnode=1,e%pnode
         aux1=acden*gpadv(jnode)         
         do inode=1,e%pnode
            tens_grvel=0.0_rp
            tens_sol=0.0_rp
         
            do itens=1, e%ndime 
               tensor(itens,:)=  aux*e%cartd(itens,inode)*aux1 + aux*e%cartd(itens,inode)*e%shape(jnode,e%igaus)*tmp 
            end do  
         
            call GetStress_ut(e%ndime,tensor,tens_sol)
            
            vector(:)=aux2*e%shape(jnode,e%igaus)*e%cartd(:,inode)
            
            if (linearization==1) then
               call GetGrVel_ut(e%ndime,vector,grvel,tens_grvel)
            end if   
               
            elmat(1:auxtens,inode,:,jnode) = tens_sol(1:auxtens,:) &  
                  + tens_grvel(1:auxtens,:) &
                  + elmat(1:auxtens,inode,:,jnode)                                                                       

         end do  
      end do

   end subroutine supm_elmbutEst1
   
   subroutine supm_elmbutEst1Lapla(e,beta,acvis,timom,dvolu,auxtens,auxoss,elmat) 
    !-----------------------------------------------------------------------
    ! tau*(Div(t),-2*beta*nu*Lapla(u))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: dvolu,acvis,timom,beta
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,itens, jtens
      real(rp)                   :: aux,auxvis    
      real(rp)                   :: tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)
      

      aux = timom*dvolu*(1.0_rp-beta)*auxoss
      
      do jnode=1,e%pnode    
         do inode=1,e%pnode 
            auxvis= (1.0_rp)*(beta*acvis)*sum(e%hessi(1:e%ndime,jnode))
            do itens=1, e%ndime
               do jtens=1, e%ndime
                  tensor(itens,jtens)=aux*e%cartd(itens,inode)*(auxvis)
               end do
            end do
            
            call GetStress_ut(e%ndime,tensor,tens_sol)
         
            elmat(1:auxtens,inode,:,jnode) = tens_sol(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode) 
         end do  
      end do     
      
   end subroutine supm_elmbutEst1Lapla
    
   subroutine supm_elmbptEst1(e,timom,dvolu,auxtens,beta,auxoss,elmat)
    !-----------------------------------------------------------------------
    !
    !    -tau1*(gra(p),div(t))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: timom,dvolu,beta
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,itens
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,e%ndime), vect_sol(auxtens)

      aux= -timom*dvolu*(1.0_rp-beta)*auxoss    
      do jnode=1,e%pnode  
         do inode=1,e%pnode  
            do itens=1, e%ndime
               tensor(itens,:)=e%cartd(itens,inode)*aux*e%cartd(:,jnode)
            end do
            
            call GetStress(e%ndime,tensor,vect_sol)
         
            elmat(:,inode,1,jnode) = vect_sol(:) + elmat(:,inode,1,jnode)
      end do  
   end do

   end subroutine supm_elmbptEst1 
   
   subroutine supm_elmrhcEst1(e,beta,acden,timom,dvolu,elext,auxtens,gpvel,grvel,auxoss,linearization,elrhs)
    !-----------------------------------------------------------------------
    !
    !    -tau1*(div(T), f) -tau1(div(T),rho*u*grad(u))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss,linearization
      real(rp),    intent(in)    :: elext(e%ndime),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,timom,beta,acden
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,itens
      real(rp)                   :: aux,aux1
      real(rp)                   :: tensor(e%ndime,e%ndime), vect_elext(auxtens)
      real(rp)                   :: tens_grvel(auxtens,e%ndime), grvel_trans(e%ndime,e%ndime)
      real(rp)                   :: vector(e%ndime), vect_grvel(auxtens), tensor2(e%ndime,e%ndime)
      real(rp)                   :: tens_grvel_tras(e%ndime,e%ndime)

      
      aux  = -timom*dvolu*(1.0_rp-beta)*auxoss
      aux1 = aux*acden

      
      do inode=1,e%pnode  
         vect_elext=0.0_rp
         vect_grvel=0.0_rp
      
         do itens=1, e%ndime
            tensor(itens,:)=e%cartd(itens,inode)*elext(:)*aux
            grvel_trans(itens,:)=grvel(:,itens)
         end do   
         call GetStress(e%ndime,tensor,vect_elext)
         
         if (linearization==1) then
            call GetGrVel_ut(e%ndime,gpvel,grvel_trans,tens_grvel)
            do itens=1, e%ndime
               tens_grvel_tras(:,itens)=tens_grvel(itens,:)
            end do   
            call GetGrVel_hu(e%ndime,tens_grvel_tras,vector)
         
            do itens=1, e%ndime
               tensor2(itens,:)=aux1*e%cartd(itens,inode)*vector(:) 
            end do
            call GetStress(e%ndime,tensor2,vect_grvel)
         end if
      
         elrhs(:,inode) = vect_elext(:) + vect_grvel(:) + elrhs(:,inode) 
      
      end do

   end subroutine supm_elmrhcEst1
   
   subroutine supm_elmbsvEst1(e,timom,dvolu,acden,gpadv,auxtens,auxoss,elmat)
    !-----------------------------------------------------------------------
    !
    !    tau1*(-grad(S) , rho*a*grad(V))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: timom,dvolu,acden,gpadv(e%pnode) 
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime), tens_sol_t(e%ndime,auxtens)

      do inode=1,e%pnode
         aux = -acden*gpadv(inode)*timom*dvolu*auxoss        
         do jnode=1,e%pnode
         
            do itens=1,e%ndime
               tensor(itens,:)=e%cartd(itens,jnode)*aux
            end do
         
            call GetStress_ut(e%ndime,tensor,tens_sol)
            tens_sol_t=transpose(tens_sol)
            elmat(1:e%ndime,inode,:,jnode) = tens_sol_t(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)                                   
      
         end do  
      end do
 
   end subroutine supm_elmbsvEst1
  
   subroutine supm_elmbsvEst1Lapla(e,beta,acvis,timom,dvolu,auxtens,auxoss,elmat)   
    !-----------------------------------------------------------------------
    !
    !    tau1*(-grad(S) , 2*beta*eta*lapla(V))
    !   
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: acvis,dvolu,timom,beta    
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens
      real(rp)                   :: aux,auxvis     
      real(rp)                   :: tensor(e%ndime,e%ndime), tens_sol(auxtens,e%ndime)
      real(rp)                   :: tens(e%ndime,auxtens)
   

      aux = (-1.0_rp)*(beta*acvis)*timom*dvolu*auxoss
      do inode=1,e%pnode
         auxvis= sum(e%hessi(1:e%ndime,inode))*aux
         do jnode=1,e%pnode
         
            do itens=1,e%ndime
               tensor(itens,:)=auxvis*e%cartd(itens,jnode)
            end do
            
            call GetStress_ut(e%ndime,tensor,tens_sol)
            tens=transpose(tens_sol)
            elmat(1:e%ndime,inode,:,jnode) = tens(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)  
          
         end do
      end do   

   end subroutine supm_elmbsvEst1Lapla  
   
   subroutine supm_elmbuvEst1(e,dvolu,beta,acden,timom,tidiv,gpadv,grvel,dtinv,auxoss,linearization,elmat)
    !-----------------------------------------------------------------------
    !
    !  tau1*(rho*a*grad(U) + rho*U*grad(a) + dU/dt, rho*a*grad(V))  + tau2(Div(U),Div(V))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxoss,linearization
      real(rp),    intent(in)    :: dvolu,acden,gpadv(e%pnode),grvel(e%ndime,e%ndime),dtinv,timom,tidiv,beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, itens
      real(rp)                   :: aux,aux1
      real(rp)                   :: value(e%ndime), tensor(e%ndime,e%ndime)

      
      aux1 = tidiv*dvolu
      do inode=1,e%pnode       
         aux  = (dvolu*timom)*(acden*gpadv(inode))
         do jnode=1,e%pnode
         
            do itens=1,e%ndime
               tensor(itens,:)= aux*acden*(e%shape(jnode,e%igaus)*grvel(itens,:))*linearization + e%cartd(itens,inode)*aux1*e%cartd(:,jnode)
            end do   
            value(:)=aux*acden*(dtinv*e%shape(jnode,e%igaus)*auxoss + gpadv(jnode))
            
            call AddSphericalComponent_uv(e%ndime,value,tensor)
            
            elmat(1:e%ndime,inode,:,jnode)=tensor(1:e%ndime,:)+ elmat(1:e%ndime,inode,:,jnode) 

         end do  
      end do
 
   end subroutine supm_elmbuvEst1
   
   subroutine supm_elmbuvEst1Lapla(e,beta,acvis,acden,timom,dvolu,auxtens,gpadv,grvel,dtinv,auxoss,elmat)   
    !-----------------------------------------------------------------------
    ! tau1*( 2*beta*eta*lapla(U) + rho*a*grad(U)+ rho*U*grad(a) + dU/dt, 2*beta*eta*lapla(V))
    ! tau1*( 2*beta*eta*lapla(U), rho*a*grad(V))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: dvolu,gpadv(e%pnode)
      real(rp),    intent(in)    :: acvis,beta,timom,dtinv,acden,grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, itens
      real(rp)                   :: aux,aux2,aux3,aux4   
       real(rp)                   :: value(e%ndime), tensor(e%ndime,e%ndime)
      
      aux=(beta*acvis)*(timom*dvolu)
      
      do inode=1,e%pnode
         aux3= aux*(sum(e%hessi(1:e%ndime,inode)))*auxoss
         aux4= (acden*timom)*gpadv(inode)*dvolu*auxoss
         do jnode=1,e%pnode         
            aux2 = -(beta*acvis)*sum(e%hessi(1:e%ndime,jnode))  
            
             do itens=1,e%ndime
                tensor(itens,:)= aux3*auxoss*acden*(e%shape(jnode,e%igaus)*grvel(itens,:))
             end do 
             value(:)=aux3*(auxoss*acden*(e%shape(jnode,e%igaus)*dtinv + gpadv(jnode))+ aux2) + aux4*aux2 
             call AddSphericalComponent_uv(e%ndime,value,tensor)
             
             elmat(1:e%ndime,inode,:,jnode)=tensor(1:e%ndime,:)+ elmat(1:e%ndime,inode,:,jnode) 
         end do  
      end do


   end subroutine supm_elmbuvEst1Lapla

   
   subroutine supm_elmbpvEst1(e,acden,timom,gpadv,dvolu,elmat) 
    !-----------------------------------------------------------------------
    !
    ! tau1*(grad(P), rho*a*grad(V))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,timom,acden,gpadv(e%pnode)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1 
      
      do inode=1,e%pnode      
         aux1 = (timom*dvolu)*(acden*gpadv(inode))      
         do jnode=1,e%pnode                    
            elmat(:,inode,1,jnode) = aux1*e%cartd(:,jnode) + elmat(:,inode,1,jnode)        
         end do
      end do

   end subroutine supm_elmbpvEst1
    
   subroutine supm_elmbpvEst1Lapla(e,beta,acvis,timom,dvolu,auxoss,elmat) 
    !-----------------------------------------------------------------------
    !
    ! tau1* (grad(P), 2*beta*eta*lapla(V))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: dvolu,timom,acvis,beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux,aux1 

      aux= (beta*acvis)*(dvolu*timom)*auxoss
      
      do inode=1,e%pnode      
         aux1= aux*sum(e%hessi(1:e%ndime,inode))     
         do jnode=1,e%pnode                       
            elmat(:,inode,1,jnode) = aux1*e%cartd(:,jnode) + elmat(:,inode,1,jnode)      
         end do
      end do

   end subroutine supm_elmbpvEst1Lapla
   
   subroutine supm_elmrhuEst1(e,timom,tidiv,acden,auxtens,dvolu,elextC,elext,gpvel,grvel,gpadv,auxoss,linearization,elrhs)
    !-----------------------------------------------------------------------
    !
    !  tau1* (f, rho*a*grad(V) + rho*U*grad(v)) + tau2* (grad(U),grad(V))
    !    
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss,linearization
      real(rp),    intent(in)    :: elextC(1),elext(e%ndime),gpadv(e%pnode),acden
      real(rp),    intent(in)    :: dvolu,tidiv,timom,gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,itens
      real(rp)                   :: aux,aux1
      real(rp)                   :: vect(e%ndime), vect2(e%ndime), tensor(e%ndime,e%ndime)

      aux = tidiv*dvolu 
      do inode=1,e%pnode
         vect2=0.0_rp
      
         aux1= (timom*acden*gpadv(inode))*dvolu 
         
         vect(:)=e%cartd(:,inode)*elextC(1)*aux + aux1*(elext(:)*auxoss)
         
         if (linearization ==1) then
            do itens=1,e%ndime
               tensor(itens,:)= aux1*acden*gpvel(:)*grvel(itens,:)
            end do
            call GetGrVel_hu(e%ndime,tensor,vect2) 
         end if   
         
         elrhs(:,inode) = vect(:) + vect2(:) + elrhs(:,inode)
      end do

   end subroutine supm_elmrhuEst1

   subroutine supm_elmrhuEst1Lapla(e,beta,timom,acden,acvis,dvolu,elext,grvel,gpvel,auxoss,elrhs)   
    !-----------------------------------------------------------------------
    !
    !   tau1* (f, beta*eta*lapla(V) + rho*U*grad(v))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: dvolu,acvis,beta,timom,acden,gpvel(e%ndime),grvel(e%ndime,e%ndime),elext(e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode, itens
      real(rp)                   :: aux,auxvis
       real(rp)                   :: vect(e%ndime), vect2(e%ndime), tensor(e%ndime,e%ndime)
                
      aux=(beta*acvis)*(timom*dvolu)*auxoss
      do inode=1,e%pnode
      
         auxvis=aux*sum(e%hessi(1:e%ndime,inode))
         
         vect(:)=auxvis*elext(:)
         
         do itens=1,e%ndime
            tensor(itens,:)= auxvis*acden*gpvel(:)*grvel(itens,:)
         end do
         
         call GetGrVel_hu(e%ndime,tensor,vect2) 
         
         elrhs(:,inode) = vect(:)+vect2(:)+elrhs(:,inode)  
      end do

   end subroutine supm_elmrhuEst1Lapla     

   subroutine supm_elmbsqEst1(e,timom,dvolu,auxtens,auxoss,elmat)
    !-----------------------------------------------------------------------
    !
    ! tau* (-Div(S), grad(q))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elmat(1,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens, jtens
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,auxtens), vect_sol(auxtens)

      aux=-timom*dvolu*auxoss    
      do jnode=1,e%pnode      
         do inode=1,e%pnode   
         
            do itens=1,e%ndime
               do jtens=1,e%ndime
                  tensor(itens,jtens)= e%cartd(itens,inode)*aux*e%cartd(jtens,jnode)
               end do   
            end do
         
            call GetStress(e%ndime,tensor,vect_sol)

            elmat(1,inode,:,jnode) = vect_sol(:) + elmat(1,inode,:,jnode)                    

         end do  
      end do

   end subroutine supm_elmbsqEst1
       
   subroutine supm_elmbuqEst1(e,timom,dvolu,acden,dtinv,grvel,gpadv,auxoss,linearization,elmat)
    !-----------------------------------------------------------------------
    !
    ! tau* (rho*a*grad(U), grad(q))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxoss,linearization
      real(rp),    intent(in)    :: dvolu,timom,acden,grvel(e%ndime,e%ndime),gpadv(e%pnode),dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,itens, Newton
      real(rp)                   :: aux 
      real(rp)                   :: tensor(e%ndime,e%ndime), vect(e%ndime), vect2(e%ndime)
      
      aux = (timom*dvolu)*auxoss      
      do inode=1,e%pnode
         do jnode=1,e%pnode
            tensor=0.0_rp
            vect2=0.0_rp
            
            do itens=1,e%ndime
               vect(itens)=aux*acden*(e%cartd(itens,inode)*(e%shape(jnode,e%igaus)*dtinv+gpadv(jnode)))
               tensor(:,itens)= aux*acden*(e%cartd(itens,inode))*(e%shape(jnode,e%igaus)*grvel(itens,:))*linearization
            end do         
            if (linearization ==1) call GetGrVel_hu(e%ndime,tensor,vect2)
  
            elmat(1,inode,:,jnode)=vect(:)+vect2(:)+elmat(1,inode,:,jnode)            
         end do
      end do

   end subroutine supm_elmbuqEst1
   

   subroutine supm_elmbuqEst1Lapla(e,beta,acvis,timom,dvolu,auxoss,elmat) 
    !-----------------------------------------------------------------------
    !
    ! tau1 * (-beta*eta*lapla(U), grad(q))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: dvolu,timom,acvis,beta
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux,aux2 
      
      aux= -(beta*acvis)*(dvolu*timom)*auxoss      
      do jnode=1,e%pnode      
         aux2= aux*sum(e%hessi(1:e%ndime,jnode))      
         do inode=1,e%pnode                  
            elmat(1,inode,:,jnode) = aux2*e%cartd(:,inode) + elmat(1,inode,:,jnode)     
         end do
      end do

   end subroutine supm_elmbuqEst1Lapla 
   
   subroutine supm_elmbpqEst1(e,timom,dvolu,elmat)
    !-----------------------------------------------------------------------
    !
    !    -tau1*(gra(P),grad(q))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,timom
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,itens
      real(rp)                   :: aux, vector(e%ndime),value

      aux  = timom*dvolu    
      
      do jnode=1,e%pnode  
         do inode=1,e%pnode    
            do itens=1,e%ndime
               vector(itens)=aux*(e%cartd(itens,inode)*e%cartd(itens,jnode)) 
            end do
           
            value=0
            do itens=1,e%ndime
               value=vector(itens)+value
            end do

            elmat(1,inode,1,jnode) = value + elmat(1,inode,1,jnode)
         
         end do  
      end do

   end subroutine supm_elmbpqEst1

   subroutine supm_elmrhpEst1(e,acden,timom,dvolu,elext,gpvel,grvel,auxoss,linearization,elrhs)
    !-----------------------------------------------------------------------
    !
    !  
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxoss,linearization
      real(rp),    intent(in)    :: acden,timom,elext(e%ndime)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,aux1, itens
      real(rp)                   :: aux, tensor(e%ndime,e%ndime),  tensor_trans(e%ndime,e%ndime)
      real(rp)                   :: value, vect2(e%ndime), vector(e%ndime)
      
      aux = timom*dvolu*auxoss
      do inode=1,e%pnode
      
         vector=0.0_rp
         if (linearization ==1) then
            do itens=1, e%ndime
               tensor(:,itens)=acden*gpvel(:)*grvel(itens,:)
               tensor_trans(itens,:)=tensor(:,itens)
            end do   
            call GetGrVel_hu(e%ndime,tensor_trans,vector)  
         end if
         
         call GetGrVel(e%ndime,elext,vector) 
         
         do itens=1,e%ndime
            vect2(itens)=aux*e%cartd(itens,inode)*vector(itens)
         end do
         call GetPlus(e%ndime,vect2,value)

         elrhs(1,inode) = value + elrhs(1,inode) 
  
      end do

   end subroutine supm_elmrhpEst1 
   
   subroutine supm_elmrhpForce(e,auxoss,timom,dvolu,elext,elrhs)
    !-----------------------------------------------------------------------
    !
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: timom,elext(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux

      aux = timom*dvolu*(1_ip-auxoss)
    
      do inode=1,e%pnode
         do idime=1,e%ndime

            elrhs(1,inode) = aux*(e%cartd(idime,inode)*elext(idime)) &
                  + elrhs(1,inode) 
         end do
      end do

   end subroutine supm_elmrhpForce    
  

   !**********************************************************************************************************************************
   !Estabilization terms of Constitutive Equation 
   
   subroutine supm_elmbstEst2(e,auxVE,auxG,auxPTT,acvis,tisig,dtinv,gpadv,gpsig,dvolu,auxtens,grvel,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !    (-(1/2eta)*tau + auxVE*a*grad(tau)+(auxVE*(2*grad(a)*tau...), R_constituitva)
    !   
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: acvis,auxVE,auxG,dvolu,gpsig(auxtens),gpadv(e%pnode),dtinv
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),tisig     
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, itens, jtens
      real(rp)                   :: aux, aux1, valuei, valuej
      real(rp)                   :: tensor(e%ndime, e%ndime), vector(auxtens) 
      real(rp)                   :: tensor2(e%ndime,e%ndime), vector2(auxtens)
      real(rp)                   :: SigmaTermsj(auxtens,auxtens), GiesekusTermj(auxtens,auxtens)
      real(rp)                   :: Matrixj(auxtens,auxtens), DefTermj(auxtens,auxtens)
      real(rp)                   :: SigmaTermsi(auxtens,auxtens),GiesekusTermi(auxtens,auxtens)
      real(rp)                   :: Matrixi(auxtens,auxtens), DefTermi(auxtens,auxtens)
      real(rp)                   :: Prod(auxtens,auxtens)
                                                                         
      aux=1.0_rp/(2.0_rp*acvis)
 
      do jnode=1,e%pnode  
         do inode=1,e%pnode
            DefTermj=0.0_rp
         
            !MATRIXj AUXTENS X AUXTENS-- jnode
            !do itens=1,e%ndime
            tensor=-auxVE*grvel*e%shape(jnode,e%igaus)
            !end do 
            call GetStress_tt(e%ndime,tensor,DefTermj) 

            vector= auxVE*(auxG/acvis)*gpsig*e%shape(jnode,e%igaus)
            call GetAdvStress_tt(e%ndime,auxPTT,vector,GiesekusTermj)
            SigmaTermsj=GiesekusTermj
          
            valuej=e%shape(jnode,e%igaus)*aux + auxVE*gpadv(jnode) + auxVE*dtinv*e%shape(jnode,e%igaus)
            call AddSphericalComponent_tt(e%ndime,valuej,SigmaTermsj)
         
            if (e%ndime==2) DefTermj(auxtens,auxtens)=0.0_rp !???
            Matrixj(1:auxtens,:)= 2.0_rp*DefTermj(1:auxtens,:) + SigmaTermsj(1:auxtens,:)
            Matrixj(e%ndime+1:auxtens,:) = Matrixj(e%ndime+1:auxtens,:)*0.5_rp
            
            
            !MATRIXi AUXTENS X AUXTENS-- inode  (Test functions)
            aux1= tisig*auxVE*e%shape(inode,e%igaus)*dvolu  
            do itens=1, e%ndime
               tensor2(itens,:)= aux1*(2.0_rp*grvel(itens,:))
            end do
            call GetStress_tt(e%ndime,tensor2,DefTermi)
            if (e%ndime==2) DefTermi(auxtens,auxtens)=0.0_rp
            
            vector2=aux1*(auxG/acvis)*gpsig
            call GetAdvStress_tt_B(e%ndime,auxPTT,vector2,GiesekusTermi)
            SigmaTermsi=GiesekusTermi
            
            valuei= tisig*auxVE*gpadv(inode)*dvolu  + e%shape(inode,e%igaus)*(-(tisig/(2.0_rp*acvis)))*dvolu
            call AddSphericalComponent_tt(e%ndime,valuei,SigmaTermsi)
            
            Matrixi=DefTermi+SigmaTermsi
            
            !Product
            Prod=matmul(Matrixi,Matrixj)       
            elmat(1:auxtens,inode,:,jnode)=Prod(1:auxtens,:) +   elmat(1:auxtens,inode,:,jnode)

         end do
      end do   

   end subroutine supm_elmbstEst2
   
   
   
   subroutine supm_elmbutEst2(e,beta,auxVE,auxG,auxPTT,acvis,tisig,gpsig,dvolu,auxtens,grvel,grsig,gpadv,elmat) 
    !-----------------------------------------------------------------------
    !
    ! tau3 *(R_const, R_const)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: auxVE,dvolu,gpsig(auxtens),grsig(auxtens,e%ndime),tisig
      real(rp),    intent(in)    :: acvis,auxG,grvel(e%ndime,e%ndime),beta,gpadv(e%pnode)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, itens
      real(rp)                   :: aux,aux3      
      real(rp)                   :: Matrixj(auxtens,e%ndime), GradVelTerm(auxtens,e%ndime)
      real(rp)                   :: tensor(auxtens,e%ndime), vector(e%ndime,e%ndime), SigmaTerm(auxtens,e%ndime)
      real(rp)                   :: ConvecTerm(auxtens, e%ndime)
      real(rp)                   :: tensor2(e%ndime,e%ndime), tens_grvel(auxtens,auxtens), vector2(auxtens)
      real(rp)                   :: tens_sig(auxtens,auxtens), Matrixi(auxtens,auxtens), c(auxtens,e%ndime)
      

      do jnode=1,e%pnode    
         do inode=1,e%pnode  
         
            !MATRIX AUXTENS X NDIME-- jnode:  -(1-beta)*grad_sym(u) + auxVE(u*grad(sigma)-g(u,sigma)+...
            do itens=1, e%ndime
               tensor(:,itens)= -2.0_rp*auxVE*gpsig(:)*e%cartd(itens,jnode)
               vector(itens,itens) = -(1.0_rp - beta)*e%cartd(itens,jnode)
            end do
            ConvecTerm(1:auxtens,:)= auxVE*e%shape(jnode,e%igaus)*grsig(1:auxtens,:)
            call GetStress_ut(e%ndime,vector,GradVelTerm)
            SigmaTerm=ConvecTerm
            call AddSigmaGradVel_ut(e%ndime,tensor,SigmaTerm)

            Matrixj(1:auxtens,:) = GradVelTerm(1:auxtens,:) + SigmaTerm(1:auxtens,:) 
            Matrixj(e%ndime+1:auxtens,:)=Matrixj(e%ndime+1:auxtens,:)*0.5_rp
            
         
            !MATRIX AUXTENS X AUXTENS -- inode:  -tau/2eta + auxVE ugrad(tau) + g(u,tau)+...
            
            aux3= ((tisig*auxVE)*e%shape(inode,e%igaus))*dvolu         
            
            do itens=1, e%ndime
               tensor2(:,itens)= aux3*(2.0_rp*grvel(:,itens))
            end do
            call GetStress_tt(e%ndime,tensor2,tens_grvel)
            if (e%ndime==2) tens_grvel(auxtens,auxtens)=0
            
            vector2(:)=aux3*((auxG/acvis)*gpsig(:))
            call GetAdvStress_tt_B(e%ndime,auxPTT,vector2,tens_sig)
            
            aux= ((tisig*auxVE)*gpadv(inode))*dvolu + (e%shape(inode,e%igaus)*(-(tisig/(2.0_rp*acvis))))*dvolu
            call AddSphericalComponent_tt(e%ndime,aux,tens_sig)
          
            Matrixi(1:auxtens,:)=tens_grvel(1:auxtens,:)+tens_sig(1:auxtens,:)
          
            c=matmul(Matrixi,Matrixj) !matrices multiplication

            elmat(1:auxtens,inode,:,jnode) = c(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode) 
          
         end do  
      end do

   end subroutine supm_elmbutEst2 
   
   

   subroutine supm_elmrhcEst2(e,auxVE,auxG,auxPTT,tisig,acvis,gpsig,grsig,gpvel,gpadv,grvel &
         ,dvolu,auxtens,elextS,elrhs)         
    !-----------------------------------------------------------------------
    !
    ! tau3*(f,...)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: gpsig(auxtens),grsig(auxtens,e%ndime),tisig,elextS(auxtens) 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpvel(e%ndime),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, itens, jtens
      real(rp)                   :: aux2,aux5, auxGrVel
      real(rp)                   :: traza, tensor(e%ndime,auxtens), vect_elext(auxtens)
      real(rp)                   :: tens_grvel(auxtens,e%ndime), grvel_trans(e%ndime,e%ndime)
      real(rp)                   :: vector(e%ndime), vect_grvel(auxtens), tensorB(auxtens,auxtens)
      real(rp)                   :: Res(auxtens), tensor2(e%ndime,e%ndime), tensGrVel(e%ndime,e%ndime)
      real(rp)                   :: vect_sig(auxtens), vect_GrSig(auxtens), value(auxtens), tens_vel2(auxtens,auxtens)
      real(rp)                   :: vector2(auxtens), tens_sig2(auxtens,auxtens), Res2(auxtens,auxtens), c(auxtens)
      
      
      do inode=1,e%pnode  
      ! VECTOR 1 * AUXTENS
      !product GrSig*Vel 
      do itens=1, e%ndime
         tensor(itens,:) = auxVE*gpvel(itens)*grsig(:,itens)
         tensGrVel(itens,:)=-2.0_rp*auxVE*grvel(itens,:)
      end do    
      call GetGrSig_hc(e%ndime,tensor,vect_GrSig) 
      
      !product GrVel*Sig
      call GetGrVel_hc(e%ndime,tensGrVel,gpsig,vect_GrVel) 
      
      !product Sig*Sig
      do jtens=1, auxtens
         tensorB(jtens,:) = auxVE*(auxG/acvis)*gpsig(jtens)*gpsig(:)
      end do  
      call GetStress_hc(e%ndime,auxPTT,tensorB,vect_sig)
      
      value(:)=elextS(:)
      call AddValues_hc(e%ndime,value,vect_sig)
        
      Res(:) = vect_GrSig(:)+ vect_GrVel(:) + vect_sig(:)
      Res(e%ndime+1:auxtens)=Res(e%ndime+1:auxtens)/2.0_rp

          
      !MATRIX AUXTENS * AUXTENS
      aux2= (tisig*auxVE)*(gpadv(inode))*dvolu  + (e%shape(inode,e%igaus)*(-(tisig/(2.0_rp*acvis))))*dvolu
      aux5= (tisig*auxVE)*(e%shape(inode,e%igaus))*dvolu
         
      do itens=1, e%ndime
         tensor2(itens,:) = aux5*2.0_rp*grvel(itens,:) 
      end do  
         
      call GetStress_tt(e%ndime,tensor2,tens_vel2)
      if (e%ndime==2) tens_vel2(auxtens,auxtens)=0.0_rp
      vector2(:)=aux5*(auxG/acvis)*gpsig(:)
      call GetAdvStress_tt_B(e%ndime,auxPTT,vector2,tens_sig2)
      call AddSphericalComponent_tt(e%ndime,aux2,tens_sig2)
         
      Res2(1:auxtens,:)=tens_vel2(1:auxtens,:) + tens_sig2(1:auxtens,:) 
         
      c=matmul(Res2,Res)
         
      elrhs(:,inode)=c(:)+ elrhs(:,inode)
 
      end do

   end subroutine supm_elmrhcEst2   

   subroutine supm_elmbsvEst2(e,beta,auxVE,auxG,auxPTT,acvis,tisig,dtinv,gpadv,gpsig,dvolu,auxtens,grvel,elmat)
    !-----------------------------------------------------------------------
    !
    !    tau3*(R_constituitva(S),R_constituitva(v))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: acvis,auxVE,beta,auxG,dvolu,gpadv(e%pnode),gpsig(auxtens)
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),tisig,dtinv     
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode,auxstab, itens
      real(rp)                   :: aux2, aux, aux1
      real(rp)                   :: traza      
      real(rp)                   :: Res(auxtens,auxtens), value, tensor(e%ndime,e%ndime), tens_grvel(auxtens,auxtens)
      real(rp)                   :: tens_sig(auxtens,auxtens), vector(auxtens), tensor2(e%ndime,e%ndime)
      real(rp)                   :: tens_grvel2(auxtens,e%ndime), Res2(e%ndime,auxtens), c(e%ndime,auxtens)
                            
      aux2 = (tisig*dvolu)
      
      traza= (gpsig(1) + gpsig(2))
      aux=1.0_rp/(2.0_rp*acvis)
      do jnode=1,e%pnode 
         do inode=1,e%pnode  
         
            ! MATRIX 1 auxtens x auxtens  -- jnode (sigma)
            do itens=1,e%ndime
               tensor(itens,:)=-auxVE*grvel(itens,:)*e%shape(jnode,e%igaus)
            end do
            call GetStress_tt(e%ndime,tensor,tens_grvel)  

            vector(:)= auxVE*(auxG/acvis)*gpsig(:)*e%shape(jnode,e%igaus)
            call GetAdvStress_tt(e%ndime,auxPTT,vector,tens_sig) 
          
            value=e%shape(jnode,e%igaus)*aux + auxVE*gpadv(jnode) + auxVE*dtinv*e%shape(jnode,e%igaus)
            call AddSphericalComponent_tt(e%ndime,value,tens_sig)
         
            if (e%ndime==2) tens_grvel(auxtens,auxtens)=0.0_rp
            Res(1:auxtens,:)= 2.0_rp*tens_grvel(1:auxtens,:) + tens_sig(1:auxtens,:)
            Res(e%ndime+1:auxtens,:) = Res(e%ndime+1:auxtens,:)/2.0

            ! MATRIX 2 ndime x auxtens  -- inode (v)
            do itens=1,e%ndime
               tensor2(itens,itens)=-aux2*e%cartd(itens,inode)
            end do   
            
            call GetStress_ut(e%ndime,tensor2,tens_grvel2)
            Res2= transpose(tens_grvel2)
            c=matmul(Res2,Res)
            
            elmat(1:e%ndime,inode,:,jnode)=c(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)
                                      
         end do
      end do   

   end subroutine supm_elmbsvEst2
   
   
  
   
   subroutine supm_elmbuvEst2(e,beta,auxVE,tisig,gpsig,dvolu,auxtens,grvel,grsig,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! tau3*(-(1-beta)*grad(u) + auxVE*(U*grad(sig)-... ), -grad(V) + auxVE*(V*grad(tau)+...) )
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: auxVE,dvolu,gpsig(auxtens),grsig(auxtens,e%ndime),tisig
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,itens
      real(rp)                   :: aux2
      real(rp)                   :: Res(auxtens,e%ndime), tens_grsig(auxtens,e%ndime), tensor(e%ndime,e%ndime)
      real(rp)                   :: tens_vel(auxtens,e%ndime), tens_gpsig(auxtens,e%ndime)
      real(rp)                   :: tensor2(e%ndime,e%ndime)
      real(rp)                   :: tens_grvel2(auxtens,e%ndime), Res2(e%ndime,auxtens), c(e%ndime,e%ndime)
      
      
      aux2=tisig*dvolu
      
      do jnode=1,e%pnode   
          do inode=1,e%pnode
         !MATRIX 1 AUXTENS X NDIME jdime (u)
         do itens=1, auxtens
            tens_grsig(itens,:) = auxVE*grsig(itens,:)*e%shape(jnode,e%igaus)
            tens_gpsig(itens,:) = - 2.0_rp*auxVE*gpsig(itens)*e%cartd(:,jnode) 
         end do    
         
         do itens=1,e%ndime
            tensor(itens,itens)=-(1.0_rp - beta)*e%cartd(itens,jnode)
         end do
         call GetStress_ut(e%ndime,tensor,tens_vel)
         
         call AddSigmaGradVel_ut(e%ndime,tens_gpsig,tens_grsig)

         Res(1:auxtens,:)= tens_vel(1:auxtens,:) + tens_grsig(1:auxtens,:)
         Res(e%ndime+1:auxtens,:)=Res(e%ndime+1:auxtens,:)/2.0_rp
         
         ! MATRIX 2 NDIME X AUXTENS idime (v)
         do itens=1,e%ndime
            tensor2(itens,itens)=-aux2*e%cartd(itens,inode)
         end do   
         call GetStress_ut(e%ndime,tensor2,tens_grvel2)
         Res2= transpose(tens_grvel2)
         
         c=matmul(Res2,Res)
            
         elmat(1:e%ndime,inode,:,jnode)=c(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)

            end do  
      end do

   end subroutine supm_elmbuvEst2 
   
   
   
   subroutine supm_elmrhuEst2(e,beta,auxVE,auxG,auxPTT,tisig,acvis,gpsig,grsig,gpvel,grvel,dvolu,auxtens,elextS,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: gpsig(auxtens),grsig(auxtens,e%ndime),tisig,beta
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpvel(e%ndime),grvel(e%ndime,e%ndime),elextS(auxtens)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,aux1,auxstab,aux3,aux,itens,jtens
      real(rp)                   :: aux2
      real(rp)                   :: Res11,Res12,Res22, Res(auxtens)
      real(rp)                   :: value(auxtens), auxGrVel
      real(rp)                   :: tensor(e%ndime,auxtens),vect_GrSig(auxtens), tensorB(auxtens,auxtens)
      real(rp)                   :: vect_GrVel(auxtens), vect_Sig(auxtens), tensGrVel(e%ndime,e%ndime)
      real(rp)                   :: tensor2(e%ndime,e%ndime)
      real(rp)                   :: tens_grvel2(auxtens,e%ndime), Res2(e%ndime,auxtens), c(e%ndime)
      
      aux2=(tisig)*dvolu      
      
      do inode=1,e%pnode 
         !VECTOR AUXTENS X 1
         !Product Vel*GrSig
         do itens=1, e%ndime
            tensor(itens,:) = auxVE*gpvel(itens)*grsig(:,itens)
            tensGrVel(itens,:)=-2.0_rp*auxVE*grvel(itens,:)
         end do    
         call GetGrSig_hc(e%ndime,tensor,vect_GrSig) 
      
         !product GrVel*Sig
         
         call GetGrVel_hc(e%ndime,tensGrVel,gpsig,vect_GrVel) 
      
         !product Sig*Sig
         do jtens=1, auxtens
            tensorB(jtens,:) = auxVE*(auxG/acvis)*gpsig(jtens)*gpsig(:)
         end do  
         call GetStress_hc(e%ndime,auxPTT,tensorB,vect_sig)
      
         value(:)=elextS(:)
         call AddValues_hc(e%ndime,value,vect_sig)
        
         Res(:) = vect_GrSig(:)+ vect_GrVel(:) + vect_sig(:)
         Res(e%ndime+1:auxtens)=Res(e%ndime+1:auxtens)/2.0_rp
      
       
         !MATRIX 2 NDIME X AUXTENS inode
         do itens=1,e%ndime
            tensor2(itens,itens)=-aux2*e%cartd(itens,inode)
         end do   
         call GetStress_ut(e%ndime,tensor2,tens_grvel2)
         Res2= transpose(tens_grvel2)
         
         c=matmul(Res2,Res)
            
         elrhs(:,inode)=c(:) + elrhs(:,inode)

      end do

   end subroutine supm_elmrhuEst2 
   
   !--------------------------------------------------------------------------------------------------------------------
   !Split OSS
   
   subroutine supm_elmrhu_splitoss(e,acden,timom,dvolu,gpconv,gpadv,auxoss,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    ! 
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: gpconv(e%ndime),timom,acden,dvolu,gpadv(e%pnode)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: aux1,inode
      real(rp)                   :: aux2,aux3
      
      aux1=(1_ip - auxoss)
      aux2=((timom)*acden*dvolu)*aux1               

      do inode=1,e%pnode
      
         aux3=gpadv(inode)*aux2
                           
         elrhs(:,inode) = aux3*acden*gpconv(:) + elrhs(:,inode)                
      end do

   end subroutine supm_elmrhu_splitoss 
   
   subroutine supm_laplatovector(e,acvis,beta,elvel,gplapl)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: acvis,elvel(e%ndime,e%pnode),beta,gplapl(e%ndime)      
      integer(ip) :: inode,idime
      real(rp) :: aux1(e%pnode)
      
      
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo   
      do idime = 1,e%ndime                              ! Contribution from the laplacian term
         gplapl(idime) = gplapl(idime) -  (1.0_rp*beta*acvis)*dot_product(aux1,elvel(idime,1:e%pnode))
      end do
   
   end subroutine      
   
   subroutine supm_elmrhu_splitossNonLinear(e,acvis,beta,timom,dvolu,gplapl,auxoss,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !+ tau1*(rho*a*grad(v), resid_momentum)+tau2(div(tau),resid_continuity)-tau3(grad_sym(v),resid_constitutive)
    !
    !-----------------------------------------------------------------------
   use typre
   implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxoss
      real(rp),    intent(in)    :: gplapl(e%ndime),acvis,beta
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime
      integer(ip)                :: aux1
      real(rp)                   :: aux,aux2,tmp2,aux3

      aux1= 0.0_rp !(1_ip - auxoss)
      
      do inode=1,e%pnode

         aux = aux1*(timom*((beta*acvis)*sum(e%hessi(1:e%ndime,inode)))*dvolu)

         elrhs(:,inode) = aux*gplapl(:) + elrhs(:,inode)
         
      end do

   end subroutine supm_elmrhu_splitossNonLinear 
   
   
   subroutine supm_elmrhc_splitoss(e,auxtens,beta,timom,dvolu,gpdivs,auxoss,elrhs)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the rhs terms in constitutive equation for OSS
   !    
   ! -tau1(div(tau),resid_momentum)- tau3/2mu(tau,resid_constitutive)
   !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: gpdivs(e%ndime)
      real(rp),    intent(in)    :: dvolu,beta,timom
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,aux2, itens
      real(rp)                   :: aux1,aux3, tensor(e%ndime,e%ndime), vect_sol(auxtens)

      aux2 = 1_ip - auxoss
      aux1 = ((1.0_rp - beta)*timom*dvolu)*aux2    
      do inode=1,e%pnode
      
         do itens=1,e%ndime
            tensor(itens,:)=e%cartd(itens,inode)*aux1*gpdivs(:)
         end do
         call GetStress(e%ndime,tensor,vect_sol)

       elrhs(:,inode) = vect_sol(:) + elrhs(:,inode)

      end do
   end subroutine supm_elmrhc_splitoss  
   
   subroutine supm_elmrhp_splitoss(e,timom,dvolu,gpgrap,auxoss,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau1*(grad q, resid_continuity)
    !
    !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxoss       
      real(rp),    intent(in)    :: gpgrap(e%ndime),timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,idime,aux1
      real(rp)                   :: tmp1

      aux1 = 1_ip - auxoss
      tmp1 = (dvolu*timom)*aux1

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = e%cartd(idime,inode)*gpgrap(idime)*tmp1 + elrhs(1,inode)
         end do
      end do

   end subroutine supm_elmrhp_splitoss 

   
   subroutine supm_elmbstDC(e,auxtens,kdisc,dvolu,elmat)
   
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,kdisc
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux, vect(e%ndime), value

      aux= kdisc*dvolu  
      do jnode=1,e%pnode  
         do inode=1,e%pnode
         
         vect(:)=aux*e%cartd(:,inode)*e%cartd(:,jnode)
         value=sum(vect(1:e%ndime))
         
         call AddSphericalComponent_tt(e%ndime,value,elmat(1:auxtens,inode,:,jnode))

         end do
      end do   

   end subroutine supm_elmbstDC
   
   !-----------------------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   
   !This suborutine computes the residual at a Gauss point
   subroutine supm_elmrfe_oto(e,dtnsi,acden,acvis,veadv, &
         vegau,gpsig,grapr,grave,grasig,auxtens,elext,resim)
      use typre
      implicit none
      class(FiniteElement)        :: e
      integer(ip) :: auxtens, itens
      real(rp) :: dtnsi,acden,veadv(e%ndime),vegau(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime),resim(auxtens+e%ndime+1)
      real(rp) :: gpsig(auxtens),grasig(auxtens,e%ndime),acvis, grasig_trans(e%ndime,auxtens)
      real(rp) :: elext(e%ndime)
      real(rp) :: aux,aux2
      real(rp) :: tensor(e%ndime,e%ndime), vect_gra(auxtens), vect_gras(e%ndime)
      
      
      integer(ip) :: idime,jdime
      
      resim = 0.0_rp
      aux=1.0_rp/(2.0_rp*acvis)
      aux2=1.0_rp
      
      do itens=1,e%ndime
         tensor(itens,:)=-aux2*grave(itens,:)
      end do
      
      call GetStress(e%ndime,tensor,vect_gra)
      vect_gra(e%ndime+1:auxtens)=vect_gra(e%ndime+1:auxtens)*0.5_rp
      
      !Residual in constitutive equation
      resim(1:auxtens)=aux*gpsig(1:auxtens)+vect_gra(1:auxtens)
      
      !Temporal derivative LHS
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*dtnsi*vegau
      
      !RHS (Temporal derivative + External forces)
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) - elext
      
      !Contribution from pressure term
      !Pressure Gradient
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + grapr     

      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do

      !Contribution from the divergence of sigma
      grasig_trans=transpose(grasig)
      call GetVel_hu(e%ndime,grasig_trans,vect_gras)
      resim(auxtens+1:auxtens+e%ndime)= resim(auxtens+1:auxtens+e%ndime) - vect_gras(1:e%ndime)

      !Contribution from the divergence term
      do idime = 1,e%ndime                              
         resim(auxtens+e%ndime+1) = resim(auxtens+e%ndime+1) + grave(idime,idime)
      end do

   end subroutine supm_elmrfe_oto 
   
   subroutine supm_elmrfeVE_oto(e,dtnsi,acden,acvis,beta,lambda,auxVE,auxG,auxPTT,veadv, &
         vegau,gpsig,grapr,grave,grasig,auxtens,elext,elextC,elextS,resim)        
         
      use typre
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in) :: auxtens,auxPTT
      real(rp), intent(in)    :: dtnsi,acden,veadv(e%ndime),vegau(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime)
      real(rp), intent(in)    :: gpsig(auxtens),grasig(auxtens,e%ndime),acvis,elextC(1),elextS(auxtens)
      real(rp), intent(in)    :: elext(e%ndime),beta,auxVE,auxG,lambda
      real(rp), intent(out)   :: resim(auxtens+e%ndime+1)
      real(rp)                :: tensor(e%ndime,auxtens), tensor2(e%ndime,e%ndime), vect_grave(auxtens)
      real(rp)                :: vect_GrSigve(auxtens), vect_GrVelsig(auxtens), tensorB(auxtens,auxtens), vect_sig(auxtens)
      real(rp)                :: aux, auxGrVel, value(auxtens) 
      real(rp)                :: grasig_trans(e%ndime,auxtens),  vect_gras(e%ndime), tensGrVel(e%ndime,e%ndime)
      integer(ip)             :: idime,jdime, itens
      
      resim = 0.0_rp
      aux=1.0_rp/(2.0_rp*acvis)
      tensGrVel=0.0_rp
      tensor=0.0_rp
      tensor2=0.0_rp
      value=0.0_rp
      
      !Residual in constitutive equation   
      do itens=1, e%ndime
         tensor(itens,:) = auxVE*veadv(itens)*grasig(:,itens)
         tensor2(itens,:) = -(1.0_rp - beta)*(grave(itens,:)) 
         tensGrVel(itens,:)=-2.0_rp*auxVE*grave(itens,:)
      end do    
      call GetStress(e%ndime,tensor2,vect_grave) 

      call GetGrSig_hc(e%ndime,tensor,vect_GrSigve) 

      call GetGrVel_hc(e%ndime,tensGrVel,gpsig,vect_GrVelsig) 
      
      do itens=1, auxtens
         tensorB(itens,:) = auxVE*(auxG/acvis)*gpsig(itens)*gpsig(:)
      end do  
      call GetStress_hc(e%ndime,auxPTT,tensorB,vect_sig)
      
      value(:)=aux*gpsig(:) + auxVE*gpsig(:)*dtnsi - elextS(:)
      call AddValues_hc(e%ndime,value,resim)

      resim(1:auxtens)=  vect_GrSigve(1:auxtens) + vect_GrVelsig(1:auxtens) + vect_grave(1:auxtens) & 
                        + vect_sig(1:auxtens) + resim(1:auxtens)
                        
      resim(e%ndime+1:auxtens)=resim(e%ndime+1:auxtens)*0.5_rp

      !Temporal derivative LHS
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*dtnsi*vegau(1:e%ndime)
      
      !RHS (Temporal derivative + External forces)
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) - elext(1:e%ndime)
      
      !Contribution from pressure term
      !Pressure Gradient
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + grapr(1:e%ndime)     

      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do
      
      !Contribution from the divergence of sigma
      grasig_trans=transpose(grasig)
      call GetVel_hu(e%ndime,grasig_trans,vect_gras)
      resim(auxtens+1:auxtens+e%ndime)= resim(auxtens+1:auxtens+e%ndime) - vect_gras(1:e%ndime)

      !Contribution from the divergence term
      do idime = 1,e%ndime                              
         resim(auxtens+e%ndime+1) = resim(auxtens+e%ndime+1) + grave(idime,idime)
      end do
      
   end subroutine supm_elmrfeVE_oto    
   
   
   subroutine supm_elmrfeVE_oto2(e,dtnsi,dtnsi_const,acden,acvis,beta,lambda,auxVE,auxG,auxPTT,veadv, &
         vegau,gpsig,grapr,grave,grasig,auxtens,elext,elextC,elextS,resim)        
         
      use typre
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in) :: auxtens,auxPTT
      real(rp), intent(in)    :: dtnsi,dtnsi_const,acden,veadv(e%ndime),vegau(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime)
      real(rp), intent(in)    :: gpsig(auxtens),grasig(auxtens,e%ndime),acvis,elextC(1),elextS(auxtens)
      real(rp), intent(in)    :: elext(e%ndime),beta,auxVE,auxG,lambda
      real(rp), intent(out)   :: resim(auxtens+e%ndime+1)
      real(rp)                :: tensor(e%ndime,auxtens), tensor2(e%ndime,e%ndime), vect_grave(auxtens)
      real(rp)                :: vect_GrSigve(auxtens), vect_GrVelsig(auxtens), tensorB(auxtens,auxtens), vect_sig(auxtens)
      real(rp)                :: aux, auxGrVel, value(auxtens) 
      real(rp)                :: grasig_trans(e%ndime,auxtens),  vect_gras(e%ndime), tensGrVel(e%ndime,e%ndime)
      integer(ip)             :: idime,jdime, itens
      
      resim = 0.0_rp
      aux=1.0_rp/(2.0_rp*acvis)
      tensGrVel=0.0_rp
      tensor=0.0_rp
      tensor2=0.0_rp
      value=0.0_rp
      
      !Residual in constitutive equation   
      do itens=1, e%ndime
         tensor(itens,:) = auxVE*veadv(itens)*grasig(:,itens)
         tensor2(itens,:) = -(1.0_rp - beta)*(grave(itens,:)) 
         tensGrVel(itens,:)=-2.0_rp*auxVE*grave(itens,:)
      end do    
      call GetStress(e%ndime,tensor2,vect_grave) 

      call GetGrSig_hc(e%ndime,tensor,vect_GrSigve) 

      call GetGrVel_hc(e%ndime,tensGrVel,gpsig,vect_GrVelsig) 
      
      do itens=1, auxtens
         tensorB(itens,:) = auxVE*(auxG/acvis)*gpsig(itens)*gpsig(:)
      end do  
      call GetStress_hc(e%ndime,auxPTT,tensorB,vect_sig)
      
      value(:)=aux*gpsig(:) + auxVE*gpsig(:)*dtnsi_const - elextS(:)
      call AddValues_hc(e%ndime,value,resim)

      resim(1:auxtens)=  vect_GrSigve(1:auxtens) + vect_GrVelsig(1:auxtens) + vect_grave(1:auxtens) & 
                        + vect_sig(1:auxtens) + resim(1:auxtens)
                        
      resim(e%ndime+1:auxtens)=resim(e%ndime+1:auxtens)*0.5_rp

      !Temporal derivative LHS
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*dtnsi*vegau(1:e%ndime)
      
      !RHS (Temporal derivative + External forces)
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) - elext(1:e%ndime)
      
      !Contribution from pressure term
      !Pressure Gradient
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + grapr(1:e%ndime)     

      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do
      
      !Contribution from the divergence of sigma
      grasig_trans=transpose(grasig)
      call GetVel_hu(e%ndime,grasig_trans,vect_gras)
      resim(auxtens+1:auxtens+e%ndime)= resim(auxtens+1:auxtens+e%ndime) - vect_gras(1:e%ndime)

      !Contribution from the divergence term
      do idime = 1,e%ndime                              
         resim(auxtens+e%ndime+1) = resim(auxtens+e%ndime+1) + grave(idime,idime)
      end do
      
   end subroutine supm_elmrfeVE_oto2   
   
   
   subroutine supm_elmrfe_oto_nonlinear(e,acvis,auxtens,beta,elvel,gprep)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: auxtens
      real(rp) :: acvis,elvel(e%ndime,e%pnode),beta,gprep(e%ndime+auxtens+1)      
      integer(ip) :: inode,idime
      real(rp) :: aux1(e%pnode)
      
      
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo   
      do idime = 1,e%ndime                              ! Contribution from the laplacian term
         gprep(auxtens + idime) = gprep(auxtens + idime) -  (1.0_rp*beta*acvis)*dot_product(aux1,elvel(idime,1:e%pnode))
      end do
   
   end subroutine  
   
   !This subroutine assemblies the contribution of the residual to the RHS for projecting it
   subroutine supm_elmrep(e,dvol,ndofn,gprep,elrep)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) :: ndofn
      real(rp)                :: dvol
      real(rp)                :: gprep(ndofn)
      real(rp)                :: elrep(ndofn,*)
      
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elrep(:,inode) = elrep(:,inode) + e%shape(inode,e%igaus)*gprep*dvol
      enddo
      
   end subroutine   
   
   subroutine supm_elmrhu_oss(e,acden,tidiv,tisig,timom,dvolu,gpadv,gprep,auxtens,auxoss,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !+ tau1*(rho*a*grad(v), resid_momentum)+tau2(div(tau),resid_continuity)-tau3(grad_sym(v),resid_constitutive)
    !
    !-----------------------------------------------------------------------
   use typre
   implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: gprep(auxtens+e%ndime+1),acden
      real(rp),    intent(in)    :: gpadv(e%pnode)
      real(rp),    intent(in)    :: tidiv,dvolu,tisig,timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode,idime, itens
      real(rp)                   :: aux,aux2,tmp2,aux3
      real (rp)                  :: tensor(e%ndime,auxtens), vect_sol(e%ndime),gprep2(auxtens)


      aux2 = tidiv*dvolu
      aux3 = tisig*dvolu
      do inode=1,e%pnode

         aux = timom*(acden*gpadv(inode)*auxoss)*dvolu
         
         gprep2(:)=gprep(1:auxtens)
         do itens=1,e%ndime
            tensor(itens,:)=- aux3*e%cartd(itens,inode)*gprep2(:)
         end do
         call GetVel_hu(e%ndime,tensor,vect_sol)

         elrhs(1:e%ndime,inode) = aux*gprep(auxtens+1:auxtens+e%ndime) & 
              + e%cartd(1:e%ndime,inode)*aux2*gprep(auxtens+e%ndime+1) & 
              + vect_sol(1:e%ndime) + elrhs(1:e%ndime,inode)

      end do

   end subroutine supm_elmrhu_oss
   
   subroutine supm_elmrhu_ossNL(e,auxtens,acvis,beta,timom,dvolu,gprep,auxoss,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !+ tau1*(rho*a*grad(v), resid_momentum)+tau2(div(tau),resid_continuity)-tau3(grad_sym(v),resid_constitutive)
    !
    !-----------------------------------------------------------------------
   use typre
   implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxoss,auxtens
      real(rp),    intent(in)    :: gprep(auxtens+e%ndime+1),acvis,beta
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,itens
      real(rp)                   :: aux,aux2,tmp2,aux3


      do inode=1,e%pnode
      
         aux=0.0_rp
         do itens=1,e%ndime
            aux = timom*((beta*acvis)*e%hessi(itens,inode)*auxoss)*dvolu +aux
         end do   

         elrhs(:,inode) = aux*gprep(auxtens+1:auxtens+e%ndime) + elrhs(:,inode)
         
      end do

   end subroutine supm_elmrhu_ossNL
   
   subroutine supm_elmrhp_oss(e,timom,dvol,gprep,auxtens,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau1*(grad q, resid_continuity)
    !
    !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens       
      real(rp),    intent(in)    :: gprep(auxtens+e%ndime+1)
      real(rp),    intent(in)    :: timom,dvol
      real(rp),    intent(inout) :: elrhs(1,e%mnode)

      integer(ip)                :: inode,idime
      real(rp)                   :: tmp1

      tmp1 = dvol*timom

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = e%cartd(idime,inode)*gprep(auxtens+idime)*tmp1 + elrhs(1,inode)
         end do
      end do

   end subroutine supm_elmrhp_oss 
   
   
   subroutine supm_elmrhc_oss(e,timom,tisig,dvolu,gprep,acvis,beta,auxtens,auxoss,elrhs)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the rhs terms in constitutive equation for OSS
   !    
   ! -tau1(div(tau),resid_momentum)- tau3/2mu(tau,resid_constitutive)
   !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: gprep(e%ndime+1+auxtens)
      real(rp),    intent(in)    :: acvis,beta
      real(rp),    intent(in)    :: timom,tisig,dvolu
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, itens
      real(rp)                   :: aux1,aux3, auxgprep(e%ndime), tensor(e%ndime,e%ndime)
      real(rp)                   :: vect_sol(auxtens), value(auxtens)

      aux1 = -(1.0_rp - beta)*timom*dvolu*auxoss
      aux3 = -1.0_rp/(2.0_rp*acvis)*tisig*dvolu      
      do inode=1,e%pnode
      
      auxgprep(1:e%ndime)=gprep(auxtens+1:auxtens+e%ndime)
      do itens=1,e%ndime
         tensor(itens,:)=e%cartd(itens,inode)*aux1*auxgprep(:)
      end do   
      call GetStress(e%ndime,tensor,vect_sol)
      
      value(1:auxtens)=aux3*e%shape(inode,e%igaus)*gprep(1:auxtens)
      call AddValues_hc(e%ndime,value,vect_sol)

      elrhs(:,inode) = vect_sol(:)+ elrhs(:,inode)

      end do
   end subroutine supm_elmrhc_oss
   
   
   subroutine supm_elmrhcVES_oss(e,auxVE,auxG,auxPTT,tisig,acvis,dvolu,auxtens,gpadv,grvel,gpsig,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: tisig,gprep(e%ndime+1+auxtens),gpadv(e%pnode),grvel(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, itens
      real(rp)                   :: aux,aux2,aux3,aux5
      real(rp)                   :: Res(auxtens)
      real(rp)                   :: tensor2(e%ndime,e%ndime)
      real(rp)                   :: vect_sig(auxtens), vect_GrSig(auxtens), value(auxtens), tens_vel2(auxtens,auxtens)
      real(rp)                   :: vector2(auxtens), tens_sig2(auxtens,auxtens), Res2(auxtens,auxtens), c(auxtens)
      
      Res(:) = gprep(1:auxtens)
      do inode=1,e%pnode
         tensor2=0.0_rp
         tens_sig2=0.0_rp
         aux2=0.0_rp

         aux2= auxVE*gpadv(inode)*tisig*dvolu      
         aux5= tisig*auxVE*e%shape(inode,e%igaus)*dvolu
         
         do itens=1, e%ndime
            tensor2(itens,:) = aux5*2.0_rp*grvel(itens,:) 
         end do  
         
         call GetStress_tt(e%ndime,tensor2,tens_vel2)
         if (e%ndime==2) tens_vel2(auxtens,auxtens)=0.0_rp
         
         vector2(:)=aux5*(auxG/acvis)*gpsig(:)
         call GetAdvStress_tt_B(e%ndime,auxPTT,vector2,tens_sig2)
         call AddSphericalComponent_tt(e%ndime,aux2,tens_sig2)
         
         Res2(1:auxtens,:)= tens_vel2(1:auxtens,:) + tens_sig2(1:auxtens,:)
         
         c=matmul(Res2,Res)
         elrhs(:,inode)=c(:)+ elrhs(:,inode)
      end do

   end subroutine supm_elmrhcVES_oss 
   
   
   subroutine supm_elmrhcVES_dss2(e,auxVE,auxG,auxPTT,dtinv,tisig,acvis,dvolu,auxtens,gpadv,grvel,gpsig,sisgs,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: tisig,sisgs(*),gpadv(e%pnode),grvel(e%ndime,e%ndime),dtinv 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, itens
      real(rp)                   :: aux,aux2,aux3,aux5
      real(rp)                   :: Res(auxtens)
      real(rp)                   :: tensor2(e%ndime,e%ndime)
      real(rp)                   :: vect_sig(auxtens), vect_GrSig(auxtens), value(auxtens), tens_vel2(auxtens,auxtens)
      real(rp)                   :: vector2(auxtens), tens_sig2(auxtens,auxtens), Res2(auxtens,auxtens), c(auxtens)
      
      Res(:) = sisgs(1:auxtens)*dtinv*auxVE
      
      do inode=1,e%pnode
         tensor2=0.0_rp
         tens_sig2=0.0_rp
         aux2=0.0_rp

         aux2= auxVE*gpadv(inode)*tisig*dvolu      
         aux5= tisig*auxVE*e%shape(inode,e%igaus)*dvolu
         
         do itens=1, e%ndime
            tensor2(itens,:) = aux5*2.0_rp*grvel(itens,:) 
         end do  
         
         call GetStress_tt(e%ndime,tensor2,tens_vel2)
         if (e%ndime==2) tens_vel2(auxtens,auxtens)=0.0_rp
         
         vector2(:)=aux5*(auxG/acvis)*gpsig(:)
         call GetAdvStress_tt_B(e%ndime,auxPTT,vector2,tens_sig2)
         call AddSphericalComponent_tt(e%ndime,aux2,tens_sig2)
         
         Res2(1:auxtens,:)= tens_vel2(1:auxtens,:) + tens_sig2(1:auxtens,:)
         
         c=matmul(Res2,Res)
         elrhs(:,inode)=c(:)+ elrhs(:,inode)

      end do

   end subroutine supm_elmrhcVES_dss2 
   
   
   
   
   subroutine supm_elmrhu_splitoss_tidiv(e,tidiv,dvolu,gprep,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! t2*(div(v),div(u))
    !
    !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: gprep(1)
      real(rp),    intent(in)    :: tidiv,dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode

      
      do inode=1,e%pnode
         elrhs(1:e%ndime,inode) =  tidiv*dvolu*e%cartd(1:e%ndime,inode)*gprep(1) & 
             + elrhs(1:e%ndime,inode)

      end do

   end subroutine supm_elmrhu_splitoss_tidiv
   
   
   subroutine supm_elmrhu_splitoss_tisig(e,tisig,beta,dvolu,auxtens,gprep,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !+ tau1*(rho*a*grad(v), resid_momentum)+tau2(div(tau),resid_continuity)-tau3(grad_sym(v),resid_constitutive)
    !
    !-----------------------------------------------------------------------
   use typre
   implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gprep(auxtens)
      real(rp),    intent(in)    :: dvolu,tisig, beta
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode, i, j
      real(rp)                   :: aux
      real (rp)                  :: tensor(e%ndime,auxtens), vector(e%ndime)
      real(rp)                   :: gprep_Mat(e%ndime,e%ndime), Term(e%ndime)

      aux = tisig*dvolu*(1.0_rp-beta)
      
      !call sup_SigmaMatrix(e%ndime,auxtens,gprep,gprep_Mat)
      
      do inode=1,e%pnode
         do i=1,e%ndime
            tensor(i,:)= aux*e%cartd(i,inode)*gprep(:) 
            !do j=1,e%ndime
             !  Term(i)= -aux*0.5_rp*gprep_Mat(j,i)*e%cartd(j,inode)+ Term(i)
             !  Term(i)= -aux*0.5_rp*gprep_Mat(i,j)*e%cartd(j,inode)+ Term(i)
            !end do
         end do
         call GetVel_hu(e%ndime,tensor,vector)
         elrhs(1:e%ndime,inode) =  vector(1:e%ndime) + elrhs(1:e%ndime,inode)
               

      end do

   end subroutine supm_elmrhu_splitoss_tisig
   
   subroutine supm_elmrhcVES_splitoss(e,auxVE,tisig,dvolu,auxtens,gpadv,grvel,gpconvsigma,gpdeform,elrhs)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: gpconvsigma(auxtens), gpdeform(auxtens)
      real(rp),    intent(in)    :: tisig, dvolu,auxVE
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k, l
      real(rp)                   :: aux,Shapi
      real(rp)                   :: gpconvsigmaMat(e%ndime,e%ndime), gpdeformMat(e%ndime,e%ndime)
      real(rp)                   :: tensor(e%ndime,e%ndime), vector(auxtens)
      real(rp)                   ::  Prod1(e%ndime,e%ndime), Prod2(e%ndime,e%ndime)
      
      call sup_SigmaMatrix(e%ndime,auxtens,gpconvsigma,gpconvsigmaMat)
      call sup_SigmaMatrix(e%ndime,auxtens,gpdeform,gpdeformMat)
      aux=tisig*dvolu
      Prod1=matmul(gpdeformMat,transpose(grvel))
      Prod2=matmul(grvel,gpdeformMat)
      

      do inode=1,e%pnode
         tensor=0.0_rp
         tensor=  gpconvsigmaMat*auxVE*gpadv(inode)
         Shapi=e%shape(inode,e%igaus)
         tensor= auxVE*Prod1*Shapi + auxVE*Prod2*Shapi +tensor
         call GetStress(e%ndime,tensor,vector)
         
         elrhs(:,inode)= aux*vector(:)+elrhs(:,inode)  

      end do

   end subroutine supm_elmrhcVES_splitoss   
   
   
   subroutine supm_elmrhcLCR_splitoss(e,auxVE,beta,lambda,lambda0,tisig,dvolu,auxtens,gpadv,grvel,gpconvsigma,gpdeform,elrhs)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: gpconvsigma(auxtens), gpdeform(auxtens)
      real(rp),    intent(in)    :: tisig, dvolu,auxVE, lambda, lambda0, beta
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k, l
      real(rp)                   :: aux, auxVEnew, Shapi
      real(rp)                   :: gpconvsigmaMat(e%ndime,e%ndime), gpdeformMat(e%ndime,e%ndime)
      real(rp)                   :: tensor(e%ndime,e%ndime)
      real(rp)                   :: vector(auxtens), Prod1(e%ndime,e%ndime), Prod2(e%ndime,e%ndime)
      
      call sup_SigmaMatrix(e%ndime,auxtens,gpconvsigma,gpconvsigmaMat)
      call sup_SigmaMatrix(e%ndime,auxtens,gpdeform,gpdeformMat)
      aux=tisig*dvolu
      gpconvsigmaMat=aux*gpconvsigmaMat*(lambda/(2.0_rp*lambda0))
      gpdeformMat=aux*gpdeformMat*(lambda/(2.0_rp*lambda0))
      Prod1=matmul(gpdeformMat,transpose(grvel))
      Prod2=matmul(grvel,gpdeformMat)

      auxVEnew=auxVE/(1.0_rp-beta)
      
      do inode=1,e%pnode
         tensor=0.0_rp
         tensor= gpconvsigmaMat*auxVEnew*gpadv(inode)
         Shapi=e%shape(inode,e%igaus)
         tensor= -Prod1*auxVEnew*Shapi -auxVEnew*Prod2*Shapi +tensor(i,j)
         call GetStress(e%ndime,tensor,vector)
         elrhs(:,inode)= vector(:)+elrhs(:,inode)  
      end do

   end subroutine supm_elmrhcLCR_splitoss   
   
    
   !*********************************************************************************
   ! Subroutines Viscoelastic Split if kfl_repro==4 is employed (Split-OSS + Dynamics)
   !
   !*********************************************************************************
   
   subroutine supm_elmbstEst2_split(e,auxVE,tisig,gpadv,gpsig,dvolu,auxtens,grvel,elmat) 
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: auxVE,dvolu,tisig
      real(rp),    intent(in)    :: gpsig(auxtens),gpadv(e%pnode)
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime)     
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i, j, k, l
      real(rp)                   :: aux, aux_def
      real(rp)                   :: term(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term1(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term2(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term3(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term4(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: tensor(auxtens,auxtens)
 
      aux=tisig*dvolu*auxVE*auxVE    

      do jnode=1,e%pnode  
         do inode=1,e%pnode
            aux_def=aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            tensor=0.0_rp
            term=0.0_rp
            term1=0.0_rp
            term2=0.0_rp
            term3=0.0_rp
            term4=0.0_rp

            do i=1,e%ndime
               do j=1,e%ndime
                  !convective term
                  term(i,j,i,j)=aux*gpadv(jnode)*gpadv(inode)
                  
                  !deformation terms:
                  do k=1,e%ndime
                     do l=1,e%ndime 
                        term1(i,k,i,j)=-2.0_rp*aux_def*grvel(l,j)*grvel(k,l) + term1(i,k,i,j)
                        term2(i,l,j,k)=-2.0_rp*aux_def*grvel(i,j)*grvel(l,k)
!                         term3(l,k,i,j)=-aux_def*grvel(k,j)*grvel(l,i)
!                         term4(i,k,j,k)=-aux_def*grvel(l,j)*grvel(i,l) + term4(i,k,j,k)
                     end do
                  end do
                  
               end do
            end do   
            
            term=term+term1+term2+term3+term4
            call PassTensor4ToTensor2(e,auxtens,term,tensor)
            
            elmat(1:auxtens,inode,:,jnode) = tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode)
         end do
      end do   

   end subroutine supm_elmbstEst2_split
   
   subroutine supm_elmbutEst2_split(e,auxVE,tisig,gpsig,dvolu,auxtens,grvel,grsig,gpadv,elmat) 
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: auxVE,dvolu,tisig
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),gpadv(e%pnode)
      real(rp),    intent(in)    :: gpsig(auxtens),grsig(auxtens,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode, i, j, k, l
      real(rp)                   :: aux, aux_def    
      real(rp)                   :: grsigMat(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: gpsigMat(e%ndime,e%ndime)
      real(rp)                   :: term(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term1(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term2(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term3(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: term4(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: Tensor(auxtens,e%ndime)
      
      aux=tisig*dvolu*auxVE
      
      call PassGradStressToTensor3(e,auxtens,grsig,grsigMat)
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,gpsigMat)
      
      do jnode=1,e%pnode    
         do inode=1,e%pnode  
            aux_def=aux*auxVE*e%shape(inode,e%igaus)
            term=0.0_rp
            term1=0.0_rp
            term2=0.0_rp
            term3=0.0_rp
            term4=0.0_rp
            
            do i=1,e%ndime
               do j=1,e%ndime
                  do k=1,e%ndime
                     !Convective term (Newton Raphson)
                     term(i,j,k)=aux*gpadv(inode)*auxVE*e%shape(jnode,e%igaus)*grsigMat(i,j,k)
                     
                     do l=1,e%ndime
                        !Deformation terms (Newton Raphson)
                        term1(i,j,k)=-aux_def*2.0_rp*gpsigMat(i,l)*e%cartd(l,jnode)*grvel(j,k) + term1(i,j,k)
                        term2(i,j,i)=-aux_def*2.0_rp*gpsigMat(l,k)*e%cartd(l,jnode)*grvel(j,k) + term2(i,j,i)
!                         term3(i,j,j)=-aux_def*gpsigMat(l,k)*e%cartd(k,jnode)*grvel(i,l) + term3(i,j,j)
!                         term4(j,k,i)=-aux_def*gpsigMat(l,k)*e%cartd(l,jnode)*grvel(j,i) + term4(j,k,i)
                     end do
                     
                  end do
               end do
            end do   
            
            term=term+term1+term2+term3+term4
            
            call PassTensor3ToTensor2left(e,auxtens,term,Tensor)
            
            elmat(1:auxtens,inode,:,jnode) = Tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode) 
            
         end do  
      end do

   end subroutine supm_elmbutEst2_split 
   
   
   subroutine supm_elmrhcEst2_split(e,auxVE,tisig,gpsig,grsig,gpvel,gpadv,grvel,dvolu,auxtens,elrhs)         
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: tisig,dvolu,auxVE
      real(rp),    intent(in)    :: gpsig(auxtens),grsig(auxtens,e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k, l
      real(rp)                   :: aux, aux_def
      real(rp)                   :: grsigMat(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: gpsigMat(e%ndime,e%ndime)
      real(rp)                   :: term(e%ndime,e%ndime),vector(auxtens)
      real(rp)                   :: term1(e%ndime,e%ndime), term2(e%ndime,e%ndime)
      real(rp)                   :: term3(e%ndime,e%ndime), term4(e%ndime,e%ndime)

      aux=dvolu*tisig*auxVE
      
      call PassGradStressToTensor3(e,auxtens,grsig,grsigMat)
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,gpsigMat)
      
      do inode=1,e%pnode  
         aux_def=aux*auxVE*e%shape(inode,e%igaus)
         term=0.0_rp
         term1=0.0_rp
         term2=0.0_rp
         term3=0.0_rp
         term4=0.0_rp
         
         do i=1,e%ndime
            do j=1,e%ndime
               do k=1,e%ndime
                  !Convective term (Newton Raphson)
                  term(i,j)=auxVE*gpvel(k)*grsigMat(i,j,k)*aux*gpadv(inode)+term(i,j)
                  
                  do l=1,e%ndime
                     !Deformation terms (Newton Raphson)
                     term1(i,j)=-aux_def*2.0_rp*gpsigMat(i,l)*grvel(k,l)*grvel(j,k) + term1(i,j)
                     term2(i,j)=-aux_def*2.0_rp*gpsigMat(l,k)*grvel(l,i)*grvel(j,k) + term2(i,j)
                  end do
                  
               end do
            end do
         end do
         
         term=term+term1+term2+term3+term4

         call GetStress(e%ndime,term,vector)
         
         elrhs(:,inode)= vector(:) + elrhs(:,inode)
      end do

   end subroutine supm_elmrhcEst2_split
   
   
   subroutine supm_elmbuvEst2_split(e,beta,tisig,gpsig,dvolu,auxtens,grvel,grsig,elmat) 
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,tisig,beta
      real(rp),    intent(in)    :: gpsig(auxtens),grsig(auxtens,e%ndime)
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode, jnode, i, j, k
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,e%ndime)
      
      aux=tisig*dvolu*(1.0_rp-beta)
      
      do jnode=1,e%pnode   
         do inode=1,e%pnode
            tensor=0.0_rp
         
            do i=1,e%ndime
               do k=1,e%ndime
                  tensor(i,i) = 0.5_rp*0.5_rp*e%cartd(k,jnode)*(e%cartd(k,inode))*aux + tensor(i,i)
                  tensor(k,i) = 0.5_rp*0.5_rp*e%cartd(k,jnode)*(e%cartd(i,inode))*aux + tensor(k,i)

                  tensor(i,k) = 0.5_rp*0.5_rp*e%cartd(i,jnode)*(e%cartd(k,inode))*aux + tensor(i,k) 
                  tensor(k,k) = 0.5_rp*0.5_rp*e%cartd(i,jnode)*(e%cartd(i,inode))*aux + tensor(k,k)
               end do   
               
            end do

            elmat(1:e%ndime,inode,:,jnode)= tensor(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)
         end do  
      end do

   end subroutine supm_elmbuvEst2_split 

   
   
   !***********************************************************************************
   ! Dynamic Subscales
   !***********************************************************************************
   
    subroutine supm_elmrhc_dss(e,dtinv,dvolu,acden,beta,timom,auxtens,vesgs,elrhs)
      implicit none
      class(FiniteElement) :: e
      !tau1*(div(v), tilde(u)*dtinv)

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: vesgs(*)
      real(rp),    intent(in)    :: dvolu, timom, dtinv, acden, beta
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,aux2, itens
      real(rp)                   :: aux, tensor(e%ndime,e%ndime), vector(auxtens)
      !Compute contributions to RHS : Block S
      aux  = -timom*dvolu*dtinv*acden*(1.0_rp-beta)
      do inode=1,e%pnode
         do itens=1, e%ndime
               tensor(itens,1:e%ndime)=e%cartd(itens,inode)*vesgs(1:e%ndime)*aux
         end do  
         
         call GetStress(e%ndime,tensor,vector)  
         elrhs(:,inode)=vector(:)+elrhs(:,inode)
      end do

   end subroutine supm_elmrhc_dss
   
   
   !***********************************************************************************
   ! Viscoelastic Dynamic Subscales
   !***********************************************************************************
    
   subroutine supm_elmrhcVES_dss(e,auxVE,beta,tisig,acvis,dvolu,auxtens,gpadv,grvel,sisgs,dtinv,auxLCR,elrhs)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxLCR
      real(rp),    intent(in)    :: sisgs(auxtens),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: tisig, dvolu,auxVE,acvis,dtinv, beta
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k
      real(rp)                   :: aux, aux1, Shapi
      real(rp)                   :: sisgsMat(e%ndime,e%ndime)
      real(rp)                   :: tensor(e%ndime,e%ndime), Prod1(e%ndime,e%ndime)
      real(rp)                   :: Prod2(e%ndime,e%ndime)
      real(rp)                   :: vector(auxtens)
      
      sisgsMat=0.0_rp
      call sup_SigmaMatrix(e%ndime,auxtens,sisgs,sisgsMat)
      aux=tisig*dvolu
      aux1=(1.0_rp/(2.0_rp*acvis*(1.0_rp-beta*auxLCR)))
      Prod1=matmul(sisgsMat,transpose(grvel))
      Prod2=matmul(grvel,sisgsMat)
      
      do inode=1,e%pnode
         tensor=0.0_rp
         vector=0.0_rp
         Shapi=e%shape(inode,e%igaus)
         
         tensor = -auxVE*dtinv*sisgsMat*aux1*Shapi&
                  +auxVE*dtinv*sisgsMat*auxVE*gpadv(inode)&
                  +auxVE*dtinv*Prod1*auxVE*Shapi&
                  +auxVE*dtinv*Prod2*auxVE*Shapi+tensor
      
         call GetStress(e%ndime,tensor,vector)
         elrhs(:,inode)= aux*vector(:)+elrhs(:,inode)  

      end do

   end subroutine supm_elmrhcVES_dss   
   
   
   
   subroutine supm_elmrhuVES_dss(e,auxVE,tisig,dvolu,sisgs,auxtens,dtinv,elrhs) 
      use typre
      implicit none

      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: sisgs(auxtens),auxVE
      real(rp),    intent(in)    :: dvolu,tisig, dtinv
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode,idime, i, j
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,auxtens), vector(e%ndime)
      real(rp)                   :: gprep_Mat(e%ndime,e%ndime)
      real(rp)                   :: Term(e%ndime)

      aux = -tisig*dvolu*dtinv*auxVE
      do inode=1,e%pnode
         vector=0.0_rp
         
         do idime=1,e%ndime
            tensor(idime,1:auxtens)= aux*e%cartd(idime,inode)*sisgs(1:auxtens)
         end do
         call GetVel_hu(e%ndime,tensor,vector)

!         call sup_SigmaMatrix(e%ndime,auxtens,gprep,gprep_Mat)
!          Term=0.0_rp
!          do j=1,e%ndime
!             do i=1,e%ndime
!                Term(j)= aux*0.5_rp*gprep_Mat(i,j)*e%cartd(i,inode)+ Term(j)
!                Term(i)= aux*0.5_rp*gprep_Mat(i,j)*e%cartd(j,inode)+ Term(i)
!             end do
!          end do
      
         elrhs(1:e%ndime,inode)= vector(1:e%ndime) + elrhs(1:e%ndime,inode)
         
      end do

   end subroutine supm_elmrhuVES_dss
   
   
   subroutine supm_elmrhcVES_dss_split2(e,auxVE,tisig,dvolu,auxtens,gpadv,sisgs,dtinv,elrhs)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: sisgs(*),gpadv(e%pnode)
      real(rp),    intent(in)    :: tisig, dvolu,auxVE,dtinv
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k
      real(rp)                   :: aux
      real(rp)                   :: sisgsMat(e%ndime,e%ndime)
      real(rp)                   :: tensor2(e%ndime,e%ndime)
      real(rp)                   :: vector(auxtens)
      
      sisgsMat=0.0_rp
      call sup_SigmaMatrix(e%ndime,auxtens,sisgs,sisgsMat)
      aux=tisig*auxVE*dtinv*dvolu
      sisgsMat=aux*sisgsMat

      do inode=1,e%pnode
         tensor2=  sisgsMat*auxVE*gpadv(inode)
         call GetStress(e%ndime,tensor2,vector)
         elrhs(:,inode)= vector(:)+elrhs(:,inode)  
      end do

   end subroutine supm_elmrhcVES_dss_split2   
   
   
   subroutine supm_elmrhcVES_dss_split3(e,auxVE,tisig,dvolu,auxtens,grvel,sisgs,dtinv,elrhs)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: sisgs(auxtens),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: tisig, dvolu,auxVE,dtinv
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k
      real(rp)                   :: aux, Shapi
      real(rp)                   :: sisgsMat(e%ndime,e%ndime)
      real(rp)                   :: tensor(e%ndime,e%ndime),tensor4(e%ndime,e%ndime)
      real(rp)                   :: tens_sum(e%ndime,e%ndime), vector(auxtens)
      real(rp)                   :: Prod1(e%ndime,e%ndime)
      real(rp)                   :: Prod2(e%ndime,e%ndime)
      
      sisgsMat=0.0_rp
      call sup_SigmaMatrix(e%ndime,auxtens,sisgs,sisgsMat)
      aux=tisig*auxVE*dtinv*dvolu
      sisgsMat=aux*sisgsMat
      
      Prod1=matmul(sisgsMat,transpose(grvel))
      Prod2=matmul(grvel,sisgsMat)

      
      do inode=1,e%pnode
         tensor=0.0_rp
         Shapi=e%shape(inode,e%igaus)
         tensor= Prod1*auxVE*Shapi+Prod2*auxVE*Shapi+tensor
         
         call GetStress(e%ndime,tensor,vector)
         
         elrhs(:,inode)= vector(:)+elrhs(:,inode)  

      end do

   end subroutine supm_elmrhcVES_dss_split3   
   



end module            


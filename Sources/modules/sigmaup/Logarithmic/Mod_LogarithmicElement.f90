module Mod_LogarithmicElement
   use typre
   use Mod_Element
   use Mod_ThreeField
   use Mod_supm_StressGenerator
   use Mod_LogOperations
   use Mod_SupOperations
   
   implicit none
   
contains

   ! Elemental components in the Three Field case
   ! Voigt notation 2d=>  s11 s22 s12  3d => s11 s22 s33 s23 s13 s12
   
   !*****************************************************************************************************************************
   !Three Field viscoelastic 
   !*****************************************************************************************************************************  

    subroutine sup_LCR2_TimeIntegrationToElext(e,auxtens,nsteps,Integrator,dtinv,lambda,lambda0,gpsig,elextS)
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      type(TimeIntegratorDt1)    :: Integrator
      integer(ip), intent(in)    :: auxtens, nsteps
      real(rp)                   :: gpsig(auxtens,*),elextS(auxtens), auxL
      real(rp)                   :: dtinv,lambda,lambda0    
      real(rp)                   :: gprhs(auxtens),expo_gprhs(auxtens,3)
      real(rp)                   :: expomatrix(e%ndime,e%ndime),expomatrix2(e%ndime,e%ndime)
      integer(ip)                :: i, j

      !Time integration     
      expo_gprhs=0.0_rp
      
      call sup_ComputeExponential(e%ndime,auxtens,gpsig(:,2),expo_gprhs(:,1),expomatrix)
      call sup_ComputeExponential(e%ndime,auxtens,gpsig(:,3),expo_gprhs(:,2),expomatrix2)
      call Integrator%GetRHS(auxtens,expo_gprhs(:,1:2),gprhs)
      
      elextS(1:auxtens) = elextS(1:auxtens) + dtinv*(lambda/(2.0_rp*lambda0))*gprhs(1:auxtens)
      
   end subroutine sup_LCR2_TimeIntegrationToElext

   !Galerkin terms 
    subroutine  supm_LCR_elmbstGal(e,lambda,lambda0,auxG,dtinv,dvolu,auxtens,gpadv,gpvel,grvel,ExpPosMatrix,ExpGiesekus,GrPsiMatrix,GrExpMatrix,&
                                    ConvExpMatrix,GrvelExp,elmat)
    !-----------------------------------------------------------------------
    ! This routine computes the Galerkin component of elmst matrix
    ! (T , 1/(2*lambda0) * exp(S) + (lambda/2lambda0)* (u*grad(exp(S)) - exp(S)*grad(u) - grad_t(u)*exp(S)) )
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,dtinv,lambda,lambda0,auxG
      real(rp),    intent(in)    :: gpvel(e%ndime),grvel(e%ndime,e%ndime),gpadv(e%pnode),ConvExpMatrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime),GrPsiMatrix(e%ndime,e%ndime,e%ndime),GrvelExp(e%ndime,e%ndime)
      real(rp),    intent(in)    :: GrExpMatrix(e%ndime,e%ndime,e%ndime),ExpGiesekus(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i,j,k,l
      real(rp)                   :: aux, aux_lamb ,lamb_div, Shapij, ConvjShapi
      real(rp)                   :: tensor(auxtens,auxtens), conv(e%pnode)
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime,e%ndime)

      aux=dvolu
      aux_lamb=1.0_rp/(2.0_rp*lambda0)
      lamb_div=lambda*aux_lamb
      
      conv=0.0_rp
      do i=1,e%ndime
         conv(:)=gpvel(i)*e%cartd(i,:)+ conv(:)
      end do

      do jnode=1,e%pnode  
         do inode=1,e%pnode
            Term=0.0_rp
            tensor=0.0_rp
            Shapij=e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)
            ConvjShapi=conv(jnode)*e%shape(inode,e%igaus)

            do concurrent (j=1:e%ndime, k=1:e%ndime)
               
               Term(:,j,k,j)=lamb_div*ExpPosMatrix(:,k)*ConvjShapi& 
                             +lamb_div*ConvExpMatrix(:,k)*Shapij& 
                             -lamb_div*GrvelExp(:,k)*Shapij&
                             +(aux_lamb*ExpPosMatrix(:,k) + lamb_div*dtinv*ExpPosMatrix(:,k)&
                             + auxG*ExpGiesekus(:,k))*Shapij + Term(:,j,k,j)
                             
                             
               do concurrent (l=1:e%ndime)       
                  Term(:,l,j,k)=-lamb_div*ExpPosMatrix(:,j)*grvel(l,k)*Shapij + Term(:,l,j,k)
               end do 
            end do

            call PassTensor4ToTensor2(e,auxtens,Term,tensor)
            elmat(1:auxtens,inode,:,jnode) = aux*(tensor(1:auxtens,:)) + elmat(1:auxtens,inode,:,jnode)
            
         end do
      end do   

   end subroutine supm_LCR_elmbstGal  

   subroutine supm_LCR_elmbutGal(e,dvolu,lambda,lambda0,auxtens,gpsig,ExpPosMatrix,GrExpPsi,linearization,elmat)
    !-----------------------------------------------------------------------
    ! This routine computes the Galerkin component of elmut matrix 
    ! (T, -grad_sym(u) + Newton-Raphson terms)
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens, linearization
      real(rp),    intent(in)    :: dvolu,lambda,lambda0
      real(rp),    intent(in)    :: gpsig(auxtens),ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: GrExpPsi(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,i,j,k
      real(rp)                   :: Shapi, auxlamb, aux
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: Tensor(auxtens,e%ndime)
      real(rp)                   :: ExpCartd(e%ndime,e%pnode), CartdExp(e%pnode,e%ndime)
      
      auxlamb=lambda/(lambda0*2.0_rp)
      aux=-1.0_rp +2.0_rp*auxlamb*linearization 
      
      ExpCartd=matmul(ExpPosMatrix,e%cartd)
      CartdExp=matmul(transpose(e%cartd),ExpPosMatrix) 
      do jnode=1,e%pnode
         do inode=1,e%pnode
            Shapi=e%shape(inode,e%igaus)*dvolu
            Tensor=0.0_rp
            Term=0.0_rp
            
            do concurrent (j=1:e%ndime)
               Term(:,j,j)=0.5_rp*e%cartd(:,jnode)*Shapi*aux&
                           !Newton-Raphson
                           -auxlamb*ExpCartd(:,jnode)*Shapi*linearization + Term(:,j,j)  
               Term(j,:,j)=0.5_rp*e%cartd(:,jnode)*Shapi*aux& 
                           !Newton-Raphson
                           -auxlamb*CartdExp(jnode,:)*Shapi*linearization + Term(j,:,j)
               do k=1,e%ndime
                  !Newton-Raphson
                  Term(:,j,k)=auxlamb*e%shape(jnode,e%igaus)*GrExpPsi(:,j,k)*Shapi*linearization + Term(:,j,k) 
               end do
            end do

            call PassTensor3ToTensor2left(e,auxtens,Term,Tensor)
            elmat(1:auxtens,inode,:,jnode) = Tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode)            
                                                
         end do  
      end do

   end subroutine supm_LCR_elmbutGal
   
   subroutine supm_LCR_elmrhcGal(e,auxG,dvolu,auxtens,elextS,ExpGiesekus,ExpPosMatrix,gpsig,RHSlin,RHSGiesekus,elrhs)
    !-----------------------------------------------------------------------
    ! This routine computes the Galerkin component of elmrhc matrix,
    ! including rhs terms from linearizations and force terms.
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu, auxG
      real(rp),    intent(in)    :: elextS(auxtens), gpsig(auxtens), ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpGiesekus(e%ndime,e%ndime),RHSlin(e%ndime,e%ndime),RHSGiesekus(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: Shapi
      real(rp)                   :: elextSMat(e%ndime,e%ndime),vector(auxtens)
      real(rp)                   :: RHS(e%ndime,e%ndime)
      
      RHS=0.0_rp
      !Compute forces
      call sup_SigmaMatrix(e%ndime,auxtens,elextS,elextSMat)
      
      !Add all RHS terms including forces
      RHS=elextSMat+RHSlin+RHSGiesekus
      
      call GetStress(e%ndime,RHS,vector)
     
      do inode=1,e%pnode
         Shapi=e%shape(inode,e%igaus)*dvolu
         elrhs(:,inode) = Shapi*vector(:) + elrhs(:,inode) 
      end do

   end subroutine supm_LCR_elmrhcGal
   
   
   subroutine supm_LCR_elmbsvGal(e,dvolu,auxtens,auxL,DivExpPsi,ExpPosMatrix,elmat)
    !-----------------------------------------------------------------------
    ! This routine computes the Galerkin component of elmbsv matrix.
    ! auxL*(exp(s)*S , grad_sym(V))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu, auxL
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime), DivExpPsi(e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i,j,k
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime), tensor(e%ndime,auxtens)
      real(rp)                   :: Shapj, Prod(e%ndime,e%pnode)
      
      Prod=0.0_rp
      do concurrent (i=1:e%ndime, j=1:e%ndime)
         Prod(j,:)=ExpPosMatrix(i,j)*e%cartd(i,:) + Prod(j,:)
      end do 

      do inode=1,e%pnode       
         do jnode=1,e%pnode
            Term=0.0_rp
            tensor=0.0_rp
            Shapj=e%shape(jnode,e%igaus)
            
            do concurrent (k=1:e%ndime)
               Term(k,:,k)=auxL*Shapj*Prod(:,inode)*0.5_rp + Term(k,:,k)
               do concurrent (j=1:e%ndime)
                  Term(:,j,k)=auxL*(ExpPosMatrix(:,j)*Shapj*e%cartd(k,inode))*0.5_rp + Term(:,j,k)
               end do      
            end do
            
            call PassTensor3ToTensor2right(e,auxtens,Term,tensor)
            elmat(1:e%ndime,inode,:,jnode) = dvolu*tensor(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)     
         end do  
      end do
 
   end subroutine supm_LCR_elmbsvGal
   
   subroutine supm_LCR_elmrhuGal(e,auxtens,auxL,acden,dvolu,elext,grvel,gpvel,gpsig,ExpGpPsi_Matrix,RHS1,RHS2,elrhs)
    !-----------------------------------------------------------------------
    ! This routine computes the rhs terms in momentum equation and force terms    
    ! (f + rho a* grad(a), V)  + (exp(S), grad_sym(V))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: acden,auxL,dvolu
      real(rp),    intent(in)    :: elext(e%ndime),grvel(e%ndime,e%ndime),gpvel(e%ndime),gpsig(auxtens)
      real(rp),    intent(in)    :: ExpGpPsi_Matrix(e%ndime,e%ndime), RHS1(e%ndime),RHS2(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: jtens,inode,Newton,i,j
      real(rp)                   :: aux 
      real(rp)                   :: Term(e%ndime)
         
      do inode=1,e%pnode  
         aux = e%shape(inode,e%igaus)*dvolu
         Term=0.0_rp
         do concurrent (i=1:e%ndime)
            Term(:)= (RHS2(i,:)+RHS2(:,i))*e%cartd(i,inode)*0.5_rp + Term(:)
         end do  

         elrhs(:,inode) = dvolu*Term(:) +  aux*RHS1(:) + elrhs(:,inode)
      end do

   end subroutine supm_LCR_elmrhuGal 
   
   
   subroutine supm_LCR_elmbpvGal(e,dvolu,auxL,elmat) 
    !-----------------------------------------------------------------------
    ! This routine computes the Galerkin component of elmbpv matrix.
    ! -(p,div v) 
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,auxL
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux
          
      do jnode=1,e%pnode
         aux=-dvolu*(e%shape(jnode,e%igaus))  
         do inode=1,e%pnode
            elmat(:,inode,1,jnode) = aux*e%cartd(:,inode) + elmat(:,inode,1,jnode)  
         end do
      end do

   end subroutine supm_LCR_elmbpvGal 
   
   
   
 !**********************************************************************************************************************************
 !Estabilization terms of Momentum and Continuity equation 
   
   subroutine supm_LCR_elmbstEst1(e,timom,dvolu,auxL,gpsig,auxtens,DivExp,ExpPosMatrix,elmat)
    !-----------------------------------------------------------------------
    ! This routine computes:
    ! tau1*( -auxL * div(Exp(S),- div( T))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu,auxL
      real(rp),    intent(in)    :: gpsig(auxtens),DivExp(e%ndime) 
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i,j,k,l
      real(rp)                   :: aux
      real(rp)                   :: tensor(auxtens,auxtens)
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp)                   :: Prod1(e%ndime,e%pnode)
      
 
      aux= timom*dvolu
      Prod1=matmul(ExpPosMatrix,e%cartd)
      
      do jnode=1,e%pnode    
         do inode=1,e%pnode
            Term=0.0_rp
            tensor=0.0_rp

            do concurrent (j=1:e%ndime, k=1:e%ndime)
               Term(:,k,j,k)= (-auxL)*(DivExp(j)*e%shape(jnode,e%igaus))*(-e%cartd(:,inode))&
                              +(-auxL)*Prod1(j,jnode)*(-e%cartd(:,inode)) + Term(:,k,j,k)                
            end do
            
            call PassTensor4ToTensor2(e,auxtens,Term,tensor)
            elmat(1:auxtens,inode,:,jnode)= aux*tensor(1:auxtens,:)+ elmat(1:auxtens,inode,:,jnode)
               
         end do
      end do   

   end subroutine supm_LCR_elmbstEst1 
   
   
   subroutine supm_LCR_elmbutEst1(e,acden,timom,dtinv,dvolu,auxtens,linearization,gpadv,grvel,elmat)
    !-----------------------------------------------------------------------
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! (-div(T) , p u*grad(u))
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens, linearization
      real(rp),    intent(in)    :: timom,dvolu,acden,dtinv
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,k,i,j
      real(rp)                   :: tmp,aux
      real(rp)                   :: tensor(auxtens,e%ndime)
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime)
      
      tmp = acden*dtinv 
      aux = timom*dvolu

      do jnode=1,e%pnode      
         do inode=1,e%pnode
            Term=0.0_rp
            tensor=0.0_rp
               do j=1,e%ndime
                  Term(:,j,j)=(acden*gpadv(jnode)+tmp*e%shape(jnode,e%igaus))*(-e%cartd(:,inode)) + Term(:,j,j)
                  do k=1,e%ndime
                     Term(k,j,:)=acden*e%shape(jnode,e%igaus)*grvel(j,:)*(-e%cartd(k,inode))*linearization + Term(k,j,:)
                  end do
               end do
            
            call PassTensor3ToTensor2Left(e,auxtens,Term,tensor) 
            
            elmat(1:auxtens,inode,:,jnode) = aux*tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode)                                                                       
         end do  
      end do
     
   end subroutine supm_LCR_elmbutEst1
   
   subroutine supm_LCR_elmbptEst1(e,timom,dvolu,auxtens,elmat)
    !-----------------------------------------------------------------------
    !  (-div(T) , grad(p))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,i,j
      real(rp)                   :: aux
      real(rp)                   :: Term(e%ndime,e%ndime), vector(auxtens)
      
      aux= timom*dvolu    
      do jnode=1,e%pnode  
         do inode=1,e%pnode
            Term=0.0_rp
            do concurrent (i=1:e%ndime,j=1:e%ndime)
               Term(j,i)=-e%cartd(i,jnode)*e%cartd(j,inode)
            end do
            
            call GetStress(e%ndime,Term,vector)
            elmat(:,inode,1,jnode) = aux*vector(:) + elmat(:,inode,1,jnode)
         
         end do  
      end do

   end subroutine supm_LCR_elmbptEst1 
   

  subroutine supm_LCR_elmrhcEst1(e,acden,timom,dvolu,elext,auxtens,gpvel,grvel,auxoss,linearization,RHSlin,elrhs)
    !-----------------------------------------------------------------------------------
    !
    !    tau1*(-div(T), f + rho*a*grad(a) - auxL*(Div( exp(s)*(s-I)))
    !
    !-----------------------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss,linearization
      real(rp),    intent(in)    :: dvolu,timom,acden
      real(rp),    intent(in)    :: elext(e%ndime),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: RHSlin(e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,i,j,Newton
      real(rp)                   :: aux
      real(rp)                   :: tensor(e%ndime,e%ndime),vector(auxtens)
      real(rp)                   :: RHS(e%ndime),convec(e%ndime)
      
      aux=timom*dvolu
      RHS=0.0_rp

      !Newton-Raphson to convective term
      convec=0.0_rp
      do i=1,e%ndime
         convec(:)=acden*gpvel(i)*grvel(:,i)*auxoss*linearization+convec(:)
      end do
      
      !Add with forces and linearization of div(exp(S))
      RHS=convec+elext*auxoss+RHSlin
      
      do inode=1,e%pnode
         tensor=0.0_rp
         do concurrent (i=1:e%ndime, j=1:e%ndime)
            tensor(i,j) = -e%cartd(i,inode)*RHS(j)
         end do
         call GetStress(e%ndime,tensor,vector)
         elrhs(:,inode) = aux*vector(:) + elrhs(:,inode)  
      end do

   end subroutine supm_LCR_elmrhcEst1
   
   
   subroutine supm_LCR_elmbsvEst1(e,timom,dvolu,acden,auxL,gpadv,DivExpPsi,ExpPosMatrix,auxtens,elmat)
    !-----------------------------------------------------------------------
    !    tau1*(-(auxL*Div(exp(s)*S), rho*a*grad(V)))
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu,acden, auxL
      real(rp),    intent(in)    :: gpadv(e%pnode), DivExpPsi(e%ndime), ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i,j,k
      real(rp)                   :: aux1,aux
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime), tensor(e%ndime,auxtens),Prod(e%ndime,e%pnode)

      aux=timom*dvolu
      
      Prod=0.0_rp 
      Prod=matmul(ExpPosMatrix,e%cartd)
      
      do inode=1,e%pnode
         aux1 = acden*gpadv(inode)       
         do jnode=1,e%pnode
            tensor=0.0_rp
            Term=0.0_rp
            do j=1,e%ndime
               Term(j,:,j)=-auxL*aux1*(DivExpPsi(:)*e%shape(jnode,e%igaus)+Prod(:,jnode)) + Term(j,:,j)
            end do
         
            call PassTensor3ToTensor2Right(e,auxtens,Term,tensor)
            elmat(1:e%ndime,inode,:,jnode) = aux*tensor(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)                                   
      
         end do  
      end do
   end subroutine supm_LCR_elmbsvEst1
   
   
   subroutine supm_LCR_elmbsvEst1Lapla(e,beta,acvis,auxL,timom,dvolu,auxtens,auxoss,DivExpPsi,ExpPosMatrix,elmat)   
    !-----------------------------------------------------------------------
    !    tau1*(-(auxL*Div(exp(s)*S), 2*beta*eta*lapla(V))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: acvis,dvolu,timom,beta,auxL,DivExpPsi(e%ndime), ExpPosMatrix(e%ndime,e%ndime)  
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i, j, k, l
      real(rp)                   :: aux ,auxvis 
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime), tensor(e%ndime,auxtens)
      real(rp)                   :: Prod(e%ndime,e%pnode)
      
      Prod=0.0_rp
      Prod=matmul(ExpPosMatrix,e%cartd)
   
      aux = timom*dvolu
      do inode=1,e%pnode
         auxvis= sum(e%hessi(1:e%ndime,inode))*(beta*acvis)
         do jnode=1,e%pnode
            Term=0.0_rp
            do j=1,e%ndime
               Term(j,:,j)=-auxvis*auxL*(DivExpPsi(:)*e%shape(jnode,e%igaus))*auxoss&
                           -auxvis*auxL*Prod(:,jnode)*auxoss + Term(j,:,j)
            end do
            call PassTensor3ToTensor2Right(e,auxtens,Term,tensor)

            elmat(1:e%ndime,inode,:,jnode) = aux*tensor(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)  
          
         end do
      end do   

   end subroutine supm_LCR_elmbsvEst1Lapla  
   
   subroutine supm_LCR_elmrhuEst1(e,timom,tidiv,acden,auxL,auxtens,dvolu,elextC,elext,gpvel,grvel,gpadv,auxoss,linearization,RHSlin,elrhs)
    !-----------------------------------------------------------------------
    !
    !  tau1* (f + rho*a*grad(a) + auxL*Div(exp(s)*(I)), rho*a*grad(V)) + tau2* (f,grad(V))

    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss,linearization
      real(rp),    intent(in)    :: auxL,elextC(1),elext(e%ndime),gpadv(e%pnode),acden
      real(rp),    intent(in)    :: dvolu,tidiv,timom,gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: RHSlin(e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,i
      real(rp)                   :: aux,aux1, Solu(e%ndime), convec(e%ndime), RHS(e%ndime)
      real(rp)                   :: vect_div(e%ndime), vect_mom(e%ndime)
    
      RHS=0.0_rp
      !Newton-Raphson to convective term
      convec=0.0_rp
      do i=1,e%ndime
         convec(:)=acden*gpvel(i)*grvel(:,i)*linearization+convec(:)
      end do
   
      !Add force terms and linearization terms of div(exp(S))
      RHS=convec+elext*auxoss+RHSlin*auxoss
      
      aux = tidiv*dvolu 
      aux1= timom*dvolu
      
      do inode=1,e%pnode
         vect_div(:)=e%cartd(:,inode)*aux*elextC(1)
         vect_mom(:)=RHS(:)*(acden*gpadv(inode))*aux1
         elrhs(:,inode) = vect_div(:)+vect_mom(:)+elrhs(:,inode)
      end do

   end subroutine supm_LCR_elmrhuEst1

   subroutine supm_LCR_elmrhuEst1Lapla(e,beta,auxtens,timom,acden,acvis,dvolu,elext,grvel,gpvel,auxoss,convection,RHSlin,elrhs)   
    !-----------------------------------------------------------------------
    !
    !   tau1* (f + rho*a*grad(a) - auxL*(Div( exp(s)*(I-s))), beta*eta*lapla(V))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxoss, auxtens, convection
      real(rp),    intent(in)    :: dvolu,acvis,beta,timom,acden,gpvel(e%ndime)
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),elext(e%ndime)
      real(rp),    intent(in)    :: RHSlin(e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,i
      real(rp)                   :: aux,auxvis
      real(rp)                   :: RHS(e%ndime), convec(e%ndime), vector(e%ndime)
           
      aux= timom*dvolu
      RHS=0.0_rp

      !Newton-Raphson to convective term
      convec=0.0_rp
      do i=1,e%ndime
         convec(:)=acden*gpvel(i)*grvel(:,i)*auxoss*convection+convec(:)
      end do
      
      !Add force terms and linearization terms of div(exp(S))
      RHS=convec+elext*auxoss+RHSlin*auxoss

      do inode=1,e%pnode
         auxvis=beta*acvis*sum(e%hessi(1:e%ndime,inode))
         vector=auxvis*RHS
    
         elrhs(:,inode) = aux*vector(:)+elrhs(:,inode)  
      end do

   end subroutine supm_LCR_elmrhuEst1Lapla     
   
   
   subroutine supm_LCR_elmbsqEst1(e,auxL,timom,dvolu,auxtens,DivExpPsi,ExpGradS,elmat)
    !-----------------------------------------------------------------------
    ! tau* (-auxL*Div(exp(s)*S), grad(q))
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: timom,dvolu,auxL
      real(rp),    intent(in)    :: DivExpPsi(e%ndime), ExpGradS(e%ndime,e%mnode)
      real(rp),    intent(inout) :: elmat(1,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i, j
      real(rp)                   :: aux
      real(rp)                   :: Term(e%ndime,e%ndime), tensor(auxtens)

      aux=timom*dvolu   
      do jnode=1,e%pnode      
         do inode=1,e%pnode   
            Term=0.0_rp

            do j=1,e%ndime
               Term(:,j)= -auxL*(DivExpPsi(:)*e%shape(jnode,e%igaus)+ExpGradS(:,jnode))*e%cartd(j,inode)
            end do

            call GetStress(e%ndime,Term,tensor)
            elmat(1,inode,:,jnode) = aux*tensor(:) + elmat(1,inode,:,jnode)                    
         end do  
      end do

   end subroutine supm_LCR_elmbsqEst1
   
   
   subroutine supm_LCR_elmrhpEst1(e,auxtens,acden,timom,dvolu,elext,gpvel,grvel,gpsig,convection,auxoss,RHSlin,elrhs)
    !-----------------------------------------------------------------------
    !   tau1* (rho*a*grad(a) - auxL*(Div( exp(s))), Div(q))
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens, convection, auxoss
      real(rp),    intent(in)    :: acden,timom,dvolu
      real(rp),    intent(in)    :: elext(e%ndime),gpsig(auxtens),gpvel(e%ndime),grvel(e%ndime,e%ndime)  
      real(rp),    intent(in)    :: RHSlin(e%ndime)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,itens, i
      real(rp)                   :: value, aux , convec(e%ndime)
      real(rp)                   :: RHS(e%ndime)

      
      aux = timom*dvolu
      RHS=0.0_rp

      !Newton-Raphson to convective term
      convec=0.0_rp
      do i=1,e%ndime
         convec(:)=acden*gpvel(i)*grvel(:,i)*convection*auxoss+convec(:)
      end do
      
      !Add force terms and linearization terms of div(exp(S))
      RHS=convec+elext*auxoss+RHSlin*auxoss
      
      do inode=1,e%pnode
         value=0.0_rp
         do itens=1,e%ndime
            value=e%cartd(itens,inode)*RHS(itens)+value
         end do
         elrhs(1,inode) = aux*value + elrhs(1,inode) 
      end do

   end subroutine supm_LCR_elmrhpEst1 
   
   
   !**********************************************************************************************************************************
   !Estabilization terms of Constitutive Equation 
   
   subroutine supm_LCR_elmbstEst2(e,acvis,beta,lambda,lambda0,auxG,tisig,dtinv,gpadv,gpvel,grvel,gpsig,&
                                    ExpPosMatrix,ExpGiesekus,GrExpMatrix,ConvExpMatrix,GrVelExp,dvolu,auxtens,auxoss,elmat) 

    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !
    ! auxLCR=1  tau3*[ 1/lambda*(exp(-s)·S), -(1/lambda)exp(-s)·T)]
    !
    !-----------------------------------------------------------------------
    
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: acvis,dvolu,dtinv,beta, auxG, tisig,lambda0,lambda
      real(rp),    intent(in)    :: gpsig(auxtens),gpadv(e%pnode),gpvel(e%ndime),grvel(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime), ExpGiesekus(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: ConvExpMatrix(e%ndime,e%ndime), GrVelExp(e%ndime,e%ndime)
      real(rp),    intent(in)    :: GrExpMatrix(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode,i,j,k,s,l,t
      real(rp)                   :: aux, aux_lamb, lamb_div, aux_visc, visc_div, Shapij, ShapjAdvi, ShapiAdvj
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime,e%ndime), tensor(auxtens,auxtens)
      real(rp)                   :: ExpminusId(e%ndime,e%ndime),conv(e%pnode), MatrixSumExp(e%ndime,e%ndime)
      real(rp)                   :: Prod1(e%ndime,e%ndime),Prod1Giesekus(e%ndime,e%ndime)
      real(rp)                   :: Prod2(e%ndime,e%ndime),Prod2Giesekus(e%ndime,e%ndime)
      real(rp)                   :: Prod3(e%ndime,e%ndime),Prod3Giesekus(e%ndime,e%ndime)
      real(rp)                   :: Prod4(e%ndime,e%ndime),Prod4Giesekus(e%ndime,e%ndime)
      real(rp)                   :: Grvel2(e%ndime,e%ndime)
             
       
      aux=dvolu*tisig
      
      aux_lamb=1.0_rp/(2.0_rp*lambda0)
      aux_visc=1.0_rp/(2.0_rp*acvis*(1.0_rp-beta))
      lamb_div=lambda*aux_lamb
      visc_div=lambda*aux_visc
      
      ExpminusId=0.0_rp
         
      conv=0.0_rp
      do i=1,e%ndime
         conv(:)= gpvel(i)*e%cartd(i,:) + conv(:)
      end do
      
      MatrixSumExp=aux_lamb*ExpPosMatrix+ lamb_div*dtinv*ExpPosMatrix+ auxG*ExpGiesekus
      Prod1=matmul(grvel,MatrixSumExp)
      Prod2=matmul(grvel,ExpPosMatrix)
      Prod3=matmul(grvel,ConvExpMatrix)
      Prod4=matmul(grvel,Prod2)
      Grvel2=matmul(grvel,grvel)
      
      Prod1Giesekus=0.0_rp
      Prod2Giesekus=0.0_rp
      Prod3Giesekus=0.0_rp
      Prod4Giesekus=0.0_rp
!       !GIESEKUSLCR  
!       if (auxG/=0.0_rp) then
!          ExpminusId=ExpPosMatrix
!          do i=1,e%ndime
!             ExpminusId(i,i)=ExpminusId(i,i)-1.0_rp
!          end do
!          ExpminusId=auxG*(lambda0/(acvis*(1.0_rp-beta)))*ExpminusId
!          Prod1Giesekus=matmul(ExpminusId,MatrixSumExp)
!          Prod2Giesekus=matmul(ExpminusId,ExpPosMatrix)
!          Prod3Giesekus=matmul(ExpminusId,ConvExpMatrix)
!          Prod4Giesekus=matmul(ExpminusId,Prod2)
!       end if   
      
      do jnode=1,e%pnode  
         do inode=1,e%pnode
            !Initializations
            Term=0.0_rp
            tensor=0.0_rp
            
            Shapij=e%shape(jnode,e%igaus)*e%shape(inode,e%igaus)
            ShapjAdvi=gpadv(inode)*e%shape(jnode,e%igaus)
            ShapiAdvj=conv(jnode)*e%shape(inode,e%igaus)
            
            do concurrent(j=1:e%ndime,k=1:e%ndime)
               
               Term(:,k,j,k)=&
               !(-T, Exp(S))
                              -(MatrixSumExp(:,j))*Shapij*aux_visc&
               !Cross (a*grad(T), Exp(S))
                              +(MatrixSumExp(:,j))*visc_div*ShapjAdvi*auxoss&
               !(a*grad(T), a*grad(Exp(S)))
                              +lamb_div*ExpPosMatrix(:,j)*conv(jnode)*visc_div*gpadv(inode)&
                              +lamb_div*ConvExpMatrix(:,j)*ShapjAdvi*visc_div&
               !Cross (-T, a*grad(Exp(S))
                              -lamb_div*ConvExpMatrix(:,j)*Shapij*aux_visc*auxoss&
                              -lamb_div*ExpPosMatrix(:,j)*aux_visc*ShapiAdvj*auxoss&
               !Cross terms: ( T*grad^t(u) + grad(u)*T , Exp(S))
                              +(Prod1(:,j)*visc_div + Prod1Giesekus(:,j))*Shapij*auxoss&
               !Cross (T*grad^t(u) + grad(u)*T , a*grad(Exp(S)))
                              +lamb_div*(visc_div*Prod2(:,j) + Prod2Giesekus(:,j))*ShapiAdvj*auxoss&
                              +Shapij*lamb_div*(visc_div*Prod3(:,j) + Prod3Giesekus(:,j))*auxoss&
               !Cross (-T, -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                              +lamb_div*Prod2(:,j)*aux_visc*Shapij*auxoss&
               !Cross  (a*grad(T), -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                              -lamb_div*Prod2(:,j)*visc_div*ShapjAdvi*auxoss&
               !( T*grad^t(u) + grad(u)*T , -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                              -visc_div*lamb_div*Shapij*(Prod4(:,j)+Prod4Giesekus(:,j)) + Term(:,k,j,k)
               
               do concurrent (s=1:e%ndime)
                  
                  Term(:,s,k,j) =&
                  !Cross terms: ( T*grad^t(u) + grad(u)*T , Exp(S))
                                 +(MatrixSumExp(:,k))*visc_div*grvel(s,j)*Shapij*auxoss&  
                  !Cross (T*grad^t(u) + grad(u)*T , a*grad(Exp(S)))
                                 +visc_div*lamb_div*Shapij*ConvExpMatrix(:,k)*grvel(s,j)*auxoss&
                                 +visc_div*lamb_div*ExpPosMatrix(:,k)*grvel(s,j)*ShapiAdvj*auxoss&
                  !Cross (-T, -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                                 +lamb_div*ExpPosMatrix(:,k)*grvel(s,j)*aux_visc*Shapij*auxoss&
                  !Cross  (a*grad(T), -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                                 -lamb_div*ExpPosMatrix(:,k)*grvel(s,j)*visc_div*ShapjAdvi*auxoss&
                  !( T*grad^t(u) + grad(u)*T , -Exp(s)*S*grad(U) - grad_sym(U)*Exp(s)*S)
                                 -visc_div*lamb_div*Shapij*grvel(s,j)*(Prod2(:,k)+Prod2Giesekus(:,k))&
                                 -visc_div*lamb_div*Shapij*Prod2(:,k)*grvel(s,j)&
                                 -visc_div*lamb_div*Shapij*ExpPosMatrix(:,k)*Grvel2(s,j) + Term(:,s,k,j)
               end do
            end do
            
            call PassTensor4ToTensor2(e,auxtens,Term,tensor)
            elmat(1:auxtens,inode,:,jnode)=aux*tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode)

         end do
      end do   
      
     
   end subroutine supm_LCR_elmbstEst2
   
   
    subroutine supm_LCR_elmbutEst2(e,lambda,lambda0,auxG,acvis,beta,tisig,gpsig,grvel,dvolu,auxtens,ExpPosMatrix,GrExpPsi,gpadv,auxoss,Newton,elmat) 

      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss, Newton
      real(rp),    intent(in)    :: tisig,lambda,lambda0,dvolu,beta, auxG
      real(rp),    intent(in)    :: gpsig(auxtens),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime), GrExpPsi(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(out)   :: elmat(auxtens,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode, jnode, i,j,k,t,l
      real(rp)                   :: aux,lamb_div, acvis, auxvisc, visc_div, aux2
      real(rp)                   :: Shai, Shaj
      real(rp)                   :: tensor(auxtens,e%ndime)
      real(rp)                   :: Term2(e%ndime,e%ndime,e%ndime),Term1(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: Term3(e%ndime,e%ndime,e%ndime),Term(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: ExpminusId(e%ndime,e%ndime)   
      real(rp)                   :: Prod1(e%ndime,e%pnode), Prod2(e%ndime,e%pnode)
      real(rp)                   :: Prod3(e%ndime,e%pnode)
      real(rp)                   :: Prod4(e%ndime,e%pnode)
      

      aux=dvolu*tisig
      
      lamb_div=lambda/(2.0_rp*lambda0)
      auxvisc=1.0_rp/(2.0_rp*acvis*(1.0_rp-beta))
      visc_div=lambda*auxvisc
      aux2=  +1.0_rp*auxoss -2.0_rp*lamb_div*Newton
      
      ExpminusId=0.0_rp
! !     GIESEKUSLCR      
!       if (auxG/=0.0_rp) then
!          ExpminusId=ExpPosMatrix
!          do i=1,e%ndime
!             ExpminusId(i,i)=ExpminusId(i,i)-1.0_rp
!          end do
!          ExpminusId=auxG*(lambda0/(acvis*(1.0_rp-beta)))*ExpminusId
!       end if
      
      Prod1=0.0_rp
      Prod2=0.0_rp
      Prod3=0.0_rp
      Prod4=0.0_rp
      Prod1=matmul(ExpminusId,e%cartd)
      Prod2=matmul(grvel,e%cartd)
      Prod3=matmul(ExpPosMatrix,e%cartd)
      Prod4=matmul(grvel,Prod3)

      do jnode=1,e%pnode    
         do inode=1,e%pnode 
            !Initializations
            Term=0.0_rp
            tensor=0.0_rp
            Shaj=e%shape(jnode,e%igaus)
            Shai=e%shape(inode,e%igaus)
            do concurrent (k=1:e%ndime)
               
               Term(:,k,k) = &
               !Cross (-T , -grad_sym(u))
                             -0.5_rp*auxvisc*Shai*(-e%cartd(:,jnode))*auxoss*aux2&
               !Cross (a*grad(T) , -grad_sym(u)) (1)
                             -0.5*visc_div*e%cartd(:,jnode)*gpadv(inode)*auxoss*aux2&
               !Cross terms: ( T*grad^t(u) + grad(u)*T , -grad_sym(u) )
                             -0.5_rp*(visc_div*Prod2(:,jnode)+Prod1(:,jnode))*Shai*aux2& 
               !(-T + a*grad(T) + T*grad^t(u) + grad(u)*T , - exp(S)*grad^t(U) - grad(U)*exp(S))               
                             -lamb_div*Prod3(:,jnode)*(-auxvisc*Shai+ visc_div*gpadv(inode))*auxoss*Newton& 
                             -lamb_div*Prod4(:,jnode)*visc_div*Shai*Newton +  Term(:,k,k) 
               
               Term(k,:,k) = &
               !Cross (-T , -grad_sym(u)) (2)
                              -0.5_rp*auxvisc*Shai*(-e%cartd(:,jnode))*auxoss*aux2&
               !Cross (a*grad(T) , -grad_sym(u)) (2)
                             -0.5*visc_div*e%cartd(:,jnode)*gpadv(inode)*auxoss*aux2&
               !Cross terms: ( T*grad^t(u) + grad(u)*T , -grad_sym(u))
                             -0.5_rp*Prod2(:,jnode)*visc_div*Shai*aux2&  
               !(-T + a*grad(T) + T*grad^t(u) + grad(u)*T , - exp(S)*grad^t(U) - grad(U)*exp(S))                
                             -lamb_div*Prod3(:,jnode)*(-auxvisc*Shai+ visc_div*gpadv(inode))*auxoss*Newton& 
                             -lamb_div*Prod4(:,jnode)*visc_div*Shai*Newton  + Term(k,:,k) 
               
               do concurrent (j=1:e%ndime)
                  Term(:,j,k) = -0.5_rp*e%cartd(:,jnode)*visc_div*Shai*grvel(j,k)*aux2 &
                                -0.5_rp*e%cartd(j,jnode)*(visc_div*grvel(:,k)+ ExpminusId(k,:))*Shai*aux2&
                                -lamb_div*Prod3(:,jnode)*visc_div*Shai*grvel(j,k)*Newton&
                                -lamb_div*Prod3(j,jnode)*visc_div*grvel(:,k)*Shai*Newton + Term(:,j,k)  
               
                  !(-T + a*grad(T) + T*grad^t(u) + grad(u)*T , U*grad(Exp(S)))  
                  Term(:,j,k)  = lamb_div*Shaj*GrExpPsi(:,j,k)*(-auxvisc*Shai*auxoss + visc_div*gpadv(inode))*Newton  + Term(:,j,k)
                  do i=1,e%ndime
                     Term(:,j,k) = lamb_div*Shaj*GrExpPsi(:,i,k)*visc_div*Shai*grvel(j,i)*auxoss*Newton + Term(:,j,k)
                     Term(:,j,k) = lamb_div*Shaj*GrExpPsi(i,j,k)*visc_div*grvel(:,i)*Shai*auxoss*Newton + Term(:,j,k)
                  end do
               
               end do
            end do   

            call PassTensor3ToTensor2left(e,auxtens,Term,tensor)
            elmat(1:auxtens,inode,:,jnode) =   aux*tensor(1:auxtens,:) + elmat(1:auxtens,inode,:,jnode)

         end do  
      end do
      
      
   end subroutine supm_LCR_elmbutEst2 
   
   
   subroutine supm_LCR_elmrhcEst2(e,lambda0,lambda,auxG,beta,acvis,tisig,gpsig,grsig,grvel,gpadv,ExpPosMatrix,ExpGiesekus, dvolu,auxtens,elextS,RHSlin,RHSGiesekus,elrhs)         
    !-----------------------------------------------------------------------
    !
    ! auxLCR=1 tau3*(-(1/lambda)*I + (1/lambda)*exp(-s)*(I+s), -(1/lambda)*exp(-s)*T)
    !                             
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: tisig,beta,lambda, acvis,dvolu, auxG,lambda0
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime), gpsig(auxtens)
      real(rp),    intent(in)    :: grsig(auxtens,e%ndime),ExpPosMatrix(e%ndime,e%ndime),ExpGiesekus(e%ndime,e%ndime)
      real(rp),    intent(in)    :: elextS(auxtens), RHSlin(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: i,j,k,inode
      real(rp)                   :: Shapi, aux,aux_visc,visc_div
      real(rp)                   :: Term1(e%ndime,e%ndime),Term2(e%ndime,e%ndime)
      real(rp)                   :: Term(e%ndime,e%ndime),vector(auxtens)
      real(rp)                   :: elextS_Mat(e%ndime,e%ndime),RHSGiesekus(e%ndime,e%ndime)
      real(rp)                   :: RHS(e%ndime,e%ndime)
      real(rp)                   :: ExpminusId(e%ndime,e%ndime)
      real(rp)                   :: Prod1(e%ndime,e%ndime), Prod1Giesekus(e%ndime,e%ndime)
      
      aux=tisig*dvolu
      RHS=0.0_rp

      aux_visc=1.0_rp/(2.0_rp*acvis*(1.0_rp-beta))
      visc_div=lambda*aux_visc
      
      ExpminusId=0.0_rp
      
      !GIESEKUSLCR         
     
      call sup_SigmaMatrix(e%ndime,auxtens,elextS,elextS_Mat)
      RHS=elextS_Mat+RHSlin+RHSGiesekus
      Prod1=matmul(grvel,RHS)
      
      Prod1Giesekus=0.0_rp
!        if (auxG/=0.0_rp) then
!          ExpminusId=ExpPosMatrix
!          do i=1,e%ndime
!             ExpminusId(i,i)=ExpminusId(i,i)-1.0_rp
!          end do
!          ExpminusId=auxG*(lambda0/(acvis*(1.0_rp-beta)))*ExpminusId
!          Prod1Giesekus=matmul(ExpminusId,RHS)
!       end if

      
      do inode=1,e%pnode 
         Shapi=e%shape(inode,e%igaus)
         Term=0.0_rp
         do concurrent (j=1:e%ndime)
            Term(:,j)= RHS(:,j)*(-aux_visc*Shapi+visc_div*gpadv(inode))&
                       +(visc_div*Prod1(:,j) + Prod1Giesekus(:,j))*Shapi + Term(:,j) 
            do k=1,e%ndime
               Term(:,j)= RHS(:,k)*visc_div*Shapi*grvel(j,k) + Term(:,j) 
            end do 
         end do
         
         call GetStress(e%ndime,Term,vector)
         elrhs(:,inode)= aux*vector(:) + elrhs(:,inode)
     end do

   end subroutine supm_LCR_elmrhcEst2   
   

!    subroutine supm_LCR_elmrhcEst2_split(e,lambda0,lambda,beta,acvis,tisig,grvel,gpadv, dvolu,auxtens,RHSlin_conv,RHSlin_deform,elrhs)         
!     !-----------------------------------------------------------------------
!     !
!     ! auxLCR=1 tau3*(-(1/lambda)*I + (1/lambda)*exp(-s)*(I+s), -(1/lambda)*exp(-s)*T)
!     !                             
!     !-----------------------------------------------------------------------
!       implicit none
!       class(FiniteElement) :: e
! 
!       integer(ip), intent(in)    :: auxtens
!       real(rp),    intent(in)    :: tisig,beta,lambda, acvis,dvolu,lambda0
!       real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime)
!       real(rp),    intent(in)    :: RHSlin_conv(e%ndime,e%ndime), RHSlin_deform(e%ndime,e%ndime)
!       real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
!       integer(ip)                :: i,j,k,inode
!       real(rp)                   :: aux, aux1,aux_visc,visc_div
!       real(rp)                   :: Term1(e%ndime,e%ndime),Term2(e%ndime,e%ndime)
!       real(rp)                   :: Term(e%ndime,e%ndime),vector(auxtens)
!       real(rp)                   :: elextS_Mat(e%ndime,e%ndime),RHSGiesekus(e%ndime,e%ndime)
!       real(rp)                   :: ExpminusId(e%ndime,e%ndime)
!       
!       aux1=tisig*dvolu
! 
!       aux_visc=1.0_rp/(2.0_rp*acvis*(1.0_rp-beta))
!       visc_div=lambda*aux_visc
!       
!       do inode=1,e%pnode 
!          aux=e%shape(inode,e%igaus)
!          Term=0.0_rp
!          do i=1,e%ndime
!             do j=1,e%ndime
!             
!                Term(i,j)= RHSlin_conv(i,j)*(visc_div*gpadv(inode)) + Term(i,j)
!                do k=1,e%ndime
!                   Term(i,j)= -RHSlin_deform(i,k)*visc_div*e%shape(inode,e%igaus)*grvel(j,k) + Term(i,j) 
!                   Term(i,j)= -RHSlin_deform(k,j)*(visc_div*grvel(i,k) + ExpminusId(k,i))*e%shape(inode,e%igaus) + Term(i,j) 
!                end do 
!                 
!             end do
!          end do
!          
!          call GetStress(e%ndime,Term,vector)
!          elrhs(:,inode)= aux1*vector(:) + elrhs(:,inode)
!      end do
! 
!    end subroutine supm_LCR_elmrhcEst2_split   
   
      
   subroutine supm_LCR_elmbsvEst2(e,auxtens,lambda,lambda0,auxG,dtinv,tisig,gpvel,grvel,gpadv,dvolu,ExpPosMatrix,ExpGiesekus,GrExpPsi,ConvExpMatrix,auxoss,elmat)
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxoss
      real(rp),    intent(in)    :: lambda0,lambda,dvolu,tisig,dtinv, auxG 
      real(rp),    intent(in)    :: gpadv(e%pnode),ExpPosMatrix(e%ndime,e%ndime),grvel(e%ndime,e%ndime), ExpGiesekus(e%ndime,e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime),GrExpPsi(e%ndime,e%ndime,e%ndime),ConvExpMatrix(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i,j,k,l
      real(rp)                   :: aux, aux_lamb, Shapj, lamb_div 
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime)
      real(rp)                   :: tensor(e%ndime,auxtens), MatrixSumExp(e%ndime,e%ndime)
      real(rp)                   :: Prod1(e%ndime,e%pnode),conv(e%pnode), Prod2(e%ndime,e%pnode)
      real(rp)                   :: Prod3(e%ndime,e%pnode), Prod4(e%ndime,e%ndime), Prod5(e%ndime,e%pnode)
      real(rp)                   :: Prod6(e%ndime,e%pnode)

      !Preliminaries
      aux = tisig*dvolu
      aux_lamb=1.0_rp/(2.0_rp*lambda0)
      lamb_div=lambda*aux_lamb
      
      MatrixSumExp=aux_lamb*ExpPosMatrix+ lamb_div*dtinv*ExpPosMatrix+ auxG*ExpGiesekus
      
      Prod1=0.0_rp
      Prod2=0.0_rp
      Prod3=0.0_rp
      Prod4=0.0_rp
      Prod5=0.0_rp
      Prod6=0.0_rp

      Prod1=matmul(MatrixSumExp,e%cartd)
      Prod2=matmul(ExpPosMatrix,e%cartd)
      Prod3=matmul(transpose(grvel),e%cartd)
      Prod4=matmul(grvel,ExpPosMatrix)
      Prod5=matmul(transpose(Prod4),e%cartd)
      Prod6=matmul(transpose(ConvExpMatrix),e%cartd)
      
      conv=0.0_rp
      do i=1,e%ndime
         conv(:)= gpvel(i)*e%cartd(i,:) + conv(:)
      end do
   
      !Loop
      do jnode=1,e%pnode 
         do inode=1,e%pnode 
            tensor=0.0_rp
            Term=0.0_rp
            
            Shapj=e%shape(jnode,e%igaus)

            do concurrent (k=1:e%ndime)
               !( -grad_sym(V), Exp(S))
               Term(k,:,k)=-0.5_rp*Prod1(:,inode)*Shapj*auxoss&
               !(-grad_sym(V), -exp(S)*grad(u) -grad^t(u)*exp(S)) (1)
                           +0.5_rp*lamb_div*Prod5(:,inode)*Shapj*auxoss&
               !(-grad_sym(V), u*grad(exp(S)) (1)
                           +0.5_rp*lamb_div*conv(jnode)*(-Prod2(:,inode))*auxoss&
                           +0.5_rp*lamb_div*(-Prod6(:,inode))*Shapj*auxoss + Term(k,:,k)
            
               do concurrent (j=1:e%ndime)
                  !( -grad_sym(V), Exp(S))
                  Term(:,k,j)=-0.5_rp*MatrixSumExp(:,k)*Shapj*e%cartd(j,inode)*auxoss&
                  !Cross terms
                  !(-grad_sym(V), u*grad(exp(S))
                              +0.5_rp*lamb_div*ConvExpMatrix(:,k)*Shapj*(-e%cartd(j,inode))*auxoss&
                              +0.5_rp*lamb_div*conv(jnode)*ExpPosMatrix(:,k)*(-e%cartd(j,inode))*auxoss&
                  !(-grad_sym(V), -exp(S)*grad(u) -grad^t(u)*exp(S))
                              +0.5_rp*lamb_div*Prod2(k,inode)*Shapj*grvel(:,j)*auxoss&
                              +0.5_rp*lamb_div*ExpPosMatrix(:,k)*Shapj*Prod3(j,inode)*auxoss&
                  !(-grad_sym(V), -exp(S)*grad(u) -grad^t(u)*exp(S)) (2)
                              +0.5_rp*lamb_div*Prod4(:,k)*Shapj*e%cartd(j,inode)*auxoss + Term(:,k,j)
               end do
            end do
         
         
            call PassTensor3ToTensor2Right(e,auxtens,Term,tensor)
         
            elmat(1:e%ndime,inode,:,jnode)=aux*tensor(1:e%ndime,:)+ elmat(1:e%ndime,inode,:,jnode)
                                    
         end do
      end do   

   end subroutine supm_LCR_elmbsvEst2
   

  subroutine supm_LCR_elmbuvEst2(e,lambda,lambda0,tisig,gpsig,dvolu,auxtens,ExpPosMatrix,GrExpPsi,auxoss,Newton,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! auxLCR=1   tau3*( -2exp(-S)*grad_sym(u) , -auxL*exp(s)*grad(V))
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,auxoss,Newton
      real(rp),    intent(in)    :: dvolu, tisig, lambda, lambda0
      real(rp),    intent(in)    :: gpsig(auxtens),ExpPosMatrix(e%ndime,e%ndime),GrExpPsi(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,i,j,k
      real(rp)                   :: aux2,aux,aux_lamb
      real(rp)                   :: Term(e%ndime,e%ndime),cartd2(e%pnode,e%pnode)
      real(rp)                   :: Prod1(e%ndime,e%pnode)
      real(rp)                   :: Prod2(e%ndime,e%ndime,e%pnode), Prod3(e%pnode,e%pnode)

      !Preliminaries
      aux_lamb=lambda/(2.0_rp*lambda0)
      aux2=tisig*dvolu
      aux=1.0_rp -2.0_rp*aux_lamb*Newton*auxoss
      
      cartd2=0.0_rp
      Prod1=0.0_rp
      Prod2=0.0_rp
      Prod3=0.0_rp
      Prod1=matmul(ExpPosMatrix,e%cartd) !Reminder: ExpS is a symmetric matrix
      do concurrent (i=1:e%pnode, k=1:e%ndime)
         cartd2(:,i) = e%cartd(k,:)*e%cartd(k,i) + cartd2(:,i)
         do j=1,e%ndime
            Prod2(:,k,i)=GrExpPsi(j,:,k)*e%cartd(j,i) + Prod2(:,k,i)
         end do
         Prod3(:,i)= Prod1(k,:)*e%cartd(k,i) + Prod3(:,i)
      end do
      
      !Loop
      do jnode=1,e%pnode   
         do inode=1,e%pnode
            Term=0.0_rp
            do concurrent (i=1:e%ndime)
               
               Term(i,i)=& 
               !(-grad_sym(V), -grad_sym(U))
                        +0.5_rp*cartd2(jnode,inode)*aux&
               !Cross (-grad_sym(V), U*grad_sym(s) -exp(s)*grad(U) -grad^t(U)*exp(s))
                        +aux_lamb*Prod3(jnode,inode)*auxoss*Newton + Term(i,i)
               
               Term(:,i) =&
               !(-grad_sym(V), -grad_sym(U))      
                          +0.5_rp*0.5_rp*e%cartd(:,jnode)*(e%cartd(i,inode))*aux&
                          +0.5_rp*0.5_rp*e%cartd(:,jnode)*(e%cartd(i,inode))*aux&
               !Cross (-grad_sym(V), U*grad_sym(s) -exp(s)*grad(U) -grad^t(U)*exp(s) )
                          +aux_lamb*Prod1(:,jnode)*e%cartd(i,inode)*auxoss*Newton&
                          -aux_lamb*e%shape(jnode,e%igaus)*Prod2(:,i,inode)*auxoss*Newton + Term(:,i)  
            end do
            elmat(1:e%ndime,inode,:,jnode)=aux2*Term(1:e%ndime,:) + elmat(1:e%ndime,inode,:,jnode)
         end do  
      end do

   end subroutine supm_LCR_elmbuvEst2 
   
   subroutine supm_LCR_elmrhuEst2(e,auxG,tisig,gpsig,ExpGiesekus,ExpPosMatrix,dvolu,auxtens,elextS,RHSlin,RHSGiesekus,elrhs) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    ! tau3*(-(1/lambda)*I + exp(-s)*(I+s) + u*grad(s), )
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: tisig,dvolu, auxG
      real(rp),    intent(in)    :: gpsig(auxtens), elextS(auxtens), ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpGiesekus(e%ndime,e%ndime), RHSlin(e%ndime,e%ndime),RHSGiesekus(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,i,j,Term_Exp,lin,lin2,lin3
      real(rp)                   :: aux,aux_lamb
      real(rp)                   :: RHS(e%ndime,e%ndime), Term(e%ndime),elextSMat(e%ndime,e%ndime)

      aux=tisig*dvolu  
      
      RHS=0.0_rp
     
      call sup_SigmaMatrix(e%ndime,auxtens,elextS,elextSMat)
      RHS=elextSMat+RHSlin
      
      do inode=1,e%pnode 
         
         Term=0.0_rp
         do concurrent (j=1:e%ndime)
            Term(:)= -0.5_rp*(RHS(j,:)+RHS(:,j))*e%cartd(j,inode)+ Term(:)
         end do
      
         elrhs(:,inode)= aux*Term(:) + elrhs(:,inode)
      end do

   end subroutine supm_LCR_elmrhuEst2 
   
   
 !-----------------------------------------------------------------------------------------------------------------------------------------
   
 
  subroutine ComputeMomentumRHS(e,auxtens,auxL,acden,gpsig,GrPsiMatrix,DivExpPsi,ExpGpPsi_Matrix,RHS)
      ! ----------------------------------------------
      !  auxL* Div(exp(s)) linearization
      !-----------------------------------------------
      
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: acden,auxL
      real(rp),    intent(in)    :: gpsig(auxtens), GrPsiMatrix(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(in)    :: DivExpPsi(e%ndime), ExpGpPsi_Matrix(e%ndime,e%ndime) 
      real(rp),    intent(out)   :: RHS(e%ndime)
      integer(ip)                :: i,j,k
      real(rp)                   :: SigmaMatrix(e%ndime,e%ndime)
      real(rp)                   :: Term(e%ndime)
      
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigmaMatrix)
      
      !Initializations
      RHS=0.0_rp
      
      Term=0.0_rp
      do concurrent (k=1:e%ndime)
         Term(:) = DivExpPsi(k)*SigmaMatrix(k,:) + Term(:)
         do concurrent (i=1:e%ndime)
            Term(:)= ExpGpPsi_Matrix(i,k)*GrPsiMatrix(k,:,i) + Term(:)
         end do
      end do
      
      RHS(:)=RHS(:) + auxL*(DivExpPsi(:)-Term(:))
  end subroutine
  
  
  subroutine ComputeMomentumGalerkinRHS(e,auxtens,auxL,acden,elext,grvel,gpvel,gpsig,ExpGpPsi_Matrix,Newton,Term1,Term2)
      ! ----------------------------------------------
      ! Compute the rhu term
      !  f + rho*a* grad(a) + auxL*(Div(exp(s)(I-s)))
      !-----------------------------------------------
      implicit none
         class(FiniteElement) :: e
         integer(ip), intent(in)    :: auxtens,Newton
         real(rp),    intent(in)    :: acden,auxL,elext(e%ndime),grvel(e%ndime,e%ndime),gpvel(e%ndime)
         real(rp),    intent(in)    :: gpsig(auxtens)
         real(rp),    intent(in)    :: ExpGpPsi_Matrix(e%ndime,e%ndime) 
         real(rp),    intent(out)   :: Term1(e%ndime),Term2(e%ndime,e%ndime)
         integer(ip)                :: i,j,k
         real(rp)                   :: convec(e%ndime),SigmaMatrix(e%ndime,e%ndime)
         real(rp)                   :: Prod(e%ndime,e%ndime)
         
         call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigmaMatrix)
         
         Term1=0.0_rp
         Term2=0.0_rp 
         
         convec=0.0_rp
         do i=1,e%ndime
            convec(:)=acden*gpvel(i)*grvel(:,i)*Newton+convec(:)
         end do
         
         Term1=elext+convec
         Prod=matmul(ExpGpPsi_Matrix,SigmaMatrix)
         Term2=auxL*(-ExpGpPsi_Matrix+Prod)
         do i=1,e%ndime
            Term2(i,i)=Term2(i,i)+auxL
         end do 
         
  end subroutine
  
  subroutine ComputeConstitutiveRHS(e,lambda,lambda0,auxtens,dtinv,gpsig,GrPsiMatrix,gpvel,grvel,ExpPosMatrix,GrExpPsi,Newton,Term) 
   
   ! ---------------------------------------------------------
   ! This subroutine computes the RHS of constitutive terms which appear because of the differents linearizations. 
   ! ----------------------------------------------------------
  
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,Newton
      real(rp),    intent(in)    :: lambda0,lambda,dtinv
      real(rp),    intent(in)    :: gpsig(auxtens),grvel(e%ndime,e%ndime),gpvel(e%ndime) 
      real(rp),    intent(in)    :: GrPsiMatrix(e%ndime,e%ndime,e%ndime),GrExpPsi(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(out)   :: Term(e%ndime,e%ndime)
      integer(ip)                :: i,j,k,l, s
      real(rp)                   :: aux_lamb,lambdiv
      real(rp)                   :: SigMatrix(e%ndime,e%ndime), ProdExpSig(e%ndime,e%ndime)
      real(rp)                   :: ConvExpPsiMatrix(e%ndime,e%ndime), ConvPsiMatrix(e%ndime,e%ndime)
      real(rp)                   :: Prod2(e%ndime,e%ndime), Prod3(e%ndime,e%ndime)
      
      
      aux_lamb=1.0_rp/(2.0_rp*lambda0)
      lambdiv=lambda/(2.0_rp*lambda0)
      
      !Initializations
      Term=0.0_rp
     
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigMatrix)
      ProdExpSig=matmul(ExpPosMatrix,SigMatrix)
      
      ConvExpPsiMatrix=0.0_rp
      ConvPsiMatrix=0.0_rp
      do concurrent (j=1:e%ndime, i=1:e%ndime)
         ConvExpPsiMatrix(:,j)=gpvel(i)*GrExpPsi(:,j,i) + ConvExpPsiMatrix(:,j)
         ConvPsiMatrix(:,j)=gpvel(i)*GrPsiMatrix(:,j,i) + ConvPsiMatrix(:,j)
      end do
      
      Prod2=matmul(matmul(ExpPosMatrix,SigMatrix),transpose(grvel))
      Prod3=matmul(matmul(grvel,ExpPosMatrix),SigMatrix)

      !Compute Term1: linearization for exp(S)
      Term= (aux_lamb+ lambdiv*dtinv)*(-ExpPosMatrix+ProdExpSig)
      
      !Compute Term2: linealization for u*grad(exp(S))
      Term=-lambdiv*ConvExpPsiMatrix*(1.0_rp-Newton) + Term
      
      !Compute Term3: linealization for - exp(s) * grad(u) - grad_t(u) * exp(s)
      Term=-lambdiv*(Prod2+Prod3)& 
                 -lambdiv*(grvel+transpose(grvel))*(1.0_rp-Newton) +Term
      
      do concurrent (i=1:e%ndime)
         Term(i,i)= aux_lamb + Term(i,i) !+ lambdiv*dtinv
         
         !Compute Term2: linealization for u*grad(exp(S))
         do concurrent (k=1:e%ndime)
            Term(i,:)=lambdiv*ConvExpPsiMatrix(i,k)*SigMatrix(k,:)&
                       +lambdiv*ExpPosMatrix(i,k)*ConvPsiMatrix(k,:) + Term(i,:)

            !Compute Term3: linealization for - exp(s) * grad(u) - grad_t(u) * exp(s)
            Term(i,:)=lambdiv*ExpPosMatrix(i,k)*grvel(:,k)*(1.0_rp-Newton)&
                       +lambdiv*grvel(i,k)*(ExpPosMatrix(k,:))*(1.0_rp-Newton) +Term(i,:)
         end do  
         
      end do

  
  end subroutine
  
  
   subroutine ComputeConstitutiveRHSconvective(e,lambda,lambda0,auxtens,gpsig,GrPsiMatrix,gpvel,ExpPosMatrix,GrExpPsi,ConvExpMatrix,Newton,Term) 
   
   ! ---------------------------------------------------------
   ! This subroutine computes the RHS of constitutive terms which appear because of the differents linearizations. 
   ! ----------------------------------------------------------
  
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,Newton
      real(rp),    intent(in)    :: lambda0,lambda
      real(rp),    intent(in)    :: gpsig(auxtens),gpvel(e%ndime),ConvExpMatrix(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: GrPsiMatrix(e%ndime,e%ndime,e%ndime),GrExpPsi(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(out)   :: Term(e%ndime,e%ndime)
      integer(ip)                :: i,j,k,l
      real(rp)                   :: aux_lamb,lambdiv
      real(rp)                   :: SigMatrix(e%ndime,e%ndime), ConvPsiMatrix(e%ndime,e%ndime)
      real(rp)                   :: Prod1(e%ndime,e%ndime), Prod2(e%ndime,e%ndime)
      
      
      aux_lamb=1.0_rp/(2.0_rp*lambda0)
      lambdiv=lambda/(2.0_rp*lambda0)
      !Initializations
      Term=0.0_rp
     
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigMatrix)
      
      do j=1,e%ndime
         do i=1,e%ndime
            ConvPsiMatrix(:,j)=gpvel(i)*GrPsiMatrix(:,i,j) + ConvPsiMatrix(:,j)
         end do   
      end do
      
      Prod1=matmul(ConvExpMatrix,SigMatrix)
      Prod2=matmul(ExpPosMatrix,ConvPsiMatrix)
      
      !Compute Term2: linealization for u*grad(exp(S))
      Term=-lambdiv*ConvExpMatrix*(1.0_rp-Newton)&
           +lambdiv*Prod1+lambdiv*Prod2 +Term

  end subroutine
  
  
  subroutine ComputeConstitutiveRHSdeformation(e,lambda,lambda0,auxtens,gpsig,grvel,ExpPosMatrix,Newton,Term) 
   
   ! ---------------------------------------------------------
   ! This subroutine computes the RHS of constitutive terms which appear because of the differents linearizations. 
   ! ----------------------------------------------------------
  
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens,Newton
      real(rp),    intent(in)    :: lambda0,lambda
      real(rp),    intent(in)    :: gpsig(auxtens),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(out)   :: Term(e%ndime,e%ndime)
      integer(ip)                :: i,j,k,l
      real(rp)                   :: lambdiv
      real(rp)                   :: SigMatrix(e%ndime,e%ndime), Prod(e%ndime,e%ndime)
      real(rp)                   :: Prod2(e%ndime,e%ndime), Prod3(e%ndime,e%ndime)
      
      lambdiv=lambda/(2.0_rp*lambda0)
      
      !Initializations
      Term=0.0_rp
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigMatrix)
      Prod=matmul(grvel,ExpPosMatrix)
      Prod2=matmul(matmul(ExpPosMatrix,SigMatrix),transpose(grvel))
      Prod3=matmul(matmul(grvel,ExpPosMatrix),SigMatrix)
      ! linealization for - exp(s) * grad(u) - grad_t(u) * exp(s)
      Term=lambdiv*(Prod*(1.0_rp-Newton) +Prod*(1.0_rp-Newton) -Prod3 -Prod2&
                -transpose(grvel)*(1.0_rp-Newton) -grvel*(1.0_rp-Newton))+ Term

  end subroutine
  
  
  
  subroutine ComputeTerm_rhsGiesekus(e,auxtens,auxG,gpsig,ExpPosGiesekus,ExpPosMatrix,Term)
      implicit none
      class(FiniteElement) :: e
      integer(ip),   intent(in)    :: auxtens
      real(rp),      intent(in)    :: auxG
      real(rp),      intent(in)    :: ExpPosMatrix(e%ndime,e%ndime)
      real(rp),      intent(in)    :: ExpPosGiesekus(e%ndime,e%ndime), gpsig(auxtens) 
      real(rp),      intent(out)   :: Term(e%ndime,e%ndime)
      integer(ip)                  :: i
      real(rp)                     :: SigmaMatrix(e%ndime,e%ndime)
      
      !GIESEKUSLCR
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigmaMatrix)
      Term=0.0_rp
      
      Term=-ExpPosGiesekus+ExpPosMatrix
      Term=matmul(ExpPosGiesekus,SigmaMatrix)+Term
      
      do i=1,e%ndime
         Term(i,i)=Term(i,i)-1.0_rp
      end do   
      
      Term=auxG*Term   
         
   end subroutine   
   
   
  !--------------------------------------------------------------------------------------------------------------------
   !Split OSS
  
   subroutine supm_LCR_elmrhc_splitoss(e,auxtens,beta,auxL,timom,dvolu,gpdivexps,elrhs)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the rhs terms in constitutive equation for OSS
   !    
   ! -tau1(div(tau),resid_momentum)- tau3/2mu(tau,resid_constitutive)
   !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gpdivexps(e%ndime)
      real(rp),    intent(in)    :: dvolu,beta,timom, auxL
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k
      real(rp)                   :: aux, tensor(e%ndime,e%ndime), vect_sol(auxtens)

      aux = timom*dvolu    
      do inode=1,e%pnode
         do j=1,e%ndime
            tensor(:,j)=e%cartd(:,inode)*auxL*gpdivexps(j)
         end do
         call GetStress(e%ndime,tensor,vect_sol)

         elrhs(:,inode) = aux*vect_sol(:) + elrhs(:,inode)

      end do
   end subroutine supm_LCR_elmrhc_splitoss 
   
   
   subroutine supm_LCR_elmrhu_splitoss(e,acden,tidiv,dvolu,gprep,auxtens,auxoss,elrhs) 
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
      real(rp),    intent(in)    :: tidiv,dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode,idime, itens
      real(rp)                   :: aux
      real (rp)                  :: tensor(e%ndime,auxtens), vect_sol(e%ndime),gprep2(auxtens)

      aux = tidiv*dvolu
      do inode=1,e%pnode
         
         elrhs(1:e%ndime,inode) =  e%cartd(1:e%ndime,inode)*aux*gprep(auxtens+e%ndime+1) & 
              + elrhs(1:e%ndime,inode)

      end do

   end subroutine supm_LCR_elmrhu_splitoss
   
   
    subroutine supm_elmrfeLCR_oto(e,dtnsi,dtnsi2,acden,acvis,lambda,lambda0,auxG,auxL,veadv, &
         vegau,gpsig,grsig,grapr,grave,ExpPosMatrix,GrExpPsi,auxtens,elext,elextC,elextS,resim)        
         
      use typre
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in) :: auxtens
      real(rp), intent(in)    :: dtnsi,dtnsi2,acden,veadv(e%ndime),vegau(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime),ExpPosMatrix(e%ndime,e%ndime)
      real(rp), intent(in)    :: gpsig(auxtens),acvis,elextC(1),elextS(auxtens),lambda0,GrExpPsi(e%ndime,e%ndime,e%ndime),grsig(auxtens,e%ndime)
      real(rp), intent(in)    :: elext(e%ndime),lambda,auxL,auxG
      real(rp), intent(out)   :: resim(auxtens+e%ndime+1)
      real(rp)                :: aux_exp(e%ndime,e%ndime), aux_gpgrav(e%ndime,e%ndime), aux_gpconvsigma(e%ndime,e%ndime)
      real(rp)                :: aux_exp1(e%ndime,e%ndime)
      real(rp)                :: aux_gpdeform(e%ndime,e%ndime),aux_gpgiesek(e%ndime,e%ndime)
      real(rp)                :: auxlamb, divlamb
      real(rp)                :: aux, value(auxtens), gpdivs(e%ndime) , aux_gpdivs(e%ndime), aux_sum(e%ndime,e%ndime)

      integer(ip)             :: idime,jdime, itens, i, j,k
      
      resim = 0.0_rp
      aux=1.0_rp/(2.0_rp*acvis)
      auxlamb=1.0_rp/(2.0_rp*lambda0)
      divlamb=lambda/(2.0_rp*lambda0)
      
      !Residual in constitutive equation  
      aux_gpconvsigma=0.0_rp
      aux_gpdeform=0.0_rp
      aux_gpgiesek=0.0_rp
      
      !aux_gpgiesek=ExpPosMatrix2-2.0_rp*ExpPosMatrix
      aux_exp=ExpPosMatrix
      aux_exp1=ExpPosMatrix
      aux_gpdeform= grave+transpose(grave) + aux_gpdeform
      aux_gpgrav=-0.5_rp*(grave+transpose(grave))
      
      
      do i=1,e%ndime
         !aux_gpgiesek(i,i)=aux_gpgiesek(i,i)+1.0_rp
         aux_exp(i,i)=ExpPosMatrix(i,i) -1.0_rp
         do k=1,e%ndime
            !u*grad(exp)
            aux_gpconvsigma(i,:)=vegau(k)*GrExpPsi(i,:,k) + aux_gpconvsigma(i,:)
            
            !Exp(S)*gradu
            aux_gpdeform(i,:)=-ExpPosMatrix(i,k)*grave(:,k) -grave(i,k)*ExpPosMatrix(k,:) + aux_gpdeform(i,:)
         end do  
      end do
      
      aux_sum=auxlamb*aux_exp + aux_gpgrav + divlamb*(aux_gpconvsigma+aux_gpdeform+dtnsi2*aux_exp1)
      call PassSymMatrixToVector(e%ndime,auxtens,aux_sum,resim(1:auxtens))

      resim(1:auxtens)= -elextS(1:auxtens) + resim(1:auxtens)
                       
      !Temporal derivative LHS
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*dtnsi*vegau(1:e%ndime)
      
      !Contribution from pressure term
      !Pressure Gradient
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + grapr(1:e%ndime)

      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do
      
      !Contribution from the divergence of sigma
      gpdivs=0.0_rp
      do i=1,e%ndime
         gpdivs(:)=GrExpPsi(:,i,i) +  gpdivs(:)
      end do
      
      resim(auxtens+1:auxtens+e%ndime)= resim(auxtens+1:auxtens+e%ndime) - auxL*gpdivs(1:e%ndime)
      
      !RHS (Temporal derivative + External forces)
      resim(auxtens+1:auxtens+e%ndime) = resim(auxtens+1:auxtens+e%ndime) - elext(1:e%ndime)

      !Contribution from the divergence term
      do idime = 1,e%ndime                              
         resim(auxtens+e%ndime+1) = resim(auxtens+e%ndime+1) + grave(idime,idime)
      end do

   end subroutine supm_elmrfeLCR_oto   
   
   
   subroutine supm_LCR_elmrhc_oss(e,lambda0,lambda,beta,acvis,tisig,grvel,gpadv,dvolu,auxtens,gprep,elrhs)         
    !-----------------------------------------------------------------------
    !
    ! auxLCR=1 tau3*(-(1/lambda)*I + (1/lambda)*exp(-s)*(I+s), -(1/lambda)*exp(-s)*T)
    !                             
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: tisig,beta,lambda, acvis,dvolu,lambda0
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: gprep(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: i,j,k,inode
      real(rp)                   :: Shapi, aux,aux_visc,visc_div
      real(rp)                   :: Term(e%ndime,e%ndime),vector(auxtens)
      real(rp)                   :: gprep_Mat(e%ndime,e%ndime)
      real(rp)                   :: RHS(e%ndime,e%ndime)
      real(rp)                   :: ExpminusId(e%ndime,e%ndime)
      
      aux=tisig*dvolu
      RHS=0.0_rp

      aux_visc=1.0_rp/(2.0_rp*acvis*(1.0_rp-beta))
      visc_div=lambda*aux_visc
      
      ExpminusId=0.0_rp
      call sup_SigmaMatrix(e%ndime,auxtens,gprep,gprep_Mat)
      RHS=gprep_Mat

      do inode=1,e%pnode 
         Shapi=e%shape(inode,e%igaus)
         Term=0.0_rp
         Term= RHS*(visc_div*gpadv(inode)) + Term
         do concurrent (i=1:e%ndime, k=1:e%ndime)   
            Term(i,:)= RHS(i,k)*visc_div*Shapi*grvel(:,k) + Term(i,:) 
            Term(i,:)= RHS(k,:)*visc_div*grvel(i,k)*Shapi + Term(i,:) 
         end do
         
         call GetStress(e%ndime,Term,vector)
         elrhs(:,inode)= aux*vector(:) + elrhs(:,inode)
      end do

   end subroutine supm_LCR_elmrhc_oss   
 
   
!----------------------------------------------------------------------------------------
! Shock Capturing
!----------------------------------------------------------------------------------------
   
   subroutine supm_elmbstDC_LCR(e,auxtens,kdisc,auxL,dvolu,ExpoMatrix,GrExpMatrix,elmat)
   
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,kdisc,auxL
      real(rp),    intent(in)    :: ExpoMatrix(e%ndime,e%ndime),GrExpMatrix(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode, i, j, k ,s
      real(rp)                   :: aux, vect(e%ndime), value, tens_sol(e%ndime,e%ndime)
      real(rp)                   :: Term(e%ndime,e%ndime,e%ndime,e%ndime),Cartd2(e%pnode,e%pnode)
      real(rp)                   :: Prod(e%ndime,e%ndime,e%pnode)

      aux= kdisc*dvolu  
      Cartd2=matmul(transpose(e%cartd),e%cartd)
      
      do concurrent (i=1:e%ndime, k=1:e%ndime, s=1:e%ndime)
         Prod(i,k,:)=GrExpMatrix(i,k,s)*e%cartd(s,:) + Prod(i,k,:)
      end do

      do jnode=1,e%pnode  
         do inode=1,e%pnode
            Term=0.0_rp
            do concurrent (i=1:e%ndime, k=1:e%ndime)
               Term(:,i,k,i)=auxL*ExpoMatrix(:,k)*Cartd2(jnode,inode) + auxL*Prod(:,k,inode)*e%shape(jnode,e%igaus)&
                              + Term(:,i,k,i)
            end do
            call PassTensor4ToTensor2(e,auxtens,Term,tens_sol)
            elmat(1:auxtens,inode,:,jnode) = aux*(tens_sol(1:auxtens,:)) + elmat(1:auxtens,inode,:,jnode)
         end do
      end do   

   end subroutine supm_elmbstDC_LCR
   
   
    subroutine supm_elmrhcDC_LCR(e,kdisc,auxL,dvolu,gpsig,auxtens,ExpPosMatrix,GrPsiMatrix,GrExpPsi,elrhs)

      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu, kdisc, auxL 
      real(rp),    intent(in)    :: GrPsiMatrix(e%ndime,e%ndime,e%ndime),ExpPosMatrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: GrExpPsi(e%ndime,e%ndime,e%ndime), gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode, i, j, k, l
      real(rp)                   :: aux1
      real(rp)                   :: tens_sol(auxtens,e%ndime), Term(e%ndime,e%ndime,e%ndime), SigMatrix(e%ndime,e%ndime)
      
      
      Term=0.0_rp
      
      call sup_SigmaMatrix(e%ndime,auxtens,gpsig,SigMatrix)
      Term=auxL*GrExpPsi
      
      do concurrent (j=1:e%ndime, k=1:e%ndime, l=1:e%ndime)
         Term(:,j,l)=auxL*(GrExpPsi(:,k,l)*SigMatrix(k,j)) + Term(:,j,l)
         Term(:,j,l)=auxL*(ExpPosMatrix(:,k)*GrPsiMatrix(k,j,l)) + Term(:,j,l)
      end do
      
      call PassTensor3ToTensor2Left(e,auxtens,Term,tens_sol)
     
      do inode=1,e%pnode
         aux1=kdisc*dvolu
         do i=1,e%ndime
            elrhs(:,inode) = aux1*tens_sol(:,i)*e%cartd(i,inode) + elrhs(:,inode) 
         end do   
      end do
      
   end subroutine supm_elmrhcDC_LCR
   
 
end module
 

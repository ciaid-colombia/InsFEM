subroutine sup_endste(a,itask)
   !This routine ends a time step of the incoma%pressible NS equations.
   use typre
   use def_parame
   use Mod_ThreeField   
   use Mod_Mesh 
   use Mod_memor
   use Mod_SupExacso      
   implicit none
   class(ThreeFieldNSProblem) :: a
   type(FemMesh)    :: Mesh  
   type(MemoryMan)  :: Memor
   type(SupExacso)  :: exacso    
   
   integer(ip) :: itask   
   integer(ip) :: ifout,ndime,npoin,ipoin,idime,ntens,auxstep,icomp,j
   integer(ip) :: ielem,nelem
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:),exsig(:),exsigr(:,:)
   real(rp)                :: expre
   !Nodal coordinates
   real(rp), pointer       :: coord(:) 
   !For log-conformation model
   real(rp), allocatable   :: Expo_Vect(:), Vect(:), Expo_Matrix(:,:), Matrix(:,:), psiNew(:), Vect2(:)
   real(rp)   ::  lambda0,acvis,auxL,beta,lamb, auxL2
   
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   ntens=(ndime-1)*(ndime-1)+2
   
   auxstep=3

   !Pre-CrankNicolson
   if (itask == 1) then
      !Compute dissipation and other endste elemental computations using u^(n+theta)
      !Also compute the dynamic subgrid scales
      call a%EnditeElmope('Endste')
   
   !Post-CrankNicolson
   elseif (itask == 2) then
   
      !When log-conformation model is employed a%sigma is the new variable psi.
      !We have to calculate initial sigma  (now called sigmaold) to save this field too.
      
      call a%Memor%alloc(ndime,ndime,Matrix,'Matrix','output')  
      call a%Memor%alloc(ntens,Vect,'Vect','output') 
      call a%Memor%alloc(ntens,Vect2,'Vect2','output') 
      call a%Memor%alloc(ntens,Expo_Vect,'ExpoVect','output') 
      call a%Memor%alloc(ntens,PsiNew,'PsiNew','output')  
      call a%Memor%alloc(ndime,ndime,Expo_Matrix,'Expo_Matrix','output') 
      call a%Mesh%GetNpoin(npoin)
      
      call a%IncrementalLambda(1,a%MatProp(1)%LawViParam(3),lamb)
      beta=a%MatProp(1)%LawViParam(2)
      acvis=a%MatProp(1)%LawViParam(1) 
      
      if (a%LogFormulation==1) then
      
         if (a%MatProp(1)%LawViParam(3)==0 .or. a%MatProp(1)%LawViParam(5)==0) then
            lambda0=0.01_rp
         else   
            if(a%kfl_exacs/=0) then
               lambda0=a%MatProp(1)%LawViParam(5) !SETLAMBDA0
            else
               if(a%kfl_timei==0) then
                  lambda0=lamb*a%MatProp(1)%LawViParam(5)
               else   
                  lambda0=a%MatProp(1)%LawViParam(5)*a%MatProp(1)%LawViParam(3)
               end if   
            end if   
         end if   
         auxL=((1.0_rp-beta)*acvis)/lambda0
         auxL2=((1.0_rp-beta)*acvis)/lamb
         
         do ipoin=1,npoin
            Vect(:)=a%sigma(:,ipoin,1)
            call sup_ComputeExponential(ndime,ntens,Vect,Expo_Vect,Expo_Matrix)
            
            call deter(Expo_Matrix,a%taudet(ipoin),ndime)!DETERMINANT OF CONFORMATION TENSOR
            
            do j=1,ndime
               Expo_Vect(j)=Expo_Vect(j)-1.0_rp
            end do 
            a%sigmaold(:,ipoin)=auxL*Expo_Vect(:)
            
            !If it is required, the variable psi "real" is calculated 
            !- in order to compare this field with other publications
            
            Vect2(:)=auxL*Expo_Vect(:)
            call sup_LogarithmicTransformation(ndime,ntens,auxL2,Vect2,psiNew)
            if(a%npp_stepi(28)>0) a%psiReal(:,ipoin)=psiNew(:)
            !Expo_Matrix = tau_matrix
         end do
      else if (a%LogFormulation==0) then     !IF THE STANDARD MODEL IS SET
         auxL=lamb/((1.0_rp-beta)*acvis)
         call a%Mesh%GetNpoin(npoin)
         do ipoin=1,npoin
            Vect(:)=a%sigma(:,ipoin,1)
            Vect(:)=(1.0_rp/auxL)*Vect(:)
            do j=1,ndime
               Vect(j)=Vect(j)+1.0_rp
            end do 
            !pasar vector a matrix.
            call sup_SigmaMatrix(ndime,ntens,Vect,Matrix)
            call deter(Matrix,a%taudet(ipoin),ndime)
         end do
      end if
      
      call a%Memor%dealloc(ndime,ndime,Matrix,'Matrix','output')  
      call a%Memor%dealloc(ntens,Vect,'Vect','output')  
      call a%Memor%dealloc(ntens,Vect2,'Vect2','output') 
      call a%Memor%dealloc(ntens,Expo_Vect,'ExpoVect','output') 
      call a%Memor%dealloc(ntens,PsiNew,'PsiNew','output')
      call a%Memor%dealloc(ndime,ndime,Expo_Matrix,'Expo_Matrix','output') 
   
      !Update a%velocity and a%pressure
      if (a%kfl_exacs==0) then      
         !Higher order components
         if (a%ncomp > 3) then
            do icomp = a%ncomp, 4,-1
               a%veloc(:,:,icomp) = a%veloc(:,:,icomp-1)
               a%press(:,icomp) = a%press(:,icomp-1)
               a%sigma(:,:,icomp) = a%sigma(:,:,icomp-1)
            enddo
         endif
         !Previous time step component
         a%veloc(:,:,3) = a%veloc(:,:,1)  ! Vn = Vn+1
         a%press(:,3) = a%press(:,1)  
         a%sigma(:,:,3) = a%sigma(:,:,1)

         !Boundary operation at the end
         call a%EndBouope('Endste')

         !Update the subgrid a%velocity DynamicSS
         if(a%kfl_tacsg/=0) then
            call a%Mesh%GetNelem(nelem)
            if (a%kfl_tacsg>0) then
               !Crank Nicolson schemes
               do ielem=1,nelem
                  a%vesgs(ielem)%a(:,2,:)=a%vesgs(ielem)%a(:,1,:)
                  if (a%MatProp(1)%lawvi<0) then
                     a%sisgs(ielem)%a(:,2,:)=a%sisgs(ielem)%a(:,1,:)
                  end if
                  if (a%kfl_repro>=2) then
                     a%vesgs2(ielem)%a(:,2,:)=a%vesgs2(ielem)%a(:,1,:)
                     a%vesgs3(ielem)%a(:,2,:)=a%vesgs3(ielem)%a(:,1,:)
                  end if
                  if (a%kfl_repro == 4 .and. a%MatProp(1)%lawvi<0) then
                     a%sisgs2(ielem)%a(:,2,:)=a%sisgs2(ielem)%a(:,1,:)
                     a%sisgs3(ielem)%a(:,2,:)=a%sisgs3(ielem)%a(:,1,:)
                  end if
               end do
            endif
         end if
         
      ! transient exact solution     
      elseif(a%kfl_exacs/=0.and.a%kfl_timei==1)then  
      
         if (a%ncomp > 3) then
            do icomp = a%ncomp, 4,-1
               a%veloc(:,:,icomp) = a%veloc(:,:,icomp-1)
               a%press(:,icomp)   = a%press(:,icomp-1)
               a%sigma(:,:,icomp) = a%sigma(:,:,icomp-1)
            enddo
         endif
      
     
         ! Allocate exact components
         call a%Memor%alloc(ndime,exvel,'exvel','sup_endste')  
         call a%Memor%alloc(ndime,ndime,exveg,'exveg','sup_endste')     
         call a%Memor%alloc(ntens,exsig,'exsig','sup_endste')  
         call a%Memor%alloc(ntens,ndime,exsigr,'exsigr','sup_endste')
         call a%Memor%alloc(ndime,exprg,'exprg','sup_endste')       
         
         do ipoin = 1,npoin
         
            a%veloc(:,:,3) = a%veloc(:,:,1)  ! Vn = Vn+1
            a%press(:,3) = a%press(:,1)
            a%sigma(:,:,3)=a%sigma(:,:,1)          
         
            call a%Mesh%GetPointCoord(ipoin,coord)

            ! Exact solution
            call exacso%sup_ComputeSolution(ndime,coord,a%ctime,a%LogFormulation,a)
            call exacso%sup_GetPressure(ndime,expre,exprg)         
            call exacso%sup_GetVelocity(ndime,exvel,exveg)
            if (a%LogFormulation/=0) then
               call exacso%sup_GetPsi(ndime,exsig,exsigr)
            else   
               call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
            end if   
            
            if(a%istep<=auxstep)then      
               do idime=1,ndime                
                  a%veloc(idime,ipoin,1) = exvel(idime)  ! Vn-2 = Vn-1               
               end do
               do idime=1,ntens
                  a%sigma(idime,ipoin,1) = exsig(idime)                   
               end do
            end if              
             
            !exact solution
            do idime=1,ndime
               a%exveloc(idime,ipoin,1) = exvel(idime)              
            end do
            do idime=1,ntens
               a%exsigma(idime,ipoin,1) = exsig(idime)               
            end do  
         end do
         ! Allocate exact components
         call a%Memor%dealloc(ndime,exvel,'exvel','sup_endste')  
         call a%Memor%dealloc(ndime,ndime,exveg,'exveg','sup_endste')     
         call a%Memor%dealloc(ntens,exsig,'exsig','sup_endste')  
         call a%Memor%dealloc(ntens,ndime,exsigr,'exsigr','sup_endste')
         call a%Memor%dealloc(ndime,exprg,'exprg','sup_endste')          
      end if      
         
   endif

end subroutine sup_endste

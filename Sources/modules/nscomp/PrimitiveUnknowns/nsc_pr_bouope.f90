subroutine nsc_pr_bouope(NSCompPrimitiveProblem)
   use typre
   use Mod_php_SetTimeIntegrator
   use Mod_TimeIntegrator
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_HangingNodes
   use Mod_nsc_pr_elmdir
   use Mod_NSCompressiblePrimitiveNBC
   implicit none
   class(NSCompressiblePrimitiveProblem), target :: NSCompPrimitiveProblem
   
   real(rp), allocatable :: wmatr(:,:,:,:)
   real(rp), allocatable :: wrhsi(:,:)

   integer(ip) :: iboun,nboun,ndime,idime
   integer(ip) :: igaub,inodb,inode,idofn
   integer(ip) :: kfl_HangingNodes

   real(rp), allocatable :: bopre(:,:),gpbp(:)
   real(rp), allocatable :: bovel(:,:,:),gpbv(:,:)
   real(rp), allocatable :: botem(:,:),gpbt(:)
   real(rp), allocatable :: tract(:)
   real(rp), allocatable :: gprhsvel(:)

   real(rp)    :: gprhspre(1),gprhstem(1)

   real(rp)    :: dsurf

   integer(ip) :: itime

   real(rp)    :: lomach,cgamma
   real(rp), allocatable :: lovel(:,:)

   !Working primitive NSCNBC Matrices
   real(rp), allocatable :: Rhmat(:)
   real(rp), allocatable :: weakcoeff(:)
   real(rp), allocatable :: norder(:)

   !Transform Matrix 
   real(rp), allocatable :: Transform(:,:)
   !Convection Matrix coefficients
   real(rp), allocatable :: At0(:,:), Ab(:,:), TAbT(:,:)
   !Reaction Matrix coefficients
   real(rp), allocatable :: Sb(:,:), SbJ(:), TSbT(:,:)



   a=>NSCompPrimitiveProblem
   

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSdtinv,nsteps)
   ReferenceDtinv = a%dtinv
   !If CrankNicolson the reference is 1/2 dt, used in dynamic subscales, taus etc
   if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = 2*a%dtinv
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','nsc_pr_bouope')
   
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,wmatr,'wmatr','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,e%mnode,wrhsi,'wrhsi','nsc_pr_bouope')

   call a%Memor%alloc(a%ndofn,Rhmat,'Rhmat','nsc_pr_bouope')

   call a%Memor%alloc(a%ndofn,a%ndofn,Ab,'Ab','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,a%ndofn,Sb,'Sb','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,a%ndofn,TAbT,'TAbT','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,a%ndofn,TSbT,'TSbT','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,SbJ,'SbJ','nsc_pr_bouope')
   call a%Memor%alloc(e%ndime+2,weakcoeff,'weakcoeff','nsc_pr_bouope')
   call a%Memor%alloc(e%mnode,norder,'norder','nsc_pr_bouope')
   call a%Memor%alloc(a%ndofn,a%ndofn,Transform,'Transform','nsc_pr_bouope')

   call a%Mesh%GetHanging(kfl_HangingNodes)

   !Variables!
   call a%Memor%alloc(a%ndofn,tract,'tract','nsc_pr_bouope')
   call a%Memor%alloc(e%pnodb,a%ncomp-1,bopre,'bopre','nsc_pr_bouope')
   call a%Memor%alloc(e%ndime,e%pnodb,a%ncomp-1,bovel,'bovel','nsc_pr_bouope')
   call a%Memor%alloc(e%pnodb,a%ncomp-1,botem,'botem','nsc_pr_bouope')
   call a%Memor%alloc(a%ncomp-1,gpbp,'gpbp','nsc_pr_bouope')
   call a%Memor%alloc(e%ndime,a%ncomp-1,gpbv,'gpbv','nsc_pr_bouope')
   call a%Memor%alloc(a%ncomp-1,gpbt,'gpbt','nsc_pr_bouope')
   call a%Memor%alloc(e%ndime,gprhsvel,'gprhsvel','nsc_pr_bouope')
   call a%Memor%alloc(e%ndime,2_ip,lovel,'lovel','nsc_pr_bouope')

   call a%Memor%alloc(a%ndofn,a%ndofn,At0,'At0','nsc_pr_bouope')

   ! Loop over boundaries
   boundaries: do iboun=1,nboun

      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      elmat=0.0_rp
      elrhs=0.0_rp

      if(    a%kfl_fixbo(iboun)==3.or.&   !Mom=Neu  Ene=Free
         & a%kfl_fixbo(iboun)==4.or.&     !Mom=Neu  Ene=Neu
         & a%kfl_fixbo(iboun)==5.or.&     !Mom=Neu  Ene=Rob
         & a%kfl_fixbo(iboun)==6.or.&     !Mom=Free  Ene=Neu
         & a%kfl_fixbo(iboun)==7.or.&     !Mom=Free  Ene=Rob
         & a%kfl_fixbo(iboun)==8.or.&     !Mom=Wall_law Ene=Free
         & a%kfl_fixbo(iboun)==9.or.&     !Mom= backlow penalty
         & a%kfl_fixbo(iboun)==11.or.&    !NSCBC weak imposition
         & a%kfl_fixbo(iboun)==13 &       !NSCBC weak imposition of non-reflecting pressure
         ) then  

         !Physical Parameters
         call a%GetPhysicalParameters(acvis,actco,accph,accvh)
         
         call e%elmdel
         call e%elmlen
         
         !Gather operations
         call e%gatherb(1_ip,bopre(:,1),a%press(:,1))
         call e%gatherb(1_ip,bopre(:,2),a%press(:,3))
         call e%gatherb(e%ndime,bovel(:,:,1),a%veloc(:,:,1))
         call e%gatherb(e%ndime,bovel(:,:,2),a%veloc(:,:,3))
         call e%gatherb(1_ip,botem(:,1),a%tempe(:,1))
         call e%gatherb(1_ip,botem(:,2),a%tempe(:,3))
         do itime = 3,nsteps ! Time bdf2 and others
            call e%gatherb(1_ip,bopre(:,itime),a%press(:,itime+1)) 
            call e%gatherb(e%ndime,bovel(:,:,itime),a%veloc(:,:,itime+1)) 
            call e%gatherb(1_ip,botem(:,itime),a%tempe(:,itime+1)) 
         enddo
         
         dsurf = 0.0_rp

         !Gauss-Point Loop
         gauss_points : do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Initialize
            wmatr=0.0_rp
            wrhsi=0.0_rp
            tract=0.0_rp
      
            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            
            call e%interpb(1_ip,bopre(:,1),gpbp(1))
            call e%interpb(1_ip,bopre(:,2),gpbp(2))
            call e%interpb(e%ndime,bovel(:,:,1),gpbv(:,1))
            call e%interpb(e%ndime,bovel(:,:,2),gpbv(:,2))
            call e%interpb(1_ip,botem(:,1),gpbt(1))
            call e%interpb(1_ip,botem(:,2),gpbt(2))
            do itime = 3,nsteps ! Time bdf2 and others
               call e%interpb(1_ip,bopre(:,itime),gpbp(itime))
               call e%interpb(e%ndime,bovel(:,:,itime),gpbv(:,itime))
               call e%interpb(1_ip,botem(:,itime),gpbt(itime))
            enddo
            
            !Derivatives at the boundary
            call e%elmderb         

            !Linearized variables with a block Picard's method at each time step
            cgamma = accph/accvh
            gpadp = gpbp(1) + a%relpre
            gpadt = gpbt(1) + a%reltem
            gpden = gpadp/((accph-accvh)*gpadt)
            aux = gpden*(accvh*gpadt+dot_product(gpbv(:,1),gpbv(:,1))/2.0_rp)+gpadp
            acalpha = 1.0_rp/gpadt
            acbeta = 1.0_rp/gpadp
            !Calculate speed of sound
            gpspd = sqrt(cgamma*(gpadp)/gpden)
            !Local Mach number
            lomach = sqrt(dot_product(gpbv(:,1),gpbv(:,1)))/gpspd
            !Local description of velocity
            lovel(:,1) = matmul(gpbv(:,1),e%baloc)
            lovel(:,2) = matmul(gpbv(:,2),e%baloc)
   
            !Matrices to zero
            Rhmat=0.0_rp
            Ab=0.0_rp
            Sb=0.0_rp
            TAbT=0.0_rp
            TSbT=0.0_rp
            SbJ=0.0_rp
            norder =0.0_rp
            Transform=0.0_rp
            weakcoeff=0.0_rp

            At0 = 0.0_rp !Atmm(i,d)

            !Transient Matrix
            !Mass equation
            At0(1,1) = At0(1,1) + gpden*acbeta !rho*beta
            At0(1,a%ndofn) = At0(1,a%ndofn)  - gpden*acalpha !-rho*alpha
            !Momentum equation
            do idime=1,e%ndime
               At0(idime+1,1)       = At0(idime+1,1)        + gpden*acbeta*gpbv(idime,1)!rho*beta*vel_i
               At0(idime+1,idime+1) = At0(idime+1,idime+1)  + gpden!rho*d_id
               At0(idime+1,a%ndofn) = At0(idime+1,a%ndofn)  - gpden*acalpha*gpbv(idime,1)!-rho*alpha*vel_i
            !Energy equation
               At0(a%ndofn,idime+1) = At0(a%ndofn,idime+1)  + gpden*gpbv(idime,1)!rho*vel_d 
            end do
            At0(a%ndofn,1) = At0(a%ndofn,1) + aux*acbeta !rho*(enthalpy+kineticEnerg)*beta
            At0(a%ndofn,1) = At0(a%ndofn,1) - acalpha*gpadt !-alpha*temperature
            At0(a%ndofn,a%ndofn) = At0(a%ndofn,a%ndofn) - aux*acalpha!- rho*(enthalpy+kineticEnerg)*alpha
            At0(a%ndofn,a%ndofn) = At0(a%ndofn,a%ndofn) + gpden*accph !rho*cp

            !Shape functions derivative in the normal direction (e%cartb(ndime,mnodb)) for convection
            norder = matmul(e%baloc(e%ndime,1:e%ndime),e%cartb)

            ! Weak imposition of NSCNBC
            if(a%kfl_fixbo(iboun)==11) then
           
               if(e%ndime == 3) then 
                    Transform(4,2:e%ndime+1) = Transform(4,2:e%ndime+1) + e%baloc(2,1:e%ndime)  
                    Ab(4,4) = lovel(e%ndime,1)  
               end if
              
               !Assembly reference transformation matrix
               Transform(1,1) = 1.0_rp
               Transform(2,2:e%ndime+1) = Transform(2,2:e%ndime+1) + e%baloc(e%ndime,1:e%ndime)  
               Transform(3,2:e%ndime+1) = Transform(3,2:e%ndime+1) + e%baloc(1,1:e%ndime)  
               Transform(e%ndime+2,e%ndime+2) = 1.0_rp

               !Assembly convective-like matrix
               Ab(1,1) = (lovel(e%ndime,1)+gpspd)/2.0_rp
               Ab(1,2) = gpden*gpspd*(lovel(e%ndime,1)+gpspd)/2.0_rp  
               Ab(2,1) = (lovel(e%ndime,1)+gpspd)/(gpden*gpspd*2.0_rp)  
               Ab(2,2) = (lovel(e%ndime,1)+gpspd)/2.0_rp  
               Ab(3,3) = lovel(e%ndime,1) 
               Ab(e%ndime+2,1) =gpadt*(cgamma-1.0_rp)*(lovel(e%ndime,1)+gpspd)/(gpden*gpspd*gpspd*2.0_rp)+gpadt*lovel(e%ndime,1)/(gpden*gpspd*gpspd)-lovel(e%ndime,1)/(gpden*(accph-accvh))
               Ab(e%ndime+2,2) = gpadt*(cgamma-1.0_rp)*(lovel(e%ndime,1)+gpspd)/(gpspd*2.0_rp)
               Ab(e%ndime+2,e%ndime+2) = lovel(e%ndime,1)


               !Transform to original reference system
               TAbT = matmul(At0,matmul(matmul(transpose(Transform),Ab),Transform))

               !Assembly reactive-like matrix
               Sb(1,1) = a%bvnat(iboun)%a(2)*gpspd*abs(1.0_rp-lomach*lomach)/(2.0_rp*a%bvnat(iboun)%a(3))
               Sb(2,1) = -a%bvnat(iboun)%a(2)*abs(1.0_rp-lomach*lomach)/(2.0_rp*gpden*a%bvnat(iboun)%a(3))   
               Sb(e%ndime+2,1) = gpadt*(cgamma-1.0_rp)*a%bvnat(iboun)%a(2)*abs(1.0_rp-lomach*lomach)/(2.0_rp*a%bvnat(iboun)%a(3)*gpden*gpspd) 

               !Transform to original reference system
               TSbT = matmul(matmul(transpose(Transform),Sb),Transform)

               !LHS Temporal contribution, in original reference system, included into TSbT matrix
               TSbT(1,1) = TSbT(1,1) + LHSdtinv
               TSbT(2,2) = TSbT(2,2) + LHSdtinv
               TSbT(3,3) = TSbT(3,3) + LHSdtinv
               TSbT(e%ndime+2,e%ndime+2) = TSbT(e%ndime+2,e%ndime+2) + LHSdtinv
               if(e%ndime == 3) then 
                  TSbT(4,4) = TSbT(4,4) + LHSdtinv
               end if
               TSbT = matmul(At0,TSbT)

               !RHS Temporal contribution
               call Integrator%GetRHS(1_ip,gpbp,gprhspre)
               call Integrator%GetRHS(e%ndime,gpbv(:,2),gprhsvel)
               call Integrator%GetRHS(1_ip,gpbt,gprhstem)

               Rhmat(1)           = Rhmat(1)           + a%dtinv*gprhspre(1)
               Rhmat(2:e%ndime+1) = Rhmat(2:e%ndime+1) + a%dtinv*gprhsvel(1:e%ndime)
               Rhmat(e%ndime+2)   = Rhmat(e%ndime+2)   + a%dtinv*gprhstem(1)

               !Reactive unperturbed contribution
               SbJ(1) = a%bvnat(iboun)%a(4)*a%bvnat(iboun)%a(2)*gpspd*abs(1.0_rp-lomach*lomach)/(2.0_rp*a%bvnat(iboun)%a(3))
               SbJ(2) = -a%bvnat(iboun)%a(4)*a%bvnat(iboun)%a(2)*abs(1.0_rp-lomach*lomach)/(2.0_rp*gpden*a%bvnat(iboun)%a(3))  
               SbJ(e%ndime+2) = gpadt*(cgamma-1.0_rp)*a%bvnat(iboun)%a(4)*a%bvnat(iboun)%a(2)*abs(1.0_rp-lomach*lomach)/(2.0_rp*gpden*a%bvnat(iboun)%a(3)*gpspd) 
               
               !Transform to original reference system and include in RHS
               Rhmat = Rhmat + matmul(transpose(Transform),SbJ)
               Rhmat = matmul(At0,Rhmat)

               !Assembly elemental matrices
               call nsc_pr_WeakElementalAssembly(e,a%bvnat(iboun)%a(1)*e%hleng(1),norder,TAbT,TSbT,Rhmat,wmatr,wrhsi)

            ! Weak imposition of NSCNBC pressure
            elseif(a%kfl_fixbo(iboun)==13) then
           
               !Assembly convective-like matrix
               TAbT(1,1) = (lovel(e%ndime,1)+gpspd)
               TAbT(1,2) = gpden*gpspd*(lovel(e%ndime,1)+gpspd)
               TAbT(1,3) = gpden*gpspd*(lovel(e%ndime,1)+gpspd)
               if(e%ndime == 3) then 
                  TAbT(1,4) = gpden*gpspd*(lovel(e%ndime,1)+gpspd)
               end if

               TAbT = matmul(At0,TAbT)

               !LHS Temporal contribution, in original reference system, included into TSbT matrix
               TSbT(1,1)= LHSdtinv

               TSbT = matmul(At0,TSbT)

               !RHS Forces contribution
               Rhmat(1)           = Rhmat(1)           + gpspd*gpden*(lovel(e%ndime,1)-lovel(e%ndime,2))*a%dtinv

               !RHS Temporal contribution
               call Integrator%GetRHS(1_ip,gpbp,gprhspre)
               Rhmat(1)           = Rhmat(1)           + a%dtinv*gprhspre(1)

               Rhmat = matmul(At0,Rhmat)

               !Assembly elemental matrices
               call nsc_pr_WeakElementalAssembly(e,a%bvnat(iboun)%a(1)*e%hleng(1),norder,TAbT,TSbT,Rhmat,wmatr,wrhsi)

            end if

            elmat = elmat + dsurf*wmatr            
           
            do inodb=1,e%pnodb
               inode = e%lboel(inodb)
               do idofn=1,a%ndofn
                     elrhs(idofn,inode)=elrhs(idofn,inode)+dsurf*(wrhsi(idofn,inode)+tract(idofn)*e%shapb(inodb,e%igaub))
               end do   
              
           end do
   
         end do gauss_points

      end if 

      if (kfl_HangingNodes == 1) call ModifyMatricesHanging
      
      !Boundary conditions
      call nsc_pr_elmdir(a,e,elmat,elrhs)
      
      ! Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)

   end do boundaries

   call a%Memor%dealloc(a%ndofn,tract,'tract','nsc_pr_bouope')
   call a%Memor%dealloc(e%pnodb,a%ncomp-1,bopre,'bopre','nsc_pr_bouope')
   call a%Memor%dealloc(e%ndime,e%pnodb,a%ncomp-1,bovel,'bovel','nsc_pr_bouope')
   call a%Memor%dealloc(e%pnodb,a%ncomp-1,botem,'botem','nsc_pr_bouope')
   call a%Memor%dealloc(a%ncomp-1,gpbp,'gpbp','nsc_pr_bouope')
   call a%Memor%dealloc(e%ndime,a%ncomp-1,gpbv,'gpbv','nsc_pr_bouope')
   call a%Memor%dealloc(a%ncomp-1,gpbt,'gpbt','nsc_pr_bouope')
   call a%Memor%dealloc(e%ndime,gprhsvel,'gprhsvel','nsc_pr_bouope')
   call a%Memor%dealloc(e%ndime,2_ip,lovel,'lovel','nsc_pr_bouope')
   
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','nsc_pr_bouope')

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,wmatr,'wmatr','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,e%mnode,wrhsi,'wrhsi','nsc_pr_bouope')

   call a%Memor%dealloc(a%ndofn,Rhmat,'Rhmat','nsc_pr_bouope')

   call a%Memor%dealloc(a%ndofn,a%ndofn,Ab,'Ab','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,a%ndofn,Sb,'Sb','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,a%ndofn,TAbT,'TAbT','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,a%ndofn,TSbT,'TSbT','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,SbJ,'SbJ','nsc_pr_bouope')   
   call a%Memor%dealloc(e%ndime+2,weakcoeff,'weakcoeff','nsc_pr_bouope')
   call a%Memor%dealloc(e%mnode,norder,'norder','nsc_pr_bouope')
   call a%Memor%dealloc(a%ndofn,a%ndofn,Transform,'Transform','nsc_pr_bouope')

   call a%Memor%dealloc(a%ndofn,a%ndofn,At0,'At0','nsc_pr_bouope')

   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsc_pr_bouope')
   
end subroutine nsc_pr_bouope

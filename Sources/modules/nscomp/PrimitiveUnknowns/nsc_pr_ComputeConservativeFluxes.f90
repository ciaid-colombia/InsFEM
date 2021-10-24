subroutine nsc_pr_ComputeConservativeFluxes(a,task)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!    The routine calculates fluxes over boundaries
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_Mesh
   use Mod_Element
   use Mod_Memor
   use Mod_iofile
   use Mod_int2str
   use Mod_NSCompressiblePrimitive
   use Mod_Postpr
   use MPI
   
   implicit none
   class(NSCompressiblePrimitiveProblem) :: a   
   character(6) :: task
   class(FiniteElement) , pointer:: e => NULL()
 

   real(rp), allocatable :: tract(:)
   real(rp), allocatable :: bopre(:,:),gpbp(:)
   real(rp), allocatable :: bovel(:,:,:),gpbv(:,:)
   real(rp), allocatable :: botem(:,:),gpbt(:)
   real(rp), allocatable :: elvel(:,:), eltem(:)
   real(rp), allocatable :: grvel(:,:), grtem(:)
   real(rp), allocatable :: qflx(:),mflx(:),mcflx(:)

   real(rp)    :: weightfactor
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: dflx,eflx,ecflx,edflx
   real(rp)    :: gpbd,gpbe,divv
   real(rp)    :: fluxes(5),fluxesR(5)
   integer(ip) :: idime,iboun,ndime,igaub,jdime,nboun,nelem,ibody
   integer(ip), pointer :: lbody => NULL()
   integer  icount,npoinLocal,npoin,inodb,ipoin,ierr,idofn

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)   
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)     
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_EndBouope')      
   
   !Allocation

   call a%Memor%alloc(ndime+1,tract,'tract','nsc_pr_EndBouope')   
   call a%Memor%alloc(e%pnodb,2,bopre,'bopre','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%pnodb,2,bovel,'bovel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%pnodb,2,botem,'botem','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%mnode,eltem,'eltem','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,grtem,'grtem','nsc_pr_EndBouope')
   call a%Memor%alloc(2,gpbp,'gpbp','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,2,gpbv,'gpbv','nsc_pr_EndBouope')
   call a%Memor%alloc(2,gpbt,'gpbt','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,qflx,'qflx','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,mflx,'mflx','nsc_pr_EndBouope')
   call a%Memor%alloc(e%ndime,mcflx,'mcflx','nsc_pr_EndBouope')
  

   !Inializations
   fluxes = 0.0_rp
   fluxesR = 0.0_rp
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)
      
      call a%Mesh%GetLbody(lbody,iboun)
       
      ibody=lbody
      
      !to delete repeated entities 
      icount = 0
      do inodb = 1,e%pnodb
         if (e%lnodb(inodb) <= npoinLocal) then
            icount = icount +1
         endif
      enddo
      weightfactor = real(icount)/real(e%pnodb)
      
      !Physical Parameters
      call a%GetPhysicalParameters(acvis,actco,accph,accvh)         
      
      call e%elmdel
      
      !Gathers in boundary elements 
      call e%gatherb(1_ip,bopre(:,1),a%press(:,1))
      call e%gatherb(e%ndime,bovel(:,:,1),a%veloc(:,:,1))
      call e%gatherb(1_ip,botem(:,1),a%tempe(:,1))
      call e%gather(ndime,elvel,a%veloc(:,:,1))      
      call e%gather(1_ip,eltem,a%tempe(:,1))       
      
      dsurf0 = 0.0_rp           

      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         !Traction at gauss point
         tract=0.0_rp
         divv=0.0_rp

         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(1_ip,bopre(:,1),gpbp(1))
         call e%interpb(e%ndime,bovel(:,:,1),gpbv(:,1))
         call e%interpb(1_ip,botem(:,1),gpbt(1))
         
         !Derivatives at the boundary
         call e%elmderb         
         call e%gradientb(ndime,elvel,grvel)
         call e%gradientb(1_ip,eltem,grtem)

         ! Ideal Gas State Law 
         gpbd = (gpbp(1)+a%relpre)/((accph - accvh) * (gpbt(1)+a%reltem))
         ! Total energy
         gpbe = gpbd*(accvh*(gpbt(1)+a%reltem)+dot_product(gpbv(:,1),gpbv(:,1))/ 2.0_rp)

         !Velocity divergence
         do idime =1,e%ndime
           divv = divv + grvel(idime,idime)
         end do

         !----------------------------------------------------------------
         !Convective (total) mass flux for conservation restrictions
         !  Fd = -n· rho*vel
         dflx = -gpbd*dot_product(gpbv(:,1),e%baloc(:,e%ndime))

         !Convective momentum flux for conservation restrictions
         !  Fmc = -n· rho*vel*vel
         mcflx(:) = -gpbd*dot_product(gpbv(:,1),e%baloc(:,e%ndime))*gpbv(:,1)
        
  
         !Traction value (also for momentum conservation restrictions)
         !  t=n·S  S=-p*I+tau
         !  tau=(2mu)*gradS(vel)-(2mu/3)*div(vel)*I

         do idime=1,ndime
            tract(idime)= tract(idime) - e%baloc(idime,ndime)*(gpbp(1)+a%relpre)&
             -(2/3)*acvis*divv*e%baloc(idime,ndime)
            do jdime=1,ndime
               tract(idime) = tract(idime) &
                  + acvis*e%baloc(jdime,ndime)*(grvel(idime,jdime)+grvel(jdime,idime))
            end do
         end do

         !Total momentum flux for conservation restrictions
         mflx(1:ndime) = mcflx(1:ndime) + tract(1:ndime)

         !Convective energy flux for conservation restrictions
         !  Fec = -n· e_tot*vel
         ecflx = -gpbe*dot_product(gpbv(:,1),e%baloc(:,e%ndime))

         !Dissipative work as energy flux for conservation restrictions
         !  Fed = n_j· S_ij*vel_i
         edflx = dot_product(gpbv(1:ndime,1),tract(1:ndime))

         !Heat flux (also for energy conservation restriction) 
         !   Q=q·n 
         !   q=-actco*grad(tempe)
      
         qflx = -actco*grtem

         tract(ndime+1) = dot_product(qflx,e%baloc(:,e%ndime))

         !Total energy flux for conservation restrictions
         eflx = ecflx + edflx - tract(ndime+1)
         
         !Total conservation fluxes
         fluxes(1) = fluxes(1) + dflx*dsurf*weightfactor 
         fluxes(2:ndime+1) = fluxes(2:ndime+1) + mflx(1:ndime)*dsurf*weightfactor 
         fluxes(ndime+2) = fluxes(ndime+2) + eflx*dsurf*weightfactor 

      end do 
      
   end do boundaries 

     
   !Collect Values
   call  MPI_REDUCE(fluxes, fluxesR, 5, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr ) 

   if (task .eq. 'BegSte' ) then 
      a%ConservativeFluxes(:,2) = fluxesR(:)
   elseif (task .eq. 'BegIte' ) then
      a%ConservativeFluxes(:,1) = fluxesR(:)
   endif

   !Dellocation

   call a%Memor%dealloc(e%pnodb,2,bopre,'bopre','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%pnodb,2,bovel,'bovel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%pnodb,2,botem,'botem','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%mnode,eltem,'eltem','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,grtem,'grtem','nsc_pr_EndBouope')
   call a%Memor%dealloc(2,gpbp,'gpbp','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,2,gpbv,'gpbv','nsc_pr_EndBouope')
   call a%Memor%dealloc(2,gpbt,'gpbt','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,qflx,'qflx','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,mflx,'mflx','nsc_pr_EndBouope')
   call a%Memor%dealloc(e%ndime,mcflx,'mcflx','nsc_pr_EndBouope')

   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_froces')
   
   101 format(4x,i9,2x,i9,20(2x,e12.6))
end subroutine


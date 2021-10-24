subroutine nsi_updbcs(a)
!-----------------------------------------------------------------------
!****f* Nstinc/nsi_updbcs
! NAME 
!    nsi_updbcs
! DESCRIPTION
!    This routine updates the velocity boundary conditions:
!    1. Before a time step begins
!    2. Before a global iteration begins
!    3. Before an inner iteration begins
! USED BY
!    nsi_begste
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_NavierStokes
   use Mod_NsiExacso     
   use Mod_BlasiusProfile
   use Mod_TurbulentInletBoundaryConditions
   implicit none
   class(NavierStokesProblem), target :: a
   class(FiniteElement), pointer :: e => NULL()
   type(NsiExacso)  :: exacso      
   
   integer(ip)             :: ibopo,idime,ndime,ipoin,npoin
   integer(ip)             :: iboun,nboun,ielem,nelem,inode,inodb
   logical                 :: isALE
   real(rp), pointer       :: exnor(:,:) => NULL()
   real(rp), pointer       :: meshve(:,:,:) => NULL()
   real(rp) :: coeff, Wavelength,Frequency
   integer(ip) :: CoordinateAxis

   integer(ip) :: funty(2), funno
   real(rp), pointer :: pointCoord(:) => NULL()
   
   type (BlasiusProfileGenerator) :: BPGenerator  
   real(rp) :: X,UFree, nu
   integer(ip) :: ifun
   
   real(rp), pointer :: r_i1(:,:) => NULL(),r_i21(:,:) => NULL()
   real(rp) :: w_r

   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   
   call a%Mesh%GetALE(isALE)
   call a%Mesh%GetMeshVeloc(meshve)

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_updbcs')
   
   if (isALE) then
      do ipoin=1,npoin
         do idime=1,ndime
            if(a%kfl_fixno(idime,ipoin)==1 .or. a%kfl_fixno(idime,ipoin)==7) then
               if(a%kfl_funno(ipoin)==0) then
                  a%bvess(idime,ipoin,1)=a%bvess(idime,ipoin,2)+meshve(idime,ipoin,1)
               end if
            end if
         end do
      end do 

      do iboun=1,nboun
         call a%Mesh%BoundaryLoad(iboun,e)
         if (a%kfl_fixbo(iboun) == 4) then
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               a%bvess(:,ipoin,1)=a%bvess(:,ipoin,2)*exnor(:,1)
            end do
         end if
      end do
   end if

! This does not work wirh AMR
   !Aitken
   if (associated(a%eveloc)) then
      do iboun=1,nboun
         call a%Mesh%BoundaryLoad(iboun,e)
         if (a%kfl_fixbo(iboun) == 51 .or. a%kfl_fixbo(iboun) == 52) then
            do inodb=1,e%pnodb
               ipoin=e%lnodb(inodb)
               !d_i1(:,:)  = d_i1(:,:)
               a%velocDD(:,ipoin,2) = a%bvess(:,ipoin,1) 
               !d_i2_guess(:,:)  = d_i2_guess(:,:)
               a%velocDD(:,ipoin,1) = a%eveloc(:,ipoin)
             end do
         end if
      end do  
   
      !r_i1(:,:)    = r_i2_past(:,:)
      a%bres(:,:,2) = a%bres(:,:,1)
      !r_i2(:,:)    = d_i2_guess(:,:)- d_i1(:,:)
      a%bres(:,:,1) = a%velocDD(:,:,1) - a%velocDD(:,:,2) 
      !r_i21(:,:)    = r_i2(:,:)    - r_i1(:,:)
      a%bres(:,:,3) = a%bres(:,:,1) - a%bres(:,:,2)
   
      if (a%kfl_doAitken) then
         r_i1   => a%bres(:,:,2)
         r_i21  => a%bres(:,:,3)
         call a%relaxVector('AITKE',a%relax,a%relax_max,r_i1,r_i21)
      endif
   
      w_r = a%relax
      do iboun=1,nboun
         call a%Mesh%BoundaryLoad(iboun,e)
         if (a%kfl_fixbo(iboun) == 51) then
            do inodb=1,e%pnodb
               ipoin=e%lnodb(inodb)
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               a%bvess(:,ipoin,1)=w_r*a%velocDD(:,ipoin,1) + (1.0_rp-w_r)*a%velocDD(:,ipoin,2)
               a%velocDD(:,ipoin,2)=w_r*a%velocDD(:,ipoin,1) + (1.0_rp-w_r)*a%velocDD(:,ipoin,2)
            end do
         elseif (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
            do inodb=1,e%pnodb
               ipoin=e%lnodb(inodb)
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               a%velocDD(:,ipoin,2)=w_r*a%velocDD(:,ipoin,1) + (1.0_rp-w_r)*a%velocDD(:,ipoin,2)
               a%tractDD(:,ipoin,2)=w_r*a%tractDD(:,ipoin,1) + (1.0_rp-w_r)*a%tractDD(:,ipoin,2)
            end do
         end if
      end do
   end if

   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_updbcs')
 
   !Spanwise forcing and other particular functions
   if (minval(a%kfl_funty) < 0.0_rp) then
      do ifun = 1,size(a%kfl_funty,1)
         !Blasius Plate initializations
         if (a%kfl_funty(ifun,1) == -2) then
            UFree = a%funpa(ifun)%a(1)
            nu = a%MatProp(1)%visco/a%MatProp(1)%densi
            X = a%funpa(ifun)%a(2)
            
            call BPGenerator%SetUFree(Ufree)
            call BPGenerator%SetNu(nu)
            call BPGenerator%SetDistanceXFromPlateStartingPoint(X)
         endif
      enddo
   
   
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         funno = a%kfl_funno(ipoin)
         if (funno .gt. 0) then 
         funty = a%kfl_funty(funno,:)
         else
         funty = 0
         end if
         !Spanwise forcing
         if (funty(1) == -1) then
            do idime = 1,ndime
               if (a%kfl_fixno(idime,ipoin) == 1) then
                  call a%Mesh%GetPointCoord(ipoin,pointCoord)
                  Wavelength = a%funpa(funno)%a(1)
                  Frequency = a%funpa(funno)%a(2)
                  CoordinateAxis = a%funpa(funno)%a(3)
                  coeff = sin(2*3.1415926535_rp/Wavelength*pointCoord(CoordinateAxis)-Frequency*a%ctime)
                  a%bvess(idime,ipoin,1) = a%bvess(idime,ipoin,2)*coeff
               endif
            enddo
            
         !Blasius boundary condition   
         elseif (funty(1) == -2) then   
            call a%Mesh%GetPointCoord(ipoin,pointCoord)
            call BPGenerator%GetU(pointCoord(2),a%bvess(1,ipoin,1))
         
         !Turbulent Inlet boundary condition
         elseif (funty(1) == -3) then
            call a%Mesh%GetPointCoord(ipoin,pointCoord)
            call a%TIBC%GetVelocity(pointCoord,a%ctime,a%bvess(1,ipoin,1))
         endif
      enddo
      
      !If there are turbulent inlet boundary conditions then we need to 
      !Ghost communicate since it might be different in different subdomains
      !(this was done in this manner for ease of implementation and expected low impact
      ! on the results)
      if ( ANY( a%kfl_funty==-3 ) ) then
         call a%Mesh%ArrayCommunicator%GhostCommunicate(size(a%bvess,1),a%bvess(:,:,1))
      endif   
   endif

end subroutine nsi_updbcs

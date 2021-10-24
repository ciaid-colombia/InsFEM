subroutine nsm_rotunk(a,itask)
!-----------------------------------------------------------------------
!****f* Nstinc/nsm_rotunk
! NAME 
!    nsm_rotunk
! DESCRIPTION
!    This routine performs rotations according to itask:
!
!    -Nodal unknowns from global to local axes to generate the initial 
!     condition for the solver (itask=1)
!
!    -Nodal unknowns computed in local axes to global ones (itask=2)
!
!    -Velocties from local values (given as initial condition) to global 
!     ones (itask=3)
!
!    -Velocties from global values (read from restart) to local 
!     ones (itask=4) (in order to prescribe bvess_nsi that is in local) 
!
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: itask
   integer(ip)             :: ibopo,iroty,itot0,itot1,itotv,npoin,ndime,ipoin
   real(rp), allocatable                :: worma(:,:), worve(:)
   real(rp), pointer       :: exnor(:,:) => NULL()


   
   
   
   if(a%kfl_local==0) then
      ! Modifications need to be done only if there exist image nodes or
      ! boundary conditions in skew systems 
      return 

   else if(itask==2) then  ! unkno - local to global
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,ndime,worma,'worma','nsm_rotunk')
      call a%Memor%alloc(ndime,worve,'worve','nsm_rotunk')
      
      do ipoin=1,npoin
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if(ibopo>0) then
            iroty=a%kfl_fixrs(ipoin)
            ! Tangent skew system
            if(iroty==-1) then
               worve = a%unkno(1:ndime,ipoin)
               call mbvab0(worma,exnor,worve,ndime,ndime)
               a%unkno(1:ndime,ipoin) = worma(1:ndime,1)

            ! Tangent skew system
            else if(iroty==-2) then
                 call runend('nsm_rotunk skews not ready')

            ! Given skew system
            else if(iroty>=1) then
               call runend('nsm_rotunk skews not ready')
            end if
         endif   
      end do
      
      call a%Memor%dealloc(ndime,ndime,worma,'worma','nsm_rotunk')
      call a%Memor%dealloc(ndime,worve,'worve','nsm_rotunk')
   
   end if
   
end subroutine nsm_rotunk


subroutine nsi_iniunk(a)
!NAME 
!   nsi_iniunk
!DESCRIPTION
!   This routine sets up the initial conditions for the a%velocity.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   use Mod_NsiExacso    
   
   implicit none
   class(NavierStokesProblem) :: a
   type(NsiExacso)  :: exacso   
   
   real(rp)    :: vmodu,vfreq,dummr
   integer(ip) :: dummi
   integer(ip) :: icomp,ipoin,idime,ndime,npoin
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:)
   !Nodal coordinates
   real(rp), pointer       :: coord(:) => NULL()
   real(rp) :: auxcoord(3) ,auxveloc(3)
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   
   do ipoin=1,npoin
      do idime=1,ndime
         if((a%kfl_fixno(idime,ipoin)==1) .or. (a%kfl_fixno(idime,ipoin) == 0) .or. (a%kfl_fixno(idime,ipoin) == 7)) then
            do icomp = 3,a%ncomp
               a%veloc(idime,ipoin,icomp) = a%bvess(idime,ipoin,1)
            enddo
         end if
      end do
   end do
   
   if (a%kfl_CoriolisForce == 1) then
      do ipoin= 1,npoin
         call a%Mesh%GetPointCoord(ipoin,coord)
         auxcoord(:) = 0.0_rp
         auxcoord(1:ndime) = coord
         
         call vecpro(-a%CoriolisW,auxcoord,auxveloc,3)
         do icomp = 1,a%ncomp
            a%veloc(:,ipoin,icomp) = auxveloc(1:ndime)
         enddo
      enddo
   endif
  

   if(a%kfl_incnd==1) then
      call runend('Nsi_iniunk: Initial conditions not implemented')
   end if

   if(a%kfl_exacs/=0) then

      ! Allocate exact components
      call a%Memor%alloc(ndime,exvel,'exvel','nsi_iniunk')   
      call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsi_iniunk')     
         do ipoin = 1,npoin 
            call a%Mesh%GetPointCoord(ipoin,coord)
            
            call exacso%nsi_ComputeSolution(ndime,coord,a%ctime,a)
            call exacso%nsi_GetVelocity(ndime,exvel,exveg)           

            do idime=1,ndime
               a%veloc(idime,ipoin,a%ncomp) = exvel(idime)
            end do
         end do
      ! Allocate exact components
      call a%Memor%dealloc(ndime,exvel,'exvel','nsi_iniunk') 
      call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsi_iniunk')     

   end if

   !Assign u(n,i,*) <-- u(n-1,*,*), initial guess after initialization (or reading restart)
   do icomp = 1,a%ncomp-1
      a%veloc(1:ndime,1:npoin,icomp) = a%veloc(1:ndime,1:npoin,a%ncomp)
      a%press(1:npoin,icomp) = a%press(1:npoin,a%ncomp) 
   enddo

   if (a%kfl_RVE == 1) call a%RVE(1)

   call a%Ifconf

end subroutine nsi_iniunk


subroutine lmn_iniunk(a)
!NAME 
!   lmn_iniunk
!DESCRIPTION
!   This routine sets up the initial conditions for the a%velocity.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   integer(ip) :: icomp,ipoin,idime,ndime,npoin

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   do ipoin=1,npoin
      do idime=1,ndime
         if((a%kfl_fixno(idime,ipoin)==1) .or. (a%kfl_fixno(idime,ipoin) == 0)) then
            do icomp = 3,a%ncomp
               a%veloc(idime,ipoin,icomp) = a%bvess(idime,ipoin,1)
               a%tempe(ipoin,a%ncomp) = a%bvess(ndime+1,ipoin,1)
            end do
         end if
      end do
   end do
  
   if (a%kfl_incnd==1) then
      call runend('Lmn_iniunk: Initial conditions not implemented')
   end if

   !Assign u(n,i,*) <-- u(n-1,*,*), initial guess after initialization (or reading restart)
   a%pther(a%ncomp) = a%itpre
   a%itemp(1:npoin) = a%tempe(1:npoin,a%ncomp)
   do icomp = 1,a%ncomp-1
      a%veloc(1:ndime,1:npoin,icomp) = a%veloc(1:ndime,1:npoin,a%ncomp)
      a%press(1:npoin,icomp) = a%press(1:npoin,a%ncomp) 
      a%tempe(1:npoin,icomp) = a%tempe(1:npoin,a%ncomp) 
      a%pther(icomp) = a%pther(a%ncomp)
   enddo

   call a%Ifconf

end subroutine lmn_iniunk


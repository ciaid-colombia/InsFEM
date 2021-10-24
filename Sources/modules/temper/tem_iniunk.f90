subroutine tem_iniunk(a)
!DESCRIPTION
!   This routine sets up the initial conditions for the a%temperature.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_Temperature
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: icomp,ipoin,npoin
   real(rp), pointer :: coord(:) => NULL()
   
   call a%Mesh%GetNpoin(npoin)
   
   do ipoin=1,npoin
      if((a%kfl_fixno(1,ipoin)==1) .or. (a%kfl_fixno(1,ipoin)==0)) then
         do icomp = 3,a%ncomp
            a%tempe(ipoin,icomp) = a%bvess(1,ipoin,1)
         enddo
      end if
   end do

   !Assign t(n,i,*) <-- t(n-1,*,*), initial guess after initialization (or reading restart)
   a%tempe(1:npoin,1) = a%tempe(1:npoin,3) 
   a%tempe(1:npoin,2) = a%tempe(1:npoin,3)

end subroutine

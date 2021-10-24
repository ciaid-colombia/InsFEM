subroutine sld_iniunk(a)
!NAME 
!   sld_iniunk
!DESCRIPTION
!   This routine sets up the initial conditions for the a%displacements.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)          :: icomp,ipoin,idime,ndime,npoin

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   do ipoin=1,npoin
      do idime=1,ndime
         if((a%kfl_fixno(idime,ipoin)==1) .or. (a%kfl_fixno(idime,ipoin) == 0)) then
            do icomp = 3,a%ncomp
                a%disp(idime,ipoin,a%ncomp) = a%bvess(idime,ipoin,1)
            end do
         end if
      end do
   end do
  
   if(a%kfl_incnd==1) then
      call runend('sld_iniunk: Initial conditions not implemented')
   end if

   !Assign u(n,i,*) <-- u(n-1,*,*), initial guess after initialization (or reading restart)
   do icomp = 1,a%ncomp-1
      a%disp(1:ndime,1:npoin,icomp) = a%disp(1:ndime,1:npoin,a%ncomp)
   enddo


end subroutine sld_iniunk


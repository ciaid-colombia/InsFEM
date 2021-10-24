subroutine plcd_iniunk(a)
!NAME 
!   plcd_iniunk
!DESCRIPTION
!   This routine sets up the initial conditions for the a%displacements.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use Mod_Mesh  
   use Mod_memor  
   implicit none
   class(PLCDProblem) :: a
   type(FemMesh)    :: Mesh  
   type(MemoryMan)  :: Memor   
   
   integer(ip) :: icomp,ipoin,idime,ndime,npoin

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   do ipoin=1,npoin
      do idime=1,ndime
         if((a%cs%kfl_fixno(idime,ipoin)==1) .or. (a%cs%kfl_fixno(idime,ipoin) == 0)) then
            do icomp = 3,a%ncomp
                a%Displacement(idime,ipoin,a%ncomp) = a%cs%bvess(idime,ipoin)*a%css%CurrentLoadFactor
            end do
         end if
      end do
   end do
  
   if(a%kfl_incnd==1) then
      call runend('plcd_iniunk: Initial conditions not implemented')
   end if

   !Assign u(n,i,*) <-- u(n-1,*,*), initial guess after initialization (or reading restart)
   do icomp = 1,a%ncomp-1
      a%Displacement(1:ndime,1:npoin,icomp) = a%Displacement(1:ndime,1:npoin,a%ncomp)
   enddo


end subroutine plcd_iniunk


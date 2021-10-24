subroutine lev_iniunk(a)
!DESCRIPTION
!   This routine sets up the initial conditions for the a%levelset.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_LevelSet
   use Mod_Element
   use Mod_CutMesh
   use def_parame
   implicit none
   class(LevelSetProblem) :: a

   integer(ip) :: icomp,ipoin,npoin,ielem,nelem,nnode,inode
   !Nodal coordinates
   real(rp), pointer       :: coord(:) => NULL()

   interface
      subroutine lev_CutElementsAndListByLayers(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine

      subroutine lev_ReinitLevel(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a

      end subroutine

   end interface


   call a%Mesh%GetNpoin(npoin)

   !Read the Gid data only to determine the fluid type and the initial interface
   do ipoin=1,npoin
      if((a%kfl_fixno(1,ipoin)==1) .or. (a%kfl_fixno(1,ipoin)==0) .or. a%kfl_fixno(1,ipoin) == 2) then
         do icomp = 1,a%ncomp
            a%level(ipoin,icomp) = a%bvess(1,ipoin,1)
         enddo
      end if
   end do

   if(a%kfl_ExactLevel==1)then
      !only for testing
      do ipoin=1,npoin
         do icomp=1,a%ncomp
            call a%Mesh%GetPointCoord(ipoin,coord)
               a%level(ipoin,icomp)=(0.5999_rp-coord(2))
         end do
      end do
   end if

   !Level Set Initialization
   call lev_CutElementsAndListByLayers(a)

   if (a%kfl_InitialRedistance == 0) then
      do icomp = 2,size(a%level,2)
         a%level(:,icomp) = a%level(:,1)
      enddo
   else
      call lev_ReinitLevel(a)
   endif

   !Assign t(n,i,*) <-- t(n-1,*,*), initial guess after initialization (or reading restart)
   a%level(1:npoin,1) = a%level(1:npoin,3)
   a%level(1:npoin,2) = a%level(1:npoin,3)

end subroutine

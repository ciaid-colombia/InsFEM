subroutine nsi_InitializeTurbulentInletBoundaryConditions(a)
   use typre
   use Mod_NavierStokes
   use MPI
   implicit none
   class(NavierStokesProblem) :: a
   
   integer(ip) :: ipoin,npoin,funno,funty,ierr
   
   real(rp) :: PlaneCoordinates(2,3), gPlaneCoordinates(2,3)
   real(rp) :: minrange(3), maxrange(3), mingrange(3),maxgrange(3),nu
   real(rp), pointer :: coord(:) => NULL()
   integer(ip) :: kfl_perio,auxcount
   logical :: PerDime(3)
   
   REAL(8), PARAMETER :: D_QNAN = TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
   
   !We require the coordinates of the inflow plane for TIBC
   PlaneCoordinates(1,:) = 1e24
   PlaneCoordinates(2,:) = -1e24
   auxcount = 0
   call a%Mesh%GetNpoin(npoin)
   do ipoin = 1,npoin
      funno = a%kfl_funno(ipoin)
      if(funno /= 0) then
         funty = a%kfl_funty(funno,1)
         !Turbulent inlet boundary condition
         if (funty == -3) then 
            call a%Mesh%GetPointCoord(ipoin,coord)
            PlaneCoordinates(1,:) = min(PlaneCoordinates(1,:),coord)
            PlaneCoordinates(2,:) = max(PlaneCoordinates(2,:),coord)
         endif
      endif 
   enddo
   
   minrange = PlaneCoordinates(1,:)
   maxrange = PlaneCoordinates(2,:)
    
   call  MPI_AllREDUCE( minrange, mingrange, 3, MPI_REAL8, MPI_MIN,a%MPIcomm, ierr )
   call  MPI_AllREDUCE( maxrange, maxgrange, 3, MPI_REAL8, MPI_MAX,a%MPIcomm, ierr )
   
   gPlaneCoordinates(1,:)= mingrange
   gPlaneCoordinates(2,:)= maxgrange    
   
   call a%TIBC%SetGlobalPlaneCoordinates(gPlaneCoordinates)
   
   call a%TIBC%SetPlaneCoordinates(PlaneCoordinates)
   
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) then
      call a%Mesh%GetArePeriodicDimensions(PerDime)
      call a%TIBC%SetArePeriodicBoundariesYZ(PerDime(2),PerDime(3))
   endif
   
   nu = a%MatProp(1)%visco/a%MatProp(1)%densi
   call a%TIBC%SetViscosity(nu)
   
   !Now we can initialize it
   call a%TIBC%Initialize
end subroutine



subroutine nsi_FinalizeTurbulentInletBoundaryConditions(a)
   use typre
   use Mod_NavierStokes
   use MPI
   implicit none
   class(NavierStokesProblem) :: a
   
   
   call a%TIBC%Finalize
   
end subroutine
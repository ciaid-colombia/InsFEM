subroutine lev_reampi(a)
   use typre
   use MPI
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_LevelSet
   implicit none
   
   class(LevelSetProblem) :: a
   integer :: ierr,igauge,ndime

   !Communicate lev_reaphy
   CALL MPI_BCAST(a%kfl_advec, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_ReinitLevel, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_ExactLevel, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   
   !Communicate lev_reanut
   CALL MPI_BCAST(a%staco, size(a%staco), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%NstepsReinitLevel, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_ForceEulerianAdvection, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Adaptivity
   CALL MPI_BCAST(a%OutLayers, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%GeneralRefinementLevels, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%InterfaceRefinementLevels, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   !Mass Correction
   CALL MPI_BCAST(a%kfl_MassCorrection, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Initial Redistancing
   CALL MPI_BCAST(a%kfl_InitialRedistance, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Communicate lev_reaous
   !HeightGauges
   CALL MPI_BCAST(a%nHeightGauges, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   if (a%MPIrank /= a%MPIroot) then
      allocate(a%HeightGauges(a%nHeightGauges))
      call a%Memor%allocObj(0,'HeightGauges','lev_reaous',a%nHeightGauges)
   endif
   call a%Mesh%GetNdime(ndime)
   do igauge = 1,a%nHeightGauges
      CALL MPI_BCAST(a%HeightGauges(igauge)%Origin, size(a%HeightGauges(igauge)%Origin), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%HeightGauges(igauge)%DirectionVector, size(a%HeightGauges(igauge)%DirectionVector), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      call a%HeightGauges(igauge)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%Mpirank)
      call a%HeightGauges(igauge)%SetNdime(ndime)
   enddo
   
   
   
   !For fixing CN
   if (a%kfl_tsche_1st_datafile == 'CNOBS') a%kfl_tsche_1st_datafile = 'CN   '
   
end subroutine

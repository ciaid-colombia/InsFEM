subroutine tem_reampi(a)
   use typre
   use MPI
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_Temperature
   implicit none
   
   class(TemperatureProblem) :: a
   integer :: ierr
   integer(ip) :: imaterial

   !Communicate tem_reaphy
   CALL MPI_BCAST(a%kfl_advec, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_sourc, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotur, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_joule, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_CouplingThreeField, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%ipsou, size(a%ipsou), MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%DT_kfl_DustTransport, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%react, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%prtur, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%sourc, size(a%sourc), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%turbu, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%rpsou, size(a%rpsou), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   if (a%ipsou(1) /= 0) then
      if (a%MPIrank /= a%MPIroot) call a%Memor%alloc(a%ipsou(2),a%tfsou,'tfsou','tem_reaphy')
      CALL MPI_BCAST(a%tfsou, size(a%tfsou), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   endif
   
   CALL MPI_BCAST(a%visco, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%react, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%DT_ParticleSize, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%DT_MediumDensity, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%DT_MediumDynamicViscosity, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%DT_GravityForce, 3, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   !Communicate tem_reanut
   CALL MPI_BCAST(a%kfl_repro, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_adapsgs, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_shock, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_wtemp, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_wlapl, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_trasg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_tacsg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_nolsg, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_stabm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%staco, size(a%staco), MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%shock, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%relax, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%relsg, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   !Communicate tem_reaous
   CALL MPI_BCAST(a%kfl_dispa, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%Avg1DIdime, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   call MPI_BCAST(a%NumberOfMaterials,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
   if (a%MPIrank /= a%MPIroot) then
      allocate(a%Materials(a%NumberOfMaterials))
      call a%Memor%allocObj(0,'Materials','tempe_memall',1*a%NumberOfMaterials)
   endif
   do imaterial = 1, a%NumberOfMaterials 
      call a%Materials(imaterial)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Materials(imaterial)%ScatterData
   end do
   
   
   !Forces and moments
   CALL MPI_BCAST(a%kfl_outfm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !For fixing CN
   if (a%kfl_tsche_1st_datafile == 'CNOBS') a%kfl_tsche_1st_datafile = 'CN   '
   
end subroutine

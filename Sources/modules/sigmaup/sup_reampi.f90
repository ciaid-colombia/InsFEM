subroutine sup_reampi(a)
   use typre
   use MPI
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: ierr
   integer(ip) :: ndime
   integer(ip) :: imat=1
   
   !Default
   call a%Mesh%Getndime(ndime)
   
   if (a%MatProp(imat)%lawvi < 0) then
      a%ndofbc = ndime + (ndime-1)*(ndime-1) + 2 
      a%ndofbcstart= 0
   end if
   
   CALL MPI_BCAST(a%erplit, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)    
   CALL MPI_BCAST(a%erulit, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%erslit, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%erph1t, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)    
   CALL MPI_BCAST(a%eruh1t, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr) 
   CALL MPI_BCAST(a%ersh1t, 2, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr) 
   
   !Communicate sup_reaphy
   CALL MPI_BCAST(a%kfl_bc_number, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%PTT_model, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%Giesekus_model, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotem, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotem_WLF, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotem_Boussinesq, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_cotem_Arrhenius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%ReferenceTemp, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%c1_WLF, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%c2_WLF, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%alpha_Arrhenius, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
  
   !LCR conformation
   CALL MPI_BCAST(a%LogFormulation, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_linearConstitutiveTerms, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_linearConvectiveTerm, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%nu0_LCR, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_ntausmooth, 1, MPI_INTEGER4, a%MPIroot, MPI_COMM_WORLD, ierr)  
      
   !Communicate sup_reanut
   CALL MPI_BCAST(a%kfl_tiacc, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr) !new
   
   CALL MPI_BCAST(a%incremental, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_splitOSSMomentum, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_splitOSSConstitutive, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   CALL MPI_BCAST(a%kfl_reproBoundZero, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   CALL MPI_BCAST(a%kfl_reproTemporalTermZero, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
end subroutine

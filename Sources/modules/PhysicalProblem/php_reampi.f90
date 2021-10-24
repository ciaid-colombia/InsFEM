subroutine php_reampi(a)
   use MPI
   use typre
   use Mod_PhysicalProblem
   use Mod_BroadCastBuffer
   implicit none
   class(PhysicalProblem), intent(inout) :: a

   integer(ip) :: ierr,i
   type(BroadCastBuffer) :: BBuffer


   a%istep = 0_ip

   call BBuffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call BBuffer%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call BBuffer%Initialize(200,800)

   !PhysicalProblem
   call BBuffer%Add(a%kfl_timei)
   call BBuffer%Add(a%kfl_stead)
   call BBuffer%Add(a%dtinv)
   call BBuffer%Add(a%dtinv2)
   call BBuffer%Add(a%kfl_SwitchOff)

   !Numerical Treatment
   call BBuffer%Add(a%kfl_linea)
   call BBuffer%Add(a%maxit)
   call BBuffer%Add(a%npica)
   call BBuffer%Add(a%kfl_exacs)
   call BBuffer%Add(a%exact)
   call BBuffer%Add(a%kfl_normc)
   call BBuffer%Add(a%kfl_ProjectionType)
   call BBuffer%Add(len(a%EndLoopQuadrature),a%EndLoopQuadrature)
   call BBuffer%Add(a%sstol)
   call BBuffer%Add(a%cotol)
   call BBuffer%Add(a%cptol)
   call BBuffer%Add(a%kfl_initRelax)
   call BBuffer%Add(a%kfl_adap)

   call BBuffer%Add(len(a%kfl_tsche_1st_datafile),a%kfl_tsche_1st_datafile)
   call BBuffer%Add(len(a%kfl_tsche_1st_current),a%kfl_tsche_1st_current)
   call BBuffer%Add(len(a%kfl_tsche_2nd_datafile),a%kfl_tsche_2nd_datafile)
   call BBuffer%Add(len(a%kfl_tsche_2nd_current),a%kfl_tsche_2nd_current)
   call BBuffer%Add(a%neule)
   call BBuffer%Add(a%neule_2nd)
   call BBuffer%Add(a%kfl_elmat_datafile)
   call BBuffer%Add(a%kfl_tiacc)
   call BBuffer%Add(a%kfl_twost)
   call BBuffer%Add(a%safet)
   call BBuffer%Add(len(a%RefinerErrorEstimator),a%RefinerErrorEstimator)
   call BBuffer%Add(len(a%RefinerErrorCriteria),a%RefinerErrorCriteria)
   call BBuffer%Add(a%RefinerErrorLimits)
   call BBuffer%Add(a%kfl_DeleteDirichletColumns)
   call BBuffer%Add(a%nrest)
   call BBuffer%Add(a%kfl_restrictions)

   !Output Strategy
   call BBuffer%Add(a%kfl_BoundaryConditionsReadStrategy)
   call BBuffer%Add(a%ManufacturedBoundaryCondition)
   call BBuffer%Add(a%npp_inits)
   call BBuffer%Add(a%npp_stepi)
   call BBuffer%Add(a%nptra)
   call BBuffer%Add(a%pos_tinit)
   call BBuffer%Add(a%pos_times)

   call BBuffer%Add(a%kfl_printBC)

   call BBuffer%Add(a%kfl_eigensolve)
   call BBuffer%Add(a%kfl_error)
   call BBuffer%Add(5,a%extramatrix)


   call BBuffer%BroadCast
   call BBuffer%Dealloc

!   if (a%nrest>0) then
!      do i=1,a%nrest
!         call MPI_BCAST(a%keyrest(i),5,MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
!      enddo
!   endif


   call a%SpecificReaMPI


end subroutine

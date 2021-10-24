module Mod_plcd_TD_StochasticTopologyOptimizationData

   use typre
   use Mod_MPIObject
   use Mod_Listen
   use Mod_Element
   implicit none

   type, extends(MPIObject) :: TD_StochasticData
      integer(ip) :: TransversalMPIcomm
      integer(ip) :: NForcePoints,NTransversalParallelProcesses,NStochasticIntegrationPoints, IntegrationPoint,ForcePoint=1,FunctionIntegrationPointIniDimension(6), FunctionForcePointInidimension(6),TransversalPoint
      real(rp), allocatable :: IntegrationPointPositions(:,:), IntegrationPointWeights(:)
      
      real(rp) :: alpha
      integer(ip) :: level = 3
      integer(ip) :: StochasticDimensions = 0
      integer(ip) :: VarianceDelaySteps = 0
      
      real(rp) :: FunctionForces(3,3,6)   !3 forces per function at most
      integer(ip) :: NFunctionForces(6) = 0
      integer(ip) :: FunctionIntegrationPointDimension(6)
      integer(ip) :: kfl_SimultaneousLoads = 1

      
contains
      procedure :: ReadData
      procedure :: ScatterData
    end type

contains

   subroutine ReadData(a,Listener)
   implicit none
      class(TD_StochasticData) :: a
      type(ListenFile) :: Listener
      
      integer(ip) :: iforce,ifun

      do while(Listener%words(1)/='ENDST')
         call Listener%listen('plcd_reaphy')
         if(Listener%words(1) == 'ALPHA')   then
            a%ALPHA = Listener%param(1)
         elseif(Listener%words(1) == 'LEVEL')   then
            a%level = Listener%param(1)
         elseif(Listener%words(1) == 'VARDE')   then
            a%VarianceDelaySteps = Listener%param(1)
         !Function forces
         elseif(Listener%words(1) == 'FUNCT') then
            call Listener%listen('plcd_reaphy')
            do while (Listener%exists('ENDFU') .eqv. .false.)
               ifun = Listener%param(1)
               a%NFunctionForces(ifun) = Listener%param(2)
               do iforce = 1,a%NFunctionForces(ifun)
                  a%FunctionForces(:,iforce,ifun) = Listener%param(3+(iforce-1)*3:2+iforce*3)
               enddo
               call Listener%listen('plcd_reaphy')
            enddo
         elseif(Listener%words(1) == 'SIMUL') then   
            if (Listener%exists('YES  ')) then
               a%kfl_SimultaneousLoads = 1
            else
               a%kfl_SimultaneousLoads = 0
            endif
         endif
       enddo
       !This is just so that in the external loop there are no problems
       Listener%words(1) = '     '
   end subroutine

   subroutine ScatterData(a)
      use MPI
      implicit none
      class(TD_StochasticData) :: a
      
      integer(ip) :: counter, ifun

      integer(ip) :: ierr
       call MPI_BCAST(a%ALPHA,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
       call MPI_BCAST(a%level,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
       call MPI_BCAST(a%VarianceDelaySteps,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
       call MPI_BCAST(a%NFunctionForces,6,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
       call MPI_BCAST(a%FunctionForces,54,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
       call MPI_BCAST(a%kfl_SimultaneousLoads,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)

       
       counter = 0
       do ifun = 1,6
         a%FunctionIntegrationPointDimension(ifun) = counter
         counter = counter + a%NFunctionForces(ifun)
       enddo

       !       call MPI_BCAST(a%chimin,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
!       call MPI_BCAST(a%chimax,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
!       call MPI_BCAST(a%ChiSpaceChoice,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
! !      call MPI_BCAST(a%kappa,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
! !      call MPI_BCAST(a%rhopenalty,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   end subroutine







end module

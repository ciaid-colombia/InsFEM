module Mod_plcd_TDData
   use typre
   use Mod_MPIObject
   use Mod_Listen
   use Mod_Element
   implicit none

   type, extends(MPIObject) :: TDDataType
      real(rp) :: VolumeFraction = 0.2_rp
      real(rp) :: chimin = 1e-5_rp, chimax = 1.0_rp
      real(rp) :: nu = 0.2_rp
      integer(ip) :: ChiSpaceChoice = 1 !0: ElementWise, 1: PointWise
      
      integer(ip) :: NVolumeIterations = 15
      real(rp)    :: VolumeTolerance = 0.01
      real(rp)    :: KappaIncreaseCoefficient(2) = (/1.5_rp,0.5_rp/), KappaPowerCoefficient = 0.5_rp, KappaMin = 0.0
      integer(ip) :: kfl_ProgressiveVolume = 0, ProgressiveVolume_inidelay, ProgressiveVolume_nsteps 
      real(rp) :: ProgressiveVolume_inivol, ProgressiveVolume_endvol
      

      real(rp), pointer :: NodalChi(:) => NULL()

      !When to perform a Topology Optimization Algorithm
      integer(ip) :: kfl_WhenPerformTopologyOptimization = 0 !0: Each iteration, 1:Each Time Step, 2:Each Substage, 3:Each Stage
      integer(ip) :: kfl_PerformTopologyIteration = 1 !0: Don't do it, 1: Do it.
      
contains
      procedure :: ReadData
      procedure :: ScatterData
    end type

contains

   subroutine ReadData(a,Listener)
   implicit none
      class(TDDataType) :: a
      type(ListenFile) :: Listener

      do while(Listener%words(1)/='ENDTO')
         call Listener%listen('plcd_reaphy')
         if(Listener%words(1) == 'VOLUM')   then
            a%VolumeFraction = Listener%param(1)
         elseif(Listener%words(1) == 'CHIMI')   then
            a%chimin = Listener%param(1)
         elseif(Listener%words(1) == 'CHIMA')   then
            a%chimax = Listener%param(1)
         elseif(Listener%words(1) == 'CHISP')   then
            a%ChiSpaceChoice = Listener%param(1)
         elseif(Listener%words(1) == 'ITEVO')   then
            a%NVolumeIterations = Listener%param(1)
         elseif(Listener%words(1) == 'TOLVO')   then
            a%VolumeTolerance = Listener%param(1)   
         elseif(Listener%words(1) == 'KPLUS') then
            a%KappaIncreaseCoefficient(1) = Listener%param(1)
         elseif(Listener%words(1) == 'KMINU') then
            a%KappaIncreaseCoefficient(2) = Listener%param(1)
         elseif(Listener%words(1) == 'KMINI') then
            a%KappaMin = Listener%param(1)
         elseif(Listener%words(1) == 'KPOWE') then
            a%KappaPowerCoefficient = Listener%param(1)   
         elseif(Listener%words(1) == 'PROVO' .and. Listener%words(2) == 'ON   ')   then
            a%kfl_ProgressiveVolume = 1
            a%ProgressiveVolume_inidelay = Listener%param(2)
            a%ProgressiveVolume_inivol = Listener%param(3)
            a%ProgressiveVolume_nsteps = Listener%param(4)
            a%ProgressiveVolume_endvol = Listener%param(5)
         elseif(Listener%words(1) == 'FREQU') then
            if (Listener%words(2) == 'ITERA') then
               a%kfl_WhenPerformTopologyOptimization = 0
            elseif (Listener%words(2) == 'TIMES') then
               a%kfl_WhenPerformTopologyOptimization = 1
               a%kfl_PerformTopologyIteration = 0
            elseif (Listener%words(2) == 'SUBST') then
               a%kfl_WhenPerformTopologyOptimization = 2
               a%kfl_PerformTopologyIteration = 0   
            elseif (Listener%words(2) == 'STAGE') then
               a%kfl_WhenPerformTopologyOptimization = 3
               a%kfl_PerformTopologyIteration = 0
            endif
!         elseif(Listener%words(1) == 'KAPPA')   then
!            a%Kappa = Listener%param(1)
!         elseif(Listener%words(1) == 'RHOPE')   then
!            a%rhopenalty = Listener%param(1)
         endif
       enddo
       !This is just so that in the external loop there are no problems
       Listener%words(1) = '     '
   end subroutine

   subroutine ScatterData(a)
      use MPI
      implicit none
      class(TDDataType) :: a

      integer(ip) :: ierr
      call MPI_BCAST(a%VolumeFraction,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%chimin,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%chimax,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ChiSpaceChoice,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%NVolumeIterations,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%VolumeTolerance,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%KappaIncreaseCoefficient,2,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%KappaPowerCoefficient,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%KappaMin,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
!      call MPI_BCAST(a%kappa,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
!      call MPI_BCAST(a%rhopenalty,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%kfl_ProgressiveVolume,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ProgressiveVolume_inidelay,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ProgressiveVolume_inivol,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ProgressiveVolume_nsteps,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ProgressiveVolume_endvol,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%kfl_WhenPerformTopologyOptimization,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%kfl_PerformTopologyIteration,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
end module

module Mod_MasterVariables
  use typre
  use Mod_memor
  use Mod_timer
  use Mod_postpr
  use Mod_ReadWrite
  use Mod_Listen
  use Mod_ParallelLibraryInterface
  use Mod_MeshInterpolator
  use Mod_DistributedContainerList
  
  type MasterVariables

     !MPI
     integer(ip) :: MPIrank 
     integer(ip) :: MPIsize
     integer(ip) :: MPIroot = 0
     integer(ip) :: MPIcomm
     
     
     !Run data
     character(150) :: title            !Problem title
     character(150) :: namda,oldnamda   !Dat file name
     character(150) :: BaseDataFolder   !Base data folder (does not have the rank if partitioned)
     character(150) :: OldDataFolder    !Old data folder (does not have the rank if partitioned)
     character(150) :: DataFolder       !Folder for reading the data (with rank number if partitioned)
     character(150) :: ResultsFolder
     character(150) :: PostProcessFolder
     character(150) :: RestartFolder,OldRestartFolder 
     character(66)  :: ReadTypeString   !Type of ReadingStrategy

     integer(ip) ::      &
          kfl_outfo,     &    ! Output format
          kfl_outco,     &    ! Output format compression
          kfl_writeType, &    ! Output write type
          kfl_iofor           ! I/O format
     
     integer(ip) :: kfl_MPIComType

     integer(ip) :: kfl_ParallelOrdering              = 0
     integer(ip) :: kfl_ParallelCommunicator          = 0
     integer(ip) :: kfl_ParallelRebalancePartitioning = 0
     class(ParallelLibraryInterface), pointer :: ParallelLibrary => NULL()
          
     type(MemoryMan) :: MasterMemo
     integer(ip)     :: kfl_memwrite = 0     
     
     type(Timer) ::  cpu_start(10)    ! CPU for starting operations
                                      ! 1: Turnon and Set Channels
                                      ! 2: Domain Operations
                                      ! 3: Output Domain
                                      ! 5: Reapro
                                      ! 6: Initialize Parallel Library
     type(Timer) :: cpu_end(10)       ! CPU for ending operations
     type(Timer) :: cpu_total         ! Total CPU time
     real(rp)    :: cpu_initi
     real(rp)    :: cpu_limit

          
     !Logical units
     integer(ip) ::   &
          lun_pdata = 0,  &
          lun_outpu = 0,  &
          lun_memor = 0,  &
          lun_conve = 0,  &
          lun_direc = 0,  &
          lun_rstar = 0

     
     type(ListenFile) :: MasterListener

     class(PostprFile), pointer :: FilePostpr => NULL()
     class(Reader),     pointer :: Readerpr   => NULL()
     class(Writer),     pointer :: Writerpr   => NULL()
     
     !for flushing
     integer(ip) :: kfl_flush = 0
          
     !for adaptive meshes output
     integer(ip) :: meshCounter = 0
  
     !Physical problem
     integer(ip) ::              &
          kfl_timco,             &    ! Time coupling strategy
          kfl_goblk,             &    ! Block coupling converged
          kfl_gocou=1,           &    ! Global problem converged
          kfl_gotim=1,           &    ! Global problem evolving in time
          mitgl,                 &    ! Maximum # of global iterations for each block
          nsmax,                 &    ! Maximum number of steps
          istep=0,               &    ! Current time step
          iiter,                 &    ! Current global iteration
          cpiter,                &    ! Current global iteration
          kfl_RemeshingStrategy, &
          kfl_multicomm = 1,     &
          MulticommColor = 0
          

     real(rp) ::           &
          timei,           &    ! Initial time
          timef,           &    ! Final time
          timfu(5),        &    ! Time step size as a function of the step
          dtime,           &    ! Time step size,   & dt
          dtin0,           &    ! Time step size,   & dt
          muldt,           &    ! for Time step size
          dtinv = 0.0_rp,  &    ! 1/dt
          ctime = 0.0_rp,  &    ! Current time
          tolgl,           &    ! Tolerance for global iterations for each block
          Interp_Tol            ! Interpolator search tolerance
          
     integer(ip) ::     times   ! Number of time samples
     
     integer(ip) :: kfl_multicase = 0
     integer(ip) :: kfl_maxCoupledIters = 100    ! Initialize to useless high number
     logical     :: kfl_doCoupledConv = .false.
     real(rp)    :: couplingBeta= 0              ! For adaptive coupling tolerance
     
     !For postprocessing
     integer(ip) :: StartPostprocessAt = 0
     integer(ip) :: PostprocessProcessorEvery = 0

     type(Interpolator) :: Int_Restart

  end type   
 
end module

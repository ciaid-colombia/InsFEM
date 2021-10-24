module Mod_caseReapro 
   use typre
   use def_parame
   use Mod_iofile
   use Mod_postpr
   use Mod_int2str
   use Mod_CaseVariables
   use Mod_GeneralCase
   use Mod_GeneralParallelLibrary
   use Mod_ParallelLibraryInterface
   use Mod_Postpr_NULL
   use Mod_Postpr_VTK
   use Mod_ReadWrite_PETSc
   use Mod_ReadWrite_F90
   use Mod_BroadCastBuffer
   use typre
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   use MPI
   use Mod_CaseInterpolator
   use Mod_DC_CaseInterpolator
   use Mod_DriverCreator
   implicit none

   interface
      subroutine CreateDriverFromKey(key,DriverCreator)
         use Mod_DriverCreator
         implicit none
         character(5) :: key
         type(DriverCreatorType) :: DriverCreator
      end subroutine
   end interface

contains

   subroutine case_openfi(a,itask)
      implicit none
      class(GeneralCase), target :: a
      integer(ip), intent(in) :: itask

      type(masterVariables), pointer :: m => NULL()
      integer(ip)             :: ioer1,ioer2,ioer3,ioer4,ioer5,ioer6
      integer(ip)             :: iargc,ilcha
      character(150)          :: fil_pdata_dom                  ! File names
      character(150)          :: fil_pdata,fil_outpu,fil_memor, fil_postp, fil_rstar, fil_postp_graf
      character(150)          :: fil_conve,fil_direc
      character(150)          :: fil_outpu_dom
      character(150)          :: extra = ''

      m  => a%caseVars%masterVars

      call cputim(m%cpu_initi)

      select case (itask)

      case (1)
         fil_pdata = trim(m%DataFolder)//'/'//adjustl(trim(m%namda))//'.dat'
         call iofile(zero,m%lun_pdata,fil_pdata,'DATA','old')

      case (2)
         fil_memor = trim(m%ResultsFolder)//'/'//adjustl(trim(m%namda))&
             &//adjustl(trim(int2str(m%MPIrank)))//'.mem'

         !Since memor is open for debugging, we open it always
         if (m%kfl_memwrite /= 0) then
            call iofile(zero,m%lun_memor,fil_memor,'MEMORY EVOLUTION')
         elseif (m%kfl_memwrite == 0) then
            m%lun_memor = 0
         endif


      case (3)

         ! Compose filenames
         fil_outpu = trim(m%ResultsFolder)//'/'//adjustl(trim(m%namda))&
             &//adjustl(trim(int2str(m%MPIrank)))//'.log'
         fil_postp_graf = trim(m%PostProcessFolder)//'/'//adjustl(trim(m%namda))&
             &//adjustl(trim(int2str(m%MPIrank)))//'.post.grf'

         ! Open files
         !Log file, only for root
         m%lun_outpu = -1
         if (m%MPIrank == m%MPIroot) then
            call iofile(zero,m%lun_outpu,fil_outpu,'RUN EVOLUTION')
            !Initialize runend
            call runend_init(m%lun_outpu)
            ! Write header
            write(m%lun_outpu,100) adjustl(trim(m%title))
         end if

         !Open postprocess file
         call m%FilePostpr%SetOutputFormat(m%kfl_outfo)
         call m%FilePostpr%SetOutputCompression(m%kfl_outco,m%kfl_writeType)
         call m%FilePostpr%SetFlush(m%kfl_flush)
         call m%FilePostpr%SetMPI(m%MPIsize,m%MPIroot,m%MPIrank,m%MasterMemo)
         call m%FilePostpr%Initialize(m%PostprocessFolder,m%namda)
         call m%FilePostpr%OpenPostFile(m%PostprocessFolder,m%namda)
         call m%FilePostpr%OpenGrafFile(fil_postp_graf)
         if (m%kfl_multicomm > 1) call m%FilePostpr%SetMulticomm(m%kfl_multicomm,m%MulticommColor)

         !Only if root, open cvg files
         if (m%MPIrank == m%MPIroot) then
            fil_conve = trim(m%ResultsFolder)//'/'//adjustl(trim(m%namda))//'.cvg'
            fil_direc = trim(m%ResultsFolder)//'/'//adjustl(trim(m%namda))//'.liv'

            ! Open files
            call iofile(zero,m%lun_conve,fil_conve,'CONVERGENCE HISTORY')
            call iofile(zero,m%lun_direc,fil_direc,'LIVE INFORMATION')
         endif


      end select


      ! Format
      100 format(///,&
               5x,'*******',/,&
               5x,'Yet another finite element software: ',a,/,&
               5x,'*******',/)
      110 format(/,5x,'>>>  a IS A TIME RESTART RUN: CONTINUING...',//)

   end subroutine

   subroutine case_rrudat(a)
      use typre
      use Mod_Listen
      use Mod_GeneralCase
      use Mod_caseVariables
      implicit none
      class(GeneralCase), target :: a


      type(masterVariables), pointer :: m => NULL()

      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      character(150) :: title
      integer(ip), pointer  :: nnpar => NULL(), nnwor => NULL()

      m  => a%caseVars%masterVars

      !Initializations.
      m%kfl_ParallelOrdering = 0                            !Deafault is My
      m%kfl_ParallelCommunicator = 0                        !Default is Petsc
      m%kfl_ParallelRebalancePartitioning = 0               !Default is Petsc

      m%kfl_outfo = 1                                       !Output format
      m%kfl_outco = 1                                       !Default is binary
      m%kfl_writeType= 0                                    !Default is parallelWriter
      m%kfl_iofor = 0                                       !I/O format
      m%kfl_memwrite = 0
      m%cpu_limit = 1.0e20                                  !Default CPU limit
      m%title = ' '                                         !Problem title
      m%kfl_MPIComType = 1                                  !Default is blocking
      m%kfl_flush = 0                                       !Default is do not flush

      !Begin.
      call m%MasterListener%SetLunits(m%lun_pdata,m%lun_outpu)
      call m%MasterListener%getarrs(words,param,nnpar,nnwor)

      call m%MasterListener%listen('RRUDAT')
      if(words(1)/='RUNDA') call runend('RRUDAT: WRONG DATA FILE TYPE')

      !Read, identify & execute appropriate command.
      call m%MasterListener%listen('RRUDAT')
      do while(words(1)/='ENDRU')
         select case (words(1))
         case ('CASE ')
               !Read and write title.
               m%title = words(2)
            case('MEMCH')
               if(m%MasterListener%exists('ON   ')) then
                  m%kfl_memwrite = 1                                                  !Write memory
               endif
            case ('CPULI')
               !Read the CPU limit.
               if(param(1)/=0.0_rp) m%cpu_limit = param(1)
            case ('OUTPU')
               !Read the output format
               if(m%MasterListener%exists('VTK  ')) then      !VTK
                  m%kfl_outfo = 1
               else if(m%MasterListener%exists('NULL ')) then      !Null
                  m%kfl_outfo = -1
               end if
               if(m%MasterListener%exists('ASCII')) then            !ASCII
                  m%kfl_outco = 0
               else if(m%MasterListener%exists('BINAR')) then       !Binary
                  m%kfl_outco = 1
               else if(m%MasterListener%exists('APPEN')) then       !Appended
                  m%kfl_outco = 2
               end if
            case ('WTYPE')
               if(m%MasterListener%exists('PARAL')) then            !Normal Parallel write
                   m%kfl_writeType= 0
               else if(m%MasterListener%exists('BLOCK')) then       !Block parallel write
                   m%kfl_writeType= 1
               end if
            case ('ORDER')
               !Read the Parallel Library to be used
               if(m%MasterListener%exists('PETSC')) then
                  m%kfl_ParallelOrdering = 0
               elseif (m%MasterListener%exists('MY   ')) then
                  m%kfl_ParallelOrdering = 1
               endif
            case ('COMMU')
               !Read the Parallel Library to be used
               if(m%MasterListener%exists('PETSC')) then
                  m%kfl_ParallelCommunicator = 0
               elseif (m%MasterListener%exists('MY   ')) then
                  m%kfl_ParallelCommunicator = 1
               endif
            case ('REBAL')
               !Read the Parallel Library to be used
               if(m%MasterListener%exists('PETSC')) then
                  m%kfl_ParallelRebalancePartitioning = 0
               elseif (m%MasterListener%exists('ZOLTA')) then
                  m%kfl_ParallelRebalancePartitioning = 1
               endif
            case ('MPICO')
               !Read the Parallel communications to be used
               if(m%MasterListener%exists('BLOCK')) then
                  m%kfl_MPIComType = 0
               elseif (m%MasterListener%exists('NONBL')) then
                  m%kfl_MPIComType = 1
               endif
            case ('FLUSH')
               !Read the output mode to be used
               if(m%MasterListener%exists('ON   ')) then
                  m%kfl_flush = 1
               elseif (m%MasterListener%exists('OFF  ')) then
                  m%kfl_flush = 0
               endif
            case ('IOFOR')
               !Read the output format
               if(m%MasterListener%exists('NATIV')) then
                  m%kfl_iofor = 0
               else if(m%MasterListener%exists('PETSC')) then
                  m%kfl_iofor = 1
               end if

            case ('POSTP')
               call m%MasterListener%listen('RRUDAT')
               do while (m%MasterListener%words(1) /= 'ENDPO')
                  if (m%MasterListener%words(1) == 'START') then
                     m%StartPostprocessAt = m%MasterListener%param(1)
                  elseif (m%MasterListener%words(1) == 'PROCE') then
                     m%PostprocessProcessorEvery = m%MasterListener%param(1)
                  endif
                  call m%MasterListener%listen('RRUDAT')
               enddo

            case ('CUSTO')

            case default
               !call runend('RRUDAT: ERROR IN RUN DATA BLOCK')
         end select

         call m%MasterListener%listen('RRUDAT')
      enddo
   end subroutine

   subroutine case_readat(a)
      implicit none
      class(GeneralCase), target :: a

      integer(ip) :: imodu,iorde,ipara

      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(), nnwor => NULL()

      type(masterVariables), pointer :: m => NULL()
      type(adaptiveVariables), pointer :: ad => NULL()
      m => a%caseVars%masterVars
      ad => a%caseVars%adaptiveVars

      !Initializations.
      m%timei         = 0.0_rp          !Initial time
      m%timef         = 0.0_rp          !Final time
      m%nsmax         = 0               !Max. number of steps
      m%dtime         = 0.0_rp          !Time step size dt

      m%kfl_timco     = 0               !Time coupling strategy
      m%mitgl         = 1               !Do at least one global iteration
      ad%NumberOfInitialUniformRefinementSteps = 0               !Levels of initial refinement
      m%Interp_Tol  = 0.0001_rp        !Interpolator tolerance

      ad%kfl_AdaptiveRefinement = 0                          !Adaptive Refinement
      ad%RefinementLeader = 'TEMPE'                          !Temperature Leads in Adaptive Refinement
      ad%RefinerRebalancingRatio = 1.5                       !Refine with a 1.5 real to optimum npoinLocal ratio
      ad%kfl_21Balancing = 0                                 !Default is do not do 21 balancing
      ad%AdaptiveDelay = -1

      m%kfl_RemeshingStrategy = 0

      call m%MasterListener%SetLunits(m%lun_pdata,m%lun_outpu)
      call m%MasterListener%getarrs(words,param,nnpar,nnwor)

      call m%MasterListener%listen('RPRODA')
      if(words(1)/='PROBL')&
            call runend('RPRODA: WRONG PROBLEM_DATA CARD')

      !Read data.
      do while(words(1)/='ENDPR')
         call m%MasterListener%listen('RPRODA')
         if(words(1)=='TIMEI') then                        !Time interval
            m%timei=param(1)
            m%timef=param(2)
         else if(words(1)=='TIMES') then                   !Time step size
            m%dtime=param(1)
            if(m%dtime>0.0) m%dtin0 = 1.0_rp/m%dtime
            if(m%dtime>0.0) m%dtinv = 1.0_rp/m%dtime
         else if(words(1)=='TIMEC') then                   !Time coupling strategy
            if(words(2)=='GLOBA') then
               m%kfl_timco=0
               if(m%MasterListener%exists('PRESC')) then
                  m%kfl_timco=0                                  !prescribed dt
               else if(m%MasterListener%exists('FROMC')) then
                  m%kfl_timco=1                                  !dt=min(dt,f_c*dt_c)
               else if(m%MasterListener%exists('FUNCT')) then
                  m%kfl_timco=3                                  !dt as a function of istep
                  m%timfu=param(3:7)
               end if
            else if(words(2)=='LOCAL') then
               m%kfl_timco=2                                     !dt=dt(ielem)
            end if
         else if(words(1)=='NUMBE') then                   !Maximum Number of steps
            m%nsmax = int(param(1))
         else if(words(1)=='MAXIM') then                   !Maximum number of iterations
            m%mitgl = int(param(1))
         else if(words(1)=='TOLER') then                   !Tolerance for block iterations
            m%tolgl = param(1)
         else if(words(1)=='CPCON') then                   !Check for coupled convergence in multicase
            if(m%MasterListener%exists('YES  ')) then
                m%kfl_doCoupledConv = .true.
            end if
         else if(words(1)=='CPMAX') then                   !Max coupled iterations in multicase
            m%kfl_maxCoupledIters = int(param(1))
         else if(words(1)=='CPBET') then                   !Adaptive stoping criteria
            m%couplingBeta = param(1)
         else if(words(1)=='ADAPT') then
            !Read the output mode to be used
            if(m%MasterListener%exists('ON   ')) then
               ad%kfl_AdaptiveRefinement = 1
               ad%RefinementLeader = words(3)
               ad%RefinerRebalancingRatio = m%MasterListener%getrea('REBAL',1.5_rp,'#Rebalancing Ratio')
               if(m%MasterListener%exists('21BAL')) then
                  ad%kfl_21Balancing = 1
               endif
               ad%RefinerMaxLevel = m%MasterListener%getint('MAXRE',0_ip,'#Maximum Refinement Level')
               ad%ProgressiveMaxLevel = m%MasterListener%getint('PROGR',0_ip,'#Maximum Refinement Level')
               ad%AdaptiveDelay = m%MasterListener%getint('DELAY',-1_ip,'#Adaptive Delay')
               ad%UpdateFreq = m%MasterListener%getint('UPDAT',-1_ip,'#Update IO')
            endif
         else if(words(1)=='MULTI' .and. words(2) == 'ADAPT') then
            !Read the output mode to be used
            if(m%MasterListener%exists('ON   ')) then
               ad%kfl_AdaptiveMulticomm = 1
            endif   

         else if(words(1)=='REMES') then
            if(m%MasterListener%exists('FIXED')) m%kfl_RemeshingStrategy = 2
            if(m%MasterListener%exists('ADVEC')) m%kfl_RemeshingStrategy = 3
         else if(words(1)=='INITI') then
            if(m%MasterListener%exists('ON   ')) then
               ad%NumberOfInitialUniformRefinementSteps = int(param(2))
            endif
         else if((words(1)=='INTER') .and. (words(2) == 'TOLER')) then
            m%Interp_Tol = param(2)
         end if
      end do
      
      !Compute DTINV and the time for the first time step, as well as the
      !maximum number of time steps allowed.
      m%dtinv = 0.0_rp
      m%times=(m%timef-m%timei)/m%dtime
      if(m%dtime>0.0)   m%dtinv = 1.0_rp/m%dtime
      if(m%timef<= m%timei) then
         m%dtime = 0.0_rp
         m%dtin0 = 0.0_rp
         m%dtinv = 0.0_rp
      end if
      m%ctime = m%timei
      if(m%nsmax==0) m%nsmax = int(m%timef*m%dtinv)+1

   end subroutine

   subroutine case_ReaproMPI(a,itask)
      class(GeneralCase), target :: a
      integer(ip) :: itask

      integer(ip) :: ierr
      type(BroadCastBuffer) :: BBuffer

      type(masterVariables), pointer :: m => NULL()
      type(adaptiveVariables), pointer :: ad => NULL()

      m  => a%caseVars%masterVars
      ad => a%caseVars%adaptiveVars

      if (itask == 1) then

         call BBuffer%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         !call BBuffer%SetLOutputFile(101,101)
         call BBuffer%Initialize(100,100)

         call BBuffer%Add(len(m%title),m%title)
         call BBuffer%Add(m%kfl_memwrite)
         call BBuffer%Add(m%kfl_outfo)
         call BBuffer%Add(m%kfl_writeType)
         call BBuffer%Add(m%kfl_outco)
         call BBuffer%Add(m%cpu_limit)
         call BBuffer%Add(m%kfl_ParallelOrdering)
         call BBuffer%Add(m%kfl_ParallelCommunicator)
         call BBuffer%Add(m%kfl_ParallelRebalancePartitioning)
         call BBuffer%Add(m%kfl_MPIComType)
         call BBuffer%Add(m%kfl_flush)
         call BBuffer%Add(m%StartPostprocessAt)
         call BBuffer%Add(m%PostprocessProcessorEvery)
         call BBuffer%Add(m%kfl_iofor)

         call BBuffer%BroadCast
         call BBuffer%Dealloc

      elseif (itask == 2) then

         call BBuffer%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         !call BBuffer%SetLOutputFile(101,101)
         call BBuffer%Initialize(100,100)

         call BBuffer%Add(m%times)
         call BBuffer%Add(m%nsmax)
         call BBuffer%Add(m%kfl_timco)
         call BBuffer%Add(m%mitgl)
         call BBuffer%Add(m%timei)
         call BBuffer%Add(m%timef)
         call BBuffer%Add(m%dtime)
         call BBuffer%Add(m%dtin0)
         call BBuffer%Add(m%timfu)
         call BBuffer%Add(m%muldt)
         call BBuffer%Add(m%ctime)
         call BBuffer%Add(m%tolgl)
         call BBuffer%Add(ad%kfl_AdaptiveRefinement)
         call BBuffer%Add(ad%kfl_AdaptiveMulticomm)
         call BBuffer%Add(5,ad%RefinementLeader)
         call BBuffer%Add(ad%RefinerRebalancingRatio)
         call BBuffer%Add(ad%RefinerMaxLevel)
         call BBuffer%Add(ad%ProgressiveMaxLevel)
         call BBuffer%Add(ad%AdaptiveDelay)
         call BBuffer%Add(ad%kfl_21Balancing)
         call BBuffer%Add(ad%UpdateIO)
         call BBuffer%Add(ad%UpdateFreq)
         call BBuffer%Add(m%kfl_RemeshingStrategy)
         call BBuffer%Add(ad%NumberOfInitialUniformRefinementSteps)
         call BBuffer%Add(m%Interp_Tol)
         call BBuffer%Add(m%kfl_doCoupledConv)
         call BBuffer%Add(m%kfl_maxCoupledIters)
         call BBuffer%Add(m%couplingBeta)

         call BBuffer%BroadCast
         call BBuffer%Dealloc
         
         !Adaptive Multicomm Initializations
         if (m%MulticommColor == 0 .and. ad%kfl_AdaptiveMulticomm == 1) then
            ad%kfl_AdaptiveRefinement = 0
            ad%NumberOfInitialUniformRefinementSteps = 0
         endif
      endif
   end subroutine

   subroutine case_AfterReadInitializations(a)
      implicit none
      class(GeneralCase), target :: a

      integer(ip) :: ierr

      type(GeneralParallelLibrary), pointer :: GeneralLibrary => NULL()

      type(masterVariables), pointer :: m => NULL()
      type(adaptiveVariables), pointer :: ad => NULL()
      
      m => a%caseVars%masterVars
      ad => a%caseVars%adaptiveVars

      call m%MasterMemo%init(m%lun_memor,m%lun_outpu)

      !-------------------------------------------------------------------
      !ParallelLibrary
      GeneralLibrary => GeneralParallelLibrary_Const()
      m%ParallelLibrary => GeneralLibrary
      call m%MasterMemo%AllocObj(0,'ParallelLibrary','CaseReapro',1)

      !Initialize Communicator type
      if (m%kfl_ParallelCommunicator == 0) then
         call GeneralLibrary%SetCommunicatorType('Petsc ')
      elseif (m%kfl_ParallelCommunicator == 1) then
         call GeneralLibrary%SetCommunicatorType('My    ')
      else
         call runend('InitializeParallelLibrary: Wrong CommunicatorType')
      endif

      !Initialize Ordering type
      if (m%kfl_ParallelOrdering == 0) then
         call GeneralLibrary%SetOrderingType('Petsc ')
      elseif (m%kfl_ParallelOrdering == 1) then
         call GeneralLibrary%SetOrderingType('My    ')
      else
         call runend('InitializeParallelLibrary: Wrong OrderingType')
      endif

      !Initialize Rebalancer type
      if (m%kfl_ParallelRebalancePartitioning == 0) then
         call GeneralLibrary%SetRebalancePartitionerType('Petsc ')
      elseif (m%kfl_ParallelRebalancePartitioning == 1) then
#ifdef ZOLTAN
      call GeneralLibrary%SetRebalancePartitionerType('Zoltan')
#endif
      else
         call runend('InitializeParallelLibrary: Wrong OrderingType')
      endif

      call GeneralLibrary%Initialize(trim(m%DataFolder)//'/'//adjustl(trim(m%namda))//'.sol ')

      !-------------------------------------------------------------
      !FilePostprocessor
      select case (m%kfl_outfo)
      case default
         m%FilePostpr => PostprFile_VTK_Const()
      case (-1)
         m%FilePostpr => PostprFile_NULL_Const()
      end select
      call m%MasterMemo%AllocObj(0,'FilePostpr','CaseReapro',1)

      select case (m%kfl_iofor)
      case (1)
         m%Readerpr => Reader_PETSc_Const()
         m%Writerpr => Writer_PETSc_Const()
         call m%Readerpr%SetReaderMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         call m%Writerpr%SetWriterMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
      case default
         m%Readerpr => Reader_F90_Const()
         m%Writerpr => Writer_F90_Const()
         call m%Readerpr%SetReaderMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         call m%Writerpr%SetWriterMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
      end select
      call m%MasterMemo%AllocObj(0,'Readerpr','CaseReapro',1)
      call m%MasterMemo%AllocObj(0,'Writerpr','CaseReapro',1)

   end subroutine

   subroutine case_ReadInterpolators(a)
      implicit none
      class(GeneralCase), target :: a

      integer(ip) :: ierr

      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(), nnwor => NULL()

      character(5) :: CaseInterpolatorKeys(25)

      type(masterVariables), pointer  :: m => NULL()
      type(caseVariables), pointer :: c => NULL()

      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      class(DistributedContainer), pointer :: myDC => NULL()

      integer(ip) :: ninterp = 0, iinterp

      type(BroadCastBuffer) :: bb

      m => a%caseVars%masterVars
      c => a%caseVars

      !Initialize the Interpolator List
      call a%CaseInterpolatorList%Initialize

      if (m%kfl_multicase == 1) then
         !Only MPIroot
         if (m%MPIrank == m%MPIroot) then
            call m%MasterListener%SetLunits(m%lun_pdata,m%lun_outpu)
            call m%MasterListener%getarrs(words,param,nnpar,nnwor)

            !Reach the section
            call m%MasterListener%rewind
            do while(words(1)/='RUNDA')
            call m%MasterListener%listen('RPRODA')
            end do
            do while(words(1)/='ENDRU')
            call m%MasterListener%listen('RPRODA')
            end do
            do while(words(1)/='PROBL')
            call m%MasterListener%listen('RPRODA')
            end do
            do while(words(1)/='ENDPR')
               call m%MasterListener%listen('RPRODA')
            enddo
            !Reach the interpolation section
            do while(words(1)/='INTER')
               call m%MasterListener%listen('RPRODA')
            enddo
            !Read the data
            ninterp = 0
            do while(words(1)/='ENDIN')
               call m%MasterListener%listen('RINTDA')

               !FROM_CASE
               if (words(1) == 'FROMC') then
                  ninterp = ninterp + 1
                  CaseInterpolatorKeys(ninterp) = words(2)

                  myCaseInterpolator => CaseInterpolator_Const()
                  call myCaseInterpolator%SetMyCase(a)
                  myDC => DC_CaseInterpolator_Const(myCaseInterpolator)
                  call myDC%SetKey(words(2))
                  call m%MasterMemo%allocObj(0,'CaseInterpolator','ReadInterpolator',1)

                  call CaseInterpolatorDataFromListener(myCaseInterpolator,m%MasterListener)

                  call a%CaseInterpolatorList%Add(myDC)
               endif
            enddo
         endif

         !Now the other processors receive the number of drivers
         call MPI_BCAST(ninterp, 1, MPI_INTEGER4, m%MPIroot, m%MPIcomm, ierr)
         call MPI_BCAST(CaseInterpolatorKeys, 5*ninterp, MPI_CHARACTER, m%MPIroot, m%MPIcomm, ierr)

         !The other processors need to add the interpolator to their list
         if (m%MPIrank /= m%MPIroot) then
            do iinterp = 1,ninterp

               !Get Particular Driver Info based on the name of the module
               !Also create the driver
               myCaseInterpolator => CaseInterpolator_Const()
               call myCaseInterpolator%SetMyCase(a)
               myDC => DC_CaseInterpolator_Const(myCaseInterpolator)
               call myDC%SetKey(CaseInterpolatorKeys(iinterp))
               call m%MasterMemo%allocObj(0,'CaseInterpolator','ReadInterpolator',1)

               call a%CaseInterpolatorList%Add(myDC)
            enddo
         endif

         !Now we loop through the interpolators, and scatter all of their info
         !MPI Communications
         call bb%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         call bb%SetLOutputFile(m%lun_memor,m%lun_outpu)
         call bb%Initialize(100,100)

         call a%CaseInterpolatorList%GetFirst(myDC)
         do while (associated(myDC))
            call ExtractCaseInterpolator(myDC,myCaseInterpolator)

            call AddCaseInterpolatorToBroadCastBuffer(myCaseInterpolator,bb)
            call a%CaseInterpolatorList%GetNext(myDC)
         enddo

         call bb%BroadCast
         call bb%Dealloc
      end if

   end subroutine

   subroutine case_ReadDrivers(a) 
      implicit none
      class(GeneralCase), target :: a
      type(DriverCreatorType) :: DriverCreator

      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(), nnwor => NULL()

      character(5) :: DriverKeys(25), DriverNicks(25)
      integer(ip)  :: nDrivers

      type(masterVariables), pointer  :: m => NULL()
      type(caseVariables), pointer :: c => NULL()

      class(DriverInterface), pointer :: MyDriver => NULL()
      class(DistributedContainer), pointer :: myDC => NULL()

      integer(ip) :: ierr
      integer(ip) :: idriver

      m => a%caseVars%masterVars
      c => a%caseVars

      !Initialize the DriverList
      call a%DriverList%Initialize

      !Only MPIroot
      if (m%MPIrank == m%MPIroot) then
         call m%MasterListener%SetLunits(m%lun_pdata,m%lun_outpu)
         call m%MasterListener%getarrs(words,param,nnpar,nnwor)

         !Reach the section
         call m%MasterListener%rewind
         do while(words(1)/='RUNDA')
           call m%MasterListener%listen('RPRODA')
         end do
         do while(words(1)/='ENDRU')
           call m%MasterListener%listen('RPRODA')
         end do
         do while(words(1)/='PROBL')
           call m%MasterListener%listen('RPRODA')
         end do

         ndrivers = 0
         !Read the data
         do while(words(1)/='ENDPR')
            call m%MasterListener%listen('RPRODA')

            !Get Particular Driver Info based on the name of the module
            !Also create the driver
            call CreateDriverFromKey(words(1),DriverCreator)
            !If the entry for a driver has been found
            if (DriverCreator%key /= 'NODRI' .and. words(2) == 'ON   ') then

               if (words(3) == 'MODNA') then
                  DriverCreator%nick = words(4)
               else
                  DriverCreator%nick = DriverCreator%key
               endif

               call a%caseVars%masterVars%MasterMemo%allocObj(0,'Driver','DriverCreator',1)

               !Point to the created Driver
               myDriver => DriverCreator%Driver
               !Initialize it
               call myDriver%Initialize(c)

               !-------------------------------------------------------------
               !ONLY the root reads the info
               !We will need to scatter it through REAPROMPI later
               call myDriver%Reapro(c,DriverCreator%endword)
               !-------------------------------------------------------------

               !Add it to the list of drivers
               myDC => DC_Driver_Const(myDriver)
               call myDC%SetKey(DriverCreator%key)
               call myDC%SetNick(DriverCreator%nick)
               call a%DriverList%Add(myDC)

               !Update the list of drivers to broadcast
               ndrivers = ndrivers+1
               DriverKeys(ndrivers) = DriverCreator%key
               DriverNicks(ndrivers) = DriverCreator%nick
            endif
         enddo
      endif

      !Now the other processors receive the number of drivers
      call MPI_BCAST(ndrivers, 1, MPI_INTEGER4, m%MPIroot, m%MPIcomm, ierr)
      call MPI_BCAST(DriverKeys, 5*ndrivers, MPI_CHARACTER, m%MPIroot, m%MPIcomm, ierr)
      call MPI_BCAST(DriverNicks, 5*ndrivers, MPI_CHARACTER, m%MPIroot, m%MPIcomm, ierr)
      c%ndrivers=ndrivers

      !The other processors need to add the drivers to their list
      if (m%MPIrank /= m%MPIroot) then
         do idriver = 1,ndrivers

            !Get Particular Driver Info based on the name of the module
            !Also create the driver
            call CreateDriverFromKey(DriverKeys(idriver),DriverCreator)
            call a%caseVars%masterVars%MasterMemo%allocObj(0,'Driver','DriverCreator',1)
            !Point to the created driver
            myDriver => DriverCreator%Driver
            !Initialize it
            call myDriver%Initialize(c)


            !Add it to the list of drivers
            myDC => DC_Driver_Const(myDriver)
            call myDC%SetKey(DriverCreator%key)
            call myDC%SetNick(DriverNicks(idriver))
            call a%DriverList%Add(myDC)
         enddo
      endif

      !call myDriver%SetDriverMPI(c)

      !All the processors
      !Now we loop through all the drivers and we scatter the info read by the root
      call a%DriverList%GetFirst(myDC)
      do while (associated(myDC))

         call ExtractDriver(myDC,myDriver)

         !Scatter the info amongst all processors
         call myDriver%ReaproMPI(c)

         !Now that everything is read we can initialize the channel lists
         call myDriver%InitializeChannelLists

         call a%DriverList%GetNext(myDC)
      enddo

   end subroutine

end module

subroutine case_Reapro(a)
   use Mod_CaseReapro
   implicit none

   class(GeneralCase), target :: a

   type(masterVariables), pointer :: m => NULL()
   m => a%caseVars%masterVars

   !This is the first call, we start the global timer
   call m%cpu_total%Tic

   call m%cpu_start(5)%Tic

   if (m%MPIrank == m%MPIroot) then
      !Get data file name and open it.
      call case_openfi(a,1)

      !Read run data.
      call case_rrudat(a)
   endif
   call case_ReaproMPI(a,1)

   !Open the memory file if necessary
   call case_openfi(a,2)

   call m%cpu_start(5)%Toc

   call m%cpu_start(6)%Tic

   !Initialize things with the just read data
   call case_AfterReadInitializations(a)

   call m%cpu_start(6)%Toc

   call m%cpu_start(5)%Tic

   !Get result file names and open them.
   call case_openfi(a,3)

   if (m%MPIrank == m%MPIroot) then

      !Read general problem data.
      call case_readat(a)

   endif
   call case_ReaproMPI(a,2)

   !Read modules existence
   call case_ReadDrivers(a)

   !Read interpolator info between cases
   call case_ReadInterpolators(a)

   call m%cpu_start(5)%Toc

end subroutine



module Mod_podromDriver
   use typre
   use MPI
   use Mod_BroadCastBuffer
   use Mod_Listen
   use Mod_podrom
   use Mod_caseVariables
   use Mod_PhysicalProblemDriver
   use Mod_PhysicalProblem
   use Mod_DistributedContainer 
   use Mod_DC_rp
   use Mod_InChannel
   use Mod_MasterVariables
   use Mod_DriverCreator
   implicit none

   private

   public :: podromDriver, podromDriver_Const

   type, extends(DriverInterface) :: podromDriver 
      
      class(PodromProblem), pointer :: Podrom => NULL()
      
      character(5) :: PhysicalType,PhysicalKey
      
      logical :: AreWeRunning = .false.     !false: build, true: run
      
      !The physical driver is MINE
      !I own it!
      class(DriverInterface), pointer :: PhysicalDriver => null()
    
contains
      
      procedure :: Lev1Reapro             => podrom_Lev1Reapro
      procedure :: Lev1ReaproMPI          => podrom_Lev1ReaproMPI
      procedure :: InitializeChannelLists => podrom_InitializeChannelLists
      procedure :: FinalizeChannelLists   => podrom_FinalizeChannelLists
      procedure :: SetInChannels          => podrom_SetInChannels
      procedure :: SetOutChannels         => podrom_SetOutChannels
      procedure :: UpdateChannels         => podrom_UpdateChannels
      procedure :: Turnon                 => podrom_Turnon
      procedure :: GetTimeStep            => podrom_GetTimeStep
      procedure :: SetTimeStep            => podrom_SetTimeStep
      procedure :: GetRefinementCriteria  => podrom_GetRefinementCriteria
      procedure :: PreRefine              => podrom_PreRefine
      procedure :: Refine                 => podrom_Refine
      procedure :: Rebalance              => podrom_Rebalance
      procedure :: Begste                 => podrom_Begste
      procedure :: MeshProjections        => podrom_MeshProjections
      procedure :: MeshAdvections         => podrom_MeshAdvections
      procedure :: Doiter                 => podrom_Doiter
      procedure :: Convergence            => podrom_Convergence
      procedure :: CoupledConvergence     => podrom_CoupledConvergence
      procedure :: Endste                 => podrom_Endste
      procedure :: Output                 => podrom_Output
      procedure :: Turnof                 => podrom_Turnof
      procedure :: WriteTimes             => podrom_WriteTimes
      
   end type podromDriver

   interface podromDriver_Const
      procedure constructor
   end interface 
   
   interface 
      subroutine CreateDriverFromKey(key,DriverCreator)
         use Mod_DriverCreator
         implicit none
         character(5) :: key
         type(DriverCreatorType) :: DriverCreator
      end subroutine
   end interface

contains

   function constructor()
       class(podromDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine podrom_Lev1Reapro(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      type(DriverCreatorType) :: DriverCreator
      
      if(a%Listener%words(1) == 'PHYSI') then
         if (a%Listener%exists('LOWMA')) then
            a%PhysicalType = 'LMACH'
            a%PhysicalKey  = 'LMACH'
         elseif (a%Listener%exists('TEMPE')) then
            a%PhysicalType = 'TEMPE'
            a%PhysicalKey  = 'TEMPE'
         elseif (a%Listener%exists('NAVIE')) then
            a%PhysicalType = 'NSTIN'
            a%PhysicalKey  = 'NSTIN'
         elseif (a%Listener%exists('SOLID')) then
            a%PhysicalType = 'SOLID'
            a%PhysicalKey  = 'SOLID'
         elseif (a%Listener%exists('IRRSP')) then
            a%PhysicalType = 'IRRSP'
            a%PhysicalKey  = 'SOLID'
         elseif (a%Listener%exists('COMPR')) then
            a%PhysicalType = 'NSCOM'
            a%PhysicalKey  = 'NSCOM'
         elseif (a%Listener%exists('BOUSS')) then
            a%PhysicalType = 'BOUSS'
            a%PhysicalKey  = 'NSTIN'
         else
            call runend('PodRom: physical Problem not ready for the ROM')
         endif
         
         !Now the interior physical driver is initialized
         !and told to read
         call CreateDriverFromKey(a%PhysicalKey,DriverCreator)
         if (DriverCreator%key == 'NODRI') call runend('wrong type of driver, rom')
         a%PhysicalDriver => DriverCreator%Driver
         call a%Memor%allocobj(0,'RomPhysicalDriver','podrom_turnon',1)
     
         a%kfl_preli = 1 ! comment this line for AMR
         a%kfl_rstar = 1
         
         !Initialize, driver!
         call a%PhysicalDriver%Initialize(c)
         
         !Read your things, driver!
         call a%PhysicalDriver%Reapro(c,DriverCreator%endword)
         
      endif
   end subroutine
   
   subroutine podrom_Lev1ReaproMPI(a,c,bb)
      class(podromDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      type(DriverCreatorType) :: DriverCreator
      
      integer(ip) :: ierr
      
      !We need to be quick with the Driver type because we need 
      !to do their reaproMPI. MPIBCAST instead of buffer.
      call MPI_BCAST(a%PhysicalType,5,MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%PhysicalKey,5,MPI_CHARACTER,a%MPIroot,a%MPIcomm,ierr)
      
      if (a%MPIrank /= a%MPIroot) then
         !Now the other processors create their drivers
         call CreateDriverFromKey(a%PhysicalKey,DriverCreator)
         a%PhysicalDriver => DriverCreator%Driver
         call a%Memor%allocobj(0,'RomPhysicalDriver','podrom_turnon',1)
      
         !Initialize, driver!
         call a%PhysicalDriver%Initialize(c)
      endif
      
      !Now we all do the reaproMPI
      call a%PhysicalDriver%ReaproMPI(c)
      
   end subroutine
   
   subroutine podrom_InitializeChannelLists(a)
      class(PodRomDriver) :: a
      
      !Initialize your channels, driver!
      call a%PhysicalDriver%InitializeChannelLists
      
      !Give me your channels, driver!
      !The Rom Channels point to the physical driver channels
      call a%PhysicalDriver%GetChannelLists(a%OutChannelList,a%InChannelList)
      
   end subroutine
   
   subroutine podrom_FinalizeChannelLists(a)
      class(PodRomDriver) :: a

      !Finalize your channels, driver!
      call a%PhysicalDriver%FinalizeChannelLists
   end subroutine
   
   subroutine podrom_SetOutChannels(a,task)
      class(podromDriver) :: a
      character(*) :: task
      
      !Set your OutChannels, driver!
      call a%PhysicalDriver%SetOutChannels(task)
   end subroutine
   
   subroutine podrom_SetInChannels(a)
      class(podromDriver) :: a
      
      !Set your inChannels, driver!
      call a%PhysicalDriver%SetInChannels
   end subroutine
   
   subroutine podrom_UpdateChannels(a)
      implicit none
      class(podromDriver) :: a
      
      !Update your channels, driver!
      call a%PhysicalDriver%UpdateChannels
   end subroutine
   
   subroutine podrom_Turnon(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      class(PhysicalProblem), pointer :: PhysicalPro => NULL()
      integer(ip) :: kfl_itask,PhysicalAlgor
         
      !Depending on the type of ROM, create a rom type
      a%PodRom => PodRomProblem_Const()

      call a%Memor%allocobj(0,'Podrom','podrom_turnon',1)
      
      !We turn on the rom
      !General data
      call a%PodRom%SetMPI(c%masterVars%MPIcomm,c%masterVars%MPIsize,c%masterVars%MPIroot,c%masterVars%MPIrank)
      call a%PodRom%SetReadType(c%masterVars%ReadTypeString)
      call a%PodRom%SetInputFolder(c%masterVars%BaseDataFolder)
      call a%PodRom%SetRestartFolder(c%masterVars%RestartFolder)
      call a%PodRom%SetOldRestartFile(c%masterVars%RestartFolder)
      call a%PodRom%SetInputFile(c%masterVars%namda)
      call a%PodRom%SetOutputFolder(c%masterVars%ResultsFolder)
      call a%PodRom%SetOutputFiles(c%masterVars%lun_memor,c%masterVars%lun_outpu)
      call a%PodRom%SetIOPointers(c%masterVars%Writerpr,c%masterVars%Readerpr)
      call a%PodRom%SetFlush(c%masterVars%kfl_flush)
      call a%PodRom%SetOldInputFile(c%masterVars%namda)
      
      !Turnon, driver!
      call a%PhysicalDriver%Turnon(c)
      
      !Give us your PhysicalProblem, driver!
      call a%PhysicalDriver%GetPhysicalProblem(PhysicalPro)
      call a%PhysicalDriver%GetPhysicalAlgorithm(PhysicalAlgor)
      
      !We will keep a pointer to it safely
      call a%PodRom%SetPhysicalProblem(PhysicalPro,a%PhysicalType,PhysicalAlgor)
     
      if (a%kfl_rstar==1) then 
         if (a%kfl_inter==1) then 
            call a%PodRom%SetInterpolation
            if (c%adaptiveVars%UpdateIO) call a%PodRom%SetUpdateFreq(c%adaptiveVars%UpdateFreq)
            call a%PodRom%SetOldRestartFile(c%masterVars%OldRestartFolder)
            call a%PodRom%SetOldInputFile(c%masterVars%oldnamda)
            call a%PodRom%SetInterp(c%masterVars%Int_Restart)
            call a%PodRom%SetOldMesh
         endif
      endif

      !Turnon operations
      call a%PodRom%Turnon

      !Adaptive Mesh Refinement
      if (c%adaptiveVars%kfl_AdaptiveRefinement == 1 .or. c%adaptiveVars%NumberOfInitialUniformRefinementSteps > 0) then
         call a%PodRom%SetAdaptiveRefiner(c%adaptiveVars%Refiner)
      endif
      
      call a%PodRom%GetItask(kfl_itask)
      if (kfl_itask == 1) then
         a%AreWeRunning = .true.
      else
         a%AreWeRunning = .false.
      endif
     
   end subroutine
      
   subroutine podrom_GetTimeStep(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
  
      if (a%AreWeRunning) return
  
      !Give us the time step, driver!
      call a%PhysicalDriver%GetTimeStep(c)
   end subroutine podrom_GetTimeStep
   
   subroutine podrom_SetTimeStep(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      call a%PodRom%Setste(c%mastervars%dtinv,c%mastervars%ctime)
      
      !Here you have the time step, driver
      call a%PhysicalDriver%SetTimeStep(c)
   end subroutine podrom_SetTimeStep
   
   subroutine podrom_GetRefinementCriteria(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      call a%PodRom%GetRefinementCriteria(c%adaptiveVars%RefinerMarkel)
   end subroutine podrom_GetRefinementCriteria
   
   subroutine podrom_PreRefine(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
  
      if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
         call a%PhysicalDriver%PreRefine(c)
      end if
   end subroutine podrom_PreRefine
   
   subroutine podrom_Refine(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
  
      if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
         call a%PhysicalDriver%Refine(c)
         call a%PodRom%Refine('Refine')
      end if
   end subroutine podrom_Refine
   
   subroutine podrom_Rebalance(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
         call a%PhysicalDriver%Rebalance(c)
         call a%PodRom%Refine('Rebalance')
      end if
   end subroutine podrom_Rebalance
   
   subroutine podrom_Begste(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c

      if (a%ndela > c%masterVars%istep - 1) return
      
      call a%PodRom%Begste
  
      if (a%AreWeRunning) return
 
      !Begste, driver!
      call a%PhysicalDriver%Begste(c)
   end subroutine podrom_Begste
   
   subroutine podrom_MeshProjections(a,c,itask)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
      
      if (a%ndela > c%masterVars%istep - 1) return
      
      call runend('I do not think the rom is ready for mesh projections')
  
      !Mesh projections, driver!
      call a%PhysicalDriver%MeshProjections(c,itask)
   end subroutine podrom_MeshProjections
   
   subroutine podrom_MeshAdvections(a,c,itask)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
      
      if (a%ndela > c%masterVars%istep - 1) return
      
      call runend('I do not think the rom is ready for mesh projections')
  
      !Mesh projections, driver!
      call a%PhysicalDriver%MeshAdvections(c,itask)
   end subroutine podrom_MeshAdvections
   
   subroutine podrom_Doiter(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c

      if (a%ndela > c%masterVars%istep - 1) return
      
      call a%PodRom%Doiter
      
      if (a%AreWeRunning) return
  
      !DoIter, driver!
      call a%PhysicalDriver%DoIter(c)
   end subroutine podrom_Doiter
   
   subroutine podrom_Convergence(a,c,glres)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres

      glres = 0.0_rp
      if (a%ndela > c%masterVars%istep - 1) return
  
      !Compute the convergence, driver!
      call a%PhysicalDriver%Convergence(c,glres)
   end subroutine podrom_Convergence

   subroutine podrom_CoupledConvergence(a,c,kfl_flag,cpres)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      logical             :: kfl_flag
      real(rp)            :: cpres

      kfl_flag = .true.
      if (a%ndela > c%masterVars%istep - 1) return

      call a%PhysicalDriver%CoupledConvergence(c,kfl_flag,cpres)

   end subroutine podrom_CoupledConvergence
   
   subroutine podrom_Endste(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      if (a%ndela > c%masterVars%istep - 1) return
      
      call a%PodRom%Endste(c%masterVars%kfl_gotim)
      
      if (a%AreWeRunning) return
  
      !Do your Enste, driver!
      call a%PhysicalDriver%Endste(c)
   end subroutine podrom_Endste

   subroutine podrom_Output(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      if (a%ndela > c%masterVars%istep - 1) return
      call a%PhysicalDriver%Output(c)

   end subroutine podrom_Output
   
   subroutine podrom_Turnof(a,c)
      implicit none
      class(podromDriver) :: a
      type(caseVariables) :: c
      
      call a%PodRom%Turnof
      
      !deallocate(a%PodRom)
      !I do not deallocate it, I need it to be alive in Cputab
      !We should do something about this
      call a%Memor%deallocobj(0,'Podrom','podrom_turnon',1)
      
      !Turnof, driver!
      call a%PhysicalDriver%Turnof(c)
      call a%Memor%deallocobj(0,'RomPhysicalDriver','podrom_turnon',1)
      
   end subroutine podrom_Turnof
   
   subroutine podrom_WriteTimes(a)
      implicit none
      class(podromDriver) :: a
      
      !Write the times, driver!
      call a%PhysicalDriver%WriteTimes
   end subroutine podrom_WriteTimes
   
end module 

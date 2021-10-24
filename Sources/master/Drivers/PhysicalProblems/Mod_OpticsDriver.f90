module Mod_OpticsDriver
    use typre
    use MPI
    use Mod_BroadCastBuffer
    use Mod_Listen
    use Mod_caseVariables
    use Mod_PhysicalProblemDriver
    use Mod_PhysicalProblem
    use Mod_DistributedContainer 
    use Mod_DC_rp
    use Mod_DC_ip
    use Mod_InChannel
    use Mod_MasterVariables
    use Mod_Optics
   implicit none

   private

   public :: OpticsDriver, OpticsDriver_Const

   type, extends(DriverInterface) :: OpticsDriver 
      
      type(OpticsProblem) :: Optics
      
      !Couplings 
      type(InChannelGroupInfo) :: InNstincInfo
      type(InChannelGroupInfo) :: InTempeInfo


contains

      procedure :: Lev1Reapro             => optics_Lev1Reapro
      procedure :: Lev1ReaproMPI          => optics_Lev1ReaproMPI
      procedure :: SetOutChannels         => optics_SetOutChannels
      procedure :: SetInChannels          => optics_SetInChannels
      procedure :: UpdateChannels         => optics_UpdateChannels
      procedure :: Turnon                 => optics_Turnon
      procedure :: GetTimeStep            => optics_GetTimeStep
      procedure :: SetTimeStep            => optics_SetTimeStep
      procedure :: GetRefinementCriteria  => optics_GetRefinementCriteria
      procedure :: Refine                 => optics_Refine
      procedure :: Rebalance              => optics_Rebalance
      procedure :: Begste                 => optics_Begste
      procedure :: MeshProjections        => optics_MeshProjections
      procedure :: MeshAdvections         => optics_MeshAdvections
      procedure :: Doiter                 => optics_Doiter
      procedure :: Convergence            => optics_Convergence
      procedure :: Endste                 => optics_Endste
      procedure :: Output                 => optics_Output
      procedure :: Turnof                 => optics_Turnof
      procedure :: WriteTimes             => optics_WriteTimes
      
   end type OpticsDriver

   interface OpticsDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(OpticsDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


    subroutine optics_Lev1Reapro(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      !Things that we want from outside
      if(a%Listener%words(1)=='EXTER') then
         !Velocity info
         if(a%Listener%words(2) == 'NSTIN') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
         endif
         !Temperature Info
         if(a%Listener%words(2) == 'TEMPE') then
            call InChannelGroupInfoFromListener(a%Listener,a%InTempeInfo)
         endif
       end if
   end subroutine
   
   subroutine optics_Lev1ReaproMPI(a,c,bb)
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InTempeInfo,bb)
   end subroutine
   
   subroutine optics_SetOutChannels(a,task)
      class(OpticsDriver) :: a
      character(*) :: task
      
      class(DistributedContainer), pointer :: myDC => NULL()
      type(InChannel), pointer :: myInChannel => NULL()
      
      !-----------------------------------------------------------
      !First the OutChannels
      !No out channels
   end subroutine
   
   subroutine optics_SetInChannels(a)
      class(OpticsDriver) :: a
      
      integer(ip) :: kfl_dissi
      
      !--------------------------------------------------------
      !Secondly the InChannels
      
      !Nstinc info to channel list
      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNstincInfo)
      
         call a%AddToInChannelList('PRESS','PRESS','INTERPOLABLE',a%InNstincInfo)
         
         call a%AddToInChannelList('VISCOSITY','VISCOSITY','DO_NOT_INTERPOLATE',a%InNstincInfo)
         
         call a%Optics%GetDissipationModel(kfl_dissi)
         if (kfl_dissi == 1) then
            !NSDissipation
            call a%AddToInChannelList('NSTINC_DISSIPATION','DISSIPATION','INTERPOLABLE',a%InNstincInfo)
         endif
         
      endif
      
      !Tempe info to channel list
      if (a%InTempeInfo%IsOn) then
         call a%AddToInChannelList('TEMPE','TEMPE','INTERPOLABLE',a%InTempeInfo)
         
         call a%AddToInChannelList('DENSITY','DENSITY','DO_NOT_INTERPOLATE',a%InTempeInfo)
         
         call a%AddToInChannelList('SPECIFIC_HEAT','SPECIFIC_HEAT','DO_NOT_INTERPOLATE',a%InTempeInfo)
         
         call a%AddToInChannelList('THERMAL_CONDUCTIVITY','THERMAL_CONDUCTIVITY','DO_NOT_INTERPOLATE',a%InTempeInfo)
      
         call a%AddToInChannelList('TEMPE','TEMPE','INTERPOLABLE',a%InTempeInfo)
      
         call a%Optics%GetDissipationModel(kfl_dissi)
         if (kfl_dissi == 1) then
            !Dissipation
            call a%AddToInChannelList('TEMPE_DISSIPATION','DISSIPATION','INTERPOLABLE',a%InTempeInfo)
            
         endif
         
      endif
      
   end subroutine
   
   subroutine optics_UpdateChannels(a)
      implicit none
      class(OpticsDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      
      real(rp), pointer :: veloc(:,:) => NULL(), press(:) => NULL(), tempe(:) => NULL()
      real(rp), pointer :: acden => NULL(), acvis => NULL(),actco => NULL(),acsph => NULL()
      real(rp), pointer :: veloc_dissipation(:) => NULL(), tempe_dissipation(:) => NULL()
      integer(ip) :: kfl_dissi
      
      if (a%InNstincInfo%IsOn) then
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         call a%Optics%SetVelocityArray(veloc)
         
         call a%GetFromInChannelList('PRESS',myDC)
         call ExtractRP(myDC,press)
         call a%Optics%SetPressureArray(press)
         
         call a%GetFromInChannelList('VISCOSITY',myDC)
         call ExtractRP(myDC,acvis)
         call a%Optics%SetViscosity(acvis)
         
         call a%Optics%GetDissipationModel(kfl_dissi)
         if (kfl_dissi == 1) then
            call a%GetFromInChannelList('NSTINC_DISSIPATION',myDC)
            call ExtractRP(myDC,veloc_dissipation)
            if (.not. associated(veloc_dissipation)) call runend('OpticsDriver: nstincDissipation not found!')
            call a%Optics%SetNSDissipationArray(veloc_dissipation)
         endif   
      endif  
      
      if (a%InTempeInfo%IsOn) then
         call a%GetFromInChannelList('TEMPE',myDC)
         call ExtractRP(myDC,tempe)
         call a%Optics%SetTemperatureArray(tempe)
         
         call a%GetFromInChannelList('DENSITY',myDC)
         call ExtractRP(myDC,acden)
         call a%Optics%SetDensity(acden)
         
         call a%GetFromInChannelList('SPECIFIC_HEAT',myDC)
         call ExtractRP(myDC,acsph)
         call a%Optics%SetSpecificHeat(acsph)
         
         call a%GetFromInChannelList('THERMAL_CONDUCTIVITY',myDC)
         call ExtractRP(myDC,actco)
         call a%Optics%SetThermalConductivity(actco)
         
         call a%Optics%GetDissipationModel(kfl_dissi)
         if (kfl_dissi == 1) then
            call a%GetFromInChannelList('TEMPE_DISSIPATION',myDC)
            call ExtractRP(myDC,tempe_dissipation)
            if (.not. associated(tempe_dissipation)) call runend('OpticsDriver: tempeDissipation not found!')
            call a%Optics%SetTempeDissipationArray(tempe_dissipation)
         endif   
      endif  
      
      
   end subroutine
      
   
   
   subroutine optics_Turnon(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%Optics)
   end subroutine
      
   subroutine optics_GetTimeStep(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%Optics)
  
   end subroutine optics_GetTimeStep
   
   subroutine optics_SetTimeStep(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%Optics)
  
   end subroutine optics_SetTimeStep
   
   subroutine optics_GetRefinementCriteria(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%Optics)
  
   end subroutine optics_GetRefinementCriteria
   
   subroutine optics_Refine(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%Optics)
  
   end subroutine optics_Refine
   
   subroutine optics_Rebalance(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%Optics)
  
   end subroutine optics_Rebalance
   
   subroutine optics_Begste(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%Optics)
  
   end subroutine optics_Begste
   
   subroutine optics_MeshProjections(a,c,itask)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%Optics,itask)
  
   end subroutine optics_MeshProjections
   
   subroutine optics_MeshAdvections(a,c,itask)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%Optics,itask)
  
   end subroutine optics_MeshAdvections
   
   subroutine optics_Doiter(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%Optics)
  
   end subroutine optics_Doiter
   
   subroutine optics_Convergence(a,c,glres)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%Optics,glres)
  
   end subroutine optics_Convergence
   
   subroutine optics_Endste(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%Optics)
  
   end subroutine optics_Endste

   subroutine optics_Output(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%Optics)
  
   end subroutine optics_Output
   
   subroutine optics_Turnof(a,c)
      implicit none
      class(OpticsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%Optics)
  
   end subroutine optics_Turnof
   
   subroutine optics_WriteTimes(a)
      implicit none
      class(OpticsDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%Optics)
  
   end subroutine optics_WriteTimes
   
   
  
  

end module 

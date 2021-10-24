module Mod_TempeDriver
    use typre
    use MPI
    use Mod_BroadCastBuffer
    use Mod_Listen
    use Mod_caseVariables
    use Mod_PhysicalProblemDriver
    use Mod_PhysicalProblem
    use Mod_DistributedContainer 
    use Mod_DC_ip
    use Mod_DC_rp
    use Mod_InChannel
    use Mod_MasterVariables
    use Mod_Temperature
   implicit none

   private

   public :: TempeDriver, TempeDriver_Const

   type, extends(DriverInterface) :: tempeDriver 
      
      type(TemperatureProblem) :: TempeProblem
      
      !Couplings 
      type(InChannelGroupInfo) :: InNstincInfo
      real(rp) :: acden,acsph,actco


contains

      procedure :: Lev1Reapro             => tempe_Lev1Reapro
      procedure :: Lev1ReaproMPI          => tempe_Lev1ReaproMPI
      procedure :: SetOutChannels         => tempe_SetOutChannels
      procedure :: SetInChannels          => tempe_SetInChannels
      procedure :: UpdateChannels         => tempe_UpdateChannels
      procedure :: Turnon                 => tempe_Turnon
      procedure :: GetTimeStep            => tempe_GetTimeStep
      procedure :: SetTimeStep            => tempe_SetTimeStep
      procedure :: GetRefinementCriteria  => tempe_GetRefinementCriteria
      procedure :: PreRefine              => tempe_PreRefine
      procedure :: Refine                 => tempe_Refine
      procedure :: Rebalance              => tempe_Rebalance
      procedure :: Begste                 => tempe_Begste
      procedure :: MeshProjections        => tempe_MeshProjections
      procedure :: MeshAdvections         => tempe_MeshAdvections
      procedure :: Doiter                 => tempe_Doiter
      procedure :: Convergence            => tempe_Convergence
      procedure :: Endste                 => tempe_Endste
      procedure :: Output                 => tempe_Output
      procedure :: Turnof                 => tempe_Turnof
      procedure :: WriteTimes             => tempe_WriteTimes
      procedure :: GetPhysicalProblem     => tempe_GetPhysicalProblem
      
   end type TempeDriver

   interface TempeDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(TempeDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine tempe_Lev1Reapro(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      !Things that we want from outside
      if(a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'NSTIN') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
         endif
      end if   
   end subroutine
   
   subroutine tempe_Lev1ReaproMPI(a,c,bb)
      class(TempeDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
   end subroutine
   
   subroutine tempe_SetOutChannels(a,task)
      class(tempeDriver) :: a
      character(*) :: task
      
      class(DistributedContainer), pointer :: myDC => NULL()
      type(InChannel), pointer :: myInChannel => NULL()
      
      real(rp), pointer :: tempe(:) => NULL()
      
      type(r2p), pointer :: tesgs(:) => NULL()
      real(rp) :: acden, acsph , actco , acrea , acsou 
      real(rp), pointer :: dissipation(:) => NULL()

      !-----------------------------------------------------------
      !First the OutChannels
         
      !Temperature
      call a%TempeProblem%GetTemperatureArray(tempe)
      call RPToList(tempe,'TEMPE',task,a%OutChannelList)
      
      !Temperature Subgrid Scales
      call a%TempeProblem%GetTemperatureSGS(tesgs)
      if(associated(tesgs)) call RPToList(tesgs,'TEMPE_SGS',task,a%OutChannelList)
      
      call a%TempeProblem%GetPhysicalParameters(1,a%acden,a%acsph,a%actco,acrea,acsou)
      call RPToList(a%acden,'DENSITY',task,a%OutChannelList)
      
      call RPToList(a%acsph,'SPECIFIC_HEAT',task,a%OutChannelList)
      
      call RPToList(a%actco,'THERMAL_CONDUCTIVITY',task,a%OutChannelList)
      
      !TempeDissipation
      call a%TempeProblem%GetDissipationArray(Dissipation)
      if(associated(Dissipation)) call RPToList(Dissipation,'DISSIPATION',task,a%OutChannelList)
   end subroutine
   
   subroutine tempe_SetInChannels(a)
      class(TempeDriver) :: a
      !--------------------------------------------------------
      !Secondly the InChannels
      
      !Veloc info to channel list
      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNstincInfo)
         if (a%TempeProblem%kfl_CouplingThreeField==1) call a%AddToInChannelList('SIGMA','SIGMA','INTERPOLABLE',a%InNstincInfo)
         
         call a%AddToInChannelList('VELOC_SGS','VELOC_SGS','SKIP_IF_INTERPOLATE',a%InNstincInfo,'OPTIONAL')
      endif
      
   end subroutine
      
   subroutine tempe_UpdateChannels(a)
      implicit none
      class(tempeDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      
      real(rp),    pointer :: veloc(:,:) => NULL()
      type(r3p),   pointer :: vesgs(:)   => NULL()
      real(rp),    pointer :: sigma(:,:) => NULL()
      integer(ip), pointer :: kfl_trasg  => NULL()
      
      if (a%InNstincInfo%IsOn) then
         
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         call a%TempeProblem%SetVelocityArray(veloc)
         
         if (a%TempeProblem%kfl_CouplingThreeField==1) then
            call a%GetFromInChannelList('SIGMA',myDC)
            call ExtractRP(myDC,sigma)
            if (associated(sigma)) call a%TempeProblem%SetSigmaNSArray(sigma)
         endif 

         call a%GetFromInChannelList('VELOC_SGS',myDC)
         call ExtractRP(myDC,vesgs)
         if (associated(vesgs)) call a%TempeProblem%SetVelocitySGS(vesgs)

      endif   
      
      
   end subroutine
   
   subroutine tempe_Turnon(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%TempeProblem)
   end subroutine
      
   subroutine tempe_GetTimeStep(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%TempeProblem)
  
   end subroutine tempe_GetTimeStep
   
   subroutine tempe_SetTimeStep(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%TempeProblem)
  
   end subroutine tempe_SetTimeStep
   
   subroutine tempe_GetRefinementCriteria(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%TempeProblem)
  
   end subroutine tempe_GetRefinementCriteria
   
   subroutine tempe_PreRefine(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%TempeProblem)
  
   end subroutine tempe_PreRefine
   
   subroutine tempe_Refine(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%TempeProblem)
  
   end subroutine tempe_Refine
   
   subroutine tempe_Rebalance(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%TempeProblem)
  
   end subroutine tempe_Rebalance
   
   subroutine tempe_Begste(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%TempeProblem)
  
   end subroutine tempe_Begste
   
   subroutine tempe_MeshAdvections(a,c,itask)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%TempeProblem,itask)
  
   end subroutine tempe_MeshAdvections
   
   subroutine tempe_MeshProjections(a,c,itask)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%TempeProblem,itask)
  
   end subroutine tempe_MeshProjections
   
   subroutine tempe_Doiter(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%TempeProblem)
  
   end subroutine tempe_Doiter
   
   subroutine tempe_Convergence(a,c,glres)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%TempeProblem,glres)
  
   end subroutine tempe_Convergence
   
   subroutine tempe_Endste(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%TempeProblem)
  
   end subroutine tempe_Endste

   subroutine tempe_Output(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%TempeProblem)
  
   end subroutine tempe_Output
   
   subroutine tempe_Turnof(a,c)
      implicit none
      class(TempeDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%TempeProblem)
  
   end subroutine tempe_Turnof
   
   subroutine tempe_WriteTimes(a)
      implicit none
      class(TempeDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%TempeProblem)
  
   end subroutine tempe_WriteTimes
   
   subroutine tempe_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(TempeDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%TempeProblem
   end subroutine
   
end module 

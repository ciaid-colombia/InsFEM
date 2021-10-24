module Mod_ALEDriver
    use typre
    use MPI
    use Mod_BroadCastBuffer
    use Mod_Listen
    use Mod_caseVariables
    use Mod_PhysicalProblemDriver
    use Mod_PhysicalProblem
    use Mod_DistributedContainer 
    use Mod_DC_rp
    use Mod_InChannel
    use Mod_MasterVariables
    use Mod_Alemov
    use Mod_sldAlemov
    use Mod_CutMesh
    use Mod_DC_CutMesh
   implicit none

   private

   public :: ALEDriver, ALEDriver_Const
   
   class(AlemovProblem), pointer :: ALEProblemDebug => NULL()
   

   type, extends(DriverInterface) :: ALEDriver 
      
      class(AlemovProblem), pointer :: ALEProblem => NULL()
      
      integer(ip) :: kfl_fixedMeshALE

      !Couplings 
      type(InChannelGroupInfo) :: InNstincInfo
      type(InChannelGroupInfo) :: InSolidsInfo
      type(InChannelGroupInfo) :: InNscompInfo
      type(InChannelGroupInfo) :: InCutMeshInfo
      type(InChannelGroupInfo) :: InPLCdInfo


contains

      procedure :: Lev1Reapro             => aledri_Lev1Reapro
      procedure :: Lev1ReaproMPI          => aledri_Lev1ReaproMPI
      procedure :: SetInChannels          => aledri_SetInChannels
      procedure :: SetOutChannels         => aledri_SetOutChannels
      procedure :: UpdateChannels         => aledri_UpdateChannels
      procedure :: Turnon                 => aledri_Turnon
      procedure :: GetTimeStep            => aledri_GetTimeStep
      procedure :: SetTimeStep            => aledri_SetTimeStep
      procedure :: GetRefinementCriteria  => aledri_GetRefinementCriteria
      procedure :: Refine                 => aledri_Refine
      procedure :: Rebalance              => aledri_Rebalance
      procedure :: Begste                 => aledri_Begste
      procedure :: MeshProjections        => aledri_MeshProjections
      procedure :: MeshAdvections         => aledri_MeshAdvections
      procedure :: Doiter                 => aledri_Doiter
      procedure :: Convergence            => aledri_Convergence
      procedure :: CoupledConvergence     => aledri_CoupledConvergence
      procedure :: Endste                 => aledri_Endste
      procedure :: Output                 => aledri_Output
      procedure :: Turnof                 => aledri_Turnof
      procedure :: WriteTimes             => aledri_WriteTimes
      procedure :: GetRemeshingCriteria   => aledri_GetRemeshingCriteria
      procedure :: GetPhysicalProblem     => aledri_GetPhysicalProblem
      
   end type ALEDriver

   interface ALEDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(ALEDriver), pointer :: constructor

       allocate(constructor)
   end function constructor
   
   subroutine ALEPointerDebug(a)
      implicit none
      class(ALEDriver), target :: a
      
      ALEProblemDebug => a%ALEProblem
   end subroutine

   subroutine aledri_Lev1Reapro(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      if(a%Listener%words(1).eq.'ALGOR') then 
         if(a%Listener%exists('SOLID')) then
             a%kfl_algor= 2
         end if
      end if   
      !Things that we want from outside
      if (a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'NSTIN') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
         elseif(a%Listener%words(2) == 'SOLID') then
            call InChannelGroupInfoFromListener(a%Listener,a%InSolidsInfo)
         elseif(a%Listener%words(2) == 'NSCOM') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNscompInfo)
         elseif(a%Listener%words(2) == 'PLCDP') then
            call InChannelGroupInfoFromListener(a%Listener,a%InPLCdInfo)
         elseif (a%Listener%words(2) == 'LEVEL' .or. a%Listener%words(2) == 'LEVSE')then
            call InChannelGroupInfoFromListener(a%Listener,a%InCutMeshInfo)
            a%InCutMeshInfo%DriverKey = 'LEVSE'
            a%InCutMeshInfo%DriverNick = 'LEVSE'
         end if
      endif
   end subroutine
   
   subroutine aledri_Lev1ReaproMPI(a,c,bb)
      class(ALEDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call bb%Add(a%kfl_algor)
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InSolidsInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InNscompInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InCutMeshInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InPLCdInfo,bb)
   end subroutine
   
   subroutine aledri_SetOutChannels(a,task)
      class(ALEDriver) :: a
      character(*) :: task
      
      real(rp), pointer :: veloc(:,:,:) => NULL()
      
      !Veloc
      call a%AleProblem%GetVelocityArray(veloc)
      call RPToList(veloc(:,:,1),'VELOC',task,a%OutChannelList)

   end subroutine
   
   subroutine aledri_SetInChannels(a)
      class(ALEDriver) :: a
      
      !--------------------------------------------------------
      !Secondly the InChannels
      
      !Veloc info to channel list
      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNstincInfo)
         
      endif

      if (a%InSolidsInfo%IsOn) then
         call a%AddToInChannelList('DISPL','DISPL','INTERPOLABLE',a%InSolidsInfo)
         
      endif

      if (a%InNscompInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNscompInfo)
         
      endif

      if (a%InCutMeshInfo%IsOn) then
         call a%AddToInChannelList('CUT_MESH','CUT_MESH','RUNEND_IF_INTERPOLATE',a%InCutMeshInfo)
      endif
      
      if (a%InPLCdInfo%IsOn) then
         call a%AddToInChannelList('DISPL','DISPLACEMENT','INTERPOLABLE',a%InPLCdInfo)
      endif
      
   end subroutine
      
   subroutine aledri_UpdateChannels(a)
      implicit none
      class(ALEDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      type(CutMesh), pointer :: CMesh => NULL()
      real(rp), pointer :: veloc(:,:) => NULL(), displ(:,:) => NULL()
      
      if (a%InNstincInfo%IsOn) then
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         call a%AleProblem%SetVelocityArray(veloc)
      endif   

      if (a%InSolidsInfo%IsOn) then
         call a%GetFromInChannelList('DISPL',myDC)
         call ExtractRP(myDC,displ)
         call a%AleProblem%SetDisplacementArray(displ)
      endif   

      if (a%InNscompInfo%IsOn) then
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         call a%AleProblem%SetVelocityArray(veloc)
      endif   
      
      if (a%InCutMeshInfo%IsOn) then
         call a%GetFromInChannelList('CUT_MESH',myDC)
         call ExtractCutMesh(myDC,CMesh)
         call a%AleProblem%SetCutMesh(CMesh)
      endif   
      
      if (a%InPLCdInfo%IsOn) then
         call a%GetFromInChannelList('DISPL',myDC)
         call ExtractRP(myDC,displ)
         call a%AleProblem%SetDisplacementArray(displ)
      endif  
      
   end subroutine
   
   subroutine aledri_Turnon(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      
      if (a%kfl_algor== 1) then
         a%ALEProblem=> AlemovProblem_Const()
      elseif (a%kfl_algor== 2) then
         a%ALEProblem=> sldAlemovProblem_Const()
      else
         call runend('ale_driver_reapro: wrong algorithm for ALE')
      endif
      
      call ALEPointerDebug(a)
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%ALEProblem)
   end subroutine
      
   subroutine aledri_GetTimeStep(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%ALEProblem)
  
   end subroutine aledri_GetTimeStep
   
   subroutine aledri_SetTimeStep(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%ALEProblem)
  
   end subroutine aledri_SetTimeStep
   
   subroutine aledri_GetRefinementCriteria(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%ALEProblem)
  
   end subroutine aledri_GetRefinementCriteria
   
   subroutine aledri_Refine(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%ALEProblem)
  
   end subroutine aledri_Refine
   
   subroutine aledri_Rebalance(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%ALEProblem)
  
   end subroutine aledri_Rebalance
   
   subroutine aledri_Begste(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%ALEProblem)

   end subroutine aledri_Begste
   
   subroutine aledri_MeshProjections(a,c,itask)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%ALEProblem,itask)
  
   end subroutine aledri_MeshProjections
   
   subroutine aledri_MeshAdvections(a,c,itask)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%ALEProblem,itask)
  
   end subroutine aledri_MeshAdvections
   
   subroutine aledri_Doiter(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%ALEProblem)
  
   end subroutine aledri_Doiter
   
   subroutine aledri_Convergence(a,c,glres)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%ALEProblem,glres)
  
   end subroutine aledri_Convergence

   subroutine aledri_CoupledConvergence(a,c,kfl_flag,cpres)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
      logical             :: kfl_flag
      real(rp):: cpres
  
      !PhysicalProblem operations
      call physical_CoupledConvergence(a,c,a%ALEProblem,kfl_flag,cpres)
   end subroutine aledri_CoupledConvergence
   
   subroutine aledri_Endste(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%ALEProblem)

      !So next steps are delayed 
      !(but the output displacements are not smoothly updated) 
      !if (a%ndela > c%masterVars%istep - 1) return
      !if (a%consec_delay > 0) a%ndela = a%ndela + a%consec_delay

   end subroutine aledri_Endste

   subroutine aledri_Output(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%ALEProblem)
   end subroutine aledri_Output
   
   subroutine aledri_Turnof(a,c)
      implicit none
      class(ALEDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%ALEProblem)
  
   end subroutine aledri_Turnof
   
   subroutine aledri_WriteTimes(a)
      implicit none
      class(ALEDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%ALEProblem)
  
   end subroutine aledri_WriteTimes
   
   subroutine aledri_GetRemeshingCriteria(a,DoRemesh)
      class(ALEDriver) :: a
      integer(ip) :: DoRemesh
      
      call a%AleProblem%GetDoRemesh(DoRemesh)
   end subroutine
   
   subroutine aledri_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(ALEDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%AleProblem
   end subroutine

end module 

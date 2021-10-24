module Mod_LevelDriver
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
    use Mod_LevelSet
    use Mod_CutMesh
    use Mod_DC_CutMesh
   implicit none

   private

   public :: LevelDriver, LevelDriver_Const

   type(LevelSetProblem), pointer :: LevelSetProblemPointer => NULL()
   
   type, extends(DriverInterface) :: LevelDriver 
      
      type(LevelSetProblem) :: Level
      
      !Couplings 
      type(InChannelGroupInfo) :: InNstincInfo
      type(InChannelGroupInfo) :: InLevseInfo

contains

      procedure :: Lev1Reapro             => Level_Lev1Reapro
      procedure :: Lev1ReaproMPI          => Level_Lev1ReaproMPI
      procedure :: SetInChannels          => Level_SetInChannels
      procedure :: SetOutChannels         => Level_SetOUtChannels
      procedure :: UpdateChannels         => Level_UpdateChannels
      procedure :: Turnon                 => Level_Turnon
      procedure :: GetTimeStep            => Level_GetTimeStep
      procedure :: SetTimeStep            => Level_SetTimeStep
      procedure :: GetRefinementCriteria  => Level_GetRefinementCriteria
      procedure :: Refine                 => Level_Refine
      procedure :: Rebalance              => Level_Rebalance
      procedure :: Begste                 => Level_Begste
      procedure :: MeshProjections        => Level_MeshProjections
      procedure :: MeshAdvections         => Level_MeshAdvections
      procedure :: Doiter                 => Level_Doiter
      procedure :: Convergence            => Level_Convergence
      procedure :: Endste                 => Level_Endste
      procedure :: Output                 => Level_Output
      procedure :: Turnof                 => Level_Turnof
      procedure :: WriteTimes             => Level_WriteTimes
      
   end type LevelDriver

   interface LevelDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(LevelDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine Level_Lev1Reapro(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      !Things that we want from outside
      if(a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'NSTIN' .or. a%Listener%words(2) == 'VELOC') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
         endif
         if(a%Listener%words(2) == 'LEVSE' .or. a%Listener%words(2) == 'LEVEL') then
            call InChannelGroupInfoFromListener(a%Listener,a%InLevseInfo)
         endif
      end if  
   end subroutine
   
   subroutine Level_Lev1ReaproMPI(a,c,bb)
      class(LevelDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InLevseInfo,bb)
   end subroutine
   
   subroutine Level_SetOutChannels(a,task)
      class(LevelDriver) :: a
      character(*) :: task
      
      class(DistributedContainer), pointer :: myDC => NULL()
      type(CutMesh), pointer :: CMesh => NULL()
      real(rp), pointer      :: SGradient(:,:)=> NULL(), level(:) => NULL()
      
      
      !-----------------------------------------------------------
      !First the OutChannels
      call a%Level%GetCutMesh(CMesh)
      if(associated(CMesh)) call CutMeshToList(CMesh,'CUT_MESH',task,a%OutChannelList)
      
      call a%Level%GetSmoothGradient(SGradient)
      if(associated(SGradient)) call RPToList(SGradient,'CUT_GRADIENT',task,a%OutChannelList)
      
      call a%Level%GetLevelsetArray(level)
      if(associated(level)) call RPToList(level,'LEVEL',task,a%OutChannelList)
   end subroutine
   
   subroutine Level_SetInChannels(a)
      class(LevelDriver) :: a
      
      !--------------------------------------------------------
      !Secondly the InChannels
      
      !Veloc info to channel list
      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNstincInfo)
      endif
      if (a%InLevseInfo%IsOn) then
         call a%AddToInChannelList('LEVEL','LEVEL','INTERPOLABLE',a%InLevseInfo)
      endif
   end subroutine
      
   subroutine Level_UpdateChannels(a)
      implicit none
      class(LevelDriver) :: a
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      real(rp), pointer :: veloc(:,:) => NULL(), levse(:) => NULL()
      
      if (a%InNstincInfo%IsOn) then
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         call a%Level%SetVelocityArray(veloc)
      endif 
      
      if (a%InLevseInfo%IsOn) then
         call a%Level%GetLevelsetArray(levse)
         call a%SetDefaultInterpolatedArray('LEVEL',levse)
         call a%GetFromInChannelList('LEVEL',myDC)
         call ExtractRP(myDC,levse)
         call a%Level%SetUnkno(levse,1)
         call a%Level%SkipSystemSolve(.true.)
         !call a%Level%FilePostpr%postpr(a%Level%level(:,1),'levse',a%Level%istep,a%Level%ctime,a%Level%Mesh)
      endif 
      
   end subroutine
   
   subroutine Level_Turnon(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
      
      call PointPointer(a) 
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%Level)
   end subroutine
   
   subroutine PointPointer(a)
      class(LevelDriver), target :: a
      
      LevelSetProblemPointer => a%Level
   end subroutine
      
   subroutine Level_GetTimeStep(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%Level)
   end subroutine Level_GetTimeStep
   
   subroutine Level_SetTimeStep(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%Level)
   end subroutine Level_SetTimeStep
   
   subroutine Level_GetRefinementCriteria(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%Level)
   end subroutine Level_GetRefinementCriteria
   
   subroutine Level_Refine(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%Level)
   end subroutine Level_Refine
   
   subroutine Level_Rebalance(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%Level)
   end subroutine Level_Rebalance
   
   subroutine Level_Begste(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%Level)
   end subroutine Level_Begste
   
   subroutine Level_MeshProjections(a,c,itask)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%Level,itask)
   end subroutine Level_MeshProjections
   
   subroutine Level_MeshAdvections(a,c,itask)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%Level,itask)
   end subroutine Level_MeshAdvections
   
   subroutine Level_Doiter(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%Level)
   end subroutine Level_Doiter
   
   subroutine Level_Convergence(a,c,glres)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%Level,glres)
   end subroutine Level_Convergence
   
   subroutine Level_Endste(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%Level)
   end subroutine Level_Endste

   subroutine Level_Output(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%Level)
   end subroutine Level_Output
   
   subroutine Level_Turnof(a,c)
      implicit none
      class(LevelDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%Level)
   end subroutine Level_Turnof
   
   subroutine Level_WriteTimes(a)
      implicit none
      class(LevelDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%Level)
   end subroutine Level_WriteTimes
   
     

end module 

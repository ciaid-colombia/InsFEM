module Mod_NstincDriver
   use typre
   use MPI
   use Mod_BroadCastBuffer
   use Mod_Listen
   use Mod_caseVariables
   use Mod_DriverInterface
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   use Mod_NSFractionalStep
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSTINC_SUP) || !defined MODULE_SELECTION_ON   
   use Mod_ThreeField
   use Mod_SUPFractionalStep
#endif
   use Mod_CutMesh
   use Mod_DistributedContainer 
   use Mod_DC_rp
   use Mod_DC_ip
   use Mod_DC_CutMesh
   use Mod_InChannel
   use Mod_MasterVariables
   use Mod_PhysicalProblemDriver 
   implicit none

   private

   public :: NstincDriver, NstincDriver_Const

   type, extends(DriverInterface) :: NstincDriver 
   
      class(NavierStokesProblem), pointer :: NavierStokes => NULL() 

      !Modify viscosity at a certain time step
      integer(ip) :: kfl_ModifyViscosity = 0
      integer(ip) :: ModifiedViscosityStep
      real(rp)    :: ModifiedViscosity

      !Couplings 
      type(InChannelGroupInfo) :: InTempeInfo, InLevelInfo, InSolidInfo, InDirichletInfo, InPLCdInfo, InNstincInfo

contains

      procedure :: Lev1Reapro             => nstinc_Lev1Reapro
      procedure :: Lev1ReaproMPI          => nstinc_Lev1ReaproMPI
      procedure :: SetOutChannels         => nstinc_SetOutChannels
      procedure :: SetInChannels          => nstinc_SetInChannels
      procedure :: UpdateChannels         => nstinc_UpdateChannels
      procedure :: Turnon                 => nstinc_Turnon
      procedure :: GetTimeStep            => nstinc_GetTimeStep
      procedure :: SetTimeStep            => nstinc_SetTimeStep
      procedure :: GetRefinementCriteria  => nstinc_GetRefinementCriteria
      procedure :: PreRefine              => nstinc_PreRefine
      procedure :: Refine                 => nstinc_Refine
      procedure :: Rebalance              => nstinc_Rebalance
      procedure :: Begste                 => nstinc_Begste
      procedure :: MeshProjections        => nstinc_MeshProjections
      procedure :: MeshAdvections         => nstinc_MeshAdvections
      procedure :: Doiter                 => nstinc_Doiter
      procedure :: Convergence            => nstinc_Convergence
      procedure :: CoupledConvergence     => nstinc_CoupledConvergence
      procedure :: Endste                 => nstinc_Endste
      procedure :: Output                 => nstinc_Output
      procedure :: Turnof                 => nstinc_Turnof
      procedure :: WriteTimes             => nstinc_WriteTimes
      procedure :: GetPhysicalProblem     => nstinc_GetPhysicalProblem
      
   end type NstincDriver

   interface NstincDriver_Const
      procedure constructor
   end interface 
   

contains

   function constructor()
       class(NstincDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine nstinc_Lev1Reapro(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      if(a%Listener%words(1).eq.'ALGOR') then 
         if(a%Listener%exists('MONOL')) then
             a%kfl_algor= 1
         end if
         if(a%Listener%exists('FRACT')) then
             a%kfl_algor= 2
         end if
         if(a%Listener%exists('THREE')) then  ! THREE_FIELD_CASE
             a%kfl_algor= 3
         end if
         if(a%Listener%exists('FVISC')) then   ! Pressure correction viscoelastic case
             a%kfl_algor= 4
         end if
      
      !Things that we want from outside
      else if(a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'TEMPE') then
            call InChannelGroupInfoFromListener(a%Listener,a%InTempeInfo)
         elseif (a%Listener%words(2) == 'LEVEL' .or. a%Listener%words(2) == 'LEVSE')then
            call InChannelGroupInfoFromListener(a%Listener,a%InLevelInfo)
            a%InLevelInfo%DriverKey = 'LEVSE'
            a%InLevelInfo%DriverNick= 'LEVSE'
         elseif (a%Listener%words(2) == 'SOLID') then
            call InChannelGroupInfoFromListener(a%Listener,a%InSolidInfo)
         elseif(a%Listener%words(2) == 'DIRVE') then
            call InChannelGroupInfoFromListener(a%Listener,a%InDirichletInfo)
         elseif(a%Listener%words(2) == 'PLCDP') then
            call InChannelGroupInfoFromListener(a%Listener,a%InPLCdInfo)
         elseif(a%Listener%words(2) == 'NSTIN') then
            call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
         endif
      else if(a%Listener%words(1)=='MODIF') then
         if(a%Listener%exists('VISCO')) then
             a%kfl_modifyViscosity = 1   !We are going to modify viscosity at some step
             a%ModifiedViscosity = a%Listener%Getrea('VALUE',1.0_rp,'*VISCOSITY VALUES')
             a%ModifiedViscosityStep = a%Listener%Getint('STEP ',1_ip,'*STEP WHEN THE MODIFICATION OCCURS')
         endif
      end if   
   end subroutine
   
   subroutine nstinc_Lev1ReaproMPI(a,c,bb)
      class(NstincDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call bb%Add(a%kfl_ModifyViscosity)
      call bb%Add(a%ModifiedViscosityStep)
      call bb%Add(a%ModifiedViscosity)
      call bb%Add(a%kfl_algor)
      
      call InChannelGroupInfoAddToBuffer(a%InTempeInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InLevelInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InSolidInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InDirichletInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InPLCdInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
   end subroutine
   
   subroutine nstinc_SetOutChannels(a, task)
      class(NstincDriver) :: a
      character(*) :: task
      
      class(DistributedContainer), pointer :: myDC => NULL()
      class(NavierStokesProblem), pointer :: NS => NULL()
      type(InChannel), pointer :: myInChannel => NULL()
      
      real(rp), pointer :: veloc(:,:) => NULL()
      real(rp), pointer :: press(:) => NULL()
      real(rp), pointer :: diss(:) => NULL()
      real(rp), pointer :: density => NULL(), viscosity => NULL()
      real(rp), pointer :: sigma(:,:) => NULL()

      integer(ip) :: kfl_cosig
      type(r3p), pointer :: vesgs(:) => NULL()
      real(rp), pointer :: trac(:,:) => NULL(), btau(:) => NULL()
      real(rp), pointer :: source(:,:) => NULL(), RHSout(:) => NULL()

      !-------------------------------------------------------------------
      !First we set the OutChannels   
      !A getter should be implemented
      call IPToList(a%NavierStokes%kfl_trasg,'KFL_TRASG',task,a%OutChannelList)
      
      !Veloc
      call a%NavierStokes%GetVelocityArray(veloc)
      if(associated(veloc)) call RPToList(veloc,'VELOC',task,a%OutChannelList)
      
      !Velocity Subgrid Scales
      call a%NavierStokes%GetVelocitySGS(vesgs)
      if(associated(vesgs)) call RPToList(vesgs,'VELOC_SGS',task,a%OutChannelList)
      
      !Press
      call a%NavierStokes%GetPressureArray(press)
      if(associated(press)) call RPToList(press,'PRESS',task,a%OutChannelList)

      !Sigma
      call a%NavierStokes%GetSigmaArray(sigma)
      if(associated(sigma)) call RPToList(sigma,'SIGMA',task,a%OutChannelList)

      !Dissipation
      call a%NavierStokes%GetDissipationArray(diss)
      if(associated(diss)) call RPToList(diss,'DISSIPATION',task,a%OutChannelList)
      
      !Density and Viscosity
      !At this moment only for the first material (only ready for one material)
      call a%NavierStokes%GetPhysicalParametersPointers(1,density,viscosity)
      
      call RPToList(density,'DENSITY',task,a%OutChannelList)
      
      call RPToList(viscosity,'VISCOSITY',task,a%OutChannelList)
      
      !Boundary coupling
      !Tractions
      call a%NavierStokes%GetTractionArray(trac)
      if(associated(trac)) call RPToList(trac,'TRACTION',task,a%OutChannelList)
      !Boundary stabilization
      call a%NavierStokes%GetBoundaryTau(btau) 
      if (associated(btau)) call RPToList(btau,'BOUTAU',task,a%OutChannelList)
      
      !Fixno
      !A getter should be implemented
      call IPToList(a%NavierStokes%kfl_fixno,'KFL_FIXNO',task,a%OutChannelList)
   end subroutine
   
   subroutine nstinc_SetInChannels(a)
      class(NstincDriver) :: a

      !-----------------------------------------------------------------------
      !Secondly we set the InChannels
      
      !Tempe info to channel list
      if (a%InTempeInfo%IsOn) then
         call a%AddToInChannelList('TEMPE','TEMPE','INTERPOLABLE',a%InTempeInfo)
         
         if (a%NavierStokes%kfl_trasg/=0) call a%AddToInChannelList('TEMPE_SGS','TEMPE_SGS','SKIP_IF_INTERPOLATE',a%InTempeInfo)
      endif
      
      !Level info to channel
      !For this one, I will need to build multiple channels
      if (a%InLevelInfo%IsOn) then
         call a%AddToInChannelList('CUT_MESH','CUT_MESH','RUNEND_IF_INTERPOLATE',a%InLevelInfo)
         
         call a%AddToInChannelList('CUT_GRADIENT','CUT_GRADIENT','RUNEND_IF_INTERPOLATE',a%InLevelInfo)
         
!          call a%AddToInChannelList('LEVEL','LEVEL','INTERPOLABLE',a%InLevelInfo)
      endif
      
      if (a%InSolidInfo%IsOn) then
         call a%AddToInChannelList('INTTRAC','INTTRAC','INTERPOLABLE',a%InSolidInfo)
         call a%NavierStokes%SetTractionFlag
      endif  
      
      if (a%InDirichletInfo%IsOn) then
         call a%AddToInChannelList('DIRVE','VELOC','INTERPOLABLE',a%InDirichletInfo)
      endif      
      
      if (a%InPLCdInfo%IsOn) then
         call a%AddToInChannelList('INTTRAC','TRACTION','INTERPOLABLE',a%InPLCdInfo)
      endif
   
      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('VELOC','VELOC','INTERPOLABLE',a%InNstincInfo)
         !call a%AddToInChannelList('PRESS','PRESS','INTERPOLABLE',a%InNstincInfo)
         call a%AddToInChannelList('TRACTION','TRACTION','INTERPOLABLE',a%InNstincInfo)
         if (a%NavierStokes%kfl_bousg == 1) then 
            call a%AddToInChannelList('BOUTAU','BOUTAU','INTERPOLABLE',a%InNstincInfo)
         end if
         call a%NavierStokes%SetTractionFlag
      endif  
      
   end subroutine
      
   subroutine nstinc_UpdateChannels(a)
      class(NstincDriver) :: a
      
      class(DistributedContainer), pointer :: myDC => NULL()
      
      real(rp), pointer :: tempe(:) => NULL()
      type(r2p), pointer :: tesgs(:) => NULL()
      
      type(CutMesh), pointer :: CMesh => NULL()
      real(rp), pointer :: press(:) => NULL(), veloc(:,:) => NULL()
      real(rp), pointer :: Sgradient(:,:) => NULL(), level(:) => NULL()
      real(rp), pointer :: trac(:,:) => NULL(), btau(:) => NULL()
   
      !Tempe info from channel list
      if (a%InTempeInfo%IsOn) then
         
         call a%GetFromInChannelList('TEMPE',myDC)
         call ExtractRP(myDC,tempe)
         call a%NavierStokes%SetTemperatureArray(tempe)

         if (a%NavierStokes%kfl_trasg/=0) then
             call a%GetFromInChannelList('TEMPE_SGS',myDC)
             call ExtractRP(myDC,tesgs)
             if (associated(tesgs)) call a%NavierStokes%SetTemperatureSGS(tesgs)
         endif
      endif
      
      !Level
      if (a%InLevelInfo%IsOn) then
         call a%GetFromInChannelList('CUT_MESH',myDC)
         call ExtractCutMesh(myDC,CMesh)
         if(associated(CMesh)) call a%NavierStokes%SetCutMesh(CMesh)
         
         call a%GetFromInChannelList('CUT_GRADIENT',myDC)
         call ExtractRP(myDC,SGradient)
         if(associated(SGradient)) call a%NavierStokes%SetCutGradient(SGradient)
      endif
      
      if (a%InDirichletInfo%IsOn) then
         call a%GetFromInChannelList('DIRVE',myDC)
         call ExtractRP(myDC,veloc)
         if(associated(veloc)) call a%NavierStokes%SetDirichletVelocity(veloc)
      endif
      
      !Solid, for now done by means of ALE
      if (a%InSolidInfo%IsOn) then
         call a%GetFromInChannelList('INTTRAC',myDC)
         call ExtractRP(myDC,trac)
         if(associated(trac)) call a%NavierStokes%SetExternalTractionArray(trac)
      endif
      
      if (a%InPLCdInfo%IsOn) then
         call a%GetFromInChannelList('INTTRAC',myDC)
         call ExtractRP(myDC,trac)
         if(associated(trac)) call a%NavierStokes%SetExternalTractionArray(trac)
      endif
      
      !Another nstinc, boundary stabilization
      if (a%InNstincInfo%IsOn) then
         call a%GetFromInChannelList('VELOC',myDC)
         call ExtractRP(myDC,veloc)
         if (associated(veloc)) call a%NavierStokes%SetVelocityArray(veloc)
         !call a%GetFromInChannelList('PRESS',myDC)
         !call ExtractRP(myDC,press)
         !if (associated(press)) call a%NavierStokes%SetPressureArray(press)
         call a%GetFromInChannelList('TRACTION',myDC)
         call ExtractRP(myDC,trac)
         if (associated(trac)) call a%NavierStokes%SetExternalTractionArray(trac)
         if (a%NavierStokes%kfl_bousg == 1) then 
            call a%GetFromInChannelList('BOUTAU',myDC)
            call ExtractRP(myDC,btau)
            if (associated(btau)) call a%NavierStokes%SetBoundaryTau(btau)
         endif
      endif
      
   end subroutine
   
   subroutine nstinc_Turnon(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      
      !After broadcast operations
      if (a%kfl_algor== 1) then
         a%NavierStokes => NavierStokesProblem_Const()
      elseif (a%kfl_algor== 2) then
         a%NavierStokes => NSFractionalStepProblem_Const()
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSTINC_SUP) || !defined MODULE_SELECTION_ON         
      elseif (a%kfl_algor== 3) then
         a%NavierStokes => ThreeFieldNSProblem_Const()
      elseif (a%kfl_algor== 4) then
         a%NavierStokes => SUPFractionalStepProblem_Const()
#endif         
      else
         call runend('nsi_reapro: wrong algorithm for NavierStokes')
      endif
      call a%Memor%allocObj(0_ip,'PhysicalProblem','nstinc_SetPhysicalProblem',1)
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%NavierStokes)
      
   end subroutine
      
   subroutine nstinc_GetTimeStep(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%NavierStokes)
  
   end subroutine nstinc_GetTimeStep
   
   subroutine nstinc_SetTimeStep(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%NavierStokes)
  
   end subroutine nstinc_SetTimeStep
   
   subroutine nstinc_GetRefinementCriteria(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%NavierStokes)
  
   end subroutine nstinc_GetRefinementCriteria
   
   subroutine nstinc_PreRefine(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%NavierStokes)
  
   end subroutine nstinc_PreRefine
   
   subroutine nstinc_Refine(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%NavierStokes)
  
   end subroutine nstinc_Refine
   
   subroutine nstinc_Rebalance(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%NavierStokes)
  
   end subroutine nstinc_Rebalance
   
   subroutine nstinc_Begste(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c

      if (a%kfl_ModifyViscosity == 1 .and. c%masterVars%istep == a%ModifiedViscosityStep) then
         call a%NavierStokes%SetParameters(1_ip,a%ModifiedViscosity)
      endif
      
      !PhysicalProblem operations
      call physical_Begste(a,c,a%NavierStokes)
  
   end subroutine nstinc_Begste
   
   subroutine nstinc_MeshProjections(a,c,itask)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%NavierStokes,itask)
  
   end subroutine nstinc_MeshProjections
   
   subroutine nstinc_MeshAdvections(a,c,itask)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%NavierStokes,itask)
  
   end subroutine nstinc_MeshAdvections
   
   subroutine nstinc_Doiter(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%NavierStokes)
  
   end subroutine nstinc_Doiter
   
   subroutine nstinc_Convergence(a,c,glres)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%NavierStokes,glres)
  
   end subroutine nstinc_Convergence

   subroutine nstinc_CoupledConvergence(a,c,kfl_flag,cpres)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
      logical             :: kfl_flag
      real(rp):: cpres
  
      !PhysicalProblem operations
      call physical_CoupledConvergence(a,c,a%NavierStokes,kfl_flag,cpres)
   end subroutine nstinc_CoupledConvergence
   
   subroutine nstinc_Endste(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%NavierStokes)
  
   end subroutine nstinc_Endste

   subroutine nstinc_Output(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%NavierStokes)
  
   end subroutine nstinc_Output
   
   subroutine nstinc_Turnof(a,c)
      implicit none
      class(NstincDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%NavierStokes)
      
      !deallocate(a%NavierStokes)
      !I do not deallocate it, I need it to be alive in Cputab
      !We should do something about this
      call a%Memor%deallocObj(0_ip,'PhysicalProblem','SetPhysicalProblem',1_ip)
     
   end subroutine nstinc_Turnof
   
   subroutine nstinc_WriteTimes(a)
      implicit none
      class(NstincDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%NavierStokes)
  
   end subroutine nstinc_WriteTimes
   
   subroutine nstinc_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(NstincDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%NavierStokes
   end subroutine
  

end module 

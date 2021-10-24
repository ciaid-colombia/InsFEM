module Mod_PLCDDriver
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
    use Mod_PLCD
   implicit none

   private

   public :: PLCDDriver, PLCDDriver_Const
   
   !This is only for debugging
   type(PLCDProblem), pointer :: PLCD => NULL()

   type, extends(DriverInterface) :: PLCDDriver 
      
   type(PLCDProblem) :: PLCD

   !Couplings 
   type(InChannelGroupInfo) :: InNstincInfo
   type(InChannelGroupInfo) :: InAlemovInfo
     

contains

      procedure :: Lev1Reapro             => PLCD_Lev1Reapro
      procedure :: Lev1ReaproMPI          => PLCD_Lev1ReaproMPI
      procedure :: SetOutChannels         => PLCD_SetOutChannels
      procedure :: SetInChannels          => PLCD_SetInChannels
      procedure :: UpdateChannels         => PLCD_UpdateChannels
      procedure :: Turnon                 => PLCD_Turnon
      procedure :: GetTimeStep            => PLCD_GetTimeStep
      procedure :: SetTimeStep            => PLCD_SetTimeStep
      procedure :: GetRefinementCriteria  => PLCD_GetRefinementCriteria
      procedure :: PreRefine              => PLCD_PreRefine
      procedure :: Refine                 => PLCD_Refine
      procedure :: Rebalance              => PLCD_Rebalance
      procedure :: Begste                 => PLCD_Begste
      procedure :: MeshProjections        => PLCD_MeshProjections
      procedure :: MeshAdvections         => PLCD_MeshAdvections
      procedure :: Doiter                 => PLCD_Doiter
      procedure :: Convergence            => PLCD_Convergence
      procedure :: Endste                 => PLCD_Endste
      procedure :: Output                 => PLCD_Output
      procedure :: Turnof                 => PLCD_Turnof
      procedure :: WriteTimes             => PLCD_WriteTimes
      procedure :: GetPhysicalProblem     => PLCD_GetPhysicalProblem
      
   end type PLCDDriver

   interface PLCDDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(PLCDDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine PLCD_Lev1Reapro(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb

      if(a%Listener%words(1)=='EXTER') then
          if(a%Listener%words(2) == 'NSTIN') then
              call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
              a%PLCD%kfl_FSI = 1
          end if   
      end if   

      if(a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'ALEMO') then
             call InChannelGroupInfoFromListener(a%Listener,a%InAlemovInfo)
         end if   
      end if    
   end subroutine
   
   subroutine PLCD_Lev1ReaproMPI(a,c,bb)
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      integer(ip) :: ierr
      
      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
      call InChannelGroupInfoAddToBuffer(a%InAlemovInfo,bb)
      call MPI_BCAST(a%PLCD%kfl_FSI,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
   
   subroutine PLCD_SetOutChannels(a,task)
      class(PLCDDriver) :: a
      character(*) :: task
      real(rp), pointer :: displ(:,:)=> NULL()
      real(rp), pointer :: veloc(:,:)=> NULL()
      real(rp), pointer :: trac(:,:)=> NULL()

      !Displacements
      call a%PLCD%GetDisplacementArray(displ)
      call RPToList(displ,'DISPLACEMENT',task,a%OutChannelList)

      !Velocity
      call a%PLCD%GetVelocityArray(veloc)
      call RPToList(veloc,'VELOCITY',task,a%OutChannelList)

      !Tractions
      call a%PLCD%GetTractionArray(trac)
      call RPToList(trac,'TRACTION',task,a%OutChannelList)
      
   end subroutine
   
   subroutine PLCD_SetInChannels(a)
      class(PLCDDriver) :: a

      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('FLUIDTRACTION','TRACTION','INTERPOLABLE',a%InNstincInfo)
      endif
   end subroutine
      
   subroutine PLCD_UpdateChannels(a)
      implicit none
      class(PLCDDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      real(rp), pointer :: trac(:,:) => NULL()

      if (a%InNstincInfo%IsOn) then
         
         call a%GetFromInChannelList('FLUIDTRACTION',myDC)
         call ExtractRP(myDC,trac)
         call a%PLCD%SetExternalTractionArray(trac)

     endif
   end subroutine
   
   subroutine PLCD_Turnon(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      
      call point_plcd(a)
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%PLCD)
      

   end subroutine
   
   subroutine point_plcd(a)
      implicit none
      class(PLCDDriver), target :: a
      
      PLCD => a%PLCD
   end subroutine 
   
   
   subroutine PLCD_GetTimeStep(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%PLCD)
   end subroutine PLCD_GetTimeStep
   
   subroutine PLCD_SetTimeStep(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%PLCD)
   end subroutine PLCD_SetTimeStep
   
   subroutine PLCD_GetRefinementCriteria(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%PLCD)
   end subroutine PLCD_GetRefinementCriteria
   
   subroutine PLCD_PreRefine(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%PLCD)
   end subroutine PLCD_PreRefine
   
   subroutine PLCD_Refine(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%PLCD)
   end subroutine PLCD_Refine
   
   subroutine PLCD_Rebalance(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%PLCD)
   end subroutine PLCD_Rebalance
   
   subroutine PLCD_Begste(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%PLCD)
   end subroutine PLCD_Begste
   
   subroutine PLCD_MeshProjections(a,c,itask)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%PLCD,itask)
   end subroutine PLCD_MeshProjections
   
   subroutine PLCD_MeshAdvections(a,c,itask)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%PLCD,itask)
   end subroutine PLCD_MeshAdvections
   
   subroutine PLCD_Doiter(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%PLCD)
   end subroutine PLCD_Doiter
   
   subroutine PLCD_Convergence(a,c,glres)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%PLCD,glres)
   end subroutine PLCD_Convergence
   
   subroutine PLCD_Endste(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%PLCD)
   end subroutine PLCD_Endste

   subroutine PLCD_Output(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%PLCD)
   end subroutine PLCD_Output
   
   subroutine PLCD_Turnof(a,c)
      implicit none
      class(PLCDDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%PLCD)
   end subroutine PLCD_Turnof
   
   subroutine PLCD_WriteTimes(a)
      implicit none
      class(PLCDDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%PLCD)
   end subroutine PLCD_WriteTimes
   
   subroutine plcd_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(PLCDDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%PLCD
   end subroutine
   
     

end module 

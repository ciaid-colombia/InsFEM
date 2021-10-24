module Mod_SolidsDriver
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
    use Mod_Solids
    use Mod_SUPSolids
    use Mod_SUPSolids_lin
    use Mod_SUPSolids_NH
    use Mod_UPSolids_NH
   implicit none

   private

   public :: SolidsDriver, SolidsDriver_Const

   type, extends(DriverInterface) :: SolidsDriver 
      
      class(SolidsProblem), pointer :: Solids => NULL()

      !Couplings 
      type(InChannelGroupInfo) :: InNstincInfo

contains

      procedure :: Lev1Reapro             => Solids_Lev1Reapro
      procedure :: Lev1ReaproMPI          => Solids_Lev1ReaproMPI
      procedure :: SetOutChannels         => solids_SetOutChannels
      procedure :: SetInChannels          => solids_SetInChannels
      procedure :: UpdateChannels         => Solids_UpdateChannels
      procedure :: Turnon                 => Solids_Turnon
      procedure :: GetTimeStep            => Solids_GetTimeStep
      procedure :: SetTimeStep            => Solids_SetTimeStep
      procedure :: GetRefinementCriteria  => Solids_GetRefinementCriteria
      procedure :: Refine                 => Solids_Refine
      procedure :: Rebalance              => Solids_Rebalance
      procedure :: Begste                 => Solids_Begste
      procedure :: MeshProjections        => Solids_MeshProjections
      procedure :: MeshAdvections         => Solids_MeshAdvections
      procedure :: Doiter                 => Solids_Doiter
      procedure :: Convergence            => Solids_Convergence
      procedure :: CoupledConvergence     => Solids_CoupledConvergence
      procedure :: Endste                 => Solids_Endste
      procedure :: Output                 => Solids_Output
      procedure :: Turnof                 => Solids_Turnof
      procedure :: WriteTimes             => Solids_WriteTimes
      procedure :: GetPhysicalProblem     => Solids_GetPhysicalProblem
      
   end type SolidsDriver

   interface SolidsDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(SolidsDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine Solids_Lev1Reapro(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb

      if(a%Listener%words(1).eq.'ALGOR') then 
          if(a%Listener%exists('DISPL')) then
              a%kfl_algor= 1
          elseif(a%Listener%exists('SIGMA')) then
              a%kfl_algor= 2
          elseif(a%Listener%exists('SUPNH')) then
              a%kfl_algor= 3
          elseif(a%Listener%exists('UPNHO')) then
              a%kfl_algor= 4
          end if
      end if

      if(a%Listener%words(1)=='EXTER') then
          if(a%Listener%words(2) == 'NSTIN') then
              call InChannelGroupInfoFromListener(a%Listener,a%InNstincInfo)
          end if   
      end if   

   end subroutine
   
   subroutine Solids_Lev1ReaproMPI(a,c,bb)
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
      call bb%Add(a%kfl_algor)

      call InChannelGroupInfoAddToBuffer(a%InNstincInfo,bb)
   end subroutine
   
   subroutine solids_SetOutChannels(a,task)
      class(solidsDriver) :: a
      character(*) :: task
      real(rp), pointer :: displ(:,:)=> NULL(),veloc(:,:)=> NULL(),trac(:,:)=> NULL()

      !Displacements
      call a%Solids%GetDisplacementArray(displ)
      if(associated(displ)) call RPToList(displ,'DISPL',task,a%OutChannelList)

      !Velocity
      call a%Solids%GetVelocityArray(veloc)
      if(associated(veloc)) call RPToList(veloc,'VELOC',task,a%OutChannelList)

      !Tractions
      call a%Solids%GetInternalTractionArray(trac)
      if(associated(trac)) call RPToList(trac,'INTTRAC',task,a%OutChannelList)

   end subroutine
   
   subroutine solids_SetInChannels(a)
      class(solidsDriver) :: a

      if (a%InNstincInfo%IsOn) then
         call a%AddToInChannelList('TRACTION','TRACTION','INTERPOLABLE',a%InNstincInfo)
      endif

   end subroutine
      
   subroutine Solids_UpdateChannels(a)
      implicit none
      class(SolidsDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      real(rp), pointer :: trac(:,:) => NULL()

      if (a%InNstincInfo%IsOn) then
         
         call a%GetFromInChannelList('TRACTION',myDC)
         call ExtractRP(myDC,trac)
         call a%Solids%SetExternalTractionArray(trac)

     endif

   end subroutine
   
   subroutine Solids_Turnon(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      
      !PhysicalProblem operations
      if (a%kfl_algor== 1) then
          a%Solids => SolidsProblem_Const()
      elseif (a%kfl_algor== 2) then
          a%Solids => SUPSolidsProblem_lin_Const()
      elseif (a%kfl_algor== 3) then
          a%Solids => SUPSolidsProblem_NH_Const()
      elseif (a%kfl_algor== 4) then
          a%Solids => UPSolidsProblem_NH_Const()
      endif

      call a%Memor%allocObj(0_ip,'PhysicalProblem','Solids_SetPhysicalProblem',1)

      call physical_Turnon(a,c,a%Solids)
   end subroutine
      
   subroutine Solids_GetTimeStep(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%Solids)
   end subroutine Solids_GetTimeStep
   
   subroutine Solids_SetTimeStep(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%Solids)
   end subroutine Solids_SetTimeStep
   
   subroutine Solids_GetRefinementCriteria(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%Solids)
   end subroutine Solids_GetRefinementCriteria
   
   subroutine Solids_Refine(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%Solids)
   end subroutine Solids_Refine
   
   subroutine Solids_Rebalance(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%Solids)
   end subroutine Solids_Rebalance
   
   subroutine Solids_Begste(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%Solids)
   end subroutine Solids_Begste
   
   subroutine Solids_MeshProjections(a,c,itask)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%Solids,itask)
   end subroutine Solids_MeshProjections
   
   subroutine Solids_MeshAdvections(a,c,itask)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%Solids,itask)
   end subroutine Solids_MeshAdvections
   
   subroutine Solids_Doiter(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      logical::kfl_folded=.false.
      integer(ip):: ierr
      integer(ip) :: kfl_gotim

      !PhysicalProblem operations
      call physical_Doiter(a,c,a%Solids)
      call a%Solids%GetFoldedElements(kfl_folded)

      if(kfl_folded) then
          if (c%masterVars%MPIrank == c%masterVars%MPIroot) then
              write(*,*) '****:Solids Folded Elements, attempting to shutdown THIS case smoothly'
          endif
          !call runend('Folded Elements')
          !End time step
          call flush
          !call MPI_Barrier(c%masterVars%MPIcomm, ierr)
          !call a%Solids%Endste(kfl_gotim)

          !Case output
          !call flush
          call a%Solids%Output(1_ip)
          call MPI_Barrier(c%masterVars%MPIcomm, ierr)

          call a%Solids%Turnof
          call MPI_Barrier(c%masterVars%MPIcomm, ierr)
          call runend('Solid Driver: Folded Elements')
      endif

   end subroutine Solids_Doiter
   
   subroutine Solids_Convergence(a,c,glres)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%Solids,glres)
   end subroutine Solids_Convergence

   subroutine Solids_CoupledConvergence(a,c,kfl_flag,cpres)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
      logical             :: kfl_flag
      real(rp):: cpres
  
      !PhysicalProblem operations
      call physical_CoupledConvergence(a,c,a%Solids,kfl_flag,cpres)
   end subroutine Solids_CoupledConvergence
   
   subroutine Solids_Endste(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%Solids)
   end subroutine Solids_Endste

   subroutine Solids_Output(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%Solids)
   end subroutine Solids_Output
   
   subroutine Solids_Turnof(a,c)
      implicit none
      class(SolidsDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%Solids)

      call a%Memor%deallocObj(0_ip,'PhysicalProblem','SetPhysicalProblem',1_ip)
   end subroutine Solids_Turnof
   
   subroutine Solids_WriteTimes(a)
      implicit none
      class(SolidsDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%Solids)
   end subroutine Solids_WriteTimes
   
   subroutine Solids_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(SolidsDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%Solids
   end subroutine
   
end module 

module Mod_LowMachDriver
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
    use Mod_LowMach
   implicit none

   private

   public :: LowMachDriver, LowMachDriver_Const

   type, extends(DriverInterface) :: LowMachDriver 
      
      type(LowMachProblem) :: LowMach
     

contains

      procedure :: Lev1Reapro             => lowmach_Lev1Reapro
      procedure :: Lev1ReaproMPI          => lowmach_Lev1ReaproMPI
      procedure :: SetOutChannels         => lowmach_SetOutChannels
      procedure :: SetInChannels          => lowmach_SetInChannels
      procedure :: UpdateChannels         => lowmach_UpdateChannels
      procedure :: Turnon                 => lowmach_Turnon
      procedure :: GetTimeStep            => lowmach_GetTimeStep
      procedure :: SetTimeStep            => lowmach_SetTimeStep
      procedure :: GetRefinementCriteria  => lowmach_GetRefinementCriteria
      procedure :: PreRefine              => lowmach_PreRefine
      procedure :: Refine                 => lowmach_Refine
      procedure :: Rebalance              => lowmach_Rebalance
      procedure :: Begste                 => lowmach_Begste
      procedure :: MeshProjections        => lowmach_MeshProjections
      procedure :: MeshAdvections         => lowmach_MeshAdvections
      procedure :: Doiter                 => lowmach_Doiter
      procedure :: Convergence            => lowmach_Convergence
      procedure :: Endste                 => lowmach_Endste
      procedure :: Output                 => lowmach_Output
      procedure :: Turnof                 => lowmach_Turnof
      procedure :: WriteTimes             => lowmach_WriteTimes
      procedure :: GetPhysicalProblem     => lowmach_GetPhysicalProblem
      
   end type LowMachDriver

   interface LowMachDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(LowMachDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine lowmach_Lev1Reapro(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
   end subroutine
   
   subroutine lowmach_Lev1ReaproMPI(a,c,bb)
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
   end subroutine
   
   subroutine lowmach_SetOutChannels(a,task)
      class(LowMachDriver) :: a
      character(*) :: task
   end subroutine
   
   subroutine lowmach_SetInChannels(a)
      class(LowMachDriver) :: a
   end subroutine
      
   subroutine lowmach_UpdateChannels(a)
      implicit none
      class(LowMachDriver) :: a
   end subroutine
   
   subroutine lowmach_Turnon(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%LowMach)
   end subroutine
      
   subroutine lowmach_GetTimeStep(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%LowMach)
   end subroutine lowmach_GetTimeStep
   
   subroutine lowmach_SetTimeStep(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%LowMach)
   end subroutine lowmach_SetTimeStep
   
   subroutine lowmach_GetRefinementCriteria(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%LowMach)
   end subroutine lowmach_GetRefinementCriteria
   
   subroutine lowmach_PreRefine(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%LowMach)
   end subroutine lowmach_PreRefine
   
   subroutine lowmach_Refine(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%LowMach)
   end subroutine lowmach_Refine
   
   subroutine lowmach_Rebalance(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%LowMach)
   end subroutine lowmach_Rebalance
   
   subroutine lowmach_Begste(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%LowMach)
   end subroutine lowmach_Begste
   
   subroutine lowmach_MeshProjections(a,c,itask)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%LowMach,itask)
   end subroutine lowmach_MeshProjections
   
   subroutine lowmach_MeshAdvections(a,c,itask)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%LowMach,itask)
   end subroutine lowmach_MeshAdvections
   
   subroutine lowmach_Doiter(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%LowMach)
   end subroutine lowmach_Doiter
   
   subroutine lowmach_Convergence(a,c,glres)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%LowMach,glres)
   end subroutine lowmach_Convergence
   
   subroutine lowmach_Endste(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%LowMach)
   end subroutine lowmach_Endste

   subroutine lowmach_Output(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%LowMach)
   end subroutine lowmach_Output
   
   subroutine lowmach_Turnof(a,c)
      implicit none
      class(LowMachDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%LowMach)
   end subroutine lowmach_Turnof
   
   subroutine lowmach_WriteTimes(a)
      implicit none
      class(LowMachDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%LowMach)
   end subroutine lowmach_WriteTimes
   
   subroutine lowmach_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(LowMachDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%LowMach
   end subroutine
   
     

end module 

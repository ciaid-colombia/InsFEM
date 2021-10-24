module Mod_emptyDriver
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
    use Mod_empty
   implicit none

   private

   public :: emptyDriver, emptyDriver_Const

   type, extends(DriverInterface) :: emptyDriver 
      
      type(emptyProblem) :: empty
     

contains

      procedure :: Lev1Reapro             => empty_Lev1Reapro
      procedure :: Lev1ReaproMPI          => empty_Lev1ReaproMPI
      procedure :: SetInChannels          => empty_SetInChannels
      procedure :: SetOutChannels         => empty_SetOutChannels
      procedure :: UpdateChannels         => empty_UpdateChannels
      procedure :: Turnon                 => empty_Turnon
      procedure :: GetTimeStep            => empty_GetTimeStep
      procedure :: SetTimeStep            => empty_SetTimeStep
      procedure :: GetRefinementCriteria  => empty_GetRefinementCriteria
      procedure :: Refine                 => empty_Refine
      procedure :: Rebalance              => empty_Rebalance
      procedure :: Begste                 => empty_Begste
      procedure :: MeshProjections        => empty_MeshProjections
      procedure :: MeshAdvections         => empty_MeshAdvections
      procedure :: Doiter                 => empty_Doiter
      procedure :: Convergence            => empty_Convergence
      procedure :: Endste                 => empty_Endste
      procedure :: Output                 => empty_Output
      procedure :: Turnof                 => empty_Turnof
      procedure :: WriteTimes             => empty_WriteTimes
      
   end type emptyDriver

   interface emptyDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(emptyDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine empty_Lev1Reapro(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
   end subroutine
   
   subroutine empty_Lev1ReaproMPI(a,c,bb)
      class(emptyDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
   end subroutine
   
   subroutine empty_SetOutChannels(a,task)
      class(emptyDriver) :: a
      character(*) :: task
   end subroutine
   
   subroutine empty_SetInChannels(a)
      class(emptyDriver) :: a
   end subroutine
   
      
   subroutine empty_UpdateChannels(a)
      implicit none
      class(emptyDriver) :: a
   end subroutine
   
   subroutine empty_Turnon(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%empty)
   end subroutine
      
   subroutine empty_GetTimeStep(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%empty)
   end subroutine empty_GetTimeStep
   
   subroutine empty_SetTimeStep(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%empty)
   end subroutine empty_SetTimeStep
   
   subroutine empty_GetRefinementCriteria(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%empty)
   end subroutine empty_GetRefinementCriteria
   
   subroutine empty_Refine(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%empty)
   end subroutine empty_Refine
   
   subroutine empty_Rebalance(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%empty)
   end subroutine empty_Rebalance
   
   subroutine empty_Begste(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%empty)
   end subroutine empty_Begste
   
   subroutine empty_MeshProjections(a,c,itask)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%empty,itask)
   end subroutine empty_MeshProjections
   
   subroutine empty_MeshAdvections(a,c,itask)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%empty,itask)
   end subroutine empty_MeshAdvections
   
   subroutine empty_Doiter(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%empty)
   end subroutine empty_Doiter
   
   subroutine empty_Convergence(a,c,glres)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%empty,glres)
   end subroutine empty_Convergence
   
   subroutine empty_Endste(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%empty)
   end subroutine empty_Endste

   subroutine empty_Output(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%empty)
   end subroutine empty_Output
   
   subroutine empty_Turnof(a,c)
      implicit none
      class(emptyDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%empty)
   end subroutine empty_Turnof
   
   subroutine empty_WriteTimes(a)
      implicit none
      class(emptyDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%empty)
   end subroutine empty_WriteTimes
   
     

end module 

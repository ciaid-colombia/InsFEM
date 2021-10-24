module Mod_NSCompDriver
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
    use Mod_NSCompressible
    use Mod_NSCompressibleExplicit
    use Mod_NSCompressibleImplicit
    use Mod_NSCompressiblePrimitive
   implicit none

   private

   public :: nscompDriver, nscompDriver_Const, CompressibleNavierStokes

   type, extends(DriverInterface) :: nscompDriver 
      
      class(NSCompressibleProblem), pointer :: CompressibleNavierStokes => NULL()
     
      !Set the ALE mesh
      integer(ip) :: kfl_ALEMesh = 0  !Default is Eulerian 

contains

      procedure :: Lev1Reapro             => nscomp_Lev1Reapro
      procedure :: Lev1ReaproMPI          => nscomp_Lev1ReaproMPI
      procedure :: SetOutChannels         => nscomp_SetOutChannels
      procedure :: SetInChannels          => nscomp_SetInChannels
      procedure :: UpdateChannels         => nscomp_UpdateChannels
      procedure :: Turnon                 => nscomp_Turnon
      procedure :: GetTimeStep            => nscomp_GetTimeStep
      procedure :: SetTimeStep            => nscomp_SetTimeStep
      procedure :: GetRefinementCriteria  => nscomp_GetRefinementCriteria
      procedure :: PreRefine              => nscomp_PreRefine
      procedure :: Refine                 => nscomp_Refine
      procedure :: Rebalance              => nscomp_Rebalance
      procedure :: Begste                 => nscomp_Begste
      procedure :: MeshProjections        => nscomp_MeshProjections
      procedure :: MeshAdvections         => nscomp_MeshAdvections
      procedure :: Doiter                 => nscomp_Doiter
      procedure :: Convergence            => nscomp_Convergence
      procedure :: Endste                 => nscomp_Endste
      procedure :: Output                 => nscomp_Output
      procedure :: Turnof                 => nscomp_Turnof
      procedure :: WriteTimes             => nscomp_WriteTimes
      procedure :: GetPhysicalProblem     => nscomp_GetPhysicalProblem
      
   end type nscompDriver

   interface nscompDriver_Const
      procedure constructor
   end interface 

   class(NSCompressibleProblem), pointer :: CompressibleNavierStokes => NULL()

contains

   function constructor()
       class(nscompDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine nscomp_Lev1Reapro(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb

      if(a%Listener%words(1).eq.'ALGOR') then 
         if(a%Listener%exists('EXPLI')) then
             a%kfl_algor= 1
         end if
         if(a%Listener%exists('IMPLI')) then
             a%kfl_algor= 2
         end if
         if(a%Listener%exists('PRIMI')) then
             a%kfl_algor= 3
         end if

      !Things that we want from outside
      else if(a%Listener%words(1)=='EXTER') then
         if(a%Listener%words(2) == 'ALEMO') then
             a%kfl_ALEMesh = 1
         end if
      end if
   end subroutine
   
   subroutine nscomp_Lev1ReaproMPI(a,c,bb)
      class(nscompDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb

      call bb%Add(a%kfl_ALEMesh)
      call bb%Add(a%kfl_algor)
      
   end subroutine
   
   subroutine nscomp_SetOutChannels(a,task)
      class(nscompDriver) :: a
      character(*) :: task
      real(rp), pointer :: veloc(:,:) => NULL()

      !Veloc
      call a%CompressibleNavierStokes%GetVelocityArray(veloc)
      call RPToList(veloc,'VELOC',task,a%OutChannelList)

   end subroutine
   
   subroutine nscomp_SetInChannels(a)
      class(nscompDriver) :: a
   end subroutine
      
   subroutine nscomp_UpdateChannels(a)
      implicit none
      class(nscompDriver) :: a
   end subroutine
   
   subroutine nscomp_Turnon(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
      
      !After broadcast operations
      if (a%kfl_algor== 1) then
         a%CompressibleNavierStokes => NSCompressibleExplicitProblem_Const()
      elseif (a%kfl_algor== 2) then
         a%CompressibleNavierStokes => NSCompressibleImplicitProblem_Const()
      elseif (a%kfl_algor== 3) then
         a%CompressibleNavierStokes => NSCompressiblePrimitiveProblem_Const()
      else
         call runend('nsc_reapro: wrong algorithm for CompressibleNavierStokes')
      endif

      !PhysicalProblem operations
      call physical_Turnon(a,c,a%CompressibleNavierStokes)
      !Set ALE mesh
      if(a%kfl_ALEMesh == 1) then
         call a%CompressibleNavierStokes%Mesh%SetALE(1)
      end if
   end subroutine
      
   subroutine nscomp_GetTimeStep(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_GetTimeStep
   
   subroutine nscomp_SetTimeStep(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_SetTimeStep
   
   subroutine nscomp_GetRefinementCriteria(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_GetRefinementCriteria

   subroutine nscomp_PreRefine(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_PreRefine
   
   subroutine nscomp_Refine(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Refine
   
   subroutine nscomp_Rebalance(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Rebalance
   
   subroutine nscomp_Begste(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Begste
   
   subroutine nscomp_MeshProjections(a,c,itask)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%CompressibleNavierStokes,itask)
   end subroutine nscomp_MeshProjections
   
   subroutine nscomp_MeshAdvections(a,c,itask)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%CompressibleNavierStokes,itask)
   end subroutine nscomp_MeshAdvections
   
   subroutine nscomp_Doiter(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Doiter
   
   subroutine nscomp_Convergence(a,c,glres)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%CompressibleNavierStokes,glres)
   end subroutine nscomp_Convergence
   
   subroutine nscomp_Endste(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Endste
   
   subroutine nscomp_Output(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Output
   
   subroutine nscomp_Turnof(a,c)
      implicit none
      class(nscompDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%CompressibleNavierStokes)
   end subroutine nscomp_Turnof
   
   subroutine nscomp_WriteTimes(a)
      implicit none
      class(nscompDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%CompressibleNavierStokes)
   end subroutine nscomp_WriteTimes
   
   subroutine nscomp_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(nscompDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%CompressibleNavierStokes
   end subroutine
     

end module 

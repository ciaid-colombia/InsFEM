module Mod_GeneralCase
   use typre
   use Mod_DistributedContainerList
   use Mod_DistributedContainer
   use Mod_caseVariables
   use Mod_GeneralParallelLibrary 
   use Mod_DriverInterface
   use Mod_Postpr_VTK
   use Mod_Postpr_NULL
   use Mod_int2str
   implicit none

   private
   public :: GeneralCase, GeneralCase_Const

   type :: GeneralCase

      !Lists of Drivers in the case
      type(DistributedContainerList) :: DriverList

      !List of interpolators between cases
      type(DistributedContainerList) :: CaseInterpolatorList
      
      !All kind of case local variables
      type(caseVariables) :: caseVars
      
contains

      procedure :: SetMPI
      procedure :: Reapro                       => case_Reapro
      procedure :: SetArgumentsAndEnvironment   => case_SetArgumentsAndEnvironment 
      procedure :: SetMulticommData             => case_SetMulticommData
      procedure :: Turnon                       => case_Turnon 
      procedure :: Domain                       => case_Domain
      procedure :: OldDomain                    => case_OldDomain
      procedure :: Adaptive                     => case_Adaptive
      procedure :: InitialUniformRefinement     => case_InitialUniformRefinement
      procedure :: SetupInterpolators           => case_SetupInterpolators
      procedure :: ReSetupInterpolators         => case_ReSetupInterpolators
      procedure :: SetupIO                      => case_SetupIO
      procedure :: DoTimeStep                   => case_DoTimeStep
      procedure :: UpdateChannels               => case_UpdateChannels
      procedure :: GetGoTimeStepping            => case_GetGoTimeStepping
      procedure :: GetDoCoupledConv             => case_GetDoCoupledConv
      procedure :: GetMaxCoupledIters           => case_GetMaxCoupledIters
      procedure :: Setste                       => case_Setste
      procedure :: Begste                       => case_Begste
      procedure :: MeshProjections              => case_MeshProjections
      procedure :: Doiter                       => case_Doiter
      procedure :: LocalCouplingConvergence     => case_LocalCouplingConvergence
      procedure :: GlobalCouplingConvergence    => case_GlobalCouplingConvergence
      procedure :: Endste                       => case_Endste 
      procedure :: Output                       => case_Output
      procedure :: Outmem                       => case_Outmem
      procedure :: Cputab                       => case_Cputab
      procedure :: Turnof                       => case_Turnof
      procedure :: FinalizeInterpolators        => case_FinalizeInterpolators




   end type 
    
interface 
        
      subroutine case_Reapro(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_Turnon(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_Domain(a,itask)
         use typre
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
         integer(ip) :: itask
      end subroutine

      subroutine case_OldDomain(a,itask)
         use typre
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
         integer(ip) :: itask
      end subroutine

      subroutine case_Adaptive(a,itask)
         use typre
         import GeneralCase
         implicit none
         integer(ip) :: itask
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_InitialUniformRefinement(a,itask)
         use typre
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
         integer(ip) :: itask
      end subroutine
      
      subroutine case_SetupInterpolators(a,bOwnCaseKey)
         use Mod_DCHashCharSize
         import GeneralCase
         implicit none
         class(GeneralCase) :: a
         character(DCHashCharSize) :: bOwnCaseKey
      end subroutine
      
      subroutine case_ReSetupInterpolators(a)
         use Mod_DCHashCharSize
         import GeneralCase
         implicit none
         class(GeneralCase) :: a
      end subroutine
      
      subroutine case_FinalizeInterpolators(a)
         import GeneralCase
         implicit none
         class(GeneralCase) :: a
      end subroutine

      subroutine case_SetupIO(a,itask)
         use typre
         import GeneralCase
         implicit none
         class(GeneralCase) :: a
         integer(ip) :: itask
      end subroutine
      
      subroutine case_DoTimeStep(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_UpdateChannels(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine
      
      subroutine case_MeshProjections(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine
      
      subroutine case_Setste(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_Begste(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine 

      subroutine case_Doiter(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine 

      subroutine case_LocalCouplingConvergence(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine 

      subroutine case_GlobalCouplingConvergence(a,cpiter,kfl_coupledconv)
         import 
         implicit none
         class(GeneralCase), target :: a
         logical :: kfl_coupledconv
         integer(ip):: cpiter
      end subroutine 

      subroutine case_Endste(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine 
      
      subroutine case_Output(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine 

      subroutine case_Turnof(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_outmem(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

      subroutine case_cputab(a)
         import GeneralCase
         implicit none
         class(GeneralCase), target :: a 
      end subroutine

   end interface
   
   interface GeneralCase_Const
      procedure constructor
   end interface

contains

   function constructor()
      type(GeneralCase), pointer :: constructor

      allocate(constructor)
   end function

   subroutine SetMPI(a,MPIcomm,MPIsize,MPIroot,MPIrank)
      class(GeneralCase), target :: a
      integer(ip) :: MPIsize,MPIroot,MPIrank,MPIcomm

      type(masterVariables), pointer :: m => NULL()

      m => a%caseVars%masterVars

      m%MPIsize = MPIsize
      m%MPIroot = MPIroot
      m%MPIrank = MPIrank
      m%MPIcomm = MPIcomm
   end subroutine
    
   subroutine case_GetGoTimeStepping(a,kfl_gotim)
      class(GeneralCase) :: a
      integer(ip) :: kfl_gotim

      kfl_gotim = a%caseVars%masterVars%kfl_gotim
   end subroutine

   subroutine case_GetDoCoupledConv(a,kfl_flag)
      class(GeneralCase) :: a
      logical :: kfl_flag

      kfl_flag = a%caseVars%masterVars%kfl_doCoupledConv
   end subroutine

   subroutine case_GetMaxCoupledIters(a,maxiters)
      class(GeneralCase) :: a
      integer(ip):: maxiters

      maxiters = a%caseVars%masterVars%kfl_maxCoupledIters
   end subroutine

   subroutine case_SetArgumentsAndEnvironment(a,nam, readType, baseData, resFol, postP, rstFol, olddataFol, oldrestarFol, kfl_multicase)
        class(GeneralCase), target :: a
        character(150) :: nam,readType, baseData, resFol, postP, rstFol, oldDataFol, oldrestarFol
        integer(ip) :: kfl_multicase
        
        type(masterVariables), pointer :: m => NULL()
        m => a%caseVars%masterVars
       
        m%kfl_multicase = kfl_multicase
        m%namda = nam
        m%ReadTypeString = readType 
        m%BaseDataFolder =  baseData
        m%ResultsFolder =   resFol
        m%PostProcessFolder = postP
        m%RestartFolder =  rstFol
        m%OldDataFolder = olddataFol
        m%OldRestartFolder = oldrestarFol

        !Modify data folder if ReadPartitioned approach
        if (m%ReadTypeString == 'SERIAL') then
            m%DataFolder = m%BaseDataFolder

        elseif (m%ReadTypeString == 'PARTITIONED') then
            m%DataFolder = trim(m%BaseDataFolder)//'/data'//int2str(m%MPIrank)
        else
            call runend('Wrong type of read strategy')
        endif
   end subroutine 
   
   subroutine case_SetMulticommData(a,multicomm,multicommColor)
      implicit none
      class(GeneralCase), target :: a
      integer(ip) :: multicomm, multicommColor
      
      type(masterVariables), pointer :: m => NULL()
      
      m => a%caseVars%masterVars
      m%kfl_multicomm = multicomm
      m%multicommColor = multicommColor
   end subroutine
   
end module Mod_GeneralCase

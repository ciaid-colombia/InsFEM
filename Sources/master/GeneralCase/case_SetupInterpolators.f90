module Mod_caseSetupInterpolators
   use Mod_GeneralCase
   use Mod_CaseInterpolator
   use Mod_DC_CaseInterpolator
   use Mod_DistributedContainer
   use Mod_CaseVariables
   use MOd_DriverInterface
   use Mod_DC_Driver
   use Mod_InChannel
   use Mod_DC_INChannel
   implicit none
   
   !character(DCHashCharSize), save :: OwnCaseKey
   type(DistributedContainerList), pointer :: CaseInterpolatorList => NULL()
   type(FemMesh), pointer :: Mesh => NULL()
   type(MemoryMan), pointer :: Memo => NULL()
   real(rp) :: InterpolatorTolerance

contains   
   
   subroutine LoopSetupInterpolators(myDC)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      
      call ExtractCaseInterpolator(myDC,myCaseInterpolator)
      
      call myCaseInterpolator%SetInterpolatorTolerance(InterpolatorTolerance)
      call myCaseInterpolator%SetupInterpolator(Mesh,Memo)
   end subroutine
   
   subroutine LoopReSetupInterpolators(myDC)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      
      call ExtractCaseInterpolator(myDC,myCaseInterpolator)
      
      call myCaseInterpolator%ReSetupInterpolator(Mesh,Memo)
   end subroutine
   
   subroutine LoopDriverSetup(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%InChannelList%LoopList(LoopInChannelSetup)
   end subroutine
   
   subroutine LoopInChannelSetup(myDC)
      class(DistributedContainer), pointer :: myDC
      type(InChannel), pointer :: myInChannel => NULL()
      
      call ExtractInChannel(myDC,myInChannel)
      
      call myInChannel%SetupInterpolation(Memo)
   end subroutine

   
   subroutine LoopFinalizeInterpolators(myDC)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      
      call ExtractCaseInterpolator(myDC,myCaseInterpolator)
      
      call myCaseInterpolator%FinalizeInterpolator(Memo)
      call Memo%deallocObj(0,'CaseInterpolator','ReadInterpolator',1)
   end subroutine
   
   subroutine LoopDriverFinalize(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%InChannelList%LoopList(LoopInChannelFinalize)
   end subroutine
   
   subroutine LoopInChannelFinalize(myDC)
      class(DistributedContainer), pointer :: myDC
      type(InChannel), pointer :: myInChannel => NULL()
      
      call ExtractInChannel(myDC,myInChannel)
      call myInChannel%FinalizeInterpolation(Memo)
   end subroutine
      
   
   
end module
   
subroutine case_SetupInterpolators(a,bOwnCaseKey)
   use Mod_caseSetupInterpolators
   implicit none
   class(GeneralCase), target :: a
   character(*) :: bOwnCaseKey
   
!   OwnCasekey = bOwnCaseKey
   Mesh => a%caseVars%domainVars%Mesh
   Memo => a%caseVars%masterVars%MasterMemo
   InterpolatorTolerance = a%caseVars%masterVars%Interp_Tol
   
   CaseInterpolatorList => a%CaseInterpolatorList
   
   !We loop through the interpolator list
   call a%CaseInterpolatorList%LoopList(LoopSetupInterpolators)
   
   !Loop through Drivers, and then through Channels to set up everything
   call a%DriverList%LoopList(LoopDriverSetup)
   

 
end subroutine

subroutine case_ReSetupInterpolators(a)
   use Mod_caseSetupInterpolators
   implicit none
   class(GeneralCase), target :: a
   
!   OwnCasekey = bOwnCaseKey
   Mesh => a%caseVars%domainVars%Mesh
   Memo => a%caseVars%masterVars%MasterMemo
   CaseInterpolatorList => a%CaseInterpolatorList
   
   !We loop through the interpolator list
   call a%CaseInterpolatorList%LoopList(LoopReSetupInterpolators)
   
   !Loop through Drivers, and then through Channels to set up everything
   call a%DriverList%LoopList(LoopDriverSetup)
   

 
end subroutine

subroutine case_FinalizeInterpolators(a)
   use Mod_caseSetupInterpolators
   implicit none
   class(GeneralCase), target :: a

   Memo => a%caseVars%masterVars%MasterMemo
   CaseInterpolatorList => a%CaseInterpolatorList
   
   !We loop through the interpolator list
   call a%CaseInterpolatorList%LoopList(LoopFinalizeInterpolators)
   
   call a%DriverList%LoopList(LoopDriverFinalize)

end subroutine

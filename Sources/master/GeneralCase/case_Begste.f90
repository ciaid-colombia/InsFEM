module Mod_Begste
   !This routine starts a time step
   use typre
   use Mod_GeneralCase
   use Mod_CaseVariables
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   use Mod_DCHashCharSize
   implicit none
   
   type(caseVariables), pointer :: c => NULL()
contains
   
   subroutine LoopBegste(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%Begste(c)
   end subroutine
end module

subroutine case_Begste(a)
   !This routine starts a time step
   use Mod_Begste
   implicit none
   class(GeneralCase), target :: a
   
   character(DCHashCharSize) :: DoFirstKeys(2)
   
   c => a%caseVars
   
   !Levset must be the first one: in case of adaptivity we need the cuts to be computed
   DoFirstKeys(1) = 'LEVSE'
   DoFirstKeys(2) = 'ALEPR'
   call a%DriverList%LoopList(LoopBegste,DoFirstKeys)

end subroutine 

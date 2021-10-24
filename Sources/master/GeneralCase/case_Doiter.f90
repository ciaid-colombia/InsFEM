module Mod_caseDoiter
   use Mod_GeneralCase
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DCHashCharSize
   use Mod_caseVariables
   use Mod_DC_Driver
   implicit none
   
   type(caseVariables), pointer :: c => NULL()
   
contains
   subroutine LoopDoIter(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%DoIter(c)
      
   end subroutine
end module   

subroutine case_Doiter(a)
   use Mod_caseDoiter
   class(GeneralCase), target :: a
   
   
   character(DCHashCharSize) :: DoFirstKeys(2)
   
   c => a%caseVars
   
   DoFirstKeys(1) = 'ALEPR'  !alefirst the first one, many modules rely on its results
   DoFirstKeys(2) = 'NSTIN'  !Navier Stokes the first one, many modules rely on its results
   call a%DriverList%LoopList(LoopDoiter,DoFirstKeys)
end subroutine

module Mod_CaseVariables
   use Mod_MasterVariables
   use Mod_DomainVariables
   use Mod_AdaptiveVariables
   
   type CaseVariables
      type (MasterVariables)   :: masterVars
      type (DomainVariables)   :: domainVars
      type (AdaptiveVariables) :: adaptiveVars

      integer(ip) :: ndrivers
   end type

end module

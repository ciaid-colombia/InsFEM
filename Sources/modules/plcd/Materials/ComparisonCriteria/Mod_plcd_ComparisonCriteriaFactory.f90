module Mod_plcd_ComparisonCriteriaFactory
   use typre
   use Mod_plcd_ComparisonCriteria
   use Mod_plcd_SimoJu
   use Mod_plcd_VonMises
   implicit none

contains   
   subroutine ComparisonCriteriaFactory(CompareCriteria,CmpCriteria)
      character(5) :: CompareCriteria
      class(ComparisonCriteria ), pointer :: CmpCriteria
      
      if (CompareCriteria == 'SIMOJ') then
         allocate(SimoJu::CmpCriteria)
      elseif (CompareCriteria == 'VONMI') then
         allocate(VonMises::CmpCriteria)
      endif
   end subroutine

 
end module
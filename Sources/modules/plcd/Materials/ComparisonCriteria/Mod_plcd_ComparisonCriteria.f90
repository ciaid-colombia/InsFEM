module Mod_plcd_ComparisonCriteria
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_Listen
   implicit none
   
   type, abstract :: ComparisonCriteria
   
   
contains
      procedure :: ReadData
      procedure :: ScatterData
      procedure :: Evaluate
      procedure :: ComputeTauValueFromUniaxialStressTest
   
   end type
      
   
   
  
   
contains
   subroutine ReadData(a,Listener)
      class(ComparisonCriteria) :: a
      type(ListenFile) :: Listener

   end subroutine
   
   subroutine ScatterData(a)
      class(ComparisonCriteria) :: a
   end subroutine
   
   subroutine Evaluate(a,Stress,Strain,CT,tau)
      class(ComparisonCriteria) :: a
      real(rp) :: Stress(:),Strain(:),tau
      class(ConstitutiveTensor) :: CT
      
      call runend('Need to implement Evaluate for this comparison Criteria')
   end subroutine  
   
   subroutine ComputeTauValueFromUniaxialStressTest(a,UniaxialStressLimit,CT,tau)
      class(ComparisonCriteria) :: a
      real(rp) :: UniaxialStressLimit,tau
      class(ConstitutiveTensor) :: CT
      
      call runend('Need to implement ComputeTauValueFromUniaxialStressTest for this comparison Criteria')
   end subroutine  
   
  
end module
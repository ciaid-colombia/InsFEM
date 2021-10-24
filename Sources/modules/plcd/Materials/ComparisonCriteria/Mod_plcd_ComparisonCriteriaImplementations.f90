module Mod_plcd_SimoJU
   use typre
   use Mod_plcd_ComparisonCriteria
   use Mod_Plcd_ConstitutiveTensor
   implicit none
   
   type, extends(ComparisonCriteria) :: SimoJu

   
      contains
      procedure :: Evaluate => SimoJuEvaluate
      procedure :: ComputeTauValueFromUniaxialStressTest => SimoJuComputeTauValueFromUniaxialStressTest
   end type
contains
   subroutine SimoJuComputeTauValueFromUniaxialStressTest(a,UniaxialStressLimit,CT,tau)
      class(SimoJu) :: a
      real(rp) :: UniaxialStressLimit,tau
      class(ConstitutiveTensor) :: CT
      
      real(rp) :: Cnorm
      real(rp) :: Direction(3), E0
      
      Direction = (/1.0_rp, 0.0_rp, 0.0_rp/)
      call CT%GetYoungModulusForDirection(Direction,E0)
!       
      call CT%GetCNorm(Cnorm)
      
      Tau = UniaxialStressLimit*sqrt(1.5_rp*Cnorm/E0)

   end subroutine 
   
   subroutine SimoJuEvaluate(a,Stress,Strain,CT,tau)
      class(SimoJU) :: a
      real(rp) :: Stress(:), Strain(:), tau
      class(ConstitutiveTensor) :: CT

       real(rp) :: Cnorm
       
       call CT%GetCNorm(Cnorm)
      
      Tau = sqrt(dot_product(Stress,strain))*sqrt(1.5_rp*Cnorm)
      
   end subroutine
   
   
end module


module Mod_plcd_VonMises
   use typre
   use Mod_plcd_ComparisonCriteria
   use Mod_Plcd_ConstitutiveTensor
   implicit none
   
   type, extends(ComparisonCriteria) :: VonMises
      contains
      procedure :: Evaluate => VonMisesEvaluate
      procedure :: ComputeTauValueFromUniaxialStressTest => VonMisesComputeTauValueFromUniaxialStressTest
   end type
contains

   subroutine VonMisesComputeTauValueFromUniaxialStressTest(a,UniaxialStressLimit,CT,tau)
      class(VonMises) :: a
      real(rp) :: UniaxialStressLimit,tau
      class(ConstitutiveTensor) :: CT
      
      Tau = UniaxialStressLimit
   end subroutine



   subroutine VonMisesEvaluate(a,Stress,Strain,CT,tau)
      class(VonMises) :: a
      real(rp) :: Stress(:),Strain(:),tau
      class(ConstitutiveTensor) :: CT
      
      real(rp) :: deviatoricStress(size(Stress,1)),DeviatoricStressNorm,meanStress
      
      call CT%ComputeMeanStress(stress,meanstress)
      
      deviatoricStress = stress
      call CT%AddSphericalComponent(-meanstress,DeviatoricStress)
      
      call CT%ComputeStressNorm(DeviatoricStress,DeviatoricStressNorm)
   
      Tau = sqrt(3.0_rp/2.0_rp)*DeviatoricStressNorm
   end subroutine

end module
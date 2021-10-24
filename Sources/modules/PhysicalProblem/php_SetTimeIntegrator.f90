module Mod_php_SetTimeIntegrator
   implicit none
   private
   public php_SetTimeIntegrator
   public php_SetNewmarkCoefficients
   public php_GetTimeSchemeDerivative

   interface php_SetTimeIntegrator
      module procedure php_SetTimeIntegrator1
      module procedure php_SetTimeIntegrator2
   end interface

contains

subroutine php_SetTimeIntegrator1(a,Integrator,LHSDtinv,nsteps)
   use typre
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   type(TimeIntegratorDt1) :: Integrator
   real(rp)               :: LHSDtinv
   integer(ip)            :: nsteps
   
   call Integrator%Init(a%kfl_tsche_1st_current)
   call Integrator%GetLHSDtinv(a%dtinv,LHSdtinv)
   call Integrator%GetNumberOfTimeSteps(nsteps)
end subroutine

subroutine php_SetTimeIntegrator2(a,Integrator,LHSDtinv2,nsteps)
   use typre
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   type(TimeIntegratorDt2)  :: Integrator
   real(rp)               :: LHSDtinv2,dtinv
   integer(ip)            :: nsteps
   
   call Integrator%SetTimeStep(a%dtinv)
   call Integrator%Init(a%kfl_tsche_2nd_current)
   call Integrator%GetLHSDtinv(a%dtinv2,LHSdtinv2)
   call Integrator%GetNumberOfTimeSteps(nsteps)
end subroutine

subroutine php_SetNewmarkCoefficients(a,Integrator,beta,gamma)
   use typre
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   type(TimeIntegratorDt2):: Integrator
   real(rp)               :: beta,gamma

   call Integrator%SetNewmarkCoefficients(beta,gamma)
end subroutine

subroutine php_GetTimeSchemeDerivative(a,Integrator,deriv)
   use typre
   use Mod_TimeIntegrator
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   type(TimeIntegratorDt2):: Integrator
   real(rp)               :: deriv

   call Integrator%GetTimeSchemeDerivative(deriv)
end subroutine

end module

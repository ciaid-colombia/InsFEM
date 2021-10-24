module Mod_TimeIntegrator
   use typre
   use Mod_Memor
   implicit none
   private
   !public TimeIntegratorBase
   public TimeIntegratorDt1  ! first  time derivative du/dt
   public TimeIntegratorDt2  ! second time derivative d2u/dt2
   public GetTimeProjection
   
   type, abstract :: TimeIntegratorBase
      real(rp)              :: Coefficients(10),dtinv
      real(rp)              :: nwmBeta    = 4.0_rp !Defaults to second order, this is 1/beta
      real(rp)              :: nwmGamma   = 0.5_rp    !Defaults to second order, this is gamma
      real(rp)              :: derivative          !df_n+1/dx_n+1, used for newton raphson
      integer(ip)           :: nsteps,accura
      contains
         procedure(Init)  , deferred :: Init
         procedure(GetRHS)              , deferred :: GetRHS
         procedure(GetLHSDtinv)         , deferred :: GetLHSDtinv
         procedure :: GetTimeSchemeCoefficients
         procedure :: SetNewmarkCoefficients
         procedure :: SetTimeStep
         procedure :: GetTimeSchemeDerivative
         procedure :: GetNumberOfTimeSteps
         procedure :: GetAccuracy
   end type
   
   abstract interface
      subroutine Init(a,TimeScheme)
         use typre
         import TimeIntegratorBase
         implicit none
         class(TimeIntegratorBase) :: a
         character(5) :: TimeScheme
      end subroutine
      
      subroutine GetRHS(a,ndime,gpval,valRHS)
         use typre
         import TimeIntegratorBase
         implicit none
         class(TimeIntegratorBase) :: a
         integer(ip) :: ndime
         real(rp)    :: gpval(ndime,a%nsteps-1), valRHS(ndime)
      end subroutine
      
      subroutine GetLHSDtinv(a,dtinv,dtinvm)
         use typre
         import TimeIntegratorBase
         implicit none
         class(TimeIntegratorBase) :: a
         real(rp) :: dtinv,dtinvm
      end subroutine

   end interface
   
   type, extends(TimeIntegratorBase) :: TimeIntegratorDt1
      contains
         procedure :: Init        => InitTimeIntegrator1st
         procedure :: GetRHS      => GetRHS1
         procedure :: GetLHSDtinv => GetLHSDtinv1
   end type
   
   type, extends(TimeIntegratorBase) :: TimeIntegratorDt2
      contains
         procedure :: Init         => InitTimeIntegrator2nd
         procedure :: GetRHS       => GetRHS2
         procedure :: GetLHSDtinv  => GetLHSDtinv2
   end type
   
contains

   subroutine NULLSUB(a,ndofn,gpval,gdpdt,valRHS)
      use typre
      implicit none
      class(TimeIntegratorDt1) :: a
      integer(ip) :: ndofn
      real(rp)    :: gpval(ndofn,a%nsteps-1), gdpdt(ndofn), valRHS(ndofn)
   end subroutine

   subroutine GetNumberOfTimeSteps(a,NumberOfTimeSteps)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      integer(ip) :: NumberOfTimeSteps
      NumberOfTimeSteps = a%nsteps
   end subroutine
   
   subroutine GetAccuracy(a,Accuracy)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      integer(ip) :: Accuracy
      Accuracy = a%accura
   end subroutine
      
   
   subroutine InitTimeIntegrator1st(a,TimeScheme)
      use typre
      implicit none
      class(TimeIntegratorDt1) :: a
      character(5) :: TimeScheme
      type(MemoryMan) :: Memor
      
      if (TimeScheme == 'NONE ') then
         a%Coefficients(1) = 0.0_rp  !no change of sign
         a%Coefficients(2) = 0.0_rp  !notice the change of sign
         a%nsteps = 2  !number of time steps needed to store
         a%accura = 0  !accuracy of the time scheme
      elseif (TimeScheme == 'BDF1 ') then
         a%Coefficients(1) = 1.0_rp  !no change of sign
         a%Coefficients(2) = 1.0_rp  !notice the change of sign
         a%derivative      = 1.0_rp
         a%nsteps = 2  !number of time steps needed to store
         a%accura = 1  !accuracy of the time scheme
      elseif (TimeScheme == 'BDF2 ') then
         a%Coefficients(1) =  1.5_rp  !no change of sign
         a%Coefficients(2) =  2.0_rp  !notice the change of sign
         a%Coefficients(3) = -0.5_rp  !notice the change of sign
         a%derivative      = 1.5_rp
         a%nsteps = 3  !number of time steps needed to store
         a%accura = 2  !accuracy of the time scheme
      elseif (TimeScheme == 'BDF3 ') then
         a%Coefficients(1) =  11.0_rp/6.0_rp  !no change of sign
         a%Coefficients(2) =  3.0_rp          !notice the change of sign
         a%Coefficients(3) = -1.5_rp          !notice the change of sign
         a%Coefficients(4) =  1.0_rp/3.0_rp   !notice the change of sign
         a%derivative      = 11.0_rp/6.0_rp
         a%nsteps = 4  !number of time steps needed to store
         a%accura = 3  !accuracy of the time scheme
      elseif (TimeScheme == 'BDF4 ') then
         a%Coefficients(1) =  25.0_rp/12.0_rp  !no change of sign
         a%Coefficients(2) =  4.0_rp           !notice the change of sign
         a%Coefficients(3) = -3.0_rp           !notice the change of sign
         a%Coefficients(4) =  4.0_rp/3.0_rp    !notice the change of sign
         a%Coefficients(5) = -1.0_rp/4.0_rp    !notice the change of sign
         a%derivative      = 25.0_rp/12.0_rp
         a%nsteps = 5  !number of time steps needed to store
         a%accura = 4  !accuracy of the time scheme
      elseif (TimeScheme == 'CNOBS') then
         a%Coefficients(1) = 1.0_rp  !no change of sign
         a%Coefficients(2) = 1.0_rp  !notice the change of sign
         a%nsteps = 2  !number of time steps needed to store
         a%accura = 2  !accuracy of the time scheme
      elseif (TimeScheme == 'CN   ') then
         a%Coefficients(1) = 2.0_rp  !no change of sign
         a%Coefficients(2) = 2.0_rp  !notice the change of sign
         a%nsteps = 2  !number of time steps needed to store
         a%accura = 2  !accuracy of the time scheme
      elseif (TimeScheme == 'TR   ') then ! BC and forces computed as 0.5*(u^n+1 + u^n) and 0.5*(f^n+1 + f^n)
         a%Coefficients(1) = 2.0_rp  !no change of sign
         a%Coefficients(2) = 2.0_rp  !notice the change of sign
         a%nsteps = 2  !number of time steps needed to store
         a%accura = 2  !accuracy of the time scheme
      elseif (TimeScheme == 'NEWMA ') then
         a%derivative = a%nwmBeta*a%nwmGamma
         a%accura = 2
      else                             ! unknown time scheme
         call runend('Mod_TimeIntegrator: InitTimeIntegrator1st: unknown time scheme: '//TimeScheme)
      endif
      
   end subroutine InitTimeIntegrator1st

   subroutine InitTimeIntegrator2nd(a,TimeScheme)
      use typre
      implicit none
      class(TimeIntegratorDt2) :: a
      character(5) :: TimeScheme
      type(MemoryMan) :: Memor
      
      a%Coefficients = 0.0
      
      if (TimeScheme == 'NONE ') then
         a%nsteps = 2
      elseif (TimeScheme == 'BDF1 ') then
         a%Coefficients(1) =  1.0_rp  !no change of sign
         a%Coefficients(2) =  2.0_rp  !notice the change of sign
         a%Coefficients(3) = -1.0_rp  !notice the change of sign
         a%derivative      =  1.0_rp
         a%nsteps = 3
         a%accura = 1
      elseif (TimeScheme == 'BDF2 ') then
         a%Coefficients(1) =  2.0_rp  !no change of sign
         a%Coefficients(2) =  5.0_rp  !notice the change of sign
         a%Coefficients(3) = -4.0_rp  !notice the change of sign
         a%Coefficients(4) =  1.0_rp  !notice the change of sign
         a%derivative      =  2.0_rp
         a%nsteps = 4
         a%accura = 2
      elseif (TimeScheme == 'BDF3 ') then
         a%Coefficients(1) =  35.0_rp/12.0_rp  !no change of sign
         a%Coefficients(2) =  26.0_rp/3.0_rp   !notice the change of sign
         a%Coefficients(3) = -19.0_rp/2.0_rp   !notice the change of sign
         a%Coefficients(4) =  14.0_rp/3.0_rp   !notice the change of sign
         a%Coefficients(5) = -11.0_rp/12.0_rp  !notice the change of sign
         a%derivative      =  35.0_rp/12.0_rp
         a%nsteps = 5
         a%accura = 3
      elseif (TimeScheme == 'BDF4 ') then
         a%Coefficients(1) =  15.0_rp/4.0_rp    !no change of sign
         a%Coefficients(2) =  77.0_rp/6.0_rp    !notice the change of sign
         a%Coefficients(3) = -107.0_rp/6.0_rp   !notice the change of sign
         a%Coefficients(4) =  13.0_rp           !notice the change of sign
         a%Coefficients(5) = -61.0_rp/12.0_rp   !notice the change of sign
         a%Coefficients(6) =  5.0_rp/6.0_rp     !notice the change of sign
         a%derivative      =  15.0_rp/4.0_rp
         a%nsteps = 6
         a%accura = 4
      elseif (TimeScheme == 'NEWMA') then
         a%Coefficients(1) =  a%nwmBeta              !no change of sign
         a%Coefficients(2) =  a%nwmBeta              !notice the change of sign
         a%Coefficients(3) =  a%nwmBeta/a%dtinv      !notice the change of sign
         a%Coefficients(4) =  ((a%nwmBeta/2.0_rp)-1.0_rp)*(1.0_rp/a%dtinv)**2_ip !notice the change of sign
         a%derivative      =  a%nwmBeta
         a%nsteps = 4
      else                             ! unknown time scheme
         call runend('Mod_TimeIntegrator: InitTimeIntegrator2nd: unknown time scheme: '//TimeScheme)
      endif
      
   end subroutine InitTimeIntegrator2nd

   subroutine SetTimeStep(a,dtinv)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      real(rp)     :: dtinv

      a%dtinv = dtinv
   end subroutine SetTimeStep

   subroutine SetNewmarkCoefficients(a,beta,gamma)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      real(rp)     :: beta, gamma

      a%nwmBeta = 1.0_rp/beta
      a%nwmGamma = gamma
   end subroutine SetNewmarkCoefficients
   
   subroutine GetTimeSchemeDerivative(a,deriv)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      real(rp)     :: deriv

      deriv = a%derivative
   end subroutine GetTimeSchemeDerivative

   subroutine GetTimeSchemeCoefficients(a,c)
      use typre
      implicit none
      class(TimeIntegratorBase) :: a
      real(rp)     :: c(10)

      c = a%Coefficients
   end subroutine GetTimeSchemeCoefficients
   
   subroutine GetRHS1(a,ndime,gpval,valRHS)
      use typre
      implicit none
      class(TimeIntegratorDt1) :: a
      integer(ip) :: ndime
      real(rp)    :: gpval(ndime,a%nsteps-1), valRHS(ndime)
      
      valRHS = matmul(gpval,a%Coefficients(2:a%nsteps))
      
   end subroutine GetRHS1
   
   subroutine GetRHS2(a,ndime,gpval,valRHS)
      use typre
      implicit none
      class(TimeIntegratorDt2) :: a
      integer(ip) :: ndime
      real(rp)    :: gpval(ndime,a%nsteps-1), valRHS(ndime)
      
      valRHS = matmul(gpval,a%Coefficients(2:a%nsteps))
      
   end subroutine
   
   subroutine GetLHSDtinv1(a,dtinv,dtinvm)
      use typre
      implicit none
      class(TimeIntegratorDt1) :: a
      real(rp) :: dtinv,dtinvm

      dtinvm = dtinv*a%Coefficients(1)

   end subroutine
   
   subroutine GetLHSDtinv2(a,dtinv,dtinvm)
      use typre
      implicit none
      class(TimeIntegratorDt2) :: a
      real(rp) :: dtinv,dtinvm
      
      dtinvm = dtinv*a%Coefficients(1)
      
   end subroutine
   
   subroutine GetTimeProjection(order,ndof1,npoin,ncomp,xunkno,unknoProy)
      use typre
      implicit none
      integer(ip), intent(in)  :: order,ndof1,npoin,ncomp
      real(rp), intent(in)     :: xunkno(ndof1,npoin,ncomp)
      real(rp), intent(out)  :: unknoProy(ndof1,npoin,1)
   
      if(order==2)then
         unknoProy(:,:,1) = 2.0_rp*xunkno(:,:,3)-xunkno(:,:,4)       
      elseif(order==3)then
        unknoProy(:,:,1) = (3.0_rp)*xunkno(:,:,3)-(3.0_rp)*xunkno(:,:,4)+(1.0_rp)*xunkno(:,:,5)        
      else                             ! unknown time scheme
         call runend('GetTimeProjection: unknown time Projection') 
      end if
   
   end subroutine
   
end module

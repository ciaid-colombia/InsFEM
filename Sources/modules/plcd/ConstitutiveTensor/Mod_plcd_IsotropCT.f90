module Mod_plcd_IsotropConstitutiveTensor
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_Listen
   implicit none
   private
   public IsotropCT

   type, extends(ConstitutiveTensor) :: IsotropCT

      real(rp) :: E
      real(rp) :: nu
      real(rp) :: invK,CNorm
      real(rp) :: C(6,6),CSph(6,6),CDev(6,6)

contains
      procedure :: ReadData => RDIsotrop
      procedure :: ScatterData => SPIsotrop
      procedure :: CheckData => CheckData
      procedure :: Build => BuildIsotrop
      procedure :: GetPointerToC => GetPointerToCIsotrop
      procedure :: GetPointerToCDeviatoric => GetPointerToCDeviatoricIsotrop
      procedure :: GetPointerToCSpherical => GetPointerToCSphericalIsotrop
      procedure :: GetYoungModulusEquiv => GetYoungModulusEquivIsotrop

      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetCNorm
      procedure :: ComputeStressNorm
      procedure :: ComputeMeanStress
      procedure :: GetYoungModulusForDirection
      procedure :: GetPoissonRatio
      procedure :: GetYoungModulus
      
   end type

contains

   subroutine RDIsotrop(a,Listener)
      implicit none
      class(IsotropCT) :: a
      type(ListenFile) :: Listener

      call Listener%Listen('RDIsotrop')
      if (Listener%exists('YOUNG')) a%E = Listener%GetRea('YOUNG',1.0_rp,'Young Modulus')
      call Listener%Listen('RDIsotrop')
      if (Listener%exists('POISS')) a%nu = Listener%GetRea('POISS',0.3_rp,'Poisson Coefficient')
   
   end subroutine

   subroutine SPIsotrop(a)
      use MPI
      implicit none
      class(IsotropCT) :: a
      integer(ip) :: ierr
      call MPI_BCAST(a%E,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nu,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
   
   subroutine CheckData(a)
      use MPI
      implicit none
      class(IsotropCT) :: a
      if(a%E <= 0.0_rp) call runend('Negative value for the Young Modulus')
      if((a%nu <= -1.0_rp).OR.(a%nu > 0.5_rp)) call runend('Poisson Coefficient out of range (-1 < v <= 0.5)')
      
   end subroutine

   subroutine BuildIsotrop(a)
      implicit none
      class(IsotropCT) :: a

      real(rp) :: raux
      real(rp) :: Eigenvalues(6)

      !Deviatoric part
      a%CDev(:,1) = a%E/(1+a%nu)*(/ 2.0_rp/3.0_rp,-1.0_rp/3.0_rp,-1.0/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
      a%CDev(:,2) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,2.0_rp/3.0_rp,-1.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
      a%CDev(:,3) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,-1.0_rp/3.0_rp,2.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
      a%CDev(:,4) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp,0.0_rp /)
      a%CDev(:,5) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp /)
      a%CDev(:,6) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp /)

      if (.not. a%OnlyDeviatoricPart) then
         raux = a%nu/(1-2*a%nu)
         a%C(:,1) = a%E/(1+a%nu)*(/ 1.0_rp+raux,raux,raux,0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,2) = a%E/(1+a%nu)*(/ raux,1.0_rp+raux,raux,0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,3) = a%E/(1+a%nu)*(/ raux,raux,1.0_rp+raux,0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,4) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp,0.0_rp /)
         a%C(:,5) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp /)
         a%C(:,6) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp /)

      elseif (a%OnlyDeviatoricPart) then
         a%C = a%CDev
      endif

      !Volumetric part
      a%invK = (3*(1-2*a%nu))/a%E
      a%CSph = a%C - a%CDev

      call MatrixEigenValues(6,a%C,Eigenvalues)

      a%CNorm = maxval(abs(Eigenvalues))


   end subroutine

   subroutine GetPointerToCIsotrop(a,C)
      implicit none
      class(IsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      C => a%C
   end subroutine

   subroutine GetPointerToCSphericalIsotrop(a,C)
      implicit none
      class(IsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      C => a%CSph
   end subroutine

   subroutine GetPointerToCDeviatoricIsotrop(a,C)
      implicit none
      class(IsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      C => a%CDev
   end subroutine

   subroutine GetInverseVolumetricDeformationModulus(a,invK)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: invK

      invK  = a%invK

   end subroutine

   subroutine GetCNorm(a,CNorm)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: CNorm

      CNorm = a%CNorm

   end subroutine

   subroutine ComputeStressNorm(a,Stress,StressNorm)
      implicit none
      class(IsotropCT), target :: a
      real(rp) :: Stress(:),StressNorm

      real(rp) :: dot1,dot2

      dot1 = dot_product(Stress(1:3),Stress(1:3))
      dot2 = dot_product(Stress(4:6),Stress(4:6))
      StressNorm = sqrt(dot1+2.0_rp*dot2)

   end subroutine

   subroutine ComputeMeanStress(a,Stress,MeanStress)
      implicit none
      class(IsotropCT), target :: a
      real(rp) :: Stress(:),MeanStress

      meanStress = sum(Stress(1:3))/3.0
   end subroutine
   
   subroutine GetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      E = a%E
   end subroutine

   subroutine GetPoissonRatio(a,nu)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: nu

      nu = a%nu
   end subroutine
   
   subroutine GetYoungModulus(a,E)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: E

      E = a%E
   end subroutine
   

   subroutine GetYoungModulusEquivIsotrop(a,E)
      implicit none
      class(IsotropCT) :: a
      real(rp) :: E

      E = a%E
   end subroutine
end module

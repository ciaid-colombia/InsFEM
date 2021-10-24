module Mod_plcd_Isotrop2DPlainStrain
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_Listen
   implicit none
   private
   public Isotrop2dPlainStrain

   type, extends(ConstitutiveTensor) :: Isotrop2dPlainStrain

      real(rp) :: E
      real(rp) :: nu
      real(rp) :: invK,CNorm
      real(rp) :: C(3,3),CDev(3,3),CSph(3,3)

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
      procedure :: SetYoungModulusForDirection
      procedure :: SetPoissonRatio
      procedure :: SetYoungModulus
   end type

contains

   subroutine RDIsotrop(a,Listener)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      type(ListenFile) :: Listener
      
      call Listener%Listen('RDIsotrop')
      if (Listener%exists('YOUNG')) a%E = Listener%GetRea('YOUNG',1.0_rp,'Young Modulus')
      call Listener%Listen('RDIsotrop')
      if (Listener%exists('POISS')) a%nu = Listener%GetRea('POISS',0.3_rp,'Poisson Coefficient')
   end subroutine

   subroutine SPIsotrop(a)
      use MPI
      implicit none
      class(Isotrop2dPlainStrain) :: a
      integer(ip) :: ierr
      call MPI_BCAST(a%E,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nu,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   end subroutine

   subroutine CheckData(a)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      if(a%E <= 0.0_rp) call runend('Negatice value for the Young Modulus')
      if((a%nu <= -1.0_rp).OR.(a%nu > 0.5_rp)) call runend('Poisson Coefficient out of range (-1 < v <= 0.5)')
      
   end subroutine
   
   
   subroutine BuildIsotrop(a)
      implicit none
      class(Isotrop2dPlainStrain) :: a

      real(rp) :: raux
      real(rp) :: Eigenvalues(3)

      !Deviatoric part
      a%CDev(:,1) = a%E/(1+a%nu)*(/ 2.0_rp/3.0_rp,-1.0_rp/3.0_rp,0.0_rp /)
      a%CDev(:,2) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,2.0_rp/3.0_rp,0.0_rp /)
      a%CDev(:,3) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,1.0_rp/2.0_rp /)


      if (.not. a%OnlyDeviatoricPart) then
         raux = a%nu/(1-2*a%nu)
         a%C(:,1) = a%E/(1+a%nu)*(/ 1.0_rp+raux,raux,0.0_rp /)
         a%C(:,2) = a%E/(1+a%nu)*(/ raux,1.0_rp+raux,0.0_rp /)
         a%C(:,3) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,1.0_rp/2.0_rp /)


      elseif (a%OnlyDeviatoricPart) then
         a%C = a%Cdev
      endif

      !Volumetric part
      a%invK = (3*(1-2*a%nu))/a%E
      a%CSph = a%C-a%CDev

      call MatrixEigenValues(3,a%C,Eigenvalues)

      a%CNorm = maxval(abs(Eigenvalues))


   end subroutine

   subroutine GetPointerToCIsotrop(a,C)
      implicit none
      class(Isotrop2dPlainStrain), target :: a
      real(rp), pointer :: C(:,:)

      C => a%C
   end subroutine

   subroutine GetPointerToCDeviatoricIsotrop(a,C)
      implicit none
      class(Isotrop2dPlainStrain), target :: a
      real(rp), pointer :: C(:,:)

      C => a%CDev
   end subroutine

   subroutine GetPointerToCSphericalIsotrop(a,C)
      implicit none
      class(Isotrop2dPlainStrain), target :: a
      real(rp), pointer :: C(:,:)

      C => a%CSph
   end subroutine

   subroutine GetInverseVolumetricDeformationModulus(a,invK)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: invK

      invk = a%invK

   end subroutine

   subroutine GetCNorm(a,CNorm)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: CNorm

      CNorm = a%CNorm

   end subroutine

    subroutine ComputeStressNorm(a,Stress,StressNorm)
      implicit none
      class(Isotrop2dPlainStrain), target :: a
      real(rp) :: Stress(:),StressNorm

      real(rp) :: dot1,dot2,dot3 ,sigma3

      dot1 = dot_product(Stress(1:2),Stress(1:2))
      dot2 = Stress(3)*Stress(3)

      sigma3 = a%nu*(Stress(1)+Stress(2))
      dot3 = sigma3*sigma3

      StressNorm = sqrt(dot1+2.0_rp*dot2+dot3)

   end subroutine

   subroutine ComputeMeanStress(a,Stress,MeanStress)
      implicit none
      class(Isotrop2dPlainStrain), target :: a
      real(rp) :: Stress(:),MeanStress

      meanStress = sum(Stress(1:2))*(1+a%nu)/3.0
   end subroutine
   
   subroutine GetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      E = a%E
   end subroutine

   subroutine GetPoissonRatio(a,nu)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: nu

      nu = a%nu
   end subroutine
   

   subroutine GetYoungModulusEquivIsotrop(a,E)
 	  implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: E

      E = a%E
   end subroutine

   subroutine GetYoungModulus(a,E)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: E

      E = a%E
   end subroutine
   
   subroutine SetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      a%E = E
   end subroutine

   subroutine SetPoissonRatio(a,nu)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: nu

      a%nu = nu
   end subroutine
   
   subroutine SetYoungModulus(a,E)
      implicit none
      class(Isotrop2dPlainStrain) :: a
      real(rp) :: E

      a%E = E
   end subroutine


end module

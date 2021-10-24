module Mod_plcd_TransIsotropConstitutiveTensor
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_Listen
   implicit none
   private
   public TransIsotropCT

   type, extends(ConstitutiveTensor) :: TransIsotropCT

      real(rp) :: Et
      real(rp) :: Ep
      real(rp) :: Gt
      real(rp) :: Gp
      real(rp) :: nutp
      real(rp) :: invK,CNorm
      real(rp) :: C(6,6),CSph(6,6),CDev(6,6)

contains
      procedure :: ReadData => RDTransIsotrop
      procedure :: ScatterData => SPTransIsotrop
      procedure :: CheckData => CheckData
      procedure :: Build => BuildTransIsotrop
      procedure :: GetPointerToC => GetPointerToCTransIsotrop
      procedure :: GetPointerToCDeviatoric => GetPointerToCDeviatoricIsotrop
      procedure :: GetPointerToCSpherical => GetPointerToCSphericalIsotrop
      procedure :: GetYoungModulusEquiv => GetYoungModulusEquivTransIsotrop
      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetCNorm
      procedure :: ComputeStressNorm
      procedure :: ComputeMeanStress
      procedure :: GetYoungModulusForDirection
      procedure :: GetPoissonRatio

   end type

contains

   subroutine RDTransIsotrop(a,Listener)
      implicit none
      class(TransIsotropCT) :: a
      type(ListenFile) :: Listener

      call Listener%Listen('RDTransIsotrop')
      if (Listener%exists('YOUNG').AND.Listener%exists('TRANS')) a%Et = Listener%GetRea('TRANS',1.0_rp,'Young Modulus Transverse (Dir X)')
      call Listener%Listen('RDTransIsotrop')
      if (Listener%exists('YOUNG').AND.Listener%exists('INPLA')) a%Ep = Listener%GetRea('INPLA',1.0_rp,'Young Modulus InPlane (Dir Y and Z)')
      
      
      call Listener%Listen('RDTransIsotrop')
      if (Listener%exists('SHEAR').AND.Listener%exists('TRANS')) a%Gt = Listener%GetRea('TRANS',1.0_rp,'Shear Modulus Trans (Gxy and Gxz )')
      call Listener%Listen('RDTransIsotrop')
      if (Listener%exists('SHEAR').AND.Listener%exists('INPLA')) a%Gp = Listener%GetRea('INPLA',1.0_rp,'Shear Modulus InPlane (Gyz)')
      
      call Listener%Listen('RDTransIsotrop')
      if (Listener%exists('POISS').AND.Listener%exists('TRANS')) a%nutp = Listener%GetRea('TRANS',0.3_rp,'Poisson Coefficient TransInPlane')
   end subroutine

   subroutine SPTransIsotrop(a)
      use MPI
      implicit none
      class(TransIsotropCT) :: a
      integer(ip) :: ierr
      call MPI_BCAST(a%Et,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%Ep,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%Gt,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%Gp,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nutp,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   end subroutine
   
   subroutine CheckData(a)
      implicit none
      class(TransIsotropCT) :: a

      real(rp) :: raux
      
      if(a%Et <= 0.0_rp) call runend('Negative value for the Young Modulus Dir Tranverse')
      if(a%Ep <= 0.0_rp) call runend('Negative value for the Young Modulus Dir InPlane')      
   
      if(a%Gt <= 0.0_rp) call runend('Negative value for the Shear Modulus Dir Transverse')
      if(a%Gp <= 0.0_rp) call runend('Negative value for the Shear Modulus Dir InPlane')
      
      if(abs(a%nutp) >= sqrt(a%Et/a%Ep)) call runend('Poisson Coefficient transInplane out of range (|vtp| < (Et/Ep)^0.5)')
      if(abs(a%Ep/(2.0_rp*a%Gp)-1.0_rp) >= 1.0_rp) call runend('Poisson Coefficient InPlane out of range (|vp| < 1')
      
      raux = 1.0_rp - (a%Ep/(2.0_rp*a%Gp)-1.0_rp)**2.0_rp - 2.0_rp*a%nutp**2.0_rp*a%Ep/a%Et - 2.0_rp*(a%Ep/(2.0_rp*a%Gp)-1.0_rp)*(a%nutp**2.0_rp*a%Ep/a%Et)
      if(raux <= 0.0_rp) call runend('Negatice value for aux value Tranverse Isotropic Material')
      
   end subroutine

   subroutine BuildTransIsotrop(a)
      implicit none
      class(TransIsotropCT) :: a

      real(rp) :: raux, vpt, vp
!      real(rp) :: Eigenvalues(6)

      !Deviatoric part
!      a%CDev(:,1) = a%E/(1+a%nu)*(/ 2.0_rp/3.0_rp,-1.0_rp/3.0_rp,-1.0/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,2) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,2.0_rp/3.0_rp,-1.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,3) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,-1.0_rp/3.0_rp,2.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,4) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,5) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp /)
!      a%CDev(:,6) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp /)

!      if (.not. a%OnlyDeviatoricPart) then
         vpt = a%Ep/a%Et*a%nutp
         vp = a%Ep/(2.0_rp*a%Gp)-1.0_rp

         raux = 1.0_rp/(1.0_rp-vp**2.0_rp-2.0_rp*a%nutp*vpt-2.0_rp*vp*vpt*a%nutp)
         
         a%C(:,1) = raux*(/ a%Et*(1.0_rp-vp**2.0_rp),a%Et*(vpt+vpt*vp),a%Et*(vpt+vpt*vp),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,2) = raux*(/ a%Et*(vpt+vpt*vp),a%Ep*(1.0_rp-a%nutp*vpt),a%Ep*(vp+a%nutp*vpt),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,3) = raux*(/ a%Et*(vpt+vpt*vp),a%Ep*(vp+a%nutp*vpt),a%Ep*(1.0_rp-a%nutp*vpt),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,4) = a%Gt*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp /)
         a%C(:,5) = a%Gt*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp /)
         a%C(:,6) = a%Gp*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp /)
         
         
         
!      elseif (a%OnlyDeviatoricPart) then
!         a%C = a%CDev
!      endif

      !Volumetric part
!      a%invK = (3*(1-2*a%nu))/a%E
!      a%CSph = a%C - a%CDev

!      call MatrixEigenValues(6,a%C,Eigenvalues)

!      a%CNorm = maxval(abs(Eigenvalues))


   end subroutine

   subroutine GetPointerToCTransIsotrop(a,C)
      implicit none
      class(TransIsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      C => a%C
   end subroutine

   subroutine GetPointerToCSphericalIsotrop(a,C)
      implicit none
      class(TransIsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetPointerToCDeviatoricIsotrop(a,C)
      implicit none
      class(TransIsotropCT), target :: a
      real(rp), pointer :: C(:,:)

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetInverseVolumetricDeformationModulus(a,invK)
      implicit none
      class(TransIsotropCT) :: a
      real(rp) :: invK

      call runend('not implmemetnd0')

   end subroutine

   subroutine GetCNorm(a,CNorm)
      implicit none
      class(TransIsotropCT) :: a
      real(rp) :: CNorm

      call runend('not implmemetnd0')

   end subroutine

   subroutine ComputeStressNorm(a,Stress,StressNorm)
      implicit none
      class(TransIsotropCT), target :: a
      real(rp) :: Stress(:),StressNorm

      real(rp) :: dot1,dot2

      dot1 = dot_product(Stress(1:3),Stress(1:3))
      dot2 = dot_product(Stress(4:6),Stress(4:6))
      StressNorm = sqrt(dot1+2.0_rp*dot2)

   end subroutine

   subroutine ComputeMeanStress(a,Stress,MeanStress)
      implicit none
      class(TransIsotropCT), target :: a
      real(rp) :: Stress(:),MeanStress

      meanStress = sum(Stress(1:3))/3.0
   end subroutine

   subroutine GetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(TransIsotropCT) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetPoissonRatio(a,nu)
      implicit none
      class(TransIsotropCT) :: a
      real(rp) :: nu

      call runend('not implmemetnd0')
  end subroutine
  
   subroutine GetYoungModulusEquivTransIsotrop(a,E)
      implicit none
      class(TransIsotropCT) :: a
      real(rp) :: E

      E = a%Et
   end subroutine

end module

module Mod_plcd_OrthotropConstitutiveTensor
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_Listen
   implicit none
   private
   public OrthotropCT

   type, extends(ConstitutiveTensor) :: OrthotropCT

      real(rp) :: E1
      real(rp) :: E2
      real(rp) :: E3
      real(rp) :: G12
      real(rp) :: G13
      real(rp) :: G23
      real(rp) :: nu12
      real(rp) :: nu13
      real(rp) :: nu23
      real(rp) :: invK,CNorm
      real(rp) :: C(6,6),CSph(6,6),CDev(6,6)

contains
      procedure :: ReadData => RDOrthotrop
      procedure :: ScatterData => SPOrthotrop
      procedure :: CheckData => CheckData
      procedure :: Build => BuildOrthotrop
      procedure :: GetPointerToC => GetPointerToCOrthotrop
      procedure :: GetPointerToCDeviatoric => GetPointerToCDeviatoricIsotrop
      procedure :: GetPointerToCSpherical => GetPointerToCSphericalIsotrop
      procedure :: GetYoungModulusEquiv => GetYoungModulusEquivOrthotrop

      procedure :: GetInverseVolumetricDeformationModulus
      procedure :: GetCNorm
      procedure :: ComputeStressNorm
      procedure :: ComputeMeanStress
      procedure :: GetYoungModulusForDirection
      procedure :: GetPoissonRatio
   end type

contains

   subroutine RDOrthotrop(a,Listener)
      implicit none
      class(OrthotropCT) :: a
      type(ListenFile) :: Listener

!      Listener%param(1)
      
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('YOUNG').AND.Listener%exists('DIR11')) a%E1 = Listener%GetRea('DIR11',1.0_rp,'Young Modulus Dir X')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('YOUNG').AND.Listener%exists('DIR22')) a%E2 = Listener%GetRea('DIR22',1.0_rp,'Young Modulus Dir Y')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('YOUNG').AND.Listener%exists('DIR33')) a%E3 = Listener%GetRea('DIR33',1.0_rp,'Young Modulus Dir Z')
      
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('SHEAR').AND.Listener%exists('DIR12')) a%G12 = Listener%GetRea('DIR12',1.0_rp,'Shear Modulus Dir XY')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('SHEAR').AND.Listener%exists('DIR13')) a%G13 = Listener%GetRea('DIR13',1.0_rp,'Shear Modulus Dir XZ')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('SHEAR').AND.Listener%exists('DIR23')) a%G23 = Listener%GetRea('DIR23',1.0_rp,'Shear Modulus Dir YZ')
            
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('POISS').AND.Listener%exists('DIR12')) a%nu12 = Listener%GetRea('DIR12',0.3_rp,'Poisson Coefficient Dir XY')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('POISS').AND.Listener%exists('DIR13')) a%nu13 = Listener%GetRea('DIR13',0.3_rp,'Poisson Coefficient Dir XZ')
      call Listener%Listen('RDOrthotrop')
      if (Listener%exists('POISS').AND.Listener%exists('DIR23')) a%nu23 = Listener%GetRea('DIR23',0.3_rp,'Poisson Coefficient Dir YZ')
      
   end subroutine

   subroutine SPOrthotrop(a)
      use MPI
      implicit none
      class(OrthotropCT) :: a
      integer(ip) :: ierr
      call MPI_BCAST(a%E1,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%E2,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%E3,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%G12,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%G13,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%G23,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nu12,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nu13,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%nu23,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   end subroutine

   subroutine CheckData(a)
      implicit none
      class(OrthotropCT) :: a
      
      real(rp) :: raux
      
      if(a%E1 <= 0.0_rp) call runend('Negative value for the Young Modulus Dir X')
      if(a%E2 <= 0.0_rp) call runend('Negative value for the Young Modulus Dir Y')      
      if(a%E3 <= 0.0_rp) call runend('Negative value for the Young Modulus Dir Z')     
   
      if(a%G12 <= 0.0_rp) call runend('Negative value for the Shear Modulus Dir XY')
      if(a%G13 <= 0.0_rp) call runend('Negative value for the Shear Modulus Dir XZ')
      if(a%G23 <= 0.0_rp) call runend('Negative value for the Shear Modulus Dir YZ')
      
      if(abs(a%nu12) >= sqrt(a%E1/a%E2)) call runend('Poisson Coefficient XY out of range (|v12| < (E1/E2)^0.5)')
      if(abs(a%nu13) >= sqrt(a%E1/a%E3)) call runend('Poisson Coefficient XZ out of range (|v13| < (E1/E3)^0.5)')
      if(abs(a%nu23) >= sqrt(a%E2/a%E3)) call runend('Poisson Coefficient YZ Out of range (|v23| < (E2/E3)^0.5)')
      
      raux = 1.0_rp - a%nu12**2.0_rp*a%E2/a%E1 - a%nu13**2.0_rp*a%E3/a%E1 - a%nu23**2.0_rp*a%E3/a%E2 - 2.0_rp*(a%nu12*a%E2/a%E1)*(a%nu23*a%E3/a%E2)*a%nu13
      if(raux <= 0.0_rp) call runend('Negatice value for aux value Orthotrop Material')
      
   end subroutine   
   

   subroutine BuildOrthotrop(a)
      implicit none
      class(OrthotropCT) :: a

      real(rp) :: raux,v21,v31,v32
!      real(rp) :: Eigenvalues(6)

      !Deviatoric part
!      a%CDev(:,1) = a%E/(1+a%nu)*(/ 2.0_rp/3.0_rp,-1.0_rp/3.0_rp,-1.0/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,2) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,2.0_rp/3.0_rp,-1.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,3) = a%E/(1+a%nu)*(/ -1.0_rp/3.0_rp,-1.0_rp/3.0_rp,2.0_rp/3.0_rp,0.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,4) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp,0.0_rp /)
!      a%CDev(:,5) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp,0.0_rp /)
!      a%CDev(:,6) = a%E/(1+a%nu)*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp/2.0_rp /)

!      if (.not. a%OnlyDeviatoricPart) then
         v21 = a%E2/a%E1*a%nu12
         v31 = a%E3/a%E1*a%nu13
         v32 = a%E3/a%E2*a%nu23

         raux = 1.0_rp/(1.0_rp-a%nu12*v21-a%nu23*v32-a%nu13*v31-2.0_rp*v21*v32*a%nu13)
         
         a%C(:,1) = raux*(/ a%E1*(1.0_rp-a%nu23*v32),a%E1*(v21+v31*a%nu23),a%E1*(v31+v21*v32),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,2) = raux*(/ a%E1*(v21+v31*a%nu23),a%E2*(1.0_rp-a%nu13*v31),a%E2*(v32+a%nu12*v31),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,3) = raux*(/ a%E1*(v31+v21*v32),a%E2*(v32+a%nu12*v31),a%E3*(1.0_rp-a%nu12*v21),0.0_rp,0.0_rp,0.0_rp /)
         a%C(:,4) = a%G12*(/ 0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp /)
         a%C(:,5) = a%G13*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp /)
         a%C(:,6) = a%G23*(/ 0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp /)

!      elseif (a%OnlyDeviatoricPart) then
!         a%C = a%CDev
!      endif

      !Volumetric part
!      a%invK = (3*(1-2*a%nu))/a%E
!      a%CSph = a%C - a%CDev

!      call MatrixEigenValues(6,a%C,Eigenvalues)

!      a%CNorm = maxval(abs(Eigenvalues))


   end subroutine

   subroutine GetPointerToCOrthotrop(a,C)
      implicit none
      class(OrthotropCT), target :: a
      real(rp), pointer :: C(:,:)

      C => a%C
   end subroutine

   subroutine GetPointerToCSphericalIsotrop(a,C)
      implicit none
      class(OrthotropCT), target :: a
      real(rp), pointer :: C(:,:)

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetPointerToCDeviatoricIsotrop(a,C)
      implicit none
      class(OrthotropCT), target :: a
      real(rp), pointer :: C(:,:)

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetInverseVolumetricDeformationModulus(a,invK)
      implicit none
      class(OrthotropCT) :: a
      real(rp) :: invK

      call runend('not implmemetnd0')

   end subroutine

   subroutine GetCNorm(a,CNorm)
      implicit none
      class(OrthotropCT) :: a
      real(rp) :: CNorm

      call runend('not implmemetnd0')

   end subroutine

   subroutine ComputeStressNorm(a,Stress,StressNorm)
      implicit none
      class(OrthotropCT), target :: a
      real(rp) :: Stress(:),StressNorm

      real(rp) :: dot1,dot2

      dot1 = dot_product(Stress(1:3),Stress(1:3))
      dot2 = dot_product(Stress(4:6),Stress(4:6))
      StressNorm = sqrt(dot1+2.0_rp*dot2)

   end subroutine

   subroutine ComputeMeanStress(a,Stress,MeanStress)
      implicit none
      class(OrthotropCT), target :: a
      real(rp) :: Stress(:),MeanStress

      meanStress = sum(Stress(1:3))/3.0
   end subroutine
  
   subroutine GetYoungModulusForDirection(a,Direction,E)
      implicit none
      class(OrthotropCT) :: a
      real(rp) :: Direction(:)
      real(rp) :: E

      call runend('not implmemetnd0')
   end subroutine

   subroutine GetPoissonRatio(a,nu)
      implicit none
      class(OrthotropCT) :: a
      real(rp) :: nu

      call runend('not implmemetnd0')
 end subroutine
 
 
    subroutine GetYoungModulusEquivOrthotrop(a,E)
      implicit none
      class(OrthotropCT) :: a
      real(rp) :: E

      E = a%E1
   end subroutine

end module

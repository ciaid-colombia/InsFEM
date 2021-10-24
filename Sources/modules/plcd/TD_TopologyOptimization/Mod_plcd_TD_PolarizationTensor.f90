module Mod_plcd_TD_PolarizationTensor
   use typre
   implicit none


   type, abstract :: PolarizationTensor


contains
      procedure(ComputePTensor0), deferred :: ComputePTensor
      procedure(ApplyPTensor0), deferred :: ApplyPTensor
   end type


   type, extends(PolarizationTensor) :: PlainStressPolarizationTensor
      real(rp) :: PTensor(3,3,2)

contains
   procedure :: ComputePTensor => ComputePTensorPlainStress
   procedure :: ApplyPTensor   => ApplyPTensorPlainStress


   end type
   
   type, extends(PolarizationTensor) :: PolarizationTensor3D
      real(rp) :: PTensor(6,6,2)

contains
   procedure :: ComputePTensor => ComputePTensor3D
   procedure :: ApplyPTensor   => ApplyPTensor3D


   end type
   
   type, extends(PolarizationTensor) :: PolarizationTensor2DLargeStrains
      real(rp) :: PTensor(3,3,2)

contains
   procedure :: ComputePTensor => ComputePTensor2DLS
   procedure :: ApplyPTensor   => ApplyPTensor2DLS


   end type
   
   type, extends(PolarizationTensor) :: PolarizationTensor3DLargeStrains
      real(rp) :: PTensor(6,6,2)

contains
   procedure :: ComputePTensor => ComputePTensor3DLS
   procedure :: ApplyPTensor   => ApplyPTensor3DLS


   end type

   abstract interface
      subroutine ComputePTensor0(a,nu,extgamma)
         import
         implicit none
         class(PolarizationTensor) :: a
         real(rp) :: nu, extgamma

      end subroutine

      subroutine ApplyPTensor0(a,relativeChi,stress,PStress)
         import
         implicit none
         class(PolarizationTensor) :: a
         real(rp) :: stress(:), PStress(:)
         real(rp) :: relativeChi

      end subroutine


   end interface

contains
   subroutine ComputePTensorPlainStress(a,nu,extgamma)
      use Mod_plcd_Isotrop2DPlainStress
      use Mod_plcd_SIMPMaterial
      implicit none
      class(PlainStressPolarizationTensor) :: a
      real(rp) :: nu, extgamma

      type(Isotrop2dPlainStress) :: CTForPolarizationTensor

      real(rp) :: gamma(2)

      real(rp) :: tau1 = 1.0_rp, tau2 = 1.0_rp, tau3 = 1.0_rp
      real(rp) :: alpha, beta

      real(rp) :: muForPolarization(2), lambdaForPolarization(2), EForPolarization(2), nuForPolarization(2)
      real(rp) :: direction(3) = 1.0_rp

      real(rp), pointer :: C(:,:)
      real(rp) :: coeff
      integer(ip) :: i


      alpha = (1.0_rp+nu)/(1.0_rp-nu)
      beta = (3_rp-nu)/(1_rp+nu)
      gamma(1) = extgamma
      gamma(2) = 1/gamma(1)

      do i = 1,2
         coeff = -0.5_rp/(beta*gamma(i)+tau1)
         muForPolarization(i) = coeff*(1+beta)*(tau1-gamma(i))
         !we divide by two so that we can enter in CT
         muForPolarization(i) = muForPolarization(i)*0.5_rp
         lambdaForPolarization(i) = coeff*0.5_rp*(alpha-beta)*(gamma(i)*(gamma(i)-2*tau3)+tau1*tau2)/(alpha*gamma(i)+tau2)

         EForPolarization(i) = muForPolarization(i)*(3*lambdaForPolarization(i)+2*muForPolarization(i))/(lambdaForPolarization(i)+muForPolarization(i))
         nuForPolarization(i) = lambdaForPolarization(i)/(2*(lambdaForPolarization(i)+muForPolarization(i)))

         
         call CTForPolarizationTensor%Initialize(2_ip,0)
         CTForPolarizationTensor%E=EForPolarization(i)
         CTForPolarizationTensor%nu=nuForPolarization(i)

         call CTForPolarizationTensor%Build
         call CTForPolarizationTensor%GetPointerToC(C)

         a%PTensor(:,:,i) = C

         !we multiply by two the \sigma_12 terms so that when we multiply by sigma everything is fine
         !The tensor usually multiplies epsilon, which is not divided by two in the cross terms
         !So we need to multiply them by 2 so that we can obtain the correct result
         a%PTensor(:,3,i) = a%PTensor(:,3,i)*2.0_rp
         
         
         
         
         
    enddo

   end subroutine

   subroutine ApplyPTensorPlainStress(a,relativeChi,stress,PStress)
      class(PlainStressPolarizationTensor) :: a
      real(rp) :: relativeChi
      real(rp) :: stress(:)
      real(rp) :: PStress(:)
      
      real(rp) :: gamma
      
      real(rp) :: dp(3,3)
      real(rp) :: stress2(3,3) = 0.0_rp, Pstress2(3,3)

      
      
      
      
      PStress = relativeChi*matmul(a%PTensor(:,:,1),stress)-(1.0_rp-relativeChi)*matmul(a%PTensor(:,:,2),stress)
      
      !PStress = stress;
      
      !if (plusminus == 1) then
      !   PStress = matmul(a%PTensor(:,:,1),stress)
      !else
      !   PStress = matmul(a%PTensor(:,:,2),stress)
      !end if
      
      
      

       
       !call Compute_dP(relativeChi,0.001_rp,1.0_rp,0.21_rp,0.21_rp,0.001_rp,1.0_rp,2,3,dp)
       
       
       !call Compute_dP(0.001_rp,0.001_rp,1.0_rp,0.21_rp,0.21_rp,0.001_rp,1.0_rp,2,3,dp)
       
       
!       PStress = matmul(dp,stress)
      
      
      
      
   end subroutine
   
   subroutine ComputePTensor3D(a,nu,extgamma)
      use Mod_plcd_IsotropConstitutiveTensor
      use Mod_plcd_SIMPMaterial
      implicit none
      class(PolarizationTensor3D) :: a
      real(rp) :: nu, extgamma

      type(IsotropCT) :: CTForPolarizationTensor

      real(rp) :: gamma(2)

      real(rp) :: alpha1, alpha2, aux1,auxlambdaNoE,auxmuNoE

      real(rp) :: muForPolarization(2), lambdaForPolarization(2), EForPolarization(2), nuForPolarization(2)

      real(rp), pointer :: C(:,:)
      integer(ip) :: i
      
     
      gamma(1) = extgamma
      gamma(2) = 1/gamma(1)

      auxlambdaNoE = nu/((1+nu)*(1-2*nu))
      auxmuNoE     = 1/(2*(1+nu))
      alpha1 = (auxlambdaNoE+auxmuNoE)/auxmuNoE
      alpha2 = (auxlambdaNoE+3*auxmuNoE)/(auxlambdaNoE+auxmuNoE)
      
      
      do i = 1,2
         aux1 = 0.5_rp*(1.0_rp - gamma(i))/(1+gamma(i)*alpha2)
         
         muForPolarization(i) = aux1*(1+alpha2)*0.5_rp
         lambdaForPolarization(i) = aux1*0.5_rp*(alpha1-alpha2)*(1-gamma(i))/(1+gamma(i)*alpha1)
         
         EForPolarization(i) = muForPolarization(i)*(3*lambdaForPolarization(i)+2*muForPolarization(i))/(lambdaForPolarization(i)+muForPolarization(i))
         nuForPolarization(i) = lambdaForPolarization(i)/(2*(lambdaForPolarization(i)+muForPolarization(i)))

         
         call CTForPolarizationTensor%Initialize(2_ip,0)
         CTForPolarizationTensor%E=EForPolarization(i)
         CTForPolarizationTensor%nu=nuForPolarization(i)

         call CTForPolarizationTensor%Build
         call CTForPolarizationTensor%GetPointerToC(C)

         a%PTensor(:,:,i) = -C

         !we multiply by two the \sigma_12 terms so that when we multiply by sigma everything is fine
         !The tensor usually multiplies epsilon, which is not divided by two in the cross terms
         !So we need to multiply them by 2 so that we can obtain the correct result
         a%PTensor(:,4:6,i) = a%PTensor(:,4:6,i)*2.0_rp
         
         
         
         
         
    enddo

   end subroutine

   subroutine ApplyPTensor3D(a,relativeChi,stress,PStress)
      class(PolarizationTensor3D) :: a
      real(rp) :: relativeChi
      real(rp) :: stress(:)
      real(rp) :: PStress(:)
      
      real(rp) :: gamma
      
      real(rp) :: dp(3,3)
      real(rp) :: stress2(3,3) = 0.0_rp, Pstress2(3,3)

      PStress = relativeChi*matmul(a%PTensor(:,:,1),stress)-(1.0_rp-relativeChi)*matmul(a%PTensor(:,:,2),stress)
   end subroutine
   

   subroutine ComputePTensor2DLS(a,nu,extgamma)
      use Mod_plcd_SIMPMaterial
      implicit none
      class(PolarizationTensor2DLargeStrains) :: a
      real(rp) :: nu, extgamma
      integer(ip) :: i,ivoigt
      
     

   end subroutine

   subroutine ApplyPTensor2DLS(a,relativeChi,stress,PStress)
      class(PolarizationTensor2DLargeStrains) :: a
      real(rp) :: relativeChi
      real(rp) :: stress(:)
      real(rp) :: PStress(:)
      
      PStress = stress
      
   end subroutine
   
   subroutine ComputePTensor3DLS(a,nu,extgamma)
      use Mod_plcd_SIMPMaterial
      implicit none
      class(PolarizationTensor3DLargeStrains) :: a
      real(rp) :: nu, extgamma
      integer(ip) :: i,ivoigt
     

   end subroutine

   subroutine ApplyPTensor3DLS(a,relativeChi,stress,PStress)
      class(PolarizationTensor3DLargeStrains) :: a
      real(rp) :: relativeChi
      real(rp) :: stress(:)
      real(rp) :: PStress(:)
      
      PStress = stress
      
   end subroutine
   
   
   
   
   
   

end module







!    SUBROUTINE Compute_dP(gamm, Eminus, Eplus, numinus, nuplus, alphaminus, alphaplus, d, nstre, dP)
!       implicit none
! !       Gamma: from 0 to 1, element density
! !       Eminus y Eplus:  Young
! !       numinus, muplus: Poisson
! !       alphaminus
! 
!       REAL(rp), INTENT(in) :: gamm
!       REAL(rp), INTENT(in) :: Eminus, Eplus, numinus, nuplus, alphaminus, alphaplus
!       INTEGER(ip), INTENT(in) :: d, nstre
!       REAL(rp) :: mu, dmu, kappa, lambda, dlambda, detinv
!       REAL(rp), DIMENSION(nstre,nstre) :: C, dC, IC, P
!       REAL(rp), INTENT(out) ,DIMENSION(nstre,nstre) :: dP
!       integer(ip) :: i
!       
!       
!       
!       mu=(Eminus**2*alphaplus**2*1.5+ &
!          Eplus**2*alphaminus**2*1.5+&
!          Eminus**2*gamm**2*1.5+&
!          Eplus**2*gamm**2*1.5+&
!          Eminus**2*alphaplus**2*nuplus+&
!          Eplus**2*alphaminus**2*numinus+&
!          Eminus**2*gamm**2*nuplus+ &
!          Eplus**2*gamm**2*numinus- &
!          Eminus**2*alphaplus**2*nuplus**2*0.5- &
!          Eplus**2*alphaminus**2*numinus**2*0.5- &
!          Eminus**2*gamm**2*nuplus**2*0.5- &
!          Eplus**2*gamm**2*numinus**2*0.5+ &
!          Eminus*Eplus*alphaplus**2*4.5+ &
!          Eminus*Eplus*alphaminus**2*4.5- &
!          Eminus*Eplus*gamm**2*3.0- &
!          Eminus**2*alphaplus*gamm*3.0- &
!          Eplus**2*alphaminus*gamm*3.0- &
!          Eminus*Eplus*alphaplus**2*nuplus*1.5- &
!          Eminus*Eplus*alphaplus**2*numinus*1.5- &
!          Eminus*Eplus*alphaminus**2*nuplus*1.5- &
!          Eminus*Eplus*alphaminus**2*numinus*1.5- &
!          Eminus*Eplus*gamm**2*nuplus- &
!          Eminus*Eplus*gamm**2*numinus- &
!          Eminus**2*alphaplus*gamm*nuplus*2.0- &
!          Eplus**2*alphaminus*gamm*numinus*2.0+ &
!          Eminus**2*alphaplus*gamm*nuplus**2+ &
!          Eplus**2*alphaminus*gamm*numinus**2- &
!          Eminus*Eplus*alphaplus*alphaminus*12.0+ &
!          Eminus*Eplus*alphaplus*gamm*3.0+ &
!          Eminus*Eplus*alphaminus*gamm*3.0+ &
!          Eminus*Eplus*alphaplus*alphaminus*nuplus*2.0+ &
!          Eminus*Eplus*alphaplus*alphaminus*numinus*2.0+ &
!          Eminus*Eplus*alphaplus*gamm*nuplus+ &
!          Eminus*Eplus*alphaplus*gamm*numinus+ &
!          Eminus*Eplus*alphaminus*gamm*nuplus+ &
!          Eminus*Eplus*alphaminus*gamm*numinus+ &
!          Eminus*Eplus*alphaplus**2*nuplus*numinus*0.5+ &
!          Eminus*Eplus*alphaminus**2*nuplus*numinus*0.5+ &
!          Eminus*Eplus*gamm**2*nuplus*numinus- &
!          Eminus*Eplus*alphaplus*gamm*nuplus*numinus- &
!          Eminus*Eplus*alphaminus*gamm*nuplus*numinus)/ &
!          ( &
!          (alphaplus-alphaminus)* &
!          (Eminus*alphaplus*3.0- &
!          Eminus*alphaminus*9.0+ &
!          Eplus*alphaplus*9.0- &
!          Eplus*alphaminus*3.0+ &
!          Eminus*gamm*6.0- &
!          Eplus*gamm*6.0+ &
!          Eminus*alphaplus*nuplus*2.0+ &
!          Eminus*alphaplus*numinus*3.0- &
!          Eminus*alphaminus*nuplus*6.0+ &
!          Eminus*alphaminus*numinus*3.0- &
!          Eplus*alphaplus*nuplus*3.0+ &
!          Eplus*alphaplus*numinus*6.0- &
!          Eplus*alphaminus*nuplus*3.0- &
!          Eplus*alphaminus*numinus*2.0+ &
!          Eminus*gamm*nuplus*4.0-Eminus*gamm*numinus*6.0+ &
!          Eplus*gamm*nuplus*6.0- &
!          Eplus*gamm*numinus*4.0- &
!          Eminus*alphaplus*nuplus**2+ &
!          Eminus*alphaminus*nuplus**2*3.0- &
!          Eplus*alphaplus*numinus**2*3.0+ &
!          Eplus*alphaminus*numinus**2- &
!          Eminus*gamm*nuplus**2*2.0+ &
!          Eplus*gamm*numinus**2*2.0- &
!          Eminus*alphaplus*nuplus**2*numinus- &
!          Eminus*alphaminus*nuplus**2*numinus+ &
!          Eplus*alphaplus*nuplus*numinus**2+ &
!          Eplus*alphaminus*nuplus*numinus**2+ &
!          Eminus*gamm*nuplus**2*numinus*2.0- &
!          Eplus*gamm*nuplus*numinus**2*2.0+ &
!          Eminus*alphaplus*nuplus*numinus*2.0+ &
!          Eminus*alphaminus*nuplus*numinus*2.0- &
!          Eplus*alphaplus*nuplus*numinus*2.0- &
!          Eplus*alphaminus*nuplus*numinus*2.0- &
!          Eminus*gamm*nuplus*numinus*4.0+ &
!          Eplus*gamm*nuplus*numinus*4.0) &
!          )
!       
!       
!       dmu= (Eminus**2*alphaplus*-3.0- &
!          Eplus**2*alphaminus*3.0+ &
!          Eminus**2*gamm*3.0+ &
!          Eplus**2*gamm*3.0+ &
!          Eminus**2*alphaplus*nuplus**2+ &
!          Eplus**2*alphaminus*numinus**2- &
!          Eminus**2*gamm*nuplus**2- &
!          Eplus**2*gamm*numinus**2+ &
!          Eminus*Eplus*alphaplus*3.0 + &
!          Eminus*Eplus*alphaminus*3.0- &
!          Eminus*Eplus*gamm*6.0- &
!          Eminus**2*alphaplus*nuplus*2.0- &
!          Eplus**2*alphaminus*numinus*2.0+ &
!          Eminus**2*gamm*nuplus*2.0+ &
!          Eplus**2*gamm*numinus*2.0+ &
!          Eminus*Eplus*alphaplus*nuplus+ &
!          Eminus*Eplus*alphaplus*numinus+ &
!          Eminus*Eplus*alphaminus*nuplus+ &
!          Eminus*Eplus*alphaminus*numinus- &
!          Eminus*Eplus*gamm*nuplus*2.0- &
!          Eminus*Eplus*gamm*numinus*2.0- &
!          Eminus*Eplus*alphaplus*nuplus*numinus- &
!          Eminus*Eplus*alphaminus*nuplus*numinus+ &
!          Eminus*Eplus*gamm*nuplus*numinus*2.0)/ &
!          ( &
!          (alphaplus-alphaminus)* &
!          (Eminus*alphaplus*3.0- &
!          Eminus*alphaminus*9.0+ &
!          Eplus*alphaplus*9.0- &
!          Eplus*alphaminus*3.0+ &
!          Eminus*gamm*6.0-Eplus*gamm*6.0+ &
!          Eminus*alphaplus*nuplus*2.0+ &
!          Eminus*alphaplus*numinus*3.0- &
!          Eminus*alphaminus*nuplus*6.0+ &
!          Eminus*alphaminus*numinus*3.0- &
!          Eplus*alphaplus*nuplus*3.0+ &
!          Eplus*alphaplus*numinus*6.0- &
!          Eplus*alphaminus*nuplus*3.0- &
!          Eplus*alphaminus*numinus*2.0+ &
!          Eminus*gamm*nuplus*4.0- &
!          Eminus*gamm*numinus*6.0+ &
!          Eplus*gamm*nuplus*6.0- &
!          Eplus*gamm*numinus*4.0- &
!          Eminus*alphaplus*nuplus**2+ &
!          Eminus*alphaminus*nuplus**2*3.0- &
!          Eplus*alphaplus*numinus**2*3.0+ &
!          Eplus*alphaminus*numinus**2- &
!          Eminus*gamm*nuplus**2*2.0+ &
!          Eplus*gamm*numinus**2*2.0- &
!          Eminus*alphaplus*nuplus**2*numinus- &
!          Eminus*alphaminus*nuplus**2*numinus+ &
!          Eplus*alphaplus*nuplus*numinus**2+ &
!          Eplus*alphaminus*nuplus*numinus**2+ &
!          Eminus*gamm*nuplus**2*numinus*2.0- &
!          Eplus*gamm*nuplus*numinus**2*2.0+ &
!          Eminus*alphaplus*nuplus*numinus*2.0+ &
!          Eminus*alphaminus*nuplus*numinus*2.0- &
!          Eplus*alphaplus*nuplus*numinus*2.0- &
!          Eplus*alphaminus*nuplus*numinus*2.0- &
!          Eminus*gamm*nuplus*numinus*4.0+ &
!          Eplus*gamm*nuplus*numinus*4.0) &
!          )- &
!          ( &
!          (Eminus*6.0-Eplus*6.0+ &
!          Eminus*nuplus*4.0- &
!          Eminus*numinus*6.0+ &
!          Eplus*nuplus*6.0- &
!          Eplus*numinus*4.0- &
!          Eminus*nuplus**2*2.0+ &
!          Eplus*numinus**2*2.0- &
!          Eminus*nuplus*numinus*4.0+ &
!          Eplus*nuplus*numinus*4.0+ &
!          Eminus*nuplus**2*numinus*2.0- &
!          Eplus*nuplus*numinus**2*2.0)* &
!          1.0/ &
!          (Eminus*alphaplus*3.0- &
!          Eminus*alphaminus*9.0+ &
!          Eplus*alphaplus*9.0- &
!          Eplus*alphaminus*3.0+ &
!          Eminus*gamm*6.0- &
!          Eplus*gamm*6.0+ &
!          Eminus*alphaplus*nuplus*2.0+ &
!          Eminus*alphaplus*numinus*3.0- &
!          Eminus*alphaminus*nuplus*6.0+ &
!          Eminus*alphaminus*numinus*3.0- &
!          Eplus*alphaplus*nuplus*3.0+ &
!          Eplus*alphaplus*numinus*6.0- &
!          Eplus*alphaminus*nuplus*3.0- &
!          Eplus*alphaminus*numinus*2.0+ &
!          Eminus*gamm*nuplus*4.0- &
!          Eminus*gamm*numinus*6.0+ &
!          Eplus*gamm*nuplus*6.0- &
!          Eplus*gamm*numinus*4.0- &
!          Eminus*alphaplus*nuplus**2+ &
!          Eminus*alphaminus*nuplus**2*3.0- &
!          Eplus*alphaplus*numinus**2*3.0+ &
!          Eplus*alphaminus*numinus**2- &
!          Eminus*gamm*nuplus**2*2.0+ &
!          Eplus*gamm*numinus**2*2.0- &
!          Eminus*alphaplus*nuplus**2*numinus- &
!          Eminus*alphaminus*nuplus**2*numinus+ &
!          Eplus*alphaplus*nuplus*numinus**2+ &
!          Eplus*alphaminus*nuplus*numinus**2+ &
!          Eminus*gamm*nuplus**2*numinus*2.0- &
!          Eplus*gamm*nuplus*numinus**2*2.0+ &
!          Eminus*alphaplus*nuplus*numinus*2.0+ &
!          Eminus*alphaminus*nuplus*numinus*2.0- &
!          Eplus*alphaplus*nuplus*numinus*2.0- &
!          Eplus*alphaminus*nuplus*numinus*2.0- &
!          Eminus*gamm*nuplus*numinus*4.0+ &
!          Eplus*gamm*nuplus*numinus*4.0)**2* &
!          (Eminus**2*alphaplus**2*3.0+ &
!          Eplus**2*alphaminus**2*3.0+ &
!          Eminus**2*gamm**2*3.0+ &
!          Eplus**2*gamm**2*3.0+ &
!          Eminus**2*alphaplus**2*nuplus*2.0+ &
!          Eplus**2*alphaminus**2*numinus*2.0+ &
!          Eminus**2*gamm**2*nuplus*2.0+Eplus**2*gamm**2*numinus*2.0- &
!          Eminus**2*alphaplus**2*nuplus**2- &
!          Eplus**2*alphaminus**2*numinus**2- &
!          Eminus**2*gamm**2*nuplus**2- &
!          Eplus**2*gamm**2*numinus**2+ &
!          Eminus*Eplus*alphaplus**2*9.0+ &
!          Eminus*Eplus*alphaminus**2*9.0- &
!          Eminus*Eplus*gamm**2*6.0- &
!          Eminus**2*alphaplus*gamm*6.0- &
!          Eplus**2*alphaminus*gamm*6.0- &
!          Eminus*Eplus*alphaplus**2*nuplus*3.0- &
!          Eminus*Eplus*alphaplus**2*numinus*3.0- &
!          Eminus*Eplus*alphaminus**2*nuplus*3.0- &
!          Eminus*Eplus*alphaminus**2*numinus*3.0- &
!          Eminus*Eplus*gamm**2*nuplus*2.0- &
!          Eminus*Eplus*gamm**2*numinus*2.0- &
!          Eminus**2*alphaplus*gamm*nuplus*4.0- &
!          Eplus**2*alphaminus*gamm*numinus*4.0+ &
!          Eminus**2*alphaplus*gamm*nuplus**2*2.0+ &
!          Eplus**2*alphaminus*gamm*numinus**2*2.0- &
!          Eminus*Eplus*alphaplus*alphaminus*2.4e1+ &
!          Eminus*Eplus*alphaplus*gamm*6.0+ &
!          Eminus*Eplus*alphaminus*gamm*6.0+ &
!          Eminus*Eplus*alphaplus*alphaminus*nuplus*4.0+ &
!          Eminus*Eplus*alphaplus*alphaminus*numinus*4.0+ &
!          Eminus*Eplus*alphaplus*gamm*nuplus*2.0+ &
!          Eminus*Eplus*alphaplus*gamm*numinus*2.0+ &
!          Eminus*Eplus*alphaminus*gamm*nuplus*2.0+ &
!          Eminus*Eplus*alphaminus*gamm*numinus*2.0+ &
!          Eminus*Eplus*alphaplus**2*nuplus*numinus+ &
!          Eminus*Eplus*alphaminus**2*nuplus*numinus+ &
!          Eminus*Eplus*gamm**2*nuplus*numinus*2.0- &
!          Eminus*Eplus*alphaplus*gamm*nuplus*numinus*2.0- &
!          Eminus*Eplus*alphaminus*gamm*nuplus*numinus*2.0)*0.5)/(alphaplus-alphaminus)
!       
!       
!       
!       kappa=( &
!       Eminus**2*alphaplus**2*0.5+ &
!       Eplus**2*alphaminus**2*0.5+ &
!       Eminus**2*gamm**2*0.5+ &
!       Eplus**2*gamm**2*0.5- &
!       Eminus**2*alphaplus**2*nuplus**2*0.5- &
!       Eplus**2*alphaminus**2*numinus**2*0.5- &
!       Eminus**2*gamm**2*nuplus**2*0.5- &
!       Eplus**2*gamm**2*numinus**2*0.5+ &
!       Eminus*Eplus*alphaplus**2*0.5+ &
!       Eminus*Eplus*alphaminus**2*0.5- &
!       Eminus*Eplus*gamm**2- &
!       Eminus**2*alphaplus*gamm- &
!       Eplus**2*alphaminus*gamm+ &
!       Eminus*Eplus*alphaplus**2*nuplus*0.5+ &
!       Eminus*Eplus*alphaplus**2*numinus*0.5+ &
!       Eminus*Eplus*alphaminus**2*nuplus*0.5+ &
!       Eminus*Eplus*alphaminus**2*numinus*0.5+ &
!       Eminus**2*alphaplus*gamm*nuplus**2+ &
!       Eplus**2*alphaminus*gamm*numinus**2- &
!       Eminus*Eplus*alphaplus*alphaminus*2.0+ &
!       Eminus*Eplus*alphaplus*gamm+Eminus*Eplus*alphaminus*gamm- &
!       Eminus*Eplus*alphaplus*alphaminus*nuplus- &
!       Eminus*Eplus*alphaplus*alphaminus*numinus+ &
!       Eminus*Eplus*alphaplus**2*nuplus*numinus*0.5+ &
!       Eminus*Eplus*alphaminus**2*nuplus*numinus*0.5+ &
!       Eminus*Eplus*gamm**2*nuplus*numinus- &
!       Eminus*Eplus*alphaplus*gamm*nuplus*numinus- &
!       Eminus*Eplus*alphaminus*gamm*nuplus*numinus)/ &
!       ( &
!       (alphaplus-alphaminus)* &
!       (Eminus*alphaplus-Eminus*alphaminus+ &
!       Eplus*alphaplus- &
!       Eplus*alphaminus- &
!       Eminus*alphaplus*numinus- &
!       Eminus*alphaminus*numinus+ &
!       Eplus*alphaplus*nuplus+ &
!       Eplus*alphaminus*nuplus+ &
!       Eminus*gamm*numinus*2.0- &
!       Eplus*gamm*nuplus*2.0- &
!       Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2- &
!       Eplus*alphaplus*numinus**2+ &
!       Eplus*alphaminus*numinus**2+ &
!       Eminus*alphaplus*nuplus**2*numinus+ &
!       Eminus*alphaminus*nuplus**2*numinus- &
!       Eplus*alphaplus*nuplus*numinus**2- &
!       Eplus*alphaminus*nuplus*numinus**2- &
!       Eminus*gamm*nuplus**2*numinus*2.0+ &
!       Eplus*gamm*nuplus*numinus**2*2.0) &
!       )
!       
!       
!       lambda=( &
!       Eminus**2*alphaplus**2*0.5+ &
!       Eplus**2*alphaminus**2*0.5+ &
!       Eminus**2*gamm**2*0.5+ &
!       Eplus**2*gamm**2*0.5- &
!       Eminus**2*alphaplus**2*nuplus**2*0.5- &
!       Eplus**2*alphaminus**2*numinus**2*0.5- &
!       Eminus**2*gamm**2*nuplus**2*0.5- &
!       Eplus**2*gamm**2*numinus**2*0.5+ &
!       Eminus*Eplus*alphaplus**2*0.5+ &
!       Eminus*Eplus*alphaminus**2*0.5- &
!       Eminus*Eplus*gamm**2- &
!       Eminus**2*alphaplus*gamm- &
!       Eplus**2*alphaminus*gamm+ &
!       Eminus*Eplus*alphaplus**2*nuplus*0.5+ &
!       Eminus*Eplus*alphaplus**2*numinus*0.5+ &
!       Eminus*Eplus*alphaminus**2*nuplus*0.5+ &
!       Eminus*Eplus*alphaminus**2*numinus*0.5+ &
!       Eminus**2*alphaplus*gamm*nuplus**2+ &
!       Eplus**2*alphaminus*gamm*numinus**2- &
!       Eminus*Eplus*alphaplus*alphaminus*2.0+ &
!       Eminus*Eplus*alphaplus*gamm+ &
!       Eminus*Eplus*alphaminus*gamm- &
!       Eminus*Eplus*alphaplus*alphaminus*nuplus- &
!       Eminus*Eplus*alphaplus*alphaminus*numinus+ &
!       Eminus*Eplus*alphaplus**2*nuplus*numinus*0.5+ &
!       Eminus*Eplus*alphaminus**2*nuplus*numinus*0.5+ &
!       Eminus*Eplus*gamm**2*nuplus*numinus- &
!       Eminus*Eplus*alphaplus*gamm*nuplus*numinus-Eminus*Eplus*alphaminus*gamm*nuplus*numinus)/ &
!       ( &
!       (alphaplus-alphaminus)* &
!       (Eminus*alphaplus- &
!       Eminus*alphaminus+ &
!       Eplus*alphaplus- &
!       Eplus*alphaminus- &
!       Eminus*alphaplus*numinus- &
!       Eminus*alphaminus*numinus+ &
!       Eplus*alphaplus*nuplus+ &
!       Eplus*alphaminus*nuplus+ &
!       Eminus*gamm*numinus*2.0- &
!       Eplus*gamm*nuplus*2.0- &
!       Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2- &
!       Eplus*alphaplus*numinus**2+ &
!       Eplus*alphaminus*numinus**2+ &
!       Eminus*alphaplus*nuplus**2*numinus+ &
!       Eminus*alphaminus*nuplus**2*numinus- &
!       Eplus*alphaplus*nuplus*numinus**2- &
!       Eplus*alphaminus*nuplus*numinus**2- &
!       Eminus*gamm*nuplus**2*numinus*2.0+ &
!       Eplus*gamm*nuplus*numinus**2*2.0) &
!       )- &
!       (Eminus**2*alphaplus**2*3.0+ &
!       Eplus**2*alphaminus**2*3.0+ &
!       Eminus**2*gamm**2*3.0+ &
!       Eplus**2*gamm**2*3.0+ &
!       Eminus**2*alphaplus**2*nuplus*2.0+ &
!       Eplus**2*alphaminus**2*numinus*2.0+ &
!       Eminus**2*gamm**2*nuplus*2.0+ &
!       Eplus**2*gamm**2*numinus*2.0- &
!       Eminus**2*alphaplus**2*nuplus**2- &
!       Eplus**2*alphaminus**2*numinus**2- &
!       Eminus**2*gamm**2*nuplus**2-Eplus**2*gamm**2*numinus**2+ &
!       Eminus*Eplus*alphaplus**2*9.0+ &
!       Eminus*Eplus*alphaminus**2*9.0- &
!       Eminus*Eplus*gamm**2*6.0-Eminus**2*alphaplus*gamm*6.0- &
!       Eplus**2*alphaminus*gamm*6.0- &
!       Eminus*Eplus*alphaplus**2*nuplus*3.0- &
!       Eminus*Eplus*alphaplus**2*numinus*3.0- &
!       Eminus*Eplus*alphaminus**2*nuplus*3.0- &
!       Eminus*Eplus*alphaminus**2*numinus*3.0- &
!       Eminus*Eplus*gamm**2*nuplus*2.0- &
!       Eminus*Eplus*gamm**2*numinus*2.0- &
!       Eminus**2*alphaplus*gamm*nuplus*4.0- &
!       Eplus**2*alphaminus*gamm*numinus*4.0+ &
!       Eminus**2*alphaplus*gamm*nuplus**2*2.0+ &
!       Eplus**2*alphaminus*gamm*numinus**2*2.0- &
!       Eminus*Eplus*alphaplus*alphaminus*2.4e1+ &
!       Eminus*Eplus*alphaplus*gamm*6.0+ &
!       Eminus*Eplus*alphaminus*gamm*6.0+ &
!       Eminus*Eplus*alphaplus*alphaminus*nuplus*4.0+ &
!       Eminus*Eplus*alphaplus*alphaminus*numinus*4.0+ &
!       Eminus*Eplus*alphaplus*gamm*nuplus*2.0+ &
!       Eminus*Eplus*alphaplus*gamm*numinus*2.0+ &
!       Eminus*Eplus*alphaminus*gamm*nuplus*2.0+ &
!       Eminus*Eplus*alphaminus*gamm*numinus*2.0+ &
!       Eminus*Eplus*alphaplus**2*nuplus*numinus+ &
!       Eminus*Eplus*alphaminus**2*nuplus*numinus+ &
!       Eminus*Eplus*gamm**2*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaplus*gamm*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaminus*gamm*nuplus*numinus*2.0)/ &
!       (d*(alphaplus-alphaminus)* &
!       (Eminus*alphaplus*3.0- &
!       Eminus*alphaminus*9.0+ &
!       Eplus*alphaplus*9.0- &
!       Eplus*alphaminus*3.0+ &
!       Eminus*gamm*6.0- &
!       Eplus*gamm*6.0+ &
!       Eminus*alphaplus*nuplus*2.0+ &
!       Eminus*alphaplus*numinus*3.0- &
!       Eminus*alphaminus*nuplus*6.0+ &
!       Eminus*alphaminus*numinus*3.0- &
!       Eplus*alphaplus*nuplus*3.0+ &
!       Eplus*alphaplus*numinus*6.0- &
!       Eplus*alphaminus*nuplus*3.0- &
!       Eplus*alphaminus*numinus*2.0+ &
!       Eminus*gamm*nuplus*4.0- &
!       Eminus*gamm*numinus*6.0+ &
!       Eplus*gamm*nuplus*6.0- &
!       Eplus*gamm*numinus*4.0- &
!       Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2*3.0- &
!       Eplus*alphaplus*numinus**2*3.0+ &
!       Eplus*alphaminus*numinus**2- &
!       Eminus*gamm*nuplus**2*2.0+ &
!       Eplus*gamm*numinus**2*2.0- &
!       Eminus*alphaplus*nuplus**2*numinus- &
!       Eminus*alphaminus*nuplus**2*numinus+ &
!       Eplus*alphaplus*nuplus*numinus**2+ &
!       Eplus*alphaminus*nuplus*numinus**2+ &
!       Eminus*gamm*nuplus**2*numinus*2.0- &
!       Eplus*gamm*nuplus*numinus**2*2.0+ &
!       Eminus*alphaplus*nuplus*numinus*2.0+ &
!       Eminus*alphaminus*nuplus*numinus*2.0- &
!       Eplus*alphaplus*nuplus*numinus*2.0- &
!       Eplus*alphaminus*nuplus*numinus*2.0- &
!       Eminus*gamm*nuplus*numinus*4.0+ &
!       Eplus*gamm*nuplus*numinus*4.0) &
!       )
!       
!       
!       dlambda=( &
!       -Eminus**2*alphaplus- &
!       Eplus**2*alphaminus+ &
!       Eminus**2*gamm+ &
!       Eplus**2*gamm+ &
!       Eminus**2*alphaplus*nuplus**2+ &
!       Eplus**2*alphaminus*numinus**2- &
!       Eminus**2*gamm*nuplus**2- &
!       Eplus**2*gamm*numinus**2+ &
!       Eminus*Eplus*alphaplus+ &
!       Eminus*Eplus*alphaminus- &
!       Eminus*Eplus*gamm*2.0- &
!       Eminus*Eplus*alphaplus*nuplus*numinus- &
!       Eminus*Eplus*alphaminus*nuplus*numinus+ &
!       Eminus*Eplus*gamm*nuplus*numinus*2.0)/ &
!       ( &
!       (alphaplus-alphaminus)* &
!       (Eminus*alphaplus- &
!       Eminus*alphaminus+ &
!       Eplus*alphaplus- &
!       Eplus*alphaminus- &
!       Eminus*alphaplus*numinus- &
!       Eminus*alphaminus*numinus+ &
!       Eplus*alphaplus*nuplus+ &
!       Eplus*alphaminus*nuplus+ &
!       Eminus*gamm*numinus*2.0- &
!       Eplus*gamm*nuplus*2.0- &
!       Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2- &
!       Eplus*alphaplus*numinus**2+ &
!       Eplus*alphaminus*numinus**2+ &
!       Eminus*alphaplus*nuplus**2*numinus+ &
!       Eminus*alphaminus*nuplus**2*numinus- &
!       Eplus*alphaplus*nuplus*numinus**2- &
!       Eplus*alphaminus*nuplus*numinus**2- &
!       Eminus*gamm*nuplus**2*numinus*2.0+ &
!       Eplus*gamm*nuplus*numinus**2*2.0) &
!       )- &
!       ( &
!       (Eminus*numinus*2.0- &
!       Eplus*nuplus*2.0- &
!       Eminus*nuplus**2*numinus*2.0+ &
!       Eplus*nuplus*numinus**2*2.0)* &
!       1.0/ &
!       (Eminus*alphaplus- &
!       Eminus*alphaminus+ &
!       Eplus*alphaplus- &
!       Eplus*alphaminus- &
!       Eminus*alphaplus*numinus- &
!       Eminus*alphaminus*numinus+ &
!       Eplus*alphaplus*nuplus+ &
!       Eplus*alphaminus*nuplus+ &
!       Eminus*gamm*numinus*2.0- &
!       Eplus*gamm*nuplus*2.0- &
!       Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2- &
!       Eplus*alphaplus*numinus**2+ &
!       Eplus*alphaminus*numinus**2+ &
!       Eminus*alphaplus*nuplus**2*numinus+ &
!       Eminus*alphaminus*nuplus**2*numinus- &
!       Eplus*alphaplus*nuplus*numinus**2- &
!       Eplus*alphaminus*nuplus*numinus**2- &
!       Eminus*gamm*nuplus**2*numinus*2.0+ &
!       Eplus*gamm*nuplus*numinus**2*2.0)**2* &
!       (Eminus**2*alphaplus**2+ &
!       Eplus**2*alphaminus**2+ &
!       Eminus**2*gamm**2+ &
!       Eplus**2*gamm**2- &
!       Eminus**2*alphaplus**2*nuplus**2- &
!       Eplus**2*alphaminus**2*numinus**2- &
!       Eminus**2*gamm**2*nuplus**2- &
!       Eplus**2*gamm**2*numinus**2+ &
!       Eminus*Eplus*alphaplus**2+ &
!       Eminus*Eplus*alphaminus**2- &
!       Eminus*Eplus*gamm**2*2.0- &
!       Eminus**2*alphaplus*gamm*2.0- &
!       Eplus**2*alphaminus*gamm*2.0+ &
!       Eminus*Eplus*alphaplus**2*nuplus+ &
!       Eminus*Eplus*alphaplus**2*numinus+ &
!       Eminus*Eplus*alphaminus**2*nuplus+ &
!       Eminus*Eplus*alphaminus**2*numinus+ &
!       Eminus**2*alphaplus*gamm*nuplus**2*2.0+ &
!       Eplus**2*alphaminus*gamm*numinus**2*2.0- &
!       Eminus*Eplus*alphaplus*alphaminus*4.0+ &
!       Eminus*Eplus*alphaplus*gamm*2.0+ &
!       Eminus*Eplus*alphaminus*gamm*2.0- &
!       Eminus*Eplus*alphaplus*alphaminus*nuplus*2.0- &
!       Eminus*Eplus*alphaplus*alphaminus*numinus*2.0+ &
!       Eminus*Eplus*alphaplus**2*nuplus*numinus+ &
!       Eminus*Eplus*alphaminus**2*nuplus*numinus+ &
!       Eminus*Eplus*gamm**2*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaplus*gamm*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaminus*gamm*nuplus*numinus*2.0)*0.5)/ &
!       (alphaplus-alphaminus)- &
!       (Eminus**2*alphaplus*-6.0- &
!       Eplus**2*alphaminus*6.0+ &
!       Eminus**2*gamm*6.0+ &
!       Eplus**2*gamm*6.0+ &
!       Eminus**2*alphaplus*nuplus**2*2.0+ &
!       Eplus**2*alphaminus*numinus**2*2.0- &
!       Eminus**2*gamm*nuplus**2*2.0- &
!       Eplus**2*gamm*numinus**2*2.0+ &
!       Eminus*Eplus*alphaplus*6.0+ &
!       Eminus*Eplus*alphaminus*6.0- &
!       Eminus*Eplus*gamm*12- &
!       Eminus**2*alphaplus*nuplus*4.0- &
!       Eplus**2*alphaminus*numinus*4.0+ &
!       Eminus**2*gamm*nuplus*4.0+ &
!       Eplus**2*gamm*numinus*4.0+ &
!       Eminus*Eplus*alphaplus*nuplus*2.0+ &
!       Eminus*Eplus*alphaplus*numinus*2.0+ &
!       Eminus*Eplus*alphaminus*nuplus*2.0+ &
!       Eminus*Eplus*alphaminus*numinus*2.0- &
!       Eminus*Eplus*gamm*nuplus*4.0- &
!       Eminus*Eplus*gamm*numinus*4.0- &
!       Eminus*Eplus*alphaplus*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaminus*nuplus*numinus*2.0+ &
!       Eminus*Eplus*gamm*nuplus*numinus*4.0)/ &
!       ( &
!       d*(alphaplus-alphaminus)*(Eminus*alphaplus*3.0- &
!       Eminus*alphaminus*9.0+Eplus*alphaplus*9.0- &
!       Eplus*alphaminus*3.0+Eminus*gamm*6.0- &
!       Eplus*gamm*6.0+Eminus*alphaplus*nuplus*2.0+ &
!       Eminus*alphaplus*numinus*3.0-Eminus*alphaminus*nuplus*6.0+ &
!       Eminus*alphaminus*numinus*3.0-Eplus*alphaplus*nuplus*3.0+ &
!       Eplus*alphaplus*numinus*6.0-Eplus*alphaminus*nuplus*3.0- &
!       Eplus*alphaminus*numinus*2.0+Eminus*gamm*nuplus*4.0- &
!       Eminus*gamm*numinus*6.0+Eplus*gamm*nuplus*6.0- &
!       Eplus*gamm*numinus*4.0-Eminus*alphaplus*nuplus**2+ &
!       Eminus*alphaminus*nuplus**2*3.0-Eplus*alphaplus*numinus**2*3.0+ &
!       Eplus*alphaminus*numinus**2-Eminus*gamm*nuplus**2*2.0+ &
!       Eplus*gamm*numinus**2*2.0-Eminus*alphaplus*nuplus**2*numinus- &
!       Eminus*alphaminus*nuplus**2*numinus+Eplus*alphaplus*nuplus*numinus**2+ &
!       Eplus*alphaminus*nuplus*numinus**2+Eminus*gamm*nuplus**2*numinus*2.0- &
!       Eplus*gamm*nuplus*numinus**2*2.0+Eminus*alphaplus*nuplus*numinus*2.0+ &
!       Eminus*alphaminus*nuplus*numinus*2.0-Eplus*alphaplus*nuplus*numinus*2.0- &
!       Eplus*alphaminus*nuplus*numinus*2.0-Eminus*gamm*nuplus*numinus*4.0+ &
!       Eplus*gamm*nuplus*numinus*4.0))+((Eminus*6.0-Eplus*6.0+Eminus*nuplus*4.0- &
!       Eminus*numinus*6.0+Eplus*nuplus*6.0-Eplus*numinus*4.0-Eminus*nuplus**2*2.0+ &
!       Eplus*numinus**2*2.0-Eminus*nuplus*numinus*4.0+Eplus*nuplus*numinus*4.0+ &
!       Eminus*nuplus**2*numinus*2.0-Eplus*nuplus*numinus**2*2.0)*1.0/ &
!       (Eminus*alphaplus*3.0-Eminus*alphaminus*9.0+Eplus*alphaplus*9.0- &
!       Eplus*alphaminus*3.0+Eminus*gamm*6.0-Eplus*gamm*6.0+Eminus*alphaplus*nuplus*2.0+ &
!       Eminus*alphaplus*numinus*3.0-Eminus*alphaminus*nuplus*6.0+ &
!       Eminus*alphaminus*numinus*3.0-Eplus*alphaplus*nuplus*3.0+ &
!       Eplus*alphaplus*numinus*6.0-Eplus*alphaminus*nuplus*3.0-Eplus*alphaminus*numinus*2.0+ &
!       Eminus*gamm*nuplus*4.0-Eminus*gamm*numinus*6.0+Eplus*gamm*nuplus*6.0- &
!       Eplus*gamm*numinus*4.0-Eminus*alphaplus*nuplus**2+Eminus*alphaminus*nuplus**2*3.0- &
!       Eplus*alphaplus*numinus**2*3.0+Eplus*alphaminus*numinus**2-Eminus*gamm*nuplus**2*2.0+ &
!       Eplus*gamm*numinus**2*2.0-Eminus*alphaplus*nuplus**2*numinus-Eminus*alphaminus*nuplus**2*numinus+ &
!       Eplus*alphaplus*nuplus*numinus**2+Eplus*alphaminus*nuplus*numinus**2+ &
!       Eminus*gamm*nuplus**2*numinus*2.0-Eplus*gamm*nuplus*numinus**2*2.0+ &
!       Eminus*alphaplus*nuplus*numinus*2.0+Eminus*alphaminus*nuplus*numinus*2.0- &
!       Eplus*alphaplus*nuplus*numinus*2.0-Eplus*alphaminus*nuplus*numinus*2.0- &
!       Eminus*gamm*nuplus*numinus*4.0+Eplus*gamm*nuplus*numinus*4.0)**2* &
!       (Eminus**2*alphaplus**2*3.0+Eplus**2*alphaminus**2*3.0+ &
!       Eminus**2*gamm**2*3.0+Eplus**2*gamm**2*3.0+Eminus**2*alphaplus**2*nuplus*2.0+ &
!       Eplus**2*alphaminus**2*numinus*2.0+Eminus**2*gamm**2*nuplus*2.0+ &
!       Eplus**2*gamm**2*numinus*2.0-Eminus**2*alphaplus**2*nuplus**2- &
!       Eplus**2*alphaminus**2*numinus**2-Eminus**2*gamm**2*nuplus**2- &
!       Eplus**2*gamm**2*numinus**2+Eminus*Eplus*alphaplus**2*9.0+ &
!       Eminus*Eplus*alphaminus**2*9.0-Eminus*Eplus*gamm**2*6.0- &
!       Eminus**2*alphaplus*gamm*6.0-Eplus**2*alphaminus*gamm*6.0- &
!       Eminus*Eplus*alphaplus**2*nuplus*3.0-Eminus*Eplus*alphaplus**2*numinus*3.0- &
!       Eminus*Eplus*alphaminus**2*nuplus*3.0-Eminus*Eplus*alphaminus**2*numinus*3.0- &
!       Eminus*Eplus*gamm**2*nuplus*2.0-Eminus*Eplus*gamm**2*numinus*2.0- &
!       Eminus**2*alphaplus*gamm*nuplus*4.0-Eplus**2*alphaminus*gamm*numinus*4.0+ &
!       Eminus**2*alphaplus*gamm*nuplus**2*2.0+Eplus**2*alphaminus*gamm*numinus**2*2.0- &
!       Eminus*Eplus*alphaplus*alphaminus*2.4e1+Eminus*Eplus*alphaplus*gamm*6.0+ &
!       Eminus*Eplus*alphaminus*gamm*6.0+Eminus*Eplus*alphaplus*alphaminus*nuplus*4.0+ &
!       Eminus*Eplus*alphaplus*alphaminus*numinus*4.0+Eminus*Eplus*alphaplus*gamm*nuplus*2.0+ &
!       Eminus*Eplus*alphaplus*gamm*numinus*2.0+Eminus*Eplus*alphaminus*gamm*nuplus*2.0+ &
!       Eminus*Eplus*alphaminus*gamm*numinus*2.0+Eminus*Eplus*alphaplus**2*nuplus*numinus+ &
!       Eminus*Eplus*alphaminus**2*nuplus*numinus+Eminus*Eplus*gamm**2*nuplus*numinus*2.0- &
!       Eminus*Eplus*alphaplus*gamm*nuplus*numinus*2.0-Eminus*Eplus*alphaminus*gamm*nuplus*numinus*2.0))/(d*(alphaplus-alphaminus))
!       
!       C(1,1) = 2*mu+lambda
!       C(1,2) = lambda
!       C(1,3) = 0
!       C(2,1) = lambda
!       C(2,2) = 2*mu+lambda
!       C(2,3) = 0
!       C(3,1) = 0
!       C(3,2) = 0
!       C(3,3) = mu
!       
!       
!       
!       dC(1,1) = 2*dmu+dlambda
!       dC(1,2) = dlambda
!       dC(1,3) = 0
!       dC(2,1) = dlambda
!       dC(2,2) = 2*dmu+dlambda
!       dC(2,3) = 0
!       dC(3,1) = 0
!       dC(3,2) = 0
!       dC(3,3) = dmu
!       
!       
!       
!       ! Calculate the inverse determinant of the matrix C
!       detinv = 1/(C(1,1)*C(2,2)*C(3,3) - C(1,1)*C(2,3)*C(3,2)&
!                   - C(1,2)*C(2,1)*C(3,3) + C(1,2)*C(2,3)*C(3,1)&
!                   + C(1,3)*C(2,1)*C(3,2) - C(1,3)*C(2,2)*C(3,1))
!       
!       ! Calculate the inverse of the matrix
!          IC(1,1) = +detinv * (C(2,2)*C(3,3) - C(2,3)*C(3,2))
!          IC(2,1) = -detinv * (C(2,1)*C(3,3) - C(2,3)*C(3,1))
!          IC(3,1) = +detinv * (C(2,1)*C(3,2) - C(2,2)*C(3,1))
!          IC(1,2) = -detinv * (C(1,2)*C(3,3) - C(1,3)*C(3,2))
!          IC(2,2) = +detinv * (C(1,1)*C(3,3) - C(1,3)*C(3,1))
!          IC(3,2) = -detinv * (C(1,1)*C(3,2) - C(1,2)*C(3,1))
!          IC(1,3) = +detinv * (C(1,2)*C(2,3) - C(1,3)*C(2,2))
!          IC(2,3) = -detinv * (C(1,1)*C(2,3) - C(1,3)*C(2,1))
!          IC(3,3) = +detinv * (C(1,1)*C(2,2) - C(1,2)*C(2,1))
!       
! 
!                   P(:,:)=matmul(IC(:,:),dC(:,:))
! 
!       
!       
!          dP(1,1)=P(3,3)+P(1,2)
!          dP(1,2)=P(1,2)
!          dP(1,3)=0
!          dP(2,1)=P(1,2)
!          dP(2,2)=P(3,3)+P(1,2)
!          dP(2,3)=0
!          dP(3,1)=0
!          dP(3,2)=0
!          dP(3,3)=P(3,3)
!       
!       
!       
!       RETURN
!       END SUBROUTINE
! ! 
! !    
!    
!    
!    
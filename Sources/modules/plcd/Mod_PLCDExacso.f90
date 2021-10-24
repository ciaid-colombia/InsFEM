module Mod_plcdExacso
   use typre
   use Mod_PLCD
  
   implicit none   
   type plcdExacso
   
   ! work variables
   real(rp)    :: p
   
   real(rp)    :: disp(3),coord(3)
   real(rp)    :: grad(3,3)=0.0_rp, gradp(3) = 0.0_rp
   real(rp)    :: hessi(3,3,3)
   real(rp)    :: nn,xvisc,xdens   
   
contains   
      procedure :: ComputeSolution
      procedure :: GetDisplacement
      procedure :: GetPressure
      procedure :: GetForce
   end type
   
contains

   subroutine ComputeSolution(a,ndime,gpcod,php)      
      implicit none
      class(plcdExacso) :: a
      Class(PLCDProblem) :: php 
      integer(ip), intent(in) :: ndime      
      real(rp),    intent(in) :: gpcod(ndime)   
      real(rp) :: pi
      real(rp) :: r,phi,alpha,omega
      real(rp) :: alpsq,alpcu
      real(rp) :: cosphisq,cosphicu
      real(rp) :: sinphisq,sinphicu, invk, phi2,sum1,sum2,sum3,sum4, phi2prima,phi2prima3
      integer(ip) :: idime,jdime,kdime
      
!       if(ndime==3)then
!          call runend('ComputeSolution: exact solution not ready for 3d case')   
!       end if
      

      !Coordinates and initializations.

      a%coord(1:ndime) = gpcod
      
      a%disp = 0.0_rp
      a%p = 0.0_rp
      a%grad = 0.0_rp
      a%gradp = 0.0_rp
      
      !Hessi is stored as follows:
      !First index : u, v, w
      !Second index: derivative with respect to (first one)
      !third index: derivative with respect to (second one)
      a%hessi = 0.0_rp
      
      
      pi= 4.*atan(1.)        
      
      ! Obtain unknowns and derivatives according to the exact solution (Kexac)

      if(php%kfl_exacs==1) then !domain [0,1]x[0,1]
      
         !u = 2xy
         !v = -y**2
         !p = xy
      
         ! velocity
         a%disp(1)     =  2.0_rp*a%coord(1)*a%coord(2)
         a%disp(2)     =  -a%coord(2)*a%coord(2)
         
         !velocity derivatives
         !u component
         a%grad(1,1)    = 2.0_rp*a%coord(2)
         a%grad(1,2)    = 2.0_rp*a%coord(1)
         a%hessi(1,1,1)   =   0.0_rp
         a%hessi(1,2,2)   =   0.0_rp
         a%hessi(1,1,2) =   2.0_rp
         !v component
         a%grad(2,1)    = 0.0_rp
         a%grad(2,2)    = -2*a%coord(2)
         a%hessi(2,1,1)   = 0.0_rp
         a%hessi(2,2,2)   = -2.0_rp
         a%hessi(2,1,2) = 0.0_rp
         !pressure and derivatives
         a%p     = a%coord(1)*a%coord(2)
         a%gradp(1)  =  a%coord(2)
         a%gradp(2)  =  a%coord(1)
         
      elseif (php%kfl_exacs==2) then  !domain [0,1]x[0,1]   
         !u = 0
         !v = 0
         !p = 2*x*x*y
         
         a%p =     2*a%coord(1)*a%coord(1)*a%coord(2)
         a%gradp(1) =  4*a%coord(2)*a%coord(1)
         a%gradp(2) =  2*a%coord(1)*a%coord(1)
       
      elseif (php%kfl_exacs == 3) then
         ! spherical coordinates
         r      = sqrt(a%coord(1)*a%coord(1) + a%coord(2)*a%coord(2))
         phi    = atan(a%coord(2)/a%coord(1))
         if(a%coord(1)<epsilon(0.0_rp)) phi = phi + pi
         if(abs(a%coord(1))<epsilon(0.0_rp)) then
            if(a%coord(2)>epsilon(0.0_rp)) phi = pi/2.0_rp
            if(a%coord(2)<epsilon(0.0_rp)) phi = 3*pi/2.0_rp
         end if  
         
         alpha = 0.544483736782464
         omega  = 4.712388980384690
         
         
         alpsq  = alpha*alpha
         alpcu  = alpha*alpha*alpha
         cosphisq = cos(phi)*cos(phi)
         cosphicu = cos(phi)*cos(phi)*cos(phi)
         sinphisq = sin(phi)*sin(phi)
       
         if(r>epsilon(0.0_rp)) then
            sum1 = sin((1+alpha)*phi)*cos(alpha*omega)/(1+alpha);
            sum2 = -cos((1+alpha)*phi);
            sum3 = sin((alpha-1)*phi)*cos(alpha*omega)/(1-alpha);
            sum4 = cos((alpha-1)*phi);

            phi2 =sum1+sum2+sum3+sum4;
            phi2prima = 2*sin(alpha*phi)*cos(phi) + 2*alpha*cos(alpha*phi)*sin(phi) - 2*cos(alpha*omega)*sin(alpha*phi)*sin(phi)
            
            
            
            phi2prima3 = sin(phi*(alpha - 1))*(alpha - 1)**3 - sin(phi*(alpha + 1))*(alpha + 1)**3 + cos(alpha*omega)*cos(phi*(alpha - 1))*(alpha - 1)**2 - cos(alpha*omega)*cos(phi*(alpha + 1))*(alpha + 1)**2
            
            
            
            a%disp(1) = (r**alpha)*(cos(phi)*phi2prima+(1+alpha)*sin(phi)*phi2);
            a%disp(2) = (r**alpha)*(sin(phi)*phi2prima-(1+alpha)*cos(phi)*phi2);
                     
            a%p = -r**(alpha-1)*((1+alpha)**2*phi2prima+phi2prima3)/(1-alpha);
            
            
            a%grad(1,1) = (alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2.0_rp
            a%grad(1,2) = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) - 6*cos(phi - alpha*phi) + sin(phi + alpha*omega - alpha*phi) + sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) + 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
            a%grad(2,1) = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) + 2*cos(phi - alpha*phi) - 3*sin(phi + alpha*omega - alpha*phi) - 3*sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) - 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
            a%grad(2,2) = -(alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2.0_rp
            a%gradp(1) = -2*alpha*r**(alpha - 2)*(cos(2*phi + alpha*omega - alpha*phi) - 2*sin(phi*(alpha - 2)) + cos(2*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)))
            a%gradp(2) = -2*alpha*r**(alpha - 2)*(sin(2*phi + alpha*omega - alpha*phi) + sin(2*phi - alpha*omega - alpha*phi) - 2*cos(phi*(alpha - 2)) + 2*alpha*cos(phi*(alpha - 2)))
            a%hessi(1,1,1) = (alpha*r**(alpha - 2)*(4*sin(4*phi - alpha*phi) - 2*cos(2*phi + alpha*omega - alpha*phi) - 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(2*phi - alpha*phi) - 6*alpha*sin(4*phi - alpha*phi) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(2*phi - alpha*phi) + 2*alpha**2*sin(4*phi - alpha*phi)))/2.0_rp
            a%hessi(1,1,2) = (alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2.0_rp
            a%hessi(1,2,2) = -(alpha*r**(alpha - 2)*(2*cos(2*phi + alpha*omega - alpha*phi) - 4*sin(phi*(alpha - 4)) - 8*sin(phi*(alpha - 2)) + 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha**2*sin(phi*(alpha - 2)) - 2*alpha**2*sin(phi*(alpha - 4)) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) + 6*alpha*sin(phi*(alpha - 2)) + 6*alpha*sin(phi*(alpha - 4))))/2.0_rp
            a%hessi(2,1,1) = -(alpha*r**(alpha - 2)*(4*sin(2*phi + alpha*omega - alpha*phi) + 4*sin(2*phi - alpha*omega - alpha*phi) - 2*sin(4*phi + alpha*omega - alpha*phi) - 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) + 4*cos(phi*(alpha - 4)) - 2*alpha**2*cos(phi*(alpha - 2)) + 2*alpha**2*cos(phi*(alpha - 4)) + 6*alpha*cos(phi*(alpha - 2)) - 6*alpha*cos(phi*(alpha - 4)) - alpha*sin(2*phi + alpha*omega - alpha*phi) - alpha*sin(2*phi - alpha*omega - alpha*phi) + alpha*sin(4*phi + alpha*omega - alpha*phi) + alpha*sin(4*phi - alpha*omega - alpha*phi)))/2.0_rp
            a%hessi(2,1,2) = (alpha*r**(alpha - 2)*(4*sin(phi*(alpha - 4)) + 2*cos(2*phi + alpha*omega - alpha*phi) + 2*cos(2*phi - alpha*omega - alpha*phi) - 2*cos(4*phi + alpha*omega - alpha*phi) - 2*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(phi*(alpha - 2)) + 2*alpha**2*sin(phi*(alpha - 4)) - alpha*cos(2*phi + alpha*omega - alpha*phi) - alpha*cos(2*phi - alpha*omega - alpha*phi) + alpha*cos(4*phi + alpha*omega - alpha*phi) + alpha*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)) - 6*alpha*sin(phi*(alpha - 4))))/2.0_rp
            a%hessi(2,2,2) = -(alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2.0_rp
         endif
         
      elseif (php%kfl_exacs == 4) then    
         a%disp(1) = +1.0_rp/2.0_rp*a%coord(1)**2.0_rp+a%coord(2)+a%coord(3)+a%coord(1)*a%coord(2)+a%coord(1)*a%coord(3)-a%coord(2)**2.0_rp/2.0_rp-a%coord(3)**2.0_rp/2.0_rp
         a%disp(2) = +1.0_rp/2.0_rp*a%coord(2)**2.0_rp-a%coord(3)-a%coord(1)-a%coord(1)**2.0_rp/2.0_rp+a%coord(2)*a%coord(1)+a%coord(2)*a%coord(3)-a%coord(3)**2.0_rp/2.0_rp 
         a%disp(3) = +1.0_rp/2.0_rp*a%coord(3)**2.0_rp-a%coord(1)+a%coord(2)-a%coord(1)**2.0_rp/2.0_rp+a%coord(3)*a%coord(2)+a%coord(3)*a%coord(1)-a%coord(2)**2.0_rp/2.0_rp 
         
         a%grad = transpose(reshape((/ a%coord(1) + a%coord(2) + a%coord(3), a%coord(1) - a%coord(2) + 1, a%coord(1) - a%coord(3) + 1, &
                                       a%coord(2) - a%coord(1) - 1, a%coord(1) + a%coord(2) + a%coord(3), a%coord(2) - a%coord(3) - 1, &
                                       a%coord(3) - a%coord(1) - 1, a%coord(3) - a%coord(2) + 1, a%coord(1) + a%coord(2) + a%coord(3)/),shape(a%grad)))
                                       
         a%hessi(:,:,1) = transpose(reshape((/ 1, 1, 1, -1, 1, 0, -1, 0, 1/),shape(a%hessi(:,:,1))))
         a%hessi(:,:,2) = transpose(reshape((/ 1, -1, 0, 1, 1, 1, 0, -1, 1/),shape(a%hessi(:,:,2))))
         a%hessi(:,:,3) = transpose(reshape((/ 1, 0, -1, 0, 1, -1, 1, 1, 1/),shape(a%hessi(:,:,3))))
 
         call php%Materials(1)%p%CT%GetInverseVolumetricDeformationModulus(invK)
         if (invk == 0) call runend('Exacso 4 for plcd is to be used with compressible material')

         !Pressure from DivU  
         call PressureFromDivU
         
         
      elseif (php%kfl_exacs == 5) then    
         call php%Materials(1)%p%CT%GetInverseVolumetricDeformationModulus(invK)
         if (invk == 0) call runend('Exacso 5 for plcd is to be used with compressible material')
      
      
         a%p = a%coord(1)**3 + a%coord(2)**3
         a%gradp(1) = 3*a%coord(1)**2
         a%gradp(2) = 3*a%coord(2)**2
      else 
         call runend('ComputeSolution exact solution number not defined')      
         !define your exact solution
      end if
      
      !Cross derivatives (second order) are symmetric
      do idime = 1,ndime
         do jdime = 1,ndime
            do kdime = jdime+1,ndime
               a%hessi(idime,kdime,jdime) = a%hessi(idime,jdime,kdime)
            enddo
         enddo
      enddo
      
      
      if (php%UseUPFormulation) then
         call php%Materials(1)%p%CT%GetInverseVolumetricDeformationModulus(invK)
         
         !Compressible material, cannot fix p
         if (invk /= 0.0_rp) then
            
!             
!             a%p = 0.0_rp
!             do idime = 1,ndime
!                a%p = a%p - (1.0_rp/invk)*(a%grad(idime,idime))
!             enddo   
!             a%gradp(1:ndime) = 0.0_rp
!             do idime = 1,ndime
!                do jdime = 1,ndime
!                   a%gradp(idime) = a%gradp(idime) - (1.0_rp/invk)*a%hessi(jdime,jdime,idime)
!                enddo
!             enddo
         else
            !The one specified by the analytical solution
         endif
      endif

contains
      subroutine PressureFromDivU
         a%p = 0.0_rp
         do idime = 1,ndime
            a%p = a%p - (1.0_rp/invk)*(a%grad(idime,idime))
         enddo   
         a%gradp(1:ndime) = 0.0_rp
         do idime = 1,ndime
            do jdime = 1,ndime
               a%gradp(idime) = a%gradp(idime) - (1.0_rp/invk)*a%hessi(jdime,jdime,idime)
            enddo
         enddo
      end subroutine
      
      
   end subroutine
   
   
   subroutine GetPressure(a,ndime,expre,exprg)
      implicit none
      class(plcdExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: expre,exprg(ndime) 
      
      expre      = a%p
      exprg(1:ndime) = a%gradp(1:ndime)
         
   end subroutine

   subroutine GetDisplacement(a,ndime,exvel,exveg)
      implicit none
      class(plcdExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exvel(ndime),exveg(ndime,ndime) 
      
      exvel(1:ndime) = a%disp(1:ndime)  
      exveg(1:ndime,1:ndime) = a%grad(1:ndime,1:ndime)
   end subroutine


   subroutine GetForce(a,ndime,elext,php)
      implicit none
      class(plcdExacso) :: a
      class(plcdProblem) :: php 
      integer(ip), intent(in)  :: ndime
      real(rp), intent(inout)  :: elext(:)
      real(rp) :: D1x,D2x,D1y,D2y,aux1,aux2,aux3
      integer(ip) :: convection
      
      real(rp), pointer :: C(:,:)
      real(rp) :: divsigma(3), invk
      integer(ip) :: idime
      

      !We assume that C is the same in all the elements
      call php%Materials(1)%p%CT%GetPointerToC(C)

      if (ndime == 2) then
         divsigma(1) = C(1,1)*a%hessi(1,1,1) + C(1,2)*a%hessi(2,1,2) + C(1,3)*(a%hessi(1,1,2)+a%hessi(2,1,1)) & 
                     + C(3,1)*a%hessi(1,1,2) + C(3,2)*a%hessi(2,2,2) + C(3,3)*(a%hessi(2,1,2)+a%hessi(1,2,2))
         divsigma(2) = C(2,1)*a%hessi(1,1,2) + C(2,2)*a%hessi(2,2,2) + C(2,3)*(a%hessi(2,1,2)+a%hessi(1,2,2)) &
                     + C(3,1)*a%hessi(1,1,1) + C(3,2)*a%hessi(2,1,2) + C(3,3)*(a%hessi(1,1,2)+a%hessi(2,1,1))
      elseif (ndime == 3) then
         divsigma(1) = C(1,1)*a%hessi(1,1,1) + C(1,4)*a%hessi(1,1,2) + C(1,5)*a%hessi(1,1,3) + C(1,2)*a%hessi(2,1,2) + C(1,4)*a%hessi(2,1,1) + C(1,6)*a%hessi(2,1,3) + C(1,3)*a%hessi(3,1,3) + C(1,5)*a%hessi(3,1,1) + C(1,6)*a%hessi(3,1,2) + C(4,1)*a%hessi(1,1,2) + C(4,4)*a%hessi(1,2,2) + C(4,5)*a%hessi(1,3,2) + C(4,2)*a%hessi(2,2,2) + C(4,4)*a%hessi(2,1,2) + C(4,6)*a%hessi(2,3,2) + C(4,3)*a%hessi(3,3,2) + C(4,5)*a%hessi(3,1,2) + C(4,6)*a%hessi(3,2,2) + C(5,1)*a%hessi(1,1,3) + C(5,4)*a%hessi(1,2,3) + C(5,5)*a%hessi(1,3,3) + C(5,2)*a%hessi(2,2,3) + C(5,4)*a%hessi(2,1,3) + C(5,6)*a%hessi(2,3,3) + C(5,3)*a%hessi(3,3,3) + C(5,5)*a%hessi(3,1,3) + C(5,6)*a%hessi(3,2,3)
         divsigma(2) = C(4,1)*a%hessi(1,1,1) + C(4,4)*a%hessi(1,1,2) + C(4,5)*a%hessi(1,1,3) + C(4,2)*a%hessi(2,1,2) + C(4,4)*a%hessi(2,1,1) + C(4,6)*a%hessi(2,1,3) + C(4,3)*a%hessi(3,1,3) + C(4,5)*a%hessi(3,1,1) + C(4,6)*a%hessi(3,1,2) + C(2,1)*a%hessi(1,1,2) + C(2,4)*a%hessi(1,2,2) + C(2,5)*a%hessi(1,3,2) + C(2,2)*a%hessi(2,2,2) + C(2,4)*a%hessi(2,1,2) + C(2,6)*a%hessi(2,3,2) + C(2,3)*a%hessi(3,3,2) + C(2,5)*a%hessi(3,1,2) + C(2,6)*a%hessi(3,2,2) + C(6,1)*a%hessi(1,1,3) + C(6,4)*a%hessi(1,2,3) + C(6,5)*a%hessi(1,3,3) + C(6,2)*a%hessi(2,2,3) + C(6,4)*a%hessi(2,1,3) + C(6,6)*a%hessi(2,3,3) + C(6,3)*a%hessi(3,3,3) + C(6,5)*a%hessi(3,1,3) + C(6,6)*a%hessi(3,2,3)
         divsigma(3) = C(5,1)*a%hessi(1,1,1) + C(5,4)*a%hessi(1,1,2) + C(5,5)*a%hessi(1,1,3) + C(5,2)*a%hessi(2,1,2) + C(5,4)*a%hessi(2,1,1) + C(5,6)*a%hessi(2,1,3) + C(5,3)*a%hessi(3,1,3) + C(5,5)*a%hessi(3,1,1) + C(5,6)*a%hessi(3,1,2) + C(6,1)*a%hessi(1,1,2) + C(6,4)*a%hessi(1,2,2) + C(6,5)*a%hessi(1,3,2) + C(6,2)*a%hessi(2,2,2) + C(6,4)*a%hessi(2,1,2) + C(6,6)*a%hessi(2,3,2) + C(6,3)*a%hessi(3,3,2) + C(6,5)*a%hessi(3,1,2) + C(6,6)*a%hessi(3,2,2) + C(3,1)*a%hessi(1,1,3) + C(3,4)*a%hessi(1,2,3) + C(3,5)*a%hessi(1,3,3) + C(3,2)*a%hessi(2,2,3) + C(3,4)*a%hessi(2,1,3) + C(3,6)*a%hessi(2,3,3) + C(3,3)*a%hessi(3,3,3) + C(3,5)*a%hessi(3,1,3) + C(3,6)*a%hessi(3,2,3)
      endif
      
      !Force term
      elext(1:ndime) = elext(1:ndime) - divsigma(1:ndime)
      
      if (php%UseUPFormulation) then
         
         elext(1:ndime) = elext(1:ndime) + a%gradp(1:ndime)
         
         !if necessary, add the source term in the pressure equation
         call php%Materials(1)%p%CT%GetInverseVolumetricDeformationModulus(invK)
         
         !Compressible material, cannot fix p
         if (invk /= 0.0_rp) then
            !Force term, volumetric part
            elext(ndime+1) = +invk*a%p
            do idime = 1,ndime
               elext(ndime+1) = elext(ndime+1) + a%grad(idime,idime)
            enddo
         endif
      endif
   
   end subroutine


end module

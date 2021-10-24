module Mod_NsiExacso
   use typre
   use Mod_NavierStokes
  
   implicit none   
   type NsiExacso
   
   ! work variables
   real(rp)    :: x,y,u,v,p
   real(rp)    :: dudx,dudy,dvdx,dvdy,dpdx,dpdy,d2udx,d2udy,d2vdx,d2vdy,d2udxdy,d2vdxdy
   real(rp)    :: dudt= 0.0_rp,dvdt=0.0_rp,dpdt= 0.0_rp
   real(rp)    :: nn,xvisc,xdens 
   real(rp)    :: rho,c2snd
   
contains   
      
      procedure :: nsi_ComputeSolution
      procedure :: nsi_GetVelocity
      procedure :: nsi_GetPressure
      procedure :: nsi_GetForce
   
   end type
   
contains

   subroutine nsi_ComputeSolution(a,ndime,gpcod,wtime,php)      
      use typre
      implicit none
      class(NsiExacso) :: a
      class(NavierStokesProblem) :: php 
      integer(ip), intent(in) :: ndime      
      real(rp),    intent(in) :: gpcod(ndime),wtime   
      real(rp) :: pi
      real(rp) :: r,phi,alpha,omega,phi2,phi2prima,phi2prima3,sum1,sum2,sum3,sum4
      real(rp) :: alpsq,alpcu
      real(rp) :: cosphisq,cosphicu
      real(rp) :: sinphisq,sinphicu
      real(rp) :: normv,machn,auxvl,vsond
      
      if(ndime==3)then
         call runend('nsi_ComputeSolution: exact solution not ready for 3d case')   
      end if
      

      !Coordinates and initializations.

      a%x = gpcod(1)
      a%y = gpcod(2)
         
      a%u =     0.0_rp
      a%v =     0.0_rp
      a%dudx =  0.0_rp
      a%dudy =  0.0_rp
      a%dvdx =  0.0_rp
      a%dvdy =  0.0_rp
      a%d2udx = 0.0_rp
      a%d2udy = 0.0_rp
      a%d2vdx = 0.0_rp
      a%d2vdy = 0.0_rp
      a%d2udxdy = 0.0_rp
      a%d2vdxdy = 0.0_rp

      a%p =     0.0_rp
      a%dpdx =  0.0_rp 
      a%dpdy =  0.0_rp

      !Time derivatives
      a%dudt = 0.0_rp
      a%dvdt = 0.0_rp

      a%xdens=php%MatProp(1)%densi
      
      if(php%MatProp(1)%lawvi==0)then
         a%xvisc = php%MatProp(1)%visco
         a%nn    = 1.0_rp      
      elseif(php%MatProp(1)%lawvi==1)then
         a%xvisc = php%MatProp(1)%LawViParam(1)
         a%nn    = php%MatProp(1)%LawViParam(2)
      elseif(php%MatProp(1)%lawvi==2)then
         call runend('nsi_ComputeSolution: exact solution not implemented for Carreau model')    
      end if      

      pi= 4.*atan(1.)        
      
      ! Obtain unknowns and derivatives according to the exact solution (Kexac)

      if(php%kfl_exacs==1) then !domain [0,1]x[0,1]
      
         ! velocity
         a%u     =  2.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)
         a%v     =  -2.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)
         !velocity derivatives
         !u component
         a%dudx    =   4.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(2.0_rp*a%y-1.0_rp)
         a%dudy    =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)
         a%d2udx   =   4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)
         a%d2udy   =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(12.0_rp*a%y-6.0_rp)
         a%d2udxdy =   4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)
         !v component
         a%dvdx    = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(12.0_rp*a%x-6.0_rp)
         a%d2vdy   = -4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)
         a%d2vdxdy = -4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)
         !pressure and derivatives
         a%p     =  sin(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdx  =  2*pi*cos(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdy  =  2*pi*sin(2*pi*a%x)*cos(2*pi*a%y)

      else if(php%kfl_exacs==2) then !domain [0,1]x[0,1] (linear field of velocity and pressure) 
      
         ! velocity
         a%u       =  4.0_rp*a%x + 6 
         a%v       = -4.0_rp*a%y - 6
         !velocity derivatives
         !u component
         a%dudx    =  4.0_rp
         a%dudy    =  0.0_rp
         a%d2udx   =  0.0_rp
         a%d2udy   =  0.0_rp
         a%d2udxdy =  0.0_rp
         !v component
         a%dvdx    =  0.0_rp
         a%dvdy    = -4.0_rp
         a%d2vdx   =  0.0_rp
         a%d2vdy   =  0.0_rp
         a%d2vdxdy =  0.0_rp
         !pressure and derivatives
         a%p     =  a%x 
         a%dpdx  =  1.0_rp
         a%dpdy  =  0.0_rp  

   else if(php%kfl_exacs==10) then !domain ([-1,1]x[-1,1]) \ ([-1,0]x[0,1]) (Stokes with singularity) 
      
         ! spherical coordinates

         r      = sqrt(a%x*a%x + a%y*a%y)
         phi    = atan(a%y/a%x)
         if(a%x<epsilon(0.0_rp)) phi = phi + pi
         if(abs(a%x)<epsilon(0.0_rp)) then
            if(a%y>epsilon(0.0_rp)) phi = pi/2
            if(a%y<epsilon(0.0_rp)) phi = 3*pi/2
         end if  
         alpha  = php%expar(1)
         omega  = php%expar(2)
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
         
      
         
         a%u = (r**alpha)*(cos(phi)*phi2prima+(1+alpha)*sin(phi)*phi2);
         a%v = (r**alpha)*(sin(phi)*phi2prima-(1+alpha)*cos(phi)*phi2);
                  
         a%p = -r**(alpha-1)*((1+alpha)**2*phi2prima+phi2prima3)/(1-alpha);
         
         
         a%dudx = (alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2
         a%dudy = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) - 6*cos(phi - alpha*phi) + sin(phi + alpha*omega - alpha*phi) + sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) + 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
         a%dvdx = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) + 2*cos(phi - alpha*phi) - 3*sin(phi + alpha*omega - alpha*phi) - 3*sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) - 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
         a%dvdy = -(alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2
         a%dpdx = -2*alpha*r**(alpha - 2)*(cos(2*phi + alpha*omega - alpha*phi) - 2*sin(phi*(alpha - 2)) + cos(2*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)))
         a%dpdy = -2*alpha*r**(alpha - 2)*(sin(2*phi + alpha*omega - alpha*phi) + sin(2*phi - alpha*omega - alpha*phi) - 2*cos(phi*(alpha - 2)) + 2*alpha*cos(phi*(alpha - 2)))
         a%d2udx = (alpha*r**(alpha - 2)*(4*sin(4*phi - alpha*phi) - 2*cos(2*phi + alpha*omega - alpha*phi) - 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(2*phi - alpha*phi) - 6*alpha*sin(4*phi - alpha*phi) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(2*phi - alpha*phi) + 2*alpha**2*sin(4*phi - alpha*phi)))/2
         a%d2udxdy = (alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2
         a%d2udy = -(alpha*r**(alpha - 2)*(2*cos(2*phi + alpha*omega - alpha*phi) - 4*sin(phi*(alpha - 4)) - 8*sin(phi*(alpha - 2)) + 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha**2*sin(phi*(alpha - 2)) - 2*alpha**2*sin(phi*(alpha - 4)) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) + 6*alpha*sin(phi*(alpha - 2)) + 6*alpha*sin(phi*(alpha - 4))))/2
         a%d2vdx = -(alpha*r**(alpha - 2)*(4*sin(2*phi + alpha*omega - alpha*phi) + 4*sin(2*phi - alpha*omega - alpha*phi) - 2*sin(4*phi + alpha*omega - alpha*phi) - 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) + 4*cos(phi*(alpha - 4)) - 2*alpha**2*cos(phi*(alpha - 2)) + 2*alpha**2*cos(phi*(alpha - 4)) + 6*alpha*cos(phi*(alpha - 2)) - 6*alpha*cos(phi*(alpha - 4)) - alpha*sin(2*phi + alpha*omega - alpha*phi) - alpha*sin(2*phi - alpha*omega - alpha*phi) + alpha*sin(4*phi + alpha*omega - alpha*phi) + alpha*sin(4*phi - alpha*omega - alpha*phi)))/2
         a%d2vdxdy = (alpha*r**(alpha - 2)*(4*sin(phi*(alpha - 4)) + 2*cos(2*phi + alpha*omega - alpha*phi) + 2*cos(2*phi - alpha*omega - alpha*phi) - 2*cos(4*phi + alpha*omega - alpha*phi) - 2*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(phi*(alpha - 2)) + 2*alpha**2*sin(phi*(alpha - 4)) - alpha*cos(2*phi + alpha*omega - alpha*phi) - alpha*cos(2*phi - alpha*omega - alpha*phi) + alpha*cos(4*phi + alpha*omega - alpha*phi) + alpha*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)) - 6*alpha*sin(phi*(alpha - 4))))/2
         a%d2vdy = -(alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2

         end if

      else if(php%kfl_exacs==21)then !Temporal case, art 33 Ramon
         !Velocity
         a%u =   20000*a%x**2*a%y*exp(-wtime)*cos(pi*wtime)*(a%x - 1)**2*(2*a%y**2 - 3*a%y + 1)
         a%v =  -20000*a%x*a%y**2*exp(-wtime)*cos(pi*wtime)*(a%y - 1)**2*(2*a%x**2 - 3*a%x + 1)
         !Velocita%y derivatives
         a%dudt = -20000*a%x**2*a%y*exp(-wtime)*(cos(pi*wtime) + pi*sin(pi*wtime))*(a%x - 1)**2*(2*a%y**2 - 3*a%y + 1) 
         a%dvdt = 20000*a%x*a%y**2*exp(-wtime)*(cos(pi*wtime) + pi*sin(pi*wtime))*(a%y - 1)**2*(2*a%x**2 - 3*a%x + 1) 
         !U component
         a%dudx      =  40000*a%x*a%y*exp(-wtime)*cos(pi*wtime)*(2*a%x**2 - 3*a%x + 1)*(2*a%y**2 - 3*a%y + 1)
         a%dudy      = 20000*a%x**2*exp(-wtime)*cos(pi*wtime)*(a%x - 1)**2*(6*a%y**2 - 6*a%y + 1)
         a%d2udx     = 40000*a%y*exp(-wtime)*cos(pi*wtime)*(6*a%x**2 - 6*a%x + 1)*(2*a%y**2 - 3*a%y + 1)      
         a%d2udy     = 20000*a%x**2*exp(-wtime)*cos(pi*wtime)*(12*a%y - 6)*(a%x - 1)**2     
         a%d2udxdy   =  40000*a%x*exp(-wtime)*cos(pi*wtime)*(2*a%x**2 - 3*a%x + 1)*(6*a%y**2 - 6*a%y + 1)    
         !V component
         a%dvdx      = -20000*a%y**2*exp(-wtime)*cos(pi*wtime)*(a%y - 1)**2*(6*a%x**2 - 6*a%x + 1)
         a%dvdy      = -40000*a%x*a%y*exp(-wtime)*cos(pi*wtime)*(2*a%x**2 - 3*a%x + 1)*(2*a%y**2 - 3*a%y + 1)
         a%d2vdx     = -20000*a%y**2*exp(-wtime)*cos(pi*wtime)*(12*a%x - 6)*(a%y - 1)**2
         a%d2vdy     = -40000*a%x*exp(-wtime)*cos(pi*wtime)*(2*a%x**2 - 3*a%x + 1)*(6*a%y**2 - 6*a%y + 1)
         a%d2vdxdy   = -40000*a%y*exp(-wtime)*cos(pi*wtime)*(6*a%x**2 - 6*a%x + 1)*(2*a%y**2 - 3*a%y + 1)

         !Pressure
         a%p         = 1e5*a%x**2
         a%dpdt      = 0
         a%dpdx      = 2e5*a%x 
         a%dpdy      = 0 
         
         !a%p = 0
         !a%dpdt = 0
         !a%dpdx = 0
         !a%dpdy = 0

      else
         call runend('nsi_ComputeSolution exact solution number not defined')      
         !define your exact solution
      end if
   end subroutine
   
   
   subroutine nsi_GetPressure(a,ndime,expre,exprg)
      use typre
      implicit none
      class(NsiExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: expre,exprg(ndime) 
      
         expre      = a%p
         exprg(1)   = a%dpdx
         exprg(2)   = a%dpdy
         
   end subroutine

   subroutine nsi_GetVelocity(a,ndime,exvel,exveg)
      use typre
      implicit none
      class(NsiExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exvel(ndime),exveg(ndime,ndime) 
      
      exvel(1)   = a%u  
      exvel(2)   = a%v
      exveg(1,1) = a%dudx
      exveg(2,1) = a%dvdx
      exveg(1,2) = a%dudy
      exveg(2,2) = a%dvdy
         
   end subroutine


   subroutine nsi_GetForce(a,ndime,elext,php)
      use typre
      implicit none
      class(NsiExacso) :: a
      class(NavierStokesProblem) :: php 
      integer(ip), intent(in)  :: ndime
      real(rp), intent(inout)  :: elext(ndime+1)
      real(rp) :: D1x,D2x,D1y,D2y,aux1,aux2,aux3
      integer(ip) :: convection
      
      aux1 =0.0_rp
      aux2 =0.0_rp
      aux3 =0.0_rp
      D1x  =0.0_rp
      D2x  =0.0_rp
      D1y  =0.0_rp
      D2y  =0.0_rp 
      convection = 1

      !Diffusive term      
      if(php%MatProp(1)%lawvi==0)then
         D1x=a%xvisc*a%d2udx
         D2x=a%xvisc*a%d2udy
         D1y=a%xvisc*a%d2vdx
         D2y=a%xvisc*a%d2vdy               
      elseif(php%MatProp(1)%lawvi > 0)then
      
         aux1=((a%dudx)**2.0_rp + 0.5_rp*(a%dudy)**2.0_rp+a%dudy*a%dvdx+0.5_rp*(a%dvdx)**2.0_rp+(a%dvdy)**2.0_rp)
         aux2=aux1**((a%nn-1.0_rp)/2.0_rp)
         aux3=aux1**((a%nn-3.0_rp)/2.0_rp)
                  
         D1x= 2.0_rp*a%xvisc*a%d2udx*aux2+2.0_rp*a%xvisc*a%dudx*((a%nn-1.0_rp)/2.0_rp)*aux3 &
            *(2.0_rp*a%dudx*a%d2udx+a%dudy*a%d2udxdy+a%dudy*a%d2vdx+a%dvdx*a%d2udxdy+a%dvdx*a%d2vdx+2.0_rp*a%dvdy*a%d2vdxdy)
      
         D2x= a%xvisc*(a%d2udy+a%d2vdxdy)*aux2+a%xvisc*(a%dudy+a%dvdx)*((a%nn-1.0_rp)/2.0_rp)*aux3 &
            *(2.0_rp*a%dudx*a%d2udxdy+a%dudy*a%d2udy+a%dudy*a%d2vdxdy+a%dvdx*a%d2udy+a%dvdx*a%d2vdxdy+2.0_rp*a%dvdy*a%d2vdy)

         D1y= a%xvisc*(a%d2udxdy+a%d2vdx)*aux2+a%xvisc*(a%dudy+a%dvdx)*((a%nn-1.0_rp)/2.0_rp)*aux3 &
            *(2.0_rp*a%dudx*a%d2udx+a%dudy*a%d2udxdy+a%dudy*a%d2vdx+a%dvdx*a%d2udxdy+a%dvdx*a%d2vdx+2.0_rp*a%dvdy*a%d2vdxdy) 

         D2y= 2.0_rp*a%xvisc*a%d2vdy*aux2 + 2.0_rp*a%xvisc*a%dvdy*((a%nn-1.0_rp)/2.0_rp)*aux3 &
            *(2.0_rp*a%dudx*a%d2udxdy+a%dudy*a%d2udy+a%dudy*a%d2vdxdy+a%dvdx*a%d2udy+a%dvdx*a%d2vdxdy+2.0_rp*a%dvdy*a%d2vdy)
      end if
      
      if(php%kfl_advec == 0)then
         convection=0
      end if
      
      !Force term
      elext(1) = elext(1) + a%xdens*a%dudt + convection*(a%xdens*(a%u*a%dudx + a%v*a%dudy)) - (D1x + D2x) +  a%dpdx 
      elext(2) = elext(2) + a%xdens*a%dvdt + convection*(a%xdens*(a%u*a%dvdx + a%v*a%dvdy)) - (D1y + D2y) +  a%dpdy 
      
   end subroutine
end module

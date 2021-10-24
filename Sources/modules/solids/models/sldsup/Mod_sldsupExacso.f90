module Mod_SldSupExacso
   use typre
   use Mod_SUPSolids
  
   implicit none   
   type SldSupExacso
   
   !work variables
   real(rp)    :: x,y,u,v,p
   real(rp)    :: dudx,dudy,dvdx,dvdy,dpdx,dpdy,d2udx,d2udy,d2vdx,d2vdy,d2udxdy,d2vdxdy
   real(rp)    :: Sxx,Syy,Sxy,dsxxdx,dsxxdy,dsyydx,dsyydy,dsxydx,dsxydy,dudt,dvdt,dsxxdt,dsyydt,dsxydt  
   real(rp)    :: Pxx, Pyy, Pxy, dPxxdx,dPxxdy,dPyydx,dPyydy,dPxydx,dPxydy,dpxxdt,dpyydt,dpxydt
   real(rp)    :: dS1dx,dS1dy,dS2dx,dS2dy,dS3dx,dS3dy 
   real(rp)    :: nn,xvisc,xdens,viscnn
   
contains   
      
      procedure :: sup_ComputeSolution
      procedure :: sup_GetVelocity
      procedure :: sup_GetPressure
      procedure :: sup_GetStress
      procedure :: sup_GetForce
   
   end type
   
contains

   subroutine sup_ComputeSolution(a,ndime,gpcod,wtime,php)      
      use typre
      implicit none
      class(SldSupExacso) :: a
      Class(SUPSolidsProblem) :: php 
      integer(ip), intent(in) :: ndime     
      real(rp),    intent(in) :: gpcod(ndime),wtime   
      real(rp) :: pi,tfunc,dtfunc,cdt,cdt2
      real(rp) :: auxdt1,auxdt2,auxdt3,auxdt4,cexac1
      integer(ip) :: auxiter
      
      if(ndime==3)then
         call runend('sldsup_ComputeSolution: exact solution not ready for 3d case')   
      end if
      
      !Coordinates and initializations.
      a%x = gpcod(1)
      a%y = gpcod(2)
      
      if (a%x==0 .or. a%y==0) then
         if (a%x==0) a%x=0.0000000001
         if (a%y==0 )a%y=0.0000000001
      end if
         
      a%u       = 0.0_rp
      a%v       = 0.0_rp
      a%dudx    = 0.0_rp
      a%dudy    = 0.0_rp
      a%dvdx    = 0.0_rp
      a%dvdy    = 0.0_rp
      a%d2udx   = 0.0_rp
      a%d2udy   = 0.0_rp
      a%d2vdx   = 0.0_rp
      a%d2vdy   = 0.0_rp
      a%d2udxdy = 0.0_rp
      a%d2vdxdy = 0.0_rp

      a%p      = 0.0_rp
      a%dpdx   = 0.0_rp 
      a%dpdy   = 0.0_rp
      
      a%Sxx    = 0.0_rp
      a%Sxy    = 0.0_rp
      a%Syy    = 0.0_rp
      a%dsxxdx = 0.0_rp
      a%dsxxdy = 0.0_rp
      a%dsyydx = 0.0_rp
      a%dsyydy = 0.0_rp
      a%dsxydx = 0.0_rp
      a%dsxydy = 0.0_rp
 
      a%dudt   = 0.0_rp
      a%dvdt   = 0.0_rp
      tfunc    = 1.0_rp
      dtfunc   = 0.0_rp
      cdt      = 1.0_rp
      cdt2     = 4.0_rp
      
      a%xdens=a%densi
    
      a%xvisc = php%mu
      a%nn    = 1.0_rp      
      
      pi= 4.*atan(1.)
      
      if(php%kfl_timei==1)then
         auxdt1 = php%ctime      
         tfunc =  cos(cdt2*pi*wtime)*exp(-cdt*wtime)
         dtfunc = -cdt2*pi*sin(cdt2*pi*wtime)*exp(-cdt*wtime) - cdt*cos(cdt2*pi*wtime)*exp(-cdt*wtime) 
      end if
  
      
      ! Obtain unknowns and derivatives according to the exact solution (Kexac)

      if(php%kfl_exacs==1) then !domain [0,1]x[0,1]
      
         ! velocity
         a%u     =  2.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%v     =  -2.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*tfunc
         !velocity derivatives
         !u component
         a%dudx    =   4.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%dudy    =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2udx   =   4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%d2udy   =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(12.0_rp*a%y-6.0_rp)*tfunc
         a%d2udxdy =   4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%dudt    =   a%u*dtfunc   
         !v component
         a%dvdx    = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(12.0_rp*a%x-6.0_rp)*tfunc
         a%d2vdy   = -4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2vdxdy = -4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%dvdt    =  a%v*dtfunc
         !pressure and derivatives
         a%p     =  sin(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdx  =  2*pi*cos(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdy  =  2*pi*sin(2*pi*a%x)*cos(2*pi*a%y)
         !stress and derivatives
         a%Sxx    = 2.0_rp*a%xvisc*a%dudx        
         a%Sxy    = a%xvisc*(a%dudy+a%dvdx)
         a%Syy    = -a%Sxx
         a%dsxxdx = 8.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1)*tfunc
         a%dsxxdy = 8.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1)*tfunc
         a%dsyydx = -a%dsxxdx*tfunc
         a%dsyydy = -a%dsxxdy*tfunc
         a%dsxydx = (-4.0_rp*(2.0_rp*a%x-1)*(6.0_rp*a%x*a%x*a%y-6.0_rp*a%x*a%x*a%y*a%y-a%x*a%x+6.0_rp*a%x*a%y*a%y-6.0_rp*a%x*a%y+a%x &
                     +3.0_rp*a%y*a%y*a%y*a%y-6.0_rp*a%y*a%y*a%y+3.0_rp*a%y*a%y))*tfunc
         a%dsxydy =(-4.0_rp*(2.0_rp*a%y-1)*(6.0_rp*a%x*a%x*a%y-6.0_rp*a%x*a%x*a%y*a%y-a%y*a%y+6.0_rp*a%x*a%y*a%y-6.0_rp*a%x*a%y+a%y &
                     +3.0_rp*a%x*a%x*a%x*a%x-6.0_rp*a%x*a%x*a%x+3.0_rp*a%x*a%x))*tfunc          

      else if(php%kfl_exacs==2) then !domain [0,1]x[0,1] (linear field of velocity and pressure) 
      
         ! velocity
         a%u       =  (4.0_rp*a%x + 6)*tfunc 
         a%v       = (-4.0_rp*a%y - 6)*tfunc
         !velocity derivatives
         a%dudt    = (4.0_rp*a%x + 6)*dtfunc
         a%dvdt    = (-4.0_rp*a%y - 6)*dtfunc
         !u component
         a%dudx    =  4.0_rp*tfunc
         a%dudy    =  0.0_rp
         a%d2udx   =  0.0_rp
         a%d2udy   =  0.0_rp
         a%d2udxdy =  0.0_rp
         !v component
         a%dvdx    =  0.0_rp
         a%dvdy    = (-4.0_rp)*tfunc
         a%d2vdx   =  0.0_rp
         a%d2vdy   =  0.0_rp
         a%d2vdxdy =  0.0_rp 

         !pressure and derivatives
         a%p     =  a%x + a%y
         a%dpdx  =  1.0_rp
         a%dpdy  =  1.0_rp

         a%Sxx    = (2.0_rp*a%x + 3)*tfunc
         a%Sxy    = (a%x + a%y)*tfunc 
         a%Syy    = (2.0_rp*a%y + 3)*tfunc
         a%dsxxdx = 2.0_rp*tfunc
         a%dsxxdy = 0.0_rp
         a%dsyydx = 0.0_rp
         a%dsyydy = 2.0_rp*tfunc
         a%dsxydx = 1.0_rp*tfunc
         a%dsxydy = 1.0_rp*tfunc 
         !temporal derivatives
         a%dsxxdt = (2.0_rp*a%x + 3)*dtfunc
         a%dsyydt = (2.0_rp*a%y + 3)*dtfunc
         a%dsxydt = (a%x + a%y)*dtfunc 

         
      else if(php%kfl_exacs==3) then !domain [0,1]x[0,1] Viscoelastic case
      
         ! velocity
         a%u     =  2.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%v     =  -2.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*tfunc
         !velocity derivatives
         !u component
         a%dudx    =   4.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%dudy    =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2udx   =   4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%d2udy   =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(12.0_rp*a%y-6.0_rp)*tfunc
         a%d2udxdy =   4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         !v component
         a%dvdx    = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(12.0_rp*a%x-6.0_rp)*tfunc
         a%d2vdy   = -4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2vdxdy = -4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         !pressure and derivatives
         a%p     =  sin(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdx  =  2*pi*cos(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdy  =  2*pi*sin(2*pi*a%x)*cos(2*pi*a%y)
         !stress and derivatives 
         a%Sxx    = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%Sxy    = sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc 
         a%Syy    = -a%Sxx
         a%dsxxdx  = 10.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxxdy  = 10.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         a%dsyydx = -a%dsxxdx
         a%dsyydy = -a%dsxxdy
         a%dsxydx  = 2.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxydy  = 2.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         !temporal derivatives
         a%dsxxdt = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsyydt = -5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsxydt = sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc 
         
      else if(php%kfl_exacs==4) then !transitory test 
      
         ! velocity
         cexac1=10.0_rp
         a%u     =  (cexac1*cexac1)*(a%x*a%x*(1.0_rp-a%x)*(1.0_rp-a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y))*tfunc
         a%v     =  -(cexac1*cexac1)*(a%y*a%y*(1.0_rp-a%y)*(1.0_rp-a%y)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x))*tfunc
         !velocity derivatives
         !u component
         a%dudx    =   (cexac1*cexac1)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y)*tfunc
         a%dudy    =   (cexac1*cexac1)*(a%x*a%x*(1.0_rp-a%x)*(1.0_rp-a%x)*(2.0_rp-12.0_rp*a%y+12.0_rp*a%y*a%y))*tfunc
         a%d2udx   =   (cexac1*cexac1)*(2.0_rp-12.0_rp*a%x+12.0_rp*a%x*a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y)*tfunc
         a%d2udy   =   (cexac1*cexac1)*(a%x*a%x*(1.0_rp-a%x)*(1.0_rp-a%x)*(-12.0_rp+24.0_rp*a%y))*tfunc
         a%d2udxdy =   (cexac1*cexac1)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x)*(2.0_rp-12.0_rp*a%y+12.0_rp*a%y*a%y)*tfunc
         !v component
         a%dvdx    = -(cexac1*cexac1)*(a%y*a%y*(1.0_rp-a%y)*(1.0_rp-a%y)*(2.0_rp-12.0_rp*a%x+12.0_rp*a%x*a%x))*tfunc
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -(cexac1*cexac1)*(a%y*a%y*(1.0_rp-a%y)*(1.0_rp-a%y)*(-12.0_rp+24.0_rp*a%x))*tfunc
         a%d2vdy   = -(cexac1*cexac1)*(2.0_rp-12.0_rp*a%y+12.0_rp*a%y*a%y)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x)*tfunc
         a%d2vdxdy = -(cexac1*cexac1)*(2.0_rp-12.0_rp*a%x+12.0_rp*a%x*a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y)*tfunc
         !temporal derivatives
         a%dudt    =  (cexac1*cexac1)*(a%x*a%x*(1.0_rp-a%x)*(1.0_rp-a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y))*dtfunc
         a%dvdt    = -(cexac1*cexac1)*(a%y*a%y*(1.0_rp-a%y)*(1.0_rp-a%y)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x))*dtfunc
         !pressure and derivatives
         a%p     =  0.0_rp !cexac1*a%x*a%x
         a%dpdx  =  0.0_rp !2.0_rp*cexac1*a%x
         a%dpdy  =  0.0_rp
         !stress and derivatives 
         a%Sxx    = 2.0_rp*a%xvisc*a%dudx
         a%Sxy    = 2.0_rp*a%xvisc*0.5_rp*(a%dudy+a%dvdx) 
         a%Syy    = 2.0_rp*a%xvisc*a%dvdy
         a%dsxxdx = 2.0_rp*a%xvisc*a%d2udx 
         a%dsxxdy = 2.0_rp*a%xvisc*a%d2udxdy
         a%dsyydx = 2.0_rp*a%xvisc*a%d2vdxdy
         a%dsyydy = 2.0_rp*a%xvisc*a%d2vdy
         a%dsxydx = 2.0_rp*a%xvisc*0.5_rp*(a%d2udxdy+a%d2vdx)
         a%dsxydy = 2.0_rp*a%xvisc*0.5_rp*(a%d2udy+a%d2vdxdy)
         !temporal derivatives
         a%dsxxdt = 2.0_rp*a%xvisc*(cexac1*cexac1)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y)*dtfunc
         a%dsyydt =-2.0_rp*a%xvisc*(cexac1*cexac1)*(2.0_rp*a%x-6.0_rp*a%x*a%x+4.0_rp*a%x*a%x*a%x)*(2.0_rp*a%y-6.0_rp*a%y*a%y+4.0_rp*a%y*a%y*a%y)*dtfunc
         a%dsxydt = 2.0_rp*a%xvisc*0.5_rp*((cexac1*cexac1)*(a%x*a%x*(1.0_rp-a%x)*(1.0_rp-a%x)*(2.0_rp-12.0_rp*a%y+12.0_rp*a%y*a%y)) &
               -(cexac1*cexac1)*(a%y*a%y*(1.0_rp-a%y)*(1.0_rp-a%y)*(2.0_rp-12.0_rp*a%x+12.0_rp*a%x*a%x)))*dtfunc        
    
      
      else if(php%kfl_exacs==5) then
      
         ! velocity
         a%u     =  200.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%v     =  -200.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*tfunc
         !velocity derivatives
         !u component
         a%dudx    =   400.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%dudy    =   200.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2udx   =   400.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%d2udy   =   200.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(12.0_rp*a%y-6.0_rp)*tfunc
         a%d2udxdy =   400.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         !v component
         a%dvdx    = -200.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -200.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(12.0_rp*a%x-6.0_rp)*tfunc
         a%d2vdy   = -400.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2vdxdy = -400.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         !pressure and derivatives
         a%p     =  sin(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdx  =  2*pi*cos(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdy  =  2*pi*sin(2*pi*a%x)*cos(2*pi*a%y)
         !stress and derivatives
         !stress and derivatives 
         a%Sxx    = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%Sxy    = sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc 
         a%Syy    = -a%Sxx
         a%dsxxdx  = 10.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxxdy  = 10.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         a%dsyydx = -a%dsxxdx
         a%dsyydy = -a%dsxxdy
         a%dsxydx  = 2.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxydy  = 2.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         !temporal derivatives
         a%dsxxdt = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsyydt = -5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsxydt = sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc      
     
     else if(php%kfl_exacs==6)then
      
         ! velocity
         a%u     =  2.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%v     =  -2.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*tfunc
         !velocity derivatives
         !u component
         a%dudx    =   4.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(2.0_rp*a%y-1.0_rp)*tfunc
         a%dudy    =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2udx   =   4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%d2udy   =   2.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x-1.0_rp)*(12.0_rp*a%y-6.0_rp)*tfunc
         a%d2udxdy =   4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         !v component
         a%dvdx    = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
         a%dvdy    = -(a%dudx)
         a%d2vdx   = -2.0_rp*a%y*a%y*(a%y-1.0_rp)*(a%y-1.0_rp)*(12.0_rp*a%x-6.0_rp)*tfunc
         a%d2vdy   = -4.0_rp*a%x*(a%x-1.0_rp)*(2.0_rp*a%x-1.0_rp)*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*tfunc
         a%d2vdxdy = -4.0_rp*a%y*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*tfunc
          
         !temporal derivatives
         a%dudt    =  2.0_rp*a%x*a%x*a%y*(a%x-1.0_rp)*(a%x-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%y-1.0_rp)*dtfunc
         a%dvdt    =  -2.0_rp*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)*dtfunc         

        !pressure and derivatives
         a%p     =  sin(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdx  =  2*pi*cos(2*pi*a%x)*sin(2*pi*a%y)
         a%dpdy  =  2*pi*sin(2*pi*a%x)*cos(2*pi*a%y)
         
         !stress and derivatives 
         a%Sxx    = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%Sxy    = sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc 
         a%Syy    = -a%Sxx
         a%dsxxdx  = 10.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxxdy  = 10.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         a%dsyydx = -a%dsxxdx
         a%dsyydy = -a%dsxxdy
         a%dsxydx  = 2.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dsxydy  = 2.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         !temporal derivatives
         a%dsxxdt = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsyydt = -5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dsxydt = sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc 
            
        end if    
   end subroutine
   
   
   subroutine sup_GetPressure(a,ndime,expre,exprg)
      use typre
      implicit none
      class(SldSupExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: expre,exprg(ndime) 
      real(rp)                 :: point1, point2
      
         point1=a%x
         point2=a%y
         expre      = a%p
         exprg(1)   = a%dpdx
         exprg(2)   = a%dpdy
         
   end subroutine

   subroutine sup_GetVelocity(a,ndime,exvel,exveg)
      use typre
      implicit none
      class(SldSupExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exvel(ndime),exveg(ndime,ndime) 
      
      exvel(1)   = a%u  
      exvel(2)   = a%v
      exveg(1,1) = a%dudx
      exveg(2,1) = a%dvdx
      exveg(1,2) = a%dudy
      exveg(2,2) = a%dvdy
         
   end subroutine

   subroutine sup_GetStress(a,ndime,tn,exsig,exsigr)
      use typre
      implicit none
      class(SupExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exsig(tn),exsigr(tn,ndime) 
      real(rp)                  :: sig(tn),tau((ndime-1)*(ndime-1)+2)
      real(rp)                  :: psir(((ndime-1)*(ndime-1)+2),ndime),aux2,taugr(((ndime-1)*(ndime-1)+2),ndime)
      real(rp)                  :: sigr(((ndime-1)*(ndime-1)+2),ndime), point1,point2
      
      
      point1=a%x
      point2=a%y
      exsig(1)    = a%Sxx
      exsig(2)    = a%Syy
      exsig(3)    = a%Sxy
      exsigr(1,1) = a%dsxxdx
      exsigr(2,1) = a%dsyydx
      exsigr(3,1) = a%dsxydx
      exsigr(1,2) = a%dsxxdy
      exsigr(2,2) = a%dsyydy
      exsigr(3,2) = a%dsxydy
            
   end subroutine   
   
   subroutine sup_GetForce(a,ndime,elext,elextC,elextS,elextEstab,elextSEstab,elextSMat,elextSEstabMat,php)
      use typre
      implicit none
      class(SupExacso) :: a
      class(ThreeFieldNSProblem) :: php 
      integer(ip), intent(in)  :: ndime
      real(rp), intent(inout)  :: elext(ndime),elextC(1),elextS(tn),elextSEstab(tn)
      real(rp), intent(inout)  :: elextEstab(ndime)
      real(rp), intent(inout)  :: elextSMat(ndime,ndime),elextSEstabMat(ndime,ndime)
      real(rp) :: D1x,D2x,D1y,D2y,aux1,aux2,aux3,aux4,dudtforce,dvdtforce,dsxxdtforce,dsyydtforce,dsxydtforce
      
      aux1 =0.0_rp
      aux2 =0.0_rp
      aux3 =0.0_rp
      D1x  =0.0_rp
      D2x  =0.0_rp
      D1y  =0.0_rp
      D2y  =0.0_rp
      !temporal derivatives
      dudtforce   = 0.0_rp
      dvdtforce   = 0.0_rp
      dsxxdtforce = 0.0_rp
      dsyydtforce = 0.0_rp
      dsxydtforce = 0.0_rp
            
      if(php%kfl_timei==1)then
         dudtforce = a%dudt*a%xdens
         dvdtforce = a%dvdt*a%xdens
      endif
      
         D1x= a%xvisc*a%d2udx 
         D2x= a%xvisc*a%d2udy
         D1y= a%xvisc*a%d2vdx 
         D2y= a%xvisc*a%d2vdy  

         !Constitutive force term  
         elextS(1)= (1.0_rp/(2.0_rp*a%xvisc))*a%Sxx- a%dudx
         elextS(2)= (1.0_rp/(2.0_rp*a%xvisc))*a%Syy- a%dvdy 
         elextS(3)= (1.0_rp/(2.0_rp*a%xvisc))*a%Sxy- (1.0_rp/(2.0_rp))*(a%dudy + a%dvdx)            
         
   end subroutine

end module

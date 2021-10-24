module Mod_SupExacso
   use typre
   use Mod_ThreeField
   use Mod_LogOperations
  
   implicit none   
   type SupExacso
   
   !work variables
   real(rp)    :: x,y,u,v,p
   real(rp)    :: dudx,dudy,dvdx,dvdy,dpdx,dpdy,d2udx,d2udy,d2vdx,d2vdy,d2udxdy,d2vdxdy
   real(rp)    :: Sxx,Syy,Sxy,dsxxdx,dsxxdy,dsyydx,dsyydy,dsxydx,dsxydy,dudt,dvdt,dsxxdt,dsyydt,dsxydt  
   real(rp)    :: Pxx, Pyy, Pxy, dPxxdx,dPxxdy,dPyydx,dPyydy,dPxydx,dPxydy,dpxxdt,dpyydt,dpxydt
   real(rp)    :: dS1dx,dS1dy,dS2dx,dS2dy,dS3dx,dS3dy 
   real(rp)    :: nn,xvisc,xdens,beta,lamb,auxNL,viscnn, auxLCR, lamb0 ,auxL_inv,auxL, auxG
   real(rp)    :: expoSyy, expoSxy, expoSxx
   
contains   
      
      procedure :: sup_ComputeSolution
      procedure :: sup_GetVelocity
      procedure :: sup_GetPressure
      procedure :: sup_GetStress
      procedure :: sup_GetPsi
      procedure :: sup_GetForce
      procedure :: sup_GetExponentialDivergence
      procedure :: sup_GetExponentialGradient
      procedure :: sup_GetExponential

   
   end type
   
contains

   subroutine sup_ComputeSolution(a,ndime,gpcod,wtime,LCR,php)      
      use typre
      implicit none
      class(SupExacso) :: a
      Class(ThreeFieldNSProblem) :: php 
      integer(ip), intent(in) :: ndime     
      real(rp),    intent(in) :: gpcod(ndime),wtime   
      real(rp) :: pi,tfunc,dtfunc,cdt,cdt2
      real(rp) :: auxdt1,auxdt2,auxdt3,auxdt4,cexac1
      integer(ip) :: auxiter,LCR
      
      if(ndime==3)then
         call runend('sup_ComputeSolution: exact solution not ready for 3d case')   
      end if
      
      !Coordinates and initializations.
      a%x = gpcod(1)
      a%y = gpcod(2)
      
      if ((a%x==0 .or. a%y==0) .and. LCR/=0) then
         if (a%x==0) a%x=0.0000000001
         if (a%y==0 )a%y=0.0000000001
      end if
         
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
      
      a%Sxx = 0.0_rp
      a%Sxy = 0.0_rp
      a%Syy = 0.0_rp
      a%dsxxdx = 0.0_rp
      a%dsxxdy = 0.0_rp
      a%dsyydx = 0.0_rp
      a%dsyydy = 0.0_rp
      a%dsxydx = 0.0_rp
      a%dsxydy = 0.0_rp
      
      if (LCR/=0) then
         a%Pxx = 0.0_rp
         a%Pyy = 0.0_rp
         a%Pxy = 0.0_rp
         
         a%dpxxdx = 0.0_rp
         a%dpxxdy = 0.0_rp
         a%dpyydx = 0.0_rp
         a%dpyydy = 0.0_rp
         a%dpxydx = 0.0_rp
         a%dpxydy = 0.0_rp
         
         a%ds1dx = 0.0_rp
         a%ds2dx = 0.0_rp
         a%ds3dx = 0.0_rp
         
         a%ds1dy = 0.0_rp
         a%ds2dy = 0.0_rp
         a%ds3dy = 0.0_rp

         a%expoSxx = 0.0_rp
         a%expoSxy = 0.0_rp
         a%expoSyy = 0.0_rp

         a%auxL_inv=0.0_rp
      end if   
 
      a%dudt   = 0.0_rp
      a%dvdt   = 0.0_rp
      tfunc = 1.0_rp
      dtfunc = 0.0_rp
      cdt = 1.0_rp
      cdt2=4.0_rp
      
      a%xdens=php%MatProp(1)%densi
    
      if(php%MatProp(1)%lawvi == 0)then
         a%xvisc = php%MatProp(1)%visco
         a%nn    = 1.0_rp      
      elseif(php%MatProp(1)%lawvi == 1)then
         a%xvisc = php%MatProp(1)%LawViParam(1)
         a%nn    = php%MatProp(1)%LawViParam(2)
      elseif(php%MatProp(1)%lawvi == 2)then
         call runend('sup_ComputeSolution: exact solution not implemented for Carreau model')
      elseif(php%MatProp(1)%lawvi == -1)then
         a%xvisc = php%MatProp(1)%LawViParam(1)
         a%Beta  = php%MatProp(1)%LawViParam(2) 
         a%lamb  = php%MatProp(1)%LawViParam(3)
         a%auxNL = 0.0_rp
      elseif(php%MatProp(1)%lawvi == -2)then
         a%xvisc = php%MatProp(1)%LawViParam(1)
         a%Beta  = php%MatProp(1)%LawViParam(2) 
         a%lamb  = php%MatProp(1)%LawViParam(3)
         a%auxNL = php%MatProp(1)%LawViParam(4)
      elseif(php%MatProp(1)%lawvi == -3)then
         call runend('sup_ComputeSolution: Phan-Thien-Tanner model not implemented')           
      end if
      
      if(php%MatProp(1)%lawvi<0)then
         auxiter=php%incremental         
         !incremental scheme for the stationary case
         if(php%itera<=auxiter .and. php%kfl_timei==0)then
            a%lamb = (php%MatProp(1)%LawViParam(3)/auxiter)*php%itera
         elseif(php%itera>auxiter .and. php%kfl_timei==0 )then
            a%lamb  = php%MatProp(1)%LawViParam(3)
         elseif(php%kfl_timei==1)then
            a%lamb =php%MatProp(1)%LawViParam(3)
         end if  
         !incremental scheme for the transitory case
         if(php%kfl_timei==1)then
            if(php%istep<=auxiter)then
               a%lamb = (php%MatProp(1)%LawViParam(3)/auxiter)*php%istep
            elseif(php%istep>auxiter)then
               a%lamb  = php%MatProp(1)%LawViParam(3) 
            end if         
         end if    
         
         if (LCR==1) then
            if (php%MatProp(1)%LawViParam(3)==0) then
               a%lamb0=0.01_rp
            else
               !a%lamb0=a%lamb*php%MatProp(1)%LawViParam(5) !SETLAMBDA0
               a%lamb0=php%MatProp(1)%LawViParam(5)
!                a%lamb=a%MatProp(1)%LawViParam(3)*a%MatProp(1)%LawviParam(5)
            end if
            
            a%auxL=(a%xvisc*(1.0_rp-a%Beta))/a%lamb0
            a%auxL_inv=1.0_rp/a%auxL
          end if
      endif
      
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
         a%viscnn = a%xvisc*(a%dudx*a%dudx+0.5_rp*(a%dudy+a%dvdx)**2+a%dvdy*a%dvdy)**((a%nn-1.0_rp)/2.0_rp)
         a%Sxx    = 2.0_rp*a%viscnn*a%dudx        
         a%Sxy    = a%viscnn*(a%dudy+a%dvdx)
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

          
         if (LCR==0) then
         
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
          
         else if (LCR/=0) then 
            !Psi and derivatives 
            a%Pxx    = (2.0_rp*a%x + 3)*tfunc
            a%Pxy    = (a%x + a%y)*tfunc 
            a%Pyy    = (2.0_rp*a%y + 3)*tfunc
            a%dpxxdx = 2.0_rp*tfunc
            a%dpxxdy = 0.0_rp
            a%dpyydx = 0.0_rp
            a%dpyydy = 2.0_rp*tfunc
            a%dpxydx = 1.0_rp*tfunc
            a%dpxydy = 1.0_rp*tfunc 
            !temporal derivatives
            a%dpxxdt = (2.0_rp*a%x + 3)*dtfunc
            a%dpyydt = (2.0_rp*a%y + 3)*dtfunc
            a%dpxydt = (a%x + a%y)*dtfunc     
            
            !Exp(Psi)
            a%expoSxx =(2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + 2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3))/(4*sqrt(a%x**2+a%y**2))
            a%expoSxy =(sqrt(2.0)*exp(a%x + a%y + 3)*sinh(sqrt(2.0)*sqrt(a%x**2+a%y**2))*(a%x + a%y))/(2*sqrt(a%x**2+a%y**2))
            a%expoSyy = (2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + 2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) - sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3))/(4*sqrt(a%x**2+a%y**2))
            
            !Grad(Exp(Psi))
            a%ds1dx= (2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*(sqrt(a%x**2 + a%y**2)**3) + 2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*(sqrt(a%x**2 + a%y**2)**3) + 3*sqrt(2.0)*a%x**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%y**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%y**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 2*a%x**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 2*a%x**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 2*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - 2*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 3*sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y))/(4*(sqrt(a%x**2 + a%y**2)**3))
            a%ds2dx= ((a%x + a%y)*(sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + 2*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
            a%ds3dx=  (exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*(2*a%x**2*exp(a%x + a%y + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%y**2*exp(a%x + a%y + 3) - sqrt(2.0)*a%y**3*exp(a%x + a%y + 3) - sqrt(2.0)*a%x**3*exp(a%x + a%y + 3) + sqrt(2.0)*a%x**3*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**3*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(a%x + a%y + 3) + 2*a%x*a%y*exp(a%x + a%y + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y**2*exp(a%x + a%y + 3) - sqrt(2.0)*a%x**2*a%y*exp(a%x + a%y + 3) + sqrt(2.0)*a%x*a%y**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%x**2*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))

            a%ds1dy= ((a%x + a%y)*(sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + 2*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
            a%ds2dy=  (2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt((a%x**2 + a%y**2)**3) + 2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt((a%x**2 + a%y**2)**3) + sqrt(2.0)*a%x**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 3*sqrt(2.0)*a%y**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%y**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 2*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 2*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 2*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - 2*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 3*sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y))/(4*sqrt((a%x**2 + a%y**2)**3))
            a%ds3dy= (sqrt(2.0)*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) + 1) + sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) - 1) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) + 1) + sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) - 1))/(4*sqrt(a%x**2 + a%y**2)) - (a%y*(sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
                  
         end if

         
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
            
      else if(php%kfl_exacs==12)then 
         !velocity
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

         !pressure and derivatives
         a%p     =  a%x+a%y
         a%dpdx  =  1.0_rp
         a%dpdy  =  1.0_rp

         !stress and derivatives 
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

         !Psi
         a%Pxx = (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2))
         a%Pyy = (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2))
         a%Pxy = ((a%x + a%y)*((sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/4 - (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/4))/sqrt(a%x**2+a%y**2)
         
         ! Grad(Psi)
         a%dPxxdx = (a%lamb0*(4*a%x**2 - 2*a%x*a%y + 2*a%y**2 - 3*sqrt(2.0_rp)*a%x*sqrt(a%x**2+a%y**2) + sqrt(2.0_rp)*a%y*sqrt(a%x**2+a%y**2)))/(4*(3*a%lamb0*a%x**2 + a%lamb0*a%x**3 + 3*a%lamb0*a%y**2 + a%lamb0*a%y**3 + a%x**2*(a%xvisc*(1-a%beta)) + (a%xvisc*(1-a%beta))*a%y**2 + a%lamb0*a%x*a%y**2 + a%lamb0*a%x**2*a%y - sqrt(2.0_rp)*a%lamb0*(sqrt(a%x**2+a%y**2))**3)) + (a%lamb0*(4*a%x**2 - 2*a%x*a%y + 2*a%y**2 + 3*sqrt(2.0_rp)*a%x*sqrt(a%x**2+a%y**2) - sqrt(2.0_rp)*a%y*sqrt(a%x**2+a%y**2)))/(4*(3*a%lamb0*a%x**2 + a%lamb0*a%x**3 + 3*a%lamb0*a%y**2 + a%lamb0*a%y**3 + a%x**2*(a%xvisc*(1-a%beta)) + (a%xvisc*(1-a%beta))*a%y**2 + a%lamb0*a%x*a%y**2 + a%lamb0*a%x**2*a%y + sqrt(2.0_rp)*a%lamb0*(sqrt(a%x**2+a%y**2))**3)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)*(a%x + a%y)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)*(a%x + a%y)) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) + (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%x - sqrt(2.0_rp)*a%y + 2*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y)) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y - sqrt(2.0_rp)*a%x + 2*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y))
         a%dPyydx = (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%x + sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%x - sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%x - sqrt(a%x**2+a%y**2))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))) + (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%x + sqrt(a%x**2+a%y**2))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2)))
         a%dPxydx = (a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/(2*(a%x**2 + a%y**2)) + (a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/(2*(a%x**2 + a%y**2)) - (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y + sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) - (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y - sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x + a%y))/(4*(sqrt(a%x**2+a%y**2))**3) + (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x + a%y))/(4*(sqrt(a%x**2+a%y**2))**3) + (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%x + sqrt(a%x**2+a%y**2))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))) + (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%x - sqrt(a%x**2+a%y**2))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2)))
 
         a%dPxxdy = (a%lamb0*(a%x + a%y)*(2*a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(3*a%lamb0*a%x**2 + a%lamb0*a%x**3 + 3*a%lamb0*a%y**2 + a%lamb0*a%y**3 + a%x**2*(a%xvisc*(1-a%beta)) + (a%xvisc*(1-a%beta))*a%y**2 + a%lamb0*a%x*a%y**2 + a%lamb0*a%x**2*a%y + sqrt(2.0_rp)*a%lamb0*(sqrt(a%x**2+a%y**2))**3)) + (a%lamb0*(a%x + a%y)*(2*a%x - sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(3*a%lamb0*a%x**2 + a%lamb0*a%x**3 + 3*a%lamb0*a%y**2 + a%lamb0*a%y**3 + a%x**2*(a%xvisc*(1-a%beta)) + (a%xvisc*(1-a%beta))*a%y**2 + a%lamb0*a%x*a%y**2 + a%lamb0*a%x**2*a%y - sqrt(2.0_rp)*a%lamb0*(sqrt(a%x**2+a%y**2))**3)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)*(a%x + a%y)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*sqrt(a%x**2+a%y**2)*(a%x + a%y)) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%x - sqrt(2.0_rp)*a%y + 2*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y)) + (sqrt(2.0_rp)*a%x*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y - sqrt(2.0_rp)*a%x + 2*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(a%x + a%y))
         a%dPyydy = (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y + sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) + (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(sqrt(2.0_rp)*a%y - sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(sqrt(a%x**2+a%y**2))**3) - (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%y - sqrt(a%x**2+a%y**2))*(a%x - a%y + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))) + (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%y + sqrt(a%x**2+a%y**2))*(a%y - a%x + sqrt(2.0_rp)*sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2)))
         a%dPxydy = (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/(4*sqrt(a%x**2+a%y**2)) - (sqrt(2.0_rp)*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta))))/(4*sqrt(a%x**2+a%y**2)) - (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x + a%y))/(4*(sqrt(a%x**2+a%y**2))**3) + (sqrt(2.0_rp)*a%y*log((3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))/(a%xvisc*(1-a%beta)))*(a%x + a%y))/(4*(sqrt(a%x**2+a%y**2))**3) + (sqrt(2.0_rp)*a%lamb0*(sqrt(2.0_rp)*a%y + sqrt(a%x**2+a%y**2))*(a%x + a%y))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y + sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2))) + (sqrt(2.0_rp)*a%lamb0*(a%x + a%y)*(sqrt(2.0_rp)*a%y - sqrt(a%x**2+a%y**2)))/(4*(a%x**2 + a%y**2)*(3*a%lamb0 + (a%xvisc*(1-a%beta)) + a%lamb0*a%x + a%lamb0*a%y - sqrt(2.0_rp)*a%lamb0*sqrt(a%x**2+a%y**2)))
         
         ! Exp(Psi)
         a%expoSxx= a%auxL_inv*a%Sxx+1.0_rp 
         a%expoSyy= a%auxL_inv*a%Syy+1.0_rp
         a%expoSxy= a%auxL_inv*a%Sxy
         
         ! Grad(Exp(Psi))
         a%ds1dx=a%auxL_inv*a%dsxxdx
         a%ds2dx=a%auxL_inv*a%dsyydx
         a%ds3dx=a%auxL_inv*a%dsxydx
         
         a%ds1dy=a%auxL_inv*a%dsxxdy
         a%ds2dy=a%auxL_inv*a%dsyydy
         a%ds3dy=a%auxL_inv*a%dsxydy
         
      else if(php%kfl_exacs==13) then
      
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
         a%p     =  a%x+a%y
         a%dpdx  =  1.0_rp
         a%dpdy  =  1.0_rp

      
         !Psi and derivatives 
         a%Pxx    = (2.0_rp*a%x + 3)
         a%Pxy    = (a%x + a%y)
         a%Pyy    = (2.0_rp*a%y + 3)
         a%dpxxdx = 2.0_rp
         a%dpxxdy = 0.0_rp
         a%dpyydx = 0.0_rp
         a%dpyydy = 2.0_rp
         a%dpxydx = 1.0_rp
         a%dpxydy = 1.0_rp
         !temporal derivatives
         a%dpxxdt = (2.0_rp*a%x + 3)*dtfunc
         a%dpyydt = (2.0_rp*a%y + 3)*dtfunc
         a%dpxydt = (a%x + a%y)*dtfunc     
         
         !Exp(Psi)
         a%expoSxx =(2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + 2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3))/(4*sqrt(a%x**2+a%y**2))
         a%expoSxy =(sqrt(2.0)*exp(a%x + a%y + 3)*sinh(sqrt(2.0)*sqrt(a%x**2+a%y**2))*(a%x + a%y))/(2*sqrt(a%x**2+a%y**2))
         a%expoSyy =(2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) + 2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3)*sqrt(a%x**2+a%y**2) - sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2+a%y**2) + 3))/(4*sqrt(a%x**2+a%y**2))
         
         !Grad(Exp(Psi))
         a%ds1dx= (2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*(sqrt(a%x**2 + a%y**2)**3) + 2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*(sqrt(a%x**2 + a%y**2)**3) + 3*sqrt(2.0)*a%x**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%y**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%y**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 2*a%x**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 2*a%x**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 2*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - 2*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 3*sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y))/(4*(sqrt(a%x**2 + a%y**2)**3))
         a%ds2dx= ((a%x + a%y)*(sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + 2*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
         a%ds3dx= (exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*(2*a%x**2*exp(a%x + a%y + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%y**2*exp(a%x + a%y + 3) - sqrt(2.0)*a%y**3*exp(a%x + a%y + 3) - sqrt(2.0)*a%x**3*exp(a%x + a%y + 3) + sqrt(2.0)*a%x**3*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**3*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(a%x + a%y + 3) + 2*a%x*a%y*exp(a%x + a%y + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y**2*exp(a%x + a%y + 3) - sqrt(2.0)*a%x**2*a%y*exp(a%x + a%y + 3) + sqrt(2.0)*a%x*a%y**2*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%x**2*a%y*exp(a%x + a%y + 2*sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))

         a%ds1dy= tfunc*((a%x + a%y)*(sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + 2*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + 2*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y**2*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y**2*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
         a%ds2dy= tfunc*(2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt((a%x**2 + a%y**2)**3) + 2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt((a%x**2 + a%y**2)**3) + sqrt(2.0)*a%x**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 3*sqrt(2.0)*a%y**3*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%y**3*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 2*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + 2*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) + sqrt(2.0)*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - sqrt(2.0)*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 2*a%x*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - 2*a%x*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y)*sqrt(a%x**2 + a%y**2) - sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + 3*sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) + sqrt(2.0)*a%x*a%y**2*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y) - 3*sqrt(2.0)*a%x**2*a%y*exp(3.0)*exp(-sqrt(2.0)*sqrt(a%x**2 + a%y**2))*exp(a%x)*exp(a%y))/(4*sqrt((a%x**2 + a%y**2)**3))
         a%ds3dy= tfunc*(sqrt(2.0)*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) + 1) + sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) - 1) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) + 1) + sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)*((sqrt(2.0)*a%y)/sqrt(a%x**2 + a%y**2) - 1))/(4*sqrt(a%x**2 + a%y**2)) - (a%y*(sqrt(2.0)*a%x*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%x*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) + sqrt(2.0)*a%y*exp(a%x + a%y + sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3) - sqrt(2.0)*a%y*exp(a%x + a%y - sqrt(2.0)*sqrt(a%x**2 + a%y**2) + 3)))/(4*sqrt((a%x**2 + a%y**2)**3))
               
         !Sigma
         a%Sxx = tfunc*a%auxL*(a%expoSxx -  1.0_rp)
         a%Syy = tfunc*a%auxL*(a%expoSyy -  1.0_rp)
         a%Sxy = tfunc*a%auxL*(a%expoSxy)
         
         a%dSxxdt = a%auxL*(a%expoSxx -  1.0_rp)*dtfunc
         a%dSyydt = a%auxL*(a%expoSyy -  1.0_rp)*dtfunc
         a%dSxydt = a%auxL*(a%expoSxy)*dtfunc

         !Grad(Sigma)
         a%dsxxdx=tfunc*a%ds1dx*a%auxL
         a%dsyydx=tfunc*a%ds2dx*a%auxL
         a%dsxydx=tfunc*a%ds3dx*a%auxL
         
         a%dsxxdy=tfunc*a%ds1dy*a%auxL
         a%dsyydy=tfunc*a%ds2dy*a%auxL
         a%dsxydy=tfunc*a%ds3dy*a%auxL
         
      else if(php%kfl_exacs==16)then 
      
         !velocity
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

         a%expoSxx= a%auxL_inv*a%Sxx+1.0_rp !Exp(Psi)
         a%expoSyy= a%auxL_inv*a%Syy+1.0_rp
         a%expoSxy= a%auxL_inv*a%Sxy
         
         a%ds1dx=a%auxL_inv*a%dsxxdx
         a%ds2dx=a%auxL_inv*a%dsyydx
         a%ds3dx=a%auxL_inv*a%dsxydx
         
         a%ds1dy=a%auxL_inv*a%dsxxdy
         a%ds2dy=a%auxL_inv*a%dsyydy
         a%ds3dy=a%auxL_inv*a%dsxydy
         
         a%Pxx=log(1 - 4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))/2 + log(4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y) + 1)/2 - (5*sqrt(26.0_rp)*log(1 - 4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)))/52 + (5*sqrt(26.0_rp)*log(4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y) + 1))/52
         a%Pyy=log(1 - 4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))/2 + log(4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y) + 1)/2 + (5*sqrt(26.0_rp)*log(1 - 4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)))/52 - (5*sqrt(26.0_rp)*log(4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y) + 1))/52
         a%Pxy=(sqrt(26.0_rp)*log(4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y) + 1))/52 - (sqrt(26.0_rp)*log(1 - 4*sqrt(26.0_rp)*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)))/52
         
         a%dPxxdx=(4*a%auxL_inv*pi*cos(pi*a%y)*(2*cos(pi*a%x)**2 - 1)*(5*sin(pi*a%x)*sin(pi*a%y) + 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)**3 + 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y) - 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y)**3 - 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)))/(sin(pi*a%x) - 416*a%auxL_inv**2*sin(pi*a%x)**3*sin(pi*a%y)**2 + 416*a%auxL_inv**2*sin(pi*a%x)**3*sin(pi*a%y)**4 + 416*a%auxL_inv**2*sin(pi*a%x)**5*sin(pi*a%y)**2 - 416*a%auxL_inv**2*sin(pi*a%x)**5*sin(pi*a%y)**4)
         a%dPyydx=-(4*a%auxL_inv*pi*cos(pi*a%y)*(2*cos(pi*a%x)**2 - 1)*(5*sin(pi*a%x)*sin(pi*a%y) - 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)**3 - 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y) + 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y)**3 + 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)))/(sin(pi*a%x) - 416*a%auxL_inv**2*sin(pi*a%x)**3*sin(pi*a%y)**2 + 416*a%auxL_inv**2*sin(pi*a%x)**3*sin(pi*a%y)**4 + 416*a%auxL_inv**2*sin(pi*a%x)**5*sin(pi*a%y)**2 - 416*a%auxL_inv**2*sin(pi*a%x)**5*sin(pi*a%y)**4)
         a%dPxydx=(4*a%auxL_inv*pi*cos(pi*a%y)*sin(pi*a%y)*(2*cos(pi*a%x)**2 - 1))/(416*a%auxL_inv**2*cos(pi*a%x)**2*cos(pi*a%y)**4 - 416*a%auxL_inv**2*cos(pi*a%x)**2*cos(pi*a%y)**2 + 416*a%auxL_inv**2*cos(pi*a%x)**4*cos(pi*a%y)**2 - 416*a%auxL_inv**2*cos(pi*a%x)**4*cos(pi*a%y)**4 + 1)
         a%dPxxdy=(4*a%auxL_inv*pi*cos(pi*a%x)*(2*cos(pi*a%y)**2 - 1)*(5*sin(pi*a%x)*sin(pi*a%y) + 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)**3 + 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y) - 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y)**3 - 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)))/(sin(pi*a%y) - 416*a%auxL_inv**2*sin(pi*a%x)**2*sin(pi*a%y)**3 + 416*a%auxL_inv**2*sin(pi*a%x)**2*sin(pi*a%y)**5 + 416*a%auxL_inv**2*sin(pi*a%x)**4*sin(pi*a%y)**3 - 416*a%auxL_inv**2*sin(pi*a%x)**4*sin(pi*a%y)**5)
         a%dPyydy=-(4*a%auxL_inv*pi*cos(pi*a%x)*(2*cos(pi*a%y)**2 - 1)*(5*sin(pi*a%x)*sin(pi*a%y) - 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)**3 - 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y) + 104*a%auxL_inv*cos(pi*a%x)**3*cos(pi*a%y)**3 + 104*a%auxL_inv*cos(pi*a%x)*cos(pi*a%y)))/(sin(pi*a%y) - 416*a%auxL_inv**2*sin(pi*a%x)**2*sin(pi*a%y)**3 + 416*a%auxL_inv**2*sin(pi*a%x)**2*sin(pi*a%y)**5 + 416*a%auxL_inv**2*sin(pi*a%x)**4*sin(pi*a%y)**3 - 416*a%auxL_inv**2*sin(pi*a%x)**4*sin(pi*a%y)**5)
         a%dPxydy=(4*a%auxL_inv*pi*cos(pi*a%x)*sin(pi*a%x)*(2*cos(pi*a%y)**2 - 1))/(416*a%auxL_inv**2*cos(pi*a%x)**2*cos(pi*a%y)**4 - 416*a%auxL_inv**2*cos(pi*a%x)**2*cos(pi*a%y)**2 + 416*a%auxL_inv**2*cos(pi*a%x)**4*cos(pi*a%y)**2 - 416*a%auxL_inv**2*cos(pi*a%x)**4*cos(pi*a%y)**4 + 1)   
            
        

        else if(php%kfl_exacs==17)then 

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
         a%Pxx    = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%Pxy    = sin(2*pi*a%x)*sin(2*pi*a%y)*tfunc 
         a%Pyy    = -a%Pxx
         a%dpxxdx  = 10.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dpxxdy  = 10.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         a%dpyydx = -a%dpxxdx
         a%dpyydy = -a%dpxxdy
         a%dpxydx  = 2.0_rp*pi*cos(2*pi*a%x)*sin(2*pi*a%y)*tfunc
         a%dpxydy  = 2.0_rp*pi*sin(2*pi*a%x)*cos(2*pi*a%y)*tfunc
         !temporal derivatives
         a%dpxxdt = 5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dpyydt = -5.0_rp*sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc
         a%dpxydt = sin(2*pi*a%x)*sin(2*pi*a%y)*dtfunc 
         
         
         a%expoSxx = (exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*(26*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + 5*sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - 5*sqrt(26.0) + 26))/52
         a%expoSxy = (sqrt(26.0)*sinh(4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)))/26
         a%expoSyy = (exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*(26*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - 5*sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + 5*sqrt(26.0) + 26))/52

         a%ds1dx = 2*pi*exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*cos(2*pi*a%x)*cos(pi*a%y)*sin(pi*a%y)*(5*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - sqrt(26.0) + 5)
         a%ds2dx = -2*pi*exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*cos(2*pi*a%x)*cos(pi*a%y)*sin(pi*a%y)*(5*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + sqrt(26.0) + 5)
         a%ds3dx = 2*pi*cosh(sqrt(26.0)*sin(2*pi*a%x)*sin(2*pi*a%y))*cos(2*pi*a%x)*sin(2*pi*a%y)        
         a%ds1dy = 2*pi*exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*cos(pi*a%x)*cos(2*pi*a%y)*sin(pi*a%x)*(5*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - sqrt(26.0) + 5)
         a%ds2dy = -2*pi*exp(-4*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y))*cos(pi*a%x)*cos(2*pi*a%y)*sin(pi*a%x)*(5*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) - sqrt(26.0)*exp(8*sqrt(26.0)*cos(pi*a%x)*cos(pi*a%y)*sin(pi*a%x)*sin(pi*a%y)) + sqrt(26.0) + 5)
         a%ds3dy = 2*pi*cosh(sqrt(26.0)*sin(2*pi*a%x)*sin(2*pi*a%y))*cos(2*pi*a%y)*sin(2*pi*a%x)
         
!          !Sigma
!          a%Sxx = a%auxL*(a%expoSxx -  1.0_rp)
!          a%Syy = a%auxL*(a%expoSyy -  1.0_rp)
!          a%Sxy = a%auxL*(a%expoSxy)
!          
!          !Grad(Sigma)
!          a%dsxxdx=a%ds1dx*a%auxL
!          a%dsyydx=a%ds2dx*a%auxL
!          a%dsxydx=a%ds3dx*a%auxL
!          
!          a%dsxxdy=a%ds1dy*a%auxL
!          a%dsyydy=a%ds2dy*a%auxL
!          a%dsxydy=a%ds3dy*a%auxL
         
         
        !Sigma
         a%Sxx = tfunc*a%auxL*(a%expoSxx -  1.0_rp)
         a%Syy = tfunc*a%auxL*(a%expoSyy -  1.0_rp)
         a%Sxy = tfunc*a%auxL*(a%expoSxy)
         
         a%dSxxdt = a%auxL*(a%expoSxx -  1.0_rp)*dtfunc
         a%dSyydt = a%auxL*(a%expoSyy -  1.0_rp)*dtfunc
         a%dSxydt = a%auxL*(a%expoSxy)*dtfunc

         !Grad(Sigma)
         a%dsxxdx=tfunc*a%ds1dx*a%auxL
         a%dsyydx=tfunc*a%ds2dx*a%auxL
         a%dsxydx=tfunc*a%ds3dx*a%auxL
         
         a%dsxxdy=tfunc*a%ds1dy*a%auxL
         a%dsyydy=tfunc*a%ds2dy*a%auxL
         a%dsxydy=tfunc*a%ds3dy*a%auxL
        
        end if    
   end subroutine
   
   
   subroutine sup_GetPressure(a,ndime,expre,exprg)
      use typre
      implicit none
      class(SupExacso) :: a
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
      class(SupExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exvel(ndime),exveg(ndime,ndime) 
      
      exvel(1)   = a%u  
      exvel(2)   = a%v
      exveg(1,1) = a%dudx
      exveg(2,1) = a%dvdx
      exveg(1,2) = a%dudy
      exveg(2,2) = a%dvdy
         
   end subroutine

   subroutine sup_GetStress(a,ndime,auxLCR,exsig,exsigr)
      use typre
      implicit none
      class(SupExacso) :: a      
      integer(ip), intent(in)   :: ndime,auxLCR
      real(rp),    intent(out)  :: exsig((ndime-1)*(ndime-1)+2),exsigr(((ndime-1)*(ndime-1)+2),ndime) 
      real(rp)                  :: sig((ndime-1)*(ndime-1)+2),psi((ndime-1)*(ndime-1)+2),tau((ndime-1)*(ndime-1)+2)
      real(rp)                  :: psir(((ndime-1)*(ndime-1)+2),ndime),auxL,aux2,taugr(((ndime-1)*(ndime-1)+2),ndime)
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
   
   
   
   subroutine sup_GetPsi(a,ndime,exsig,exsigr)
      use typre
      implicit none
      class(SupExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exsig((ndime-1)*(ndime-1)+2),exsigr(((ndime-1)*(ndime-1)+2),ndime) 
      real(rp)                  :: sig((ndime-1)*(ndime-1)+2),psi((ndime-1)*(ndime-1)+2),tau((ndime-1)*(ndime-1)+2)
      real(rp)                  :: psir(((ndime-1)*(ndime-1)+2),ndime),auxL,aux2,taugr(((ndime-1)*(ndime-1)+2),ndime)
      real(rp)                  :: sigr(((ndime-1)*(ndime-1)+2),ndime), point1,point2
      
      point1=a%x
      point2=a%y
      exsig(1)    = a%Pxx
      exsig(2)    = a%Pyy
      exsig(3)    = a%Pxy
      exsigr(1,1) = a%dpxxdx
      exsigr(2,1) = a%dpyydx
      exsigr(3,1) = a%dpxydx
      exsigr(1,2) = a%dpxxdy
      exsigr(2,2) = a%dpyydy
      exsigr(3,2) = a%dpxydy

   end subroutine   
   
   subroutine sup_GetExponential(a,ndime,Expo)
      use typre
      implicit none
      class(SupExacso) :: a   
      real(rp) :: Expo((ndime-1)*(ndime-1)+2),point1,point2
      integer(ip) :: ndime
      
      
      point1=a%x
      point2=a%y
      Expo(1)= a%expoSxx
      Expo(2)= a%expoSyy
      Expo(3)= a%expoSxy
      
   end subroutine
    
   
    subroutine sup_GetExponentialDivergence(a,ndime,Div)
      use typre
      implicit none
      class(SupExacso) :: a   
      real(rp) :: Div(ndime),point1,point2,exsigr(((ndime-1)*(ndime-1)+2),ndime)
      integer(ip) :: ndime
      
      point1=a%x
      point2=a%y
      Div(1)=a%ds1dx + a%ds3dy 
      Div(2)=a%ds3dx + a%ds2dy
   end subroutine
   
 
   
   subroutine sup_GetExponentialGradient(a,ndime,Grad)
      use typre
      implicit none
      class(SupExacso) :: a   
      real(rp) :: Grad((ndime-1)*(ndime-1)+2,ndime), point1, point2
      integer(ip) :: ndime
      
      point1=a%x
      point2=a%y
      Grad(1,1)=a%ds1dx
      Grad(2,1)=a%ds2dx
      Grad(3,1)=a%ds3dx 
      Grad(1,2)=a%ds1dy
      Grad(2,2)=a%ds2dy
      Grad(3,2)=a%ds3dy

   end subroutine
   

   subroutine sup_GetForce(a,ndime,auxLCR,kfl_LogFormulation,elext,elextC,elextS,elextEstab,elextEstab2,elextSEstab,elextEstab3,&
                              elextSEstab2,elextEstab4,elextSEstab3,elextSEstab4,elextEstab5,elextSMat,elextSEstabMat,php)
      use typre
      implicit none
      class(SupExacso) :: a
      class(ThreeFieldNSProblem) :: php 
      integer(ip), intent(in)  :: ndime,auxLCR,kfl_LogFormulation
      real(rp), intent(inout)  :: elext(ndime),elextC(1),elextS((ndime-1)*(ndime-1)+2),elextSEstab((ndime-1)*(ndime-1)+2)
      real(rp), intent(inout)  :: elextEstab(ndime),elextEstab2(ndime), elextEstab3(ndime),elextSEstab2((ndime-1)*(ndime-1)+2)
      real(rp), intent(inout)  :: elextSEstab3((ndime-1)*(ndime-1)+2), elextSEstab4((ndime-1)*(ndime-1)+2) 
      real(rp), intent(inout)  :: elextEstab4(ndime), elextSMat(ndime,ndime),elextSEstabMat(ndime,ndime),elextEstab5((ndime-1)*(ndime-1)+2)
      real(rp) :: D1x,D2x,D1y,D2y,aux1,aux2,aux3,aux4,dudtforce,dvdtforce,dsxxdtforce,dsyydtforce,dsxydtforce,auxL
      integer(ip) :: convection,auxconv,auxtrac
      real(rp) :: exsig((ndime-1)*(ndime-1)+2),exsigr(((ndime-1)*(ndime-1)+2),ndime) 
      real(rp) :: tau((ndime-1)*(ndime-1)+2), tau_matrix(ndime,ndime), R(ndime,ndime), D(ndime)
      real(rp) :: psi((ndime-1)*(ndime-1)+2),Exp_psi((ndime-1)*(ndime-1)+2)
      real(rp) :: CTens(ndime,ndime,ndime,ndime), DTens(ndime,ndime,ndime,ndime)
      real(rp) :: veg(ndime,ndime),exvel(ndime), B(ndime,ndime), Omega(ndime,ndime)
      real(rp) :: taugr(((ndime-1)*(ndime-1)+2),ndime)
      real(rp) :: minus_psi((ndime-1)*(ndime-1)+2),minus_exp((ndime-1)*(ndime-1)+2),minus_exp_matrix(ndime,ndime)
      real(rp) :: U_descomp(ndime,ndime), Prod(ndime,ndime), res((ndime-1)*(ndime-1)+2), psi_matrix(ndime,ndime)
      real(rp) :: lambda,aux_lamb, lambda0, Div1x, Div2x, Div1y, Div2y, gpcod(ndime), pres1, pres2, point1,point2, aux3b
      real(rp) :: ExpGradU11, ExpGradU22, ExpGradU12, ExpGradU21, ExpoCuad((ndime-1)*(ndime-1)+2)
      real(rp) :: M(ndime,ndime)
      
      M(1,2)=0.0_rp
      M(2,1)=0.0_rp
      M(1,1)=1.0_rp
      M(2,2)=0.0_rp
      
      aux1 =0.0_rp
      aux2 =0.0_rp
      aux3 =0.0_rp
      D1x  =0.0_rp
      D2x  =0.0_rp
      D1y  =0.0_rp
      D2y  =0.0_rp
      !temporal derivatives
      dudtforce = 0.0_rp
      dvdtforce = 0.0_rp
      dsxxdtforce = 0.0_rp
      dsyydtforce = 0.0_rp
      dsxydtforce =0.0_rp
      convection=1
            
      
      if(php%kfl_timei==1)then
         dudtforce = a%dudt*a%xdens
         dvdtforce = a%dvdt*a%xdens
      endif
      
      if(php%kfl_timei==1 .and. php%MatProp(1)%lawvi<0)then
         if(auxLCR==0) then
            dsxxdtforce = a%dsxxdt*(a%lamb/(2.0_rp*a%xvisc))
            dsyydtforce = a%dsyydt*(a%lamb/(2.0_rp*a%xvisc))            
            dsxydtforce = a%dsxydt*(a%lamb/(2.0_rp*a%xvisc))
         elseif(auxLCR/=0) then
            dsxxdtforce = a%dsxxdt*(a%lamb/(2.0_rp*a%xvisc))
            dsyydtforce = a%dsyydt*(a%lamb/(2.0_rp*a%xvisc))         
            dsxydtforce = a%dsxydt*(a%lamb/(2.0_rp*a%xvisc))
         end if
      end if

      !Diffusive term      
      if(php%MatProp(1)%lawvi==0)then
      
         D1x= a%xvisc*a%d2udx 
         D2x= a%xvisc*a%d2udy
         D1y= a%xvisc*a%d2vdx 
         D2y= a%xvisc*a%d2vdy  

         !Constitutive force term  
         elextS(1)= (1.0_rp/(2.0_rp*a%xvisc))*a%Sxx- a%dudx
         elextS(2)= (1.0_rp/(2.0_rp*a%xvisc))*a%Syy- a%dvdy 
         elextS(3)= (1.0_rp/(2.0_rp*a%xvisc))*a%Sxy- (1.0_rp/(2.0_rp))*(a%dudy + a%dvdx)            
         
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

         !Constitutive force term  
         elextS(1)= (1.0_rp/(2.0_rp*a%viscnn))*a%Sxx- a%dudx
         elextS(2)= (1.0_rp/(2.0_rp*a%viscnn))*a%Syy- a%dvdy 
         elextS(3)= (1.0_rp/(2.0_rp*a%viscnn))*a%Sxy- (1.0_rp/(2.0_rp))*(a%dudy + a%dvdx) 
         
       elseif(php%MatProp(1)%lawvi < 0)then 
       
         aux2 = (1.0_rp - a%beta)
         auxconv=1.0_rp
         aux3 = 1.0_rp*a%beta*a%xvisc
         aux1 = 1.0_rp/(2.0_rp*a%xvisc)
         
         if (auxLCR==0) then
            
            aux4 = a%auxNL*a%lamb/((2.0_rp*aux2)*(a%xvisc**2))
            auxtrac=1
                
            D1x= aux3*a%d2udx + a%dsxxdx
            D2x= aux3*a%d2udy + a%dsxydy  

            D1y= aux3*a%d2vdx + a%dsxydx
            D2y= aux3*a%d2vdy + a%dsyydy  
            
            !Constitutive force term  
            elextS(1)= elextS(1) + aux1*a%Sxx - aux2*a%dudx + a%lamb*aux1*(auxconv*(a%u*a%dsxxdx + a%v*a%dsxxdy) - auxtrac*2.0_rp*(a%Sxx*a%dudx + a%Sxy*a%dudy)) &
               + aux4*(a%Sxx*a%Sxx + a%Sxy*a%Sxy) + dsxxdtforce
            elextS(2)= elextS(2) + aux1*a%Syy - aux2*a%dvdy + a%lamb*aux1*(auxconv*(a%u*a%dsyydx + a%v*a%dsyydy) - auxtrac*2.0_rp*(a%Sxy*a%dvdx + a%Syy*a%dvdy)) &
               + aux4*(a%Sxy*a%Sxy + a%Syy*a%Syy) + dsyydtforce
            elextS(3)= elextS(3) + aux1*a%Sxy - (aux2/2.0_rp)*(a%dudy + a%dvdx) + a%lamb*aux1*(auxconv*(a%u*a%dsxydx + a%v*a%dsxydy) &
               - auxtrac*(a%Sxx*a%dvdx + a%Syy*a%dudy)) + aux4*a%Sxy*(a%Sxx + a%Syy)  + dsxydtforce
                             
         else if (auxLCR==1) then
    
            lambda0=a%lamb0
            aux_lamb=a%lamb/(2.0_rp*lambda0)
            auxtrac=1
            
            a%auxG=(php%MatProp(1)%LawViParam(4)*a%lamb)/(2.0_rp*(a%lamb0)**2_ip)

            D1x= aux3*a%d2udx + a%auxL*a%DS1dx
            D2x= aux3*a%d2udy + a%auxL*a%DS3dy
            D1y= aux3*a%d2vdx + a%auxL*a%DS3dx
            D2y= aux3*a%d2vdy + a%auxL*a%DS2dy  
           
            elextS(1)= elextS(1)  + a%auxG  - 2.0_rp*a%auxG*a%expoSxx + a%auxG*((a%expoSxx**2)+a%expoSxy**2)
            elextS(2)= elextS(2)  + a%auxG  - 2.0_rp*a%auxG*a%expoSyy + a%auxG*((a%expoSyy**2)+a%expoSxy**2)
            elextS(3)= elextS(3)            - 2.0_rp*a%auxG*a%expoSxy + a%auxG*a%expoSxy*(a%expoSxx+a%expoSyy)
            
            elextS(1)= elextS(1)  -a%dudx 
            elextS(2)= elextS(2)  -a%dvdy
            elextS(3)= elextS(3)   -0.5_rp*(a%dudy + a%dvdx) 

            elextS(1)= elextS(1) - 1.0_rp/(2.0_rp*lambda0)  +(1.0_rp/(2.0_rp*lambda0))*a%expoSxx   + aux_lamb*auxconv*(a%u*a%ds1dx + a%v*a%ds1dy) 
            elextS(2)= elextS(2) - 1.0_rp/(2.0_rp*lambda0)  +(1.0_rp/(2.0_rp*lambda0))*a%expoSyy   + aux_lamb*auxconv*(a%u*a%ds2dx + a%v*a%ds2dy)
            elextS(3)= elextS(3)                            +(1.0_rp/(2.0_rp*lambda0))*a%expoSxy   + aux_lamb*auxconv*(a%u*a%ds3dx + a%v*a%ds3dy)
                
            elextS(1)= elextS(1) + dsxxdtforce
            elextS(2)= elextS(2) + dsyydtforce
            elextS(3)= elextS(3) + dsxydtforce     
                            
            elextS(1)=  elextS(1) + auxtrac*aux_lamb*2.0_rp*a%dudx      - auxtrac*aux_lamb*2.0_rp*((a%expoSxx*a%dudx + a%expoSxy*a%dudy)) 
            elextS(2)=  elextS(2) + auxtrac*aux_lamb*2.0_rp*a%dvdy      - auxtrac*aux_lamb*2.0_rp*((a%expoSxy*a%dvdx + a%expoSyy*a%dvdy)) 
            elextS(3)=  elextS(3) + auxtrac*aux_lamb*(a%dudy + a%dvdx)  - auxtrac*aux_lamb*(a%expoSxx*a%dvdx + a%expoSyy*a%dudy) 
            
            
            elextSEstab(1)= elextSEstab(1)   - 1.0_rp/(2.0_rp*lambda0)  +(1.0_rp/(2.0_rp*lambda0))*a%expoSxx 
            elextSEstab(2)= elextSEstab(2)   - 1.0_rp/(2.0_rp*lambda0)  +(1.0_rp/(2.0_rp*lambda0))*a%expoSyy 
            elextSEstab(3)= elextSEstab(3)                              +(1.0_rp/(2.0_rp*lambda0))*a%expoSxy 


            elextSEstab2(1)= elextSEstab2(1)  -a%dudx 
            elextSEstab2(2)= elextSEstab2(2)  -a%dvdy
            elextSEstab2(3)= elextSEstab2(3)  -0.5_rp*(a%dudy + a%dvdx) 
  
            elextSEstab3(1)= elextSEstab3(1)  + aux_lamb*auxconv*(a%u*a%ds1dx + a%v*a%ds1dy) 
            elextSEstab3(2)= elextSEstab3(2)  + aux_lamb*auxconv*(a%u*a%ds2dx + a%v*a%ds2dy) 
            elextSEstab3(3)= elextSEstab3(3)  + aux_lamb*auxconv*(a%u*a%ds3dx + a%v*a%ds3dy)

            
            elextSEstab4(1) = elextSEstab4(1) - auxtrac*aux_lamb*((a%expoSxx*a%dudx + a%expoSxy*a%dudy)- a%dudx )*2.0_rp
            elextSEstab4(2) = elextSEstab4(2) - auxtrac*aux_lamb*((a%expoSxy*a%dvdx + a%expoSyy*a%dvdy)- a%dvdy )*2.0_rp 
            elextSEstab4(3) = elextSEstab4(3) - auxtrac*aux_lamb*((a%expoSxx*a%dvdx + a%expoSyy*a%dudy)- (a%dudy + a%dvdx))*0.5_rp*2.0_rp

          end if
                   
        end if 
      
      if(php%kfl_advec == 0)then
         convection=0
      end if
      
      !Continuity force term
      elextC   = elextC + a%dudx + a%dvdy 
    
      !Force term 
      elext(1) = elext(1)   + dudtforce  + a%dpdx  + convection*(a%xdens*(a%u*a%dudx + a%v*a%dudy)) - (D1x + D2x)    
      elext(2) = elext(2)   + dvdtforce  + a%dpdy  + convection*(a%xdens*(a%u*a%dvdx + a%v*a%dvdy)) - (D1y + D2y)    
    

      elextEstab(1) = elextEstab(1) + convection*(a%xdens*(a%u*a%dudx + a%v*a%dudy))  
      elextEstab(2) = elextEstab(2) + convection*(a%xdens*(a%u*a%dvdx + a%v*a%dvdy))

      elextEstab2(1) = elextEstab2(1) +  a%dpdx
      elextEstab2(2) = elextEstab2(2) +  a%dpdy
      
      elextEstab3(1) = elextEstab3(1) - (a%DS1dx + a%DS3dy)*a%auxL
      elextEstab3(2) = elextEstab3(2) - (a%DS3dx + a%DS2dy)*a%auxL
      
      elextEstab4(1) = elextEstab4(1) - aux3*(a%d2udx + a%d2udy)
      elextEstab4(2) = elextEstab4(2) - aux3*(a%d2vdx + a%d2vdy)
      
      elextEstab5(1)= elextEstab5(1)   + a%auxL*a%expoSxx -a%auxL
      elextEstab5(2)= elextEstab5(2)   + a%auxL*a%expoSyy -a%auxL
      elextEstab5(3)= elextEstab5(3)   + a%auxL*a%expoSxy
       
     

   end subroutine
     

end module
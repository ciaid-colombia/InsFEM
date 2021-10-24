module Mod_NscExacso
   use typre
   use Mod_NSCompressible
  
   implicit none   
   type NscExacso
   
   ! work variables
   real(rp)    :: x,y,z,d,u,v,w,p,t
   real(rp)    :: dddx,dddy,dddz,d2ddx,d2ddy,d2ddz,d2ddxdy,d2ddxdz,d2ddydz
   real(rp)    :: dudx,dudy,dudz,d2udx,d2udy,d2udz,d2udxdy,d2udxdz,d2udydz
   real(rp)    :: dvdx,dvdy,dvdz,d2vdx,d2vdy,d2vdz,d2vdxdy,d2vdxdz,d2vdydz
   real(rp)    :: dwdx,dwdy,dwdz,d2wdx,d2wdy,d2wdz,d2wdxdy,d2wdxdz,d2wdydz
   real(rp)    :: dpdx,dpdy,dpdz,d2pdx,d2pdy,d2pdz,d2pdxdy,d2pdxdz,d2pdydz
   real(rp)    :: dtdx,dtdy,dtdz,d2tdx,d2tdy,d2tdz,d2tdxdy,d2tdxdz,d2tdydz
   real(rp)    :: acvis,actco,accph,accvh
   
contains   
      
      procedure :: nsc_ComputeSolution
      procedure :: nsc_GetDensity
      procedure :: nsc_GetMomentum
      procedure :: nsc_GetEnergy
      procedure :: nsc_GetConservativeForce
      procedure :: nsc_GetPrimitiveForce
      procedure :: nsc_GetTemperature
      procedure :: nsc_GetVelocity
      procedure :: nsc_GetPressure
   
   end type
   
contains

   subroutine nsc_ComputeSolution(a,ndime,gpcod,php)      
      use typre
      implicit none
      class(NscExacso) :: a
      Class(NSCompressibleProblem) :: php 
      integer(ip), intent(in) :: ndime      
      real(rp),    intent(in) :: gpcod(ndime)   
      real(rp) :: vis,tco,cp,cv,pi
      real(rp)    :: r,phi,alpha,omega
      real(rp)    :: alpsq,alpcu
      real(rp)    :: cosphisq,cosphicu
      real(rp)    :: sinphisq,sinphicu
      real(rp)    :: phi2,sum1,sum2,sum3,sum4, phi2prima,phi2prima3
      
      !Coordinates and initializations.

      a%x = gpcod(1)
      a%y = gpcod(2)
        
      a%d =     0.0_rp 
      a%dddx =  0.0_rp 
      a%dddy =  0.0_rp
      a%dddz =  0.0_rp
      a%d2ddx=  0.0_rp
      a%d2ddy=  0.0_rp
      a%d2ddz=  0.0_rp
      a%d2ddxdy=0.0_rp
      a%d2ddxdz=0.0_rp
      a%d2ddydz=0.0_rp

      a%u =     0.0_rp
      a%v =     0.0_rp
      a%w =     0.0_rp
      a%dudx =  0.0_rp
      a%dudy =  0.0_rp
      a%dudz =  0.0_rp
      a%dvdx =  0.0_rp
      a%dvdy =  0.0_rp
      a%dvdz =  0.0_rp
      a%dwdx =  0.0_rp
      a%dwdy =  0.0_rp
      a%dwdz =  0.0_rp
      a%d2udx = 0.0_rp
      a%d2udy = 0.0_rp
      a%d2udz = 0.0_rp
      a%d2vdx = 0.0_rp
      a%d2vdy = 0.0_rp
      a%d2vdz = 0.0_rp
      a%d2wdx = 0.0_rp
      a%d2wdy = 0.0_rp
      a%d2wdz = 0.0_rp
      a%d2udxdy = 0.0_rp
      a%d2udxdz = 0.0_rp
      a%d2udydz = 0.0_rp
      a%d2vdxdy = 0.0_rp
      a%d2vdxdz = 0.0_rp
      a%d2vdydz = 0.0_rp
      a%d2wdxdy = 0.0_rp
      a%d2wdxdz = 0.0_rp
      a%d2wdydz = 0.0_rp

      a%p =     0.0_rp
      a%dpdx =  0.0_rp 
      a%dpdy =  0.0_rp
      a%dpdz =  0.0_rp
      a%d2pdx=  0.0_rp
      a%d2pdy=  0.0_rp
      a%d2pdz=  0.0_rp
      a%d2pdxdy=0.0_rp
      a%d2pdxdz=0.0_rp
      a%d2pdydz=0.0_rp
      
      a%t =     0.0_rp
      a%dtdx =  0.0_rp 
      a%dtdy =  0.0_rp
      a%dtdz =  0.0_rp
      a%d2tdx=  0.0_rp
      a%d2tdy=  0.0_rp
      a%d2tdz=  0.0_rp
      a%d2tdxdy=0.0_rp
      a%d2tdxdz=0.0_rp
      a%d2tdydz=0.0_rp

      call php%GetPhysicalParameters(vis,tco,cp,cv)
      a%acvis=vis
      a%actco=tco
      a%accph=cp
      a%accvh=cv

      pi= 4.0_rp*atan(1.)        
      
      ! Obtain unknowns and derivatives according to the exact solution (Kexac)

      if(php%kfl_exacs==1) then !domain [0,1]x[0,1]
      
         ! density and derivatives
         a%d =     2.0_rp*a%x*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)+1 
         a%dddx =  2.0_rp*a%y*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)
         a%dddy =  2.0_rp*a%x*(2.0_rp*a%y-1.0_rp)*(a%x-1.0_rp)
         a%d2ddx=  4.0_rp*a%y*(a%y-1.0_rp)
         a%d2ddy=  4.0_rp*a%x*(a%x-1.0_rp)
         a%d2ddxdy=2.0_rp*(2.0_rp*a%y-1.0_rp)*(2.0_rp*a%x-1.0_rp)

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
         a%p     =  a%x*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)+10_rp
         a%dpdx  =  a%x*a%y*a%y*(a%y-1.0_rp)*(3.0_rp*a%x-2.0_rp)
         a%dpdy  =  a%x*a%x*a%y*(3.0_rp*a%y-2.0_rp)*(a%x-1.0_rp)
         a%d2pdx =  2.0_rp*a%y*a%y*(a%y-1.0_rp)*(3.0_rp*a%x-2.0_rp)                 
         a%d2pdy =  2.0_rp*a%x*a%x*(3.0_rp*a%y-2.0_rp)*(a%x-1.0_rp)
         a%d2pdxdy= a%x*a%y*(3.0_rp*a%y-2.0_rp)*(3.0_rp*a%x-2.0_rp)

         ! Temperature and derivatives
         a%t      = a%p/(a%d*(a%accph-a%accvh))
         a%dtdx   = (a%dpdx-a%p*a%dddx/a%d)/(a%d*(a%accph-a%accvh))
         a%dtdy   = (a%dpdy-a%p*a%dddy/a%d)/(a%d*(a%accph-a%accvh))

      else if(php%kfl_exacs==2) then !domain [0,1]x[0,1]

         !pressure and derivatives
         a%p     =  a%x*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)
         a%dpdx  =  a%x*a%y*a%y*(a%y-1.0_rp)*(3.0_rp*a%x-2.0_rp)
         a%dpdy  =  a%x*a%x*a%y*(3.0_rp*a%y-2.0_rp)*(a%x-1.0_rp)
         a%d2pdx =  2.0_rp*a%y*a%y*(a%y-1.0_rp)*(3.0_rp*a%x-1.0_rp)                 
         a%d2pdy =  2.0_rp*a%x*a%x*(3.0_rp*a%y-1.0_rp)*(a%x-1.0_rp)
         a%d2pdxdy= a%x*a%y*(3.0_rp*a%y-2.0_rp)*(3.0_rp*a%x-2.0_rp)

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

         ! Temperature and derivatives
         a%t =     2.0_rp*a%x*a%x*a%y*(a%x*a%x-1.0_rp)*(a%y-1.0_rp)
         a%dtdx =  4.0_rp*a%x*a%y*(2.0_rp*a%x*a%x-1.0_rp)*(a%y-1.0_rp)
         a%dtdy =  2.0_rp*a%x*a%x*(2.0_rp*a%y-1.0_rp)*(a%x-1.0_rp)*(a%x+1.0_rp)
         a%d2tdx=  4.0_rp*a%y*(a%y-1.0_rp)*(6.0_rp*a%x*a%x-1.0_rp)
         a%d2tdy=  4.0_rp*a%x*a%x*(a%x-1.0_rp)*(a%x+1.0_rp)
         a%d2tdxdy=4.0_rp*a%x*(2.0_rp*a%y-1.0_rp)*(2.0_rp*a%x*a%x-1.0_rp)

         if(ndime==3) then !domain [0,1]x[0,1]x[0,1]       

            a%z = gpcod(3)
   
            !pressure and derivatives
            a%p     =  a%p*a%z*a%z*(1.0_rp-a%z)
            a%dpdx  =  a%dpdx*a%z*a%z*(1.0_rp-a%z)
            a%dpdy  =  a%dpdy*a%z*a%z*(1.0_rp-a%z)
            a%d2pdx =  a%d2pdx*a%z*a%z*(1.0_rp-a%z)                  
            a%d2pdy =  a%d2pdy*a%z*a%z*(1.0_rp-a%z)
            a%d2pdxdy= a%d2pdxdy*a%z*a%z*(1.0_rp-a%z)
   
            a%dpdz  =  -a%x*a%x*a%y*a%y*a%z*(a%x-1.0_rp)*(a%y-1.0_rp)*(3.0_rp*a%z-2.0_rp)
            a%d2pdz =  -2.0_rp*a%x*a%x*a%y*a%y*(a%x-1.0_rp)*(a%y-1.0_rp)*(3.0_rp*a%z-1.0_rp)                  
            a%d2pdxdz= -a%x*a%y*a%y*a%z*(3.0_rp*a%z-2.0_rp)*(3.0_rp*a%x-2.0_rp)*(a%y-1.0_rp)
            a%d2pdydz= -a%x*a%x*a%y*a%z*(3.0_rp*a%z-2.0_rp)*(3.0_rp*a%y-2.0_rp)*(a%x-1.0_rp)

            ! velocity
            a%u     = a%u*a%z*a%z*(1.0_rp-a%z)  
            a%v     = a%v*a%z*a%z*(1.0_rp-a%z)  
            a%w     = 2.0_rp*a%x*a%z*a%z*(a%x-1.0_rp)*(1.0_rp-a%z)*(1.0_rp-a%z)*(2.0_rp*a%x-1.0_rp)*a%y*a%y*(a%y-1.0_rp)

            !velocity derivatives
            !u component
            a%dudx    =  a%dudx*a%z*a%z*(1.0_rp-a%z)    
            a%dudy    =  a%dudy*a%z*a%z*(1.0_rp-a%z)    
            a%d2udx   =  a%d2udx*a%z*a%z*(1.0_rp-a%z)   
            a%d2udy   =  a%d2udy*a%z*a%z*(1.0_rp-a%z)   
            a%d2udxdy =  a%d2udxdy*a%z*a%z*(1.0_rp-a%z) 

            a%dudz  =  -2.0_rp*(a%x-1.0_rp)*(a%x-1.0_rp)*a%x*a%x*(a%y-1.0_rp)*a%y*(2.0_rp*a%y-1.0_rp)*a%z*(3.0_rp*a%z-2.0_rp)
            a%d2udz =  -4.0_rp*(a%x-1.0_rp)*(a%x-1.0_rp)*a%x*a%x*(a%y-1.0_rp)*a%y*(2.0_rp*a%y-1.0_rp)*(3.0_rp*a%z-1.0_rp)
            a%d2udxdz= -4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*(2.0_rp*a%y-1.0_rp)*a%z*(3.0_rp*a%z-2.0_rp)
            a%d2udydz= -2.0_rp*(a%x-1.0_rp)*(a%x-1.0_rp)*a%x*a%x*(6.0_rp*a%y*a%y-6.0_rp*a%y+1.0_rp)*a%z*(3.0_rp*a%z-2.0_rp)

            !v component
            a%dvdx    = a%dvdx*a%z*a%z*(1.0_rp-a%z)   
            a%dvdy    = a%dvdy*a%z*a%z*(1.0_rp-a%z)   
            a%d2vdx   = a%d2vdx*a%z*a%z*(1.0_rp-a%z)  
            a%d2vdy   = a%d2vdy*a%z*a%z*(1.0_rp-a%z)  
            a%d2vdxdy = a%d2vdxdy*a%z*a%z*(1.0_rp-a%z)
   
            a%dvdz  =  2.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*a%z*(3.0_rp*a%z-2.0_rp) 
            a%d2vdz =  4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(3.0_rp*a%z-1.0_rp)
            a%d2vdxdz= 2.0_rp*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*(a%y-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*a%z*(3.0_rp*a%z-2.0_rp)      
            a%d2vdydz= 4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*(2.0_rp*a%y-1.0_rp)*a%z*(3.0_rp*a%z-2.0_rp)

            !w component
            a%dwdx    = 2.0_rp*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(a%z-1.0_rp)*(a%z-1.0_rp)*a%z*a%z
            a%dwdy    = 2.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*a%y*(3.0_rp*a%y-2.0_rp)*(a%z-1.0_rp)*(a%z-1.0_rp)*a%z*a%z
            a%dwdz    = 4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(a%z-1.0_rp)*a%z*(2.0_rp*a%z-1.0_rp)
            a%d2wdx   = 12.0_rp*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(a%z-1.0_rp)*(a%z-1.0_rp)*a%z*a%z
            a%d2wdy   = 4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(3.0_rp*a%y-1.0_rp)*(a%z-1.0_rp)*(a%z-1.0_rp)*a%z*a%z
            a%d2wdz   = 4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(6.0_rp*a%z*a%z-6.0_rp*a%z+1.0_rp)
            a%d2wdxdy = 2.0_rp*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*a%y*(3.0_rp*a%y-2.0_rp)*(a%z-1.0_rp)*(a%z-1.0_rp)*a%z*a%z
            a%d2wdxdz = 4.0_rp*(6.0_rp*a%x*a%x-6.0_rp*a%x+1.0_rp)*(a%y-1.0_rp)*a%y*a%y*(a%z-1.0_rp)*a%z*(2.0_rp*a%z-1.0_rp)
            a%d2wdydz = 4.0_rp*(a%x-1.0_rp)*a%x*(2.0_rp*a%x-1.0_rp)*a%y*(3.0_rp*a%y-2.0_rp)*(a%z-1.0_rp)*a%z*(2.0_rp*a%z-1.0_rp)

            ! Temperature and derivatives
            a%t      = a%t*a%z*a%z*(1.0_rp-a%z)      
            a%dtdx   = a%dtdx*a%z*a%z*(1.0_rp-a%z)     
            a%dtdy   = a%dtdy*a%z*a%z*(1.0_rp-a%z)    
            a%d2tdx  = a%d2tdx*a%z*a%z*(1.0_rp-a%z)    
            a%d2tdy  = a%d2tdy*a%z*a%z*(1.0_rp-a%z)    
            a%d2tdxdy= a%d2tdxdy*a%z*a%z*(1.0_rp-a%z)
   
            a%dtdz  =  -2.0_rp*(a%x-1.0_rp)*a%x*a%x*(a%x+1.0_rp)*(a%y-1.0_rp)*a%y*a%z*(3.0_rp*a%z-2.0_rp)
            a%d2tdz =  -4.0_rp*(a%x-1.0_rp)*a%x*a%x*(a%x+1.0_rp)*(a%y-1.0_rp)*a%y*(3.0_rp*a%z-1.0_rp)
            a%d2tdxdz= -4.0_rp*a%x*(2.0_rp*a%x*a%x-1.0_rp)*(a%y-1.0_rp)*a%y*a%z*(3.0_rp*a%z-2.0_rp)               
            a%d2tdydz= -2.0_rp*(a%x-1.0_rp)*a%x*a%x*(a%x+1.0_rp)*(2.0_rp*a%y-1.0_rp)*a%z*(3.0_rp*a%z-2.0_rp)

      end if 

      !Density and derivatives
      a%d =  (a%p+php%relpre)/((a%accph-a%accvh)*(a%t+php%reltem))
      a%dddx   = (a%dpdx-(a%p+php%relpre)*a%dtdx/(a%t+php%reltem))/((a%t+php%reltem)*(a%accph-a%accvh))
      a%dddy   = (a%dpdy-(a%p+php%relpre)*a%dtdy/(a%t+php%reltem))/((a%t+php%reltem)*(a%accph-a%accvh))
      a%d2ddx   = (a%d2pdx-(a%p+php%relpre)*a%d2tdx/(a%t+php%reltem))/((a%t+php%reltem)*(a%accph-a%accvh))&
                 - 2.0_rp*a%dpdx*a%dtdx/((a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh)) &
                 + 2.0_rp*a%dtdx*a%dtdx*(a%p+php%relpre)/((a%t+php%reltem)*(a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh))
      a%d2ddx   = (a%d2pdy-(a%p+php%relpre)*a%d2tdy/(a%t+php%reltem))/((a%t+php%reltem)*(a%accph-a%accvh))&
                 - 2.0_rp*a%dpdy*a%dtdy/((a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh)) &
                 + 2.0_rp*a%dtdy*a%dtdy*(a%p+php%relpre)/((a%t+php%reltem)*(a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh))
      a%d2ddxdy = a%d2pdxdy/((a%t+php%reltem)*(a%accph-a%accvh)) - (a%p+php%relpre)*a%d2tdxdy/((a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh))&
                 - (a%dpdx*a%dtdy+a%dtdx*a%dpdy)/((a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh)) &
                 + 2.0_rp*a%dtdx*a%dtdy*(a%p+php%relpre)/((a%t+php%reltem)*(a%t+php%reltem)*(a%t+php%reltem)*(a%accph-a%accvh))

      else if(php%kfl_exacs==3) then !domain [0,1]x[0,1]

         !pressure and derivatives
         a%p     =  a%x 
         a%dpdx  =  1.0_rp
         a%dpdy  =  0.0_rp  
  
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

         !Temperature and derivatives
         a%t     =  a%y 
         a%dtdx  =  0.0_rp
         a%dtdy  =  1.0_rp  

         !Density and derivatives
         a%d =  (a%p+php%relpre)/((a%accph-a%accvh)*(a%t+php%reltem))

      else if(php%kfl_exacs==4) then !domain [0,1]x[0,1]

         !pressure and derivatives
         a%p     =  a%x*a%x 
         a%dpdx  =  2.0_rp*a%x
         a%dpdy  =  0.0_rp  
         a%d2pdx  =  2.0_rp
          
         ! velocity
         a%u       =  4.0_rp*a%x*a%x + 6 
         a%v       =  -4.0_rp*a%y*a%y + 6 
         !velocity derivatives
         !u component
         a%dudx    =  8.0_rp*a%x
         a%dudy    =  0.0_rp
         a%d2udx   =  8.0_rp
         a%d2udy   =  0.0_rp
         a%d2udxdy =  0.0_rp
         !v component
         a%dvdx    = -8.0_rp*a%x
         a%dvdy    = -4.0_rp
         a%d2vdx   =  0.0_rp
         a%d2vdy   =  -8.0_rp
         a%d2vdxdy =  0.0_rp

         !Temperature and derivatives
         a%t     =  a%y*a%y 
         a%dtdx  =  0.0_rp  
         a%dtdy  =  2.0_rp*a%y
         a%d2tdy  =  2.0_rp

         !Density and derivatives
         a%d =  (a%p+php%relpre)/((a%accph-a%accvh)*(a%t+php%reltem))

      else if(php%kfl_exacs==5) then !domain ([-1,1]x[-1,1]) \ ([-1,0]x[0,1]) (Stokes with singularity)

         ! spherical coordinates
         r      = sqrt(a%x*a%x + a%y*a%y)
         phi    = atan(a%y/a%x)
         if(a%x<epsilon(0.0_rp)) phi = phi + pi
         if(abs(a%x)<epsilon(0.0_rp)) then
            if(a%y>epsilon(0.0_rp)) phi = pi/2
            if(a%y<epsilon(0.0_rp)) phi = 3*pi/2
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
            
            
            
            a%u = (r**alpha)*(cos(phi)*phi2prima+(1+alpha)*sin(phi)*phi2);
            a%v = (r**alpha)*(sin(phi)*phi2prima-(1+alpha)*cos(phi)*phi2);
                     
            a%p = -r**(alpha-1)*((1+alpha)**2*phi2prima+phi2prima3)/(1-alpha);
            a%dpdx = -2*alpha*r**(alpha - 2)*(cos(2*phi + alpha*omega - alpha*phi) - 2*sin(phi*(alpha - 2)) + cos(2*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)))
            a%dpdy = -2*alpha*r**(alpha - 2)*(sin(2*phi + alpha*omega - alpha*phi) + sin(2*phi - alpha*omega - alpha*phi) - 2*cos(phi*(alpha - 2)) + 2*alpha*cos(phi*(alpha - 2)))
            
            
            a%dudx = (alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2
            a%dudy = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) - 6*cos(phi - alpha*phi) + sin(phi + alpha*omega - alpha*phi) + sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) + 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
            a%dvdx = (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) + 2*cos(phi - alpha*phi) - 3*sin(phi + alpha*omega - alpha*phi) - 3*sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) - 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
            a%dvdy = -(alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2
            a%d2udx = (alpha*r**(alpha - 2)*(4*sin(4*phi - alpha*phi) - 2*cos(2*phi + alpha*omega - alpha*phi) - 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(2*phi - alpha*phi) - 6*alpha*sin(4*phi - alpha*phi) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(2*phi - alpha*phi) + 2*alpha**2*sin(4*phi - alpha*phi)))/2
            a%d2udxdy = (alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2
            a%d2udy = -(alpha*r**(alpha - 2)*(2*cos(2*phi + alpha*omega - alpha*phi) - 4*sin(phi*(alpha - 4)) - 8*sin(phi*(alpha - 2)) + 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha**2*sin(phi*(alpha - 2)) - 2*alpha**2*sin(phi*(alpha - 4)) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) + 6*alpha*sin(phi*(alpha - 2)) + 6*alpha*sin(phi*(alpha - 4))))/2
            a%d2vdx = -(alpha*r**(alpha - 2)*(4*sin(2*phi + alpha*omega - alpha*phi) + 4*sin(2*phi - alpha*omega - alpha*phi) - 2*sin(4*phi + alpha*omega - alpha*phi) - 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) + 4*cos(phi*(alpha - 4)) - 2*alpha**2*cos(phi*(alpha - 2)) + 2*alpha**2*cos(phi*(alpha - 4)) + 6*alpha*cos(phi*(alpha - 2)) - 6*alpha*cos(phi*(alpha - 4)) - alpha*sin(2*phi + alpha*omega - alpha*phi) - alpha*sin(2*phi - alpha*omega - alpha*phi) + alpha*sin(4*phi + alpha*omega - alpha*phi) + alpha*sin(4*phi - alpha*omega - alpha*phi)))/2
            a%d2vdxdy = (alpha*r**(alpha - 2)*(4*sin(phi*(alpha - 4)) + 2*cos(2*phi + alpha*omega - alpha*phi) + 2*cos(2*phi - alpha*omega - alpha*phi) - 2*cos(4*phi + alpha*omega - alpha*phi) - 2*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(phi*(alpha - 2)) + 2*alpha**2*sin(phi*(alpha - 4)) - alpha*cos(2*phi + alpha*omega - alpha*phi) - alpha*cos(2*phi - alpha*omega - alpha*phi) + alpha*cos(4*phi + alpha*omega - alpha*phi) + alpha*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 2)) - 6*alpha*sin(phi*(alpha - 4))))/2
            a%d2vdy = -(alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2
         endif

         !Temperature and derivatives
         a%t = (r**alpha)*(cos(phi)*phi2prima+(1+alpha)*sin(phi)*phi2);
         a%dtdx = (alpha*r**(alpha - 1)*(2*sin(phi*(alpha - 1)) + 2*sin(phi*(alpha - 3)) + cos(phi + alpha*omega - alpha*phi) + cos(phi - alpha*omega - alpha*phi) - cos(3*phi + alpha*omega - alpha*phi) - cos(3*phi - alpha*omega - alpha*phi) + 2*alpha*sin(phi*(alpha - 1)) - 2*alpha*sin(phi*(alpha - 3))))/2
         a%dtdy= (alpha*r**(alpha - 1)*(sin(3*phi + alpha*omega - alpha*phi) + sin(3*phi - alpha*omega - alpha*phi) - 6*cos(phi - alpha*phi) + sin(phi + alpha*omega - alpha*phi) + sin(phi - alpha*omega - alpha*phi) - 2*cos(3*phi - alpha*phi) + 4*alpha*cos(phi - alpha*phi) + alpha*sin(phi + alpha*omega - alpha*phi) + alpha*sin(phi - alpha*omega - alpha*phi) + 4*alpha*cos(3*phi - alpha*phi) + 2*alpha**2*cos(phi - alpha*phi) - alpha*sin(3*phi + alpha*omega - alpha*phi) - alpha*sin(3*phi - alpha*omega - alpha*phi) - 2*alpha**2*cos(3*phi - alpha*phi)))/(2*(alpha - 1))
         a%d2tdx= (alpha*r**(alpha - 2)*(4*sin(4*phi - alpha*phi) - 2*cos(2*phi + alpha*omega - alpha*phi) - 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha*sin(2*phi - alpha*phi) - 6*alpha*sin(4*phi - alpha*phi) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) - 2*alpha**2*sin(2*phi - alpha*phi) + 2*alpha**2*sin(4*phi - alpha*phi)))/2
         a%d2tdy= -(alpha*r**(alpha - 2)*(2*cos(2*phi + alpha*omega - alpha*phi) - 4*sin(phi*(alpha - 4)) - 8*sin(phi*(alpha - 2)) + 2*cos(2*phi - alpha*omega - alpha*phi) + 2*cos(4*phi + alpha*omega - alpha*phi) + 2*cos(4*phi - alpha*omega - alpha*phi) + 2*alpha**2*sin(phi*(alpha - 2)) - 2*alpha**2*sin(phi*(alpha - 4)) + alpha*cos(2*phi + alpha*omega - alpha*phi) + alpha*cos(2*phi - alpha*omega - alpha*phi) - alpha*cos(4*phi + alpha*omega - alpha*phi) - alpha*cos(4*phi - alpha*omega - alpha*phi) + 6*alpha*sin(phi*(alpha - 2)) + 6*alpha*sin(phi*(alpha - 4))))/2
         a%d2tdxdy= (alpha*r**(alpha - 2)*(2*sin(4*phi + alpha*omega - alpha*phi) + 2*sin(4*phi - alpha*omega - alpha*phi) - 4*cos(phi*(alpha - 2)) - 4*cos(phi*(alpha - 4)) + 2*alpha**2*cos(phi*(alpha - 2)) - 2*alpha**2*cos(phi*(alpha - 4)) + 2*alpha*cos(phi*(alpha - 2)) + 6*alpha*cos(phi*(alpha - 4)) + alpha*sin(2*phi + alpha*omega - alpha*phi) + alpha*sin(2*phi - alpha*omega - alpha*phi) - alpha*sin(4*phi + alpha*omega - alpha*phi) - alpha*sin(4*phi - alpha*omega - alpha*phi)))/2

         !Density and derivatives
         a%d =  (a%p+php%relpre)/((a%accph-a%accvh)*(a%t+php%reltem))

      else if(php%kfl_exacs>5) then
         call runend('nsc_ComputeSolution exact solution number not defined')      
         !define your exact solution
      end if
   end subroutine
   
   
   subroutine nsc_GetDensity(a,ndime,exden,exdeg)
      use typre
      implicit none
      class(NscExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: exden,exdeg(ndime) 
      
      exden      = a%d
      exdeg(1)   = a%dddx
      exdeg(2)   = a%dddy
      if(ndime==3) exdeg(3) = a%dddz
         
   end subroutine

   subroutine nsc_GetVelocity(a,ndime,exvel,exveg)
      use typre
      implicit none
      class(NscExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exvel(ndime),exveg(ndime,ndime) 
      
      exvel(1)   = a%u
      exvel(2)   = a%v
      exveg(1,1) = a%dudx 
      exveg(2,1) = a%dvdx
      exveg(1,2) = a%dudy
      exveg(2,2) = a%dvdy
         
      if(ndime==3) then
         exvel(3) = a%w
         exveg(1,3) = a%dudz
         exveg(2,3) = a%dvdz
         exveg(3,1) = a%dwdx 
         exveg(3,2) = a%dwdy
         exveg(3,3) = a%dwdz
      end if

   end subroutine

   subroutine nsc_GetMomentum(a,ndime,exmom,exmog)
      use typre
      implicit none
      class(NscExacso) :: a      
      integer(ip), intent(in)   :: ndime
      real(rp),    intent(out)  :: exmom(ndime),exmog(ndime,ndime) 
      
      exmom(1)   = a%u*a%d 
      exmom(2)   = a%v*a%d
      exmog(1,1) = a%dudx*a%d + a%u*a%dddx  
      exmog(2,1) = a%dvdx*a%d + a%v*a%dddx
      exmog(1,2) = a%dudy*a%d + a%u*a%dddy
      exmog(2,2) = a%dvdy*a%d + a%v*a%dddy
         
   end subroutine

   subroutine nsc_GetEnergy(a,ndime,exene,exeng)
      use typre
      implicit none
      class(NscExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: exene,exeng(ndime) 
      
      exene      = a%accvh*a%p/(a%accph-a%accvh)+0.5_rp*a%d*(a%u*a%u+a%v*a%v)
      exeng(1)   = a%accvh*a%dpdx/(a%accph-a%accvh)+0.5_rp*a%dddx*(a%u*a%u+a%v*a%v)+a%d*(a%u*a%dudx+a%v*a%dvdx)
      exeng(2)   = a%accvh*a%dpdy/(a%accph-a%accvh)+0.5_rp*a%dddy*(a%u*a%u+a%v*a%v)+a%d*(a%u*a%dudy+a%v*a%dvdy)
         
   end subroutine

   subroutine nsc_GetTemperature(a,ndime,extem,exteg)
      use typre
      implicit none
      class(NscExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: extem,exteg(ndime) 
      
      extem      = a%t
      exteg(1)   = a%dtdx
      exteg(2)   = a%dtdy
      if(ndime==3) exteg(3) = a%dtdz
         
   end subroutine

   subroutine nsc_GetPressure(a,ndime,expre,exprg)
      use typre
      implicit none
      class(NscExacso) :: a
      integer(ip), intent(in)  :: ndime
      real(rp),    intent(out) :: expre,exprg(ndime) 
      
      expre      = a%p
      exprg(1)   = a%dpdx
      exprg(2)   = a%dpdy
      if(ndime==3) exprg(3) = a%dpdz
         
   end subroutine

   subroutine nsc_GetConservativeForce(a,ndime,elexd,elexm,elexe,php)
      use typre
      implicit none
      class(NscExacso) :: a
      class(NSCompressibleProblem) :: php 
      integer(ip), intent(in)  :: ndime
      real(rp), intent(inout)  :: elexd(1),elexm(ndime),elexe(1)
      real(rp) :: Cd,Cu,Cv,Du,Dv,Ce1,Ce2,De1,De2,Deq
      real(rp) :: aux1,aux2,aux3,aux4
      
      Cd =0.0_rp
      Cu =0.0_rp
      Cv =0.0_rp
      Du  =0.0_rp
      Dv  =0.0_rp
      Ce1  =0.0_rp
      Ce2  =0.0_rp 
      De1  =0.0_rp
      De2  =0.0_rp
      Deq  =0.0_rp
      aux1 = 0.0_rp
      aux2 = 0.0_rp
      aux3 = 0.0_rp
      aux4 = 0.0_rp


      if (php%kfl_advec == 1) then
         !Mass convective term      
         Cd= a%d*a%dudx + a%u*a%dddx + a%d*a%dvdy + a%v*a%dddy
   
         !Momentum convective term      
         Cu= a%u*(a%u*a%dddx+a%v*a%dddy) + a%d*a%u*(a%dudx+a%dvdy) + a%d*(a%u*a%dudx+a%v*a%dudy) + a%dpdx
         Cv= a%v*(a%u*a%dddx+a%v*a%dddy) + a%d*a%v*(a%dudx+a%dvdy) + a%d*(a%u*a%dvdx+a%v*a%dvdy) + a%dpdy
   
         !Energy convective term      
         !dj rho u_j (cv*thet + P + 0.5*u_i*u_i) = ()*(uj*dj rho + rho*dj uj) + rho*uj*(cv*dj thet + dj P + ui dj ui)
         aux1= a%accvh*a%p/(a%d*(a%accph-a%accvh)) + a%p + (a%u*a%u+a%v*a%v)/2.0_rp 
         aux2= a%d*(a%dudx+a%dvdy) + a%u*a%dddx+a%v*a%dddy 
         Ce1 = aux1*aux2
         aux3= (a%accvh/(a%accph-a%accvh))*(a%u*a%dpdx+a%v*a%dpdy)&
              -(a%accvh*a%p/(a%d*(a%accph-a%accvh)))*(a%u*a%dddx+a%v*a%dddy)
         aux4= a%d*(a%u*(a%u*a%dudx+a%v*a%dudy)+a%v*(a%u*a%dvdx+a%v*a%dvdy))
         Ce2 = aux3 + a%d*(a%u*a%dpdx+a%v*a%dpdy) + aux4

      end if
   
      if (php%kfl_visco == 1) then
         !Momentum diffusive term      
         Du= a%acvis*(a%d2udx+a%d2vdxdy + a%d2udx+a%d2udy)&
             - 2.0_rp*a%acvis*(a%d2udx+a%d2vdxdy)/3.0_rp 
         Dv= a%acvis*(a%d2udxdy+a%d2vdy + a%d2vdx+a%d2vdy)&
             - 2.0_rp*a%acvis*(a%d2udxdy+a%d2vdy)/3.0_rp
   
         !Energy viscous dissipation term      
         !dj (-ui*tau_ij) = -tau_ij*dj ui - ui*dj tau_ij = -De1 -De2 */
         De1 =  2.0_rp*a%acvis*(a%dudx-(a%dudx+a%dvdy)/3.0_rp)*a%dudx &
                + a%acvis*(a%dudy+a%dvdx)*a%dudy + a%acvis*(a%dudy+a%dvdx)*a%dvdx&
               +2.0_rp*a%acvis*(a%dvdy-(a%dudx+a%dvdy)/3.0_rp)*a%dvdy 
   
         De2 = a%u*(a%acvis*(a%d2udx+a%d2vdxdy+a%d2udx+a%d2udy)-2.0_rp*a%acvis*(a%d2udx+a%d2vdxdy)/3.0_rp)& 
              +a%v*(a%acvis*(a%d2udxdy+a%d2vdy+a%d2vdx+a%d2vdy)-2.0_rp*a%acvis*(a%d2udxdy+a%d2vdy)/3.0_rp) 
    
         !Energy heat flux term      
         !di qi =  - lambd/R * (-(2/rho^2)*di P*di rho 
         !                      +(2P/rho^3)*di rho*di rho + (1/rho)*di di P - (P/rho^2)*di di rho)
          Deq =  - (a%actco/(a%accph-a%accvh))*(-2.0_rp*(a%dpdx*a%dddx+a%dpdy*a%dddy)/(a%d*a%d) &
                                          +2.0_rp*a%p*(a%dddx*a%dddx+a%dddy*a%dddy)/(a%d*a%d*a%d) &
                                          +(a%d2pdx+a%d2pdy)/a%d -a%p*(a%d2ddx+a%d2ddy)/(a%d*a%d))
      end if

      !Continuity force term
      elexd(1) = elexd(1) + Cd

      !Momentum force term
      elexm(1) = elexm(1) + Cu - Du  
      elexm(2) = elexm(2) + Cv - Dv  

      !Energy force term
      !dj (rho u_j (cv*thet + P + 0.5*u_i*u_i) - ui*tau_ij + qj) 

      elexe(1) =  elexe(1) + Ce1 + Ce2 - De1 - De2 + Deq 

   end subroutine

   subroutine nsc_GetPrimitiveForce(a,ndime,elexc,elexm,elexe,php)
      use typre
      implicit none
      class(NscExacso) :: a
      class(NSCompressibleProblem) :: php 
      integer(ip), intent(in)  :: ndime
      real(rp), intent(inout)  :: elexc(1),elexm(ndime),elexe(1)
      real(rp) :: Cpp,Cpu,Cpt,Ctp,Ctu,Ctt
      real(rp) :: Cup,Cuu,Cut,Cvp,Cvu,Cvt,Cwp,Cwu,Cwt
      real(rp) :: Du,Dv,Dw,De1,De2,Deq
      real(rp) :: auxterm,acalpha,acbeta
      
      Cpp =0.0_rp
      Cpu =0.0_rp
      Cpt =0.0_rp
      Cup =0.0_rp
      Cuu =0.0_rp
      Cut =0.0_rp
      Cvp =0.0_rp
      Cvu =0.0_rp
      Cvt =0.0_rp
      Cwp =0.0_rp
      Cwu =0.0_rp
      Cwt =0.0_rp
      Ctp =0.0_rp
      Ctu =0.0_rp
      Ctt =0.0_rp
      Du  =0.0_rp
      Dv  =0.0_rp
      Dw  =0.0_rp
      De1  =0.0_rp
      De2  =0.0_rp
      Deq  =0.0_rp

      auxterm = 0.0_rp
      acalpha = 0.0_rp
      acbeta =  0.0_rp

      !State law
      if (php%lawde /= 0) then
         if (php%lawde == 1) then
            auxterm = a%d*(a%accvh*(a%t+php%reltem)+(a%u*a%u+a%v*a%v+a%w*a%w)/2.0_rp)+(a%p+php%relpre)
            acalpha = 1.0_rp/(a%t+php%reltem)
            acbeta = 1.0_rp/(a%p+php%relpre)
         else if (php%lawde /= 1) then
            call runend('Nsc_exacsol: Non-ideal state law not ready')
         endif
      endif


      if (php%kfl_advec == 1) then
         !Continuity convective term      
         Cpp = a%d*acbeta*(a%u*a%dpdx + a%v*a%dpdy + a%w*a%dpdz)
         Cpu = a%d*(a%dudx + a%dvdy + a%dwdz)
         Cpt = -a%d*acalpha*(a%u*a%dtdx + a%v*a%dtdy + a%w*a%dtdz)
   
         !Momentum convective term      
         Cup = a%d*acbeta*a%u*(a%u*a%dpdx + a%v*a%dpdy + a%w*a%dpdz) + a%dpdx 
         Cuu = a%d*a%u*(a%dudx + a%dvdy + a%dwdz) + a%d*(a%u*a%dudx + a%v*a%dudy + a%w*a%dudz)
         Cut = -a%d*acalpha*a%u*(a%u*a%dtdx + a%v*a%dtdy + a%w*a%dtdz)
         Cvp = a%d*acbeta*a%v*(a%u*a%dpdx + a%v*a%dpdy + a%w*a%dpdz) + a%dpdy 
         Cvu = a%d*a%v*(a%dudx + a%dvdy + a%dwdz) + a%d*(a%u*a%dvdx + a%v*a%dvdy + a%w*a%dvdz)
         Cvt = -a%d*acalpha*a%v*(a%u*a%dtdx + a%v*a%dtdy + a%w*a%dtdz)
         Cwp = a%d*acbeta*a%w*(a%u*a%dpdx + a%v*a%dpdy + a%w*a%dpdz) + a%dpdz 
         Cwu = a%d*a%w*(a%dudx + a%dvdy + a%dwdz) + a%d*(a%u*a%dwdx + a%v*a%dwdy + a%w*a%dwdz)
         Cwt = -a%d*acalpha*a%w*(a%u*a%dtdx + a%v*a%dtdy + a%w*a%dtdz)
   
         !Energy convective term      
         Ctp = (auxterm*acbeta-acalpha*(a%t+php%reltem)+1.0_rp)*(a%u*a%dpdx + a%v*a%dpdy + a%w*a%dpdz)
         Ctu =  auxterm*(a%dudx + a%dvdy + a%dwdz) &
               + a%d*a%u*(a%u*a%dudx + a%v*a%dudy + a%w*a%dudz) &
               + a%d*a%v*(a%u*a%dvdx + a%v*a%dvdy + a%w*a%dvdz) &
               + a%d*a%w*(a%u*a%dwdx + a%v*a%dwdy + a%w*a%dwdz) 
         Ctt = (a%d*a%accph-auxterm*acalpha)*(a%u*a%dtdx + a%v*a%dtdy + a%w*a%dtdz)

      end if
   
      if (php%kfl_visco == 1) then
         !Momentum diffusive term second order     
         Du= a%d2udx+a%d2vdxdy+a%d2wdxdz + a%d2udx+a%d2udy+a%d2udz&
             - 2.0_rp*(a%d2udx+a%d2vdxdy+a%d2wdxdz)/3.0_rp 
         Dv= a%d2udxdy+a%d2vdy+a%d2wdydz + a%d2vdx+a%d2vdy+a%d2vdz&
             - 2.0_rp*(a%d2udxdy+a%d2vdy+a%d2wdydz)/3.0_rp
         Dw= a%d2udxdz+a%d2vdydz+a%d2wdz + a%d2wdx+a%d2wdy+a%d2wdz&
             - 2.0_rp*(a%d2udxdz+a%d2vdydz+a%d2wdz)/3.0_rp
   
         !Energy viscous dissipation term first order     
         De1 = (a%dudx+a%dudy+a%dudz)*(a%dudx+a%dudy+a%dudz) &
             + (a%dvdx+a%dvdy+a%dvdz)*(a%dvdx+a%dvdy+a%dvdz) &
             + (a%dwdx+a%dwdy+a%dwdz)*(a%dwdx+a%dwdy+a%dwdz) &
             + (a%dudx*a%dudx+a%dudy*a%dvdx+a%dudz*a%dwdx) &
             + (a%dvdx*a%dudy+a%dvdy*a%dvdy+a%dvdz*a%dwdy) &
             + (a%dwdx*a%dudz+a%dwdy*a%dvdz+a%dwdz*a%dwdz) &
             - 2.0_rp*(a%dudx+a%dvdy+a%dwdz)*(a%dudx+a%dvdy+a%dwdz)/3.0_rp
       
         !Energy viscous dissipation term second order 
         De2 = a%u*(a%d2udx+a%d2vdxdy+a%d2wdxdz+a%d2udx+a%d2udy+a%d2udz)&
             - a%u*2.0_rp*(a%d2udx+a%d2vdxdy+a%d2wdxdz)/3.0_rp& 
             + a%v*(a%d2udxdy+a%d2vdy+a%d2wdydz+a%d2vdx+a%d2vdy+a%d2vdz)&
             - a%v*2.0_rp*(a%d2udxdy+a%d2vdy+a%d2wdydz)/3.0_rp& 
             + a%w*(a%d2udxdz+a%d2vdydz+a%d2wdz+a%d2wdx+a%d2wdy+a%d2wdz)&
             - a%w*2.0_rp*(a%d2udxdz+a%d2vdydz+a%d2wdz)/3.0_rp 

         !Energy heat flux term second order
         Deq =  a%d2tdx+a%d2tdy+a%d2tdz

      end if

      !Continuity force term
      elexc(1) = elexc(1) + Cpp + Cpu + Cpt

      !Momentum force term
      elexm(1) = elexm(1) + Cup + Cuu + Cut - a%acvis*Du  
      elexm(2) = elexm(2) + Cvp + Cvu + Cvt - a%acvis*Dv  
      if(ndime==3) elexm(3) = elexm(3) + Cwp + Cwu + Cwt - a%acvis*Dw 

      !Energy force term
      elexe(1) =  elexe(1) + Ctp + Ctu + Ctt - a%acvis*De1 - a%acvis*De2 - a%actco*Deq 

   end subroutine


end module

module Mod_NscDamping
   use typre
   use Mod_NSCompressible
  
   implicit none   
   type NscDamping
   
   ! work variables
   real(rp)    :: x,y,z
   real(rp)    :: xo,yo,zo,xf,yf,zf
   real(rp)    :: cx,cy,cz
   real(rp)    :: p,u,v,w,t
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: den,mom(3),vel(3),ene
   real(rp)    :: coefficient
   
contains   
      
      procedure :: nsc_ComputeRectangularDamping
      procedure :: nsc_ComputeRadialDamping
      procedure :: nsc_GetDampingForce
      procedure :: nsc_pr_GetDampingForce
      procedure :: nsc_GetDampingCoefficient
   
   end type
   
contains

   subroutine nsc_ComputeRectangularDamping(a,ndime,gpcod,php)      
      use typre
      implicit none
      class(NscDamping) :: a
      Class(NSCompressibleProblem) :: php 
      integer(ip), intent(in) :: ndime      
      real(rp),    intent(in) :: gpcod(ndime)   
      integer(ip) :: idamp,jdamp,idime,indamp,iamin      
      real(rp) :: vis,tco,cp,cv
      real(rp) :: Cx,Cy,Cz
      real(rp) :: Xf(3),i(3),X(3),Xi(3)
      
      a%coefficient = 0.0_rp 

      iamin = 0   
      do idamp = 1,php%ndamp
         indamp = 0
         !Translate the mesh reference to Xo 
         Xf(1) = php%dampxf(1,idamp)-php%dampxo(1,idamp)
         Xf(2) = php%dampxf(2,idamp)-php%dampxo(2,idamp)
         Xf(3) = php%dampxf(3,idamp)-php%dampxo(3,idamp)
         i(1) = Xf(1)/abs(Xf(1))
         i(2) = Xf(2)/abs(Xf(2))
         i(3) = Xf(3)/abs(Xf(3))
         X(1) = gpcod(1)-php%dampxo(1,idamp)
         X(2) = gpcod(2)-php%dampxo(2,idamp)
         Xi(1) = X(1)*i(1) 
         Xi(2) = X(2)*i(2) 
         if(ndime == 3) then
            X(3) = gpcod(3)-php%dampxo(3,idamp)
            Xi(3) = X(3)*i(3)
         end if 
         !Measure if the translated position is inside the Xf-Xo range 
         do idime = 1,ndime
            if (Xi(idime)>epsilon(0.0_rp) .and. abs(X(idime))<abs(Xf(idime))) indamp = indamp + 1
         end do
         if(indamp == ndime) then
            iamin = iamin + 1
            jdamp = idamp
         end if
      end do

      if(iamin == 0) then
          return
      else if (iamin > 1) then
           call runend('nsc_ComputeDamping: Overlaping Damping zones')   
      end if
 
      !Coordinates and initializations.
      a%x = gpcod(1)
      a%y = gpcod(2)
      if(ndime == 3) a%z = gpcod(3)

      a%cx = php%dampco(1,jdamp) 
      a%cy = php%dampco(2,jdamp) 
      a%cz = php%dampco(3,jdamp) 

      a%xo = php%dampxo(1,jdamp) 
      a%yo = php%dampxo(2,jdamp) 
      a%zo = php%dampxo(3,jdamp) 

      a%xf = php%dampxf(1,jdamp) 
      a%yf = php%dampxf(2,jdamp) 
      a%zf = php%dampxf(3,jdamp) 

      a%p  = php%dampff(1,jdamp)   
      a%u  = php%dampff(2,jdamp)   
      a%v  = php%dampff(3,jdamp)   
      a%w  = php%dampff(4,jdamp)   
      a%t  = php%dampff(5,jdamp)   

      Cx = a%cx*(a%x-a%xo)/(a%xf-a%xo) 
      Cy = a%cy*(a%y-a%yo)/(a%yf-a%yo) 
      Cz = 0.0_rp
      if(ndime == 3) Cz = a%cz*(a%z-a%zo)/(a%zf-a%zo) 

      call php%GetPhysicalParameters(vis,tco,cp,cv)

      a%acvis = vis
      a%actco = tco
      a%accph = cp
      a%accvh = cv
      
      a%den    = a%p/((a%accph-a%accvh)*a%t)
      a%vel(1) = a%u 
      a%vel(2) = a%v
      a%vel(3) = a%w
      a%mom(1) = a%u*a%den 
      a%mom(2) = a%v*a%den
      a%mom(3) = a%w*a%den
      a%ene    = a%den*(a%accvh*a%t+0.5_rp*(a%u*a%u+a%v*a%v+a%w*a%w))

      a%coefficient = Cx + Cy + Cz  

   end subroutine
   
   subroutine nsc_ComputeRadialDamping(a,ndime,gpcod,php)      
      use typre
      implicit none
      class(NscDamping) :: a
      Class(NSCompressibleProblem) :: php 
      integer(ip), intent(in) :: ndime      
      real(rp),    intent(in) :: gpcod(ndime)   
      real(rp) :: vis,tco,cp,cv
      real(rp) :: C
      real(rp) :: X(3) = 0.0_rp
      real(rp) :: rad = 0.0_rp, xrad
      
      a%coefficient = 0.0_rp 

      !Calculate distance to Xo 
      X(1) = gpcod(1)-php%rdampxo(1)
      X(2) = gpcod(2)-php%rdampxo(2)
      if(ndime == 3) then
         X(3) = gpcod(3)-php%rdampxo(3)
      end if 
      rad = sqrt(dot_product(X,X)) 
      !Measure if the translated position is inside the radius 
      xrad = rad - php%rdampro
      if (xrad < epsilon(0.0_rp)) return

      !Coordinates and initializations.

      C = xrad**int(php%rdampex)

      call php%GetPhysicalParameters(vis,tco,cp,cv)

      a%acvis = vis
      a%actco = tco
      a%accph = cp
      a%accvh = cv
      
      a%p  = php%rdampff(1)   
      a%u  = php%rdampff(2)   
      a%v  = php%rdampff(3)   
      a%w  = php%rdampff(4)   
      a%t  = php%rdampff(5)   

      a%den    = a%p/((a%accph-a%accvh)*a%t)
      a%vel(1) = a%u 
      a%vel(2) = a%v
      a%vel(3) = a%w
      a%mom(1) = a%u*a%den 
      a%mom(2) = a%v*a%den
      a%mom(3) = a%w*a%den
      a%ene    = a%den*(a%accvh*a%t+0.5_rp*(a%u*a%u+a%v*a%v+a%w*a%w))

      a%coefficient = C*php%rdampco 

   end subroutine
   
   subroutine nsc_GetDampingCoefficient(a,dampc)
      use typre
      implicit none
      class(NscDamping) :: a
      real(rp),    intent(out) :: dampc 
      
         dampc = a%coefficient
         
   end subroutine

   subroutine nsc_GetDampingForce(a,ndime,dtinv,elexd,elexm,elexe)
      use typre
      implicit none
      class(NscDamping) :: a
      integer(ip), intent(in)  :: ndime
      real(rp), intent(in)  :: dtinv
      real(rp), intent(inout)  :: elexd(1),elexm(ndime),elexe(1)

      !Continuity force term
      elexd(1) = elexd(1) + a%den*a%coefficient*dtinv
                               
      !Momentum force term     
      elexm(1:ndime) = elexm(1:ndime) + a%mom(1:ndime)*a%coefficient*dtinv

      !Energy force term
      elexe(1) =  elexe(1) + a%ene*a%coefficient*dtinv

   end subroutine

   subroutine nsc_pr_GetDampingForce(a,ndime,dtinv,gpden,elexd,elexm,elexe)
      use typre
      implicit none
      class(NscDamping) :: a
      integer(ip), intent(in)  :: ndime
      real(rp), intent(in)  :: gpden,dtinv
      real(rp), intent(inout)  :: elexd(1),elexm(ndime),elexe(1)

      !Continuity force term
      elexd(1) = elexd(1) + gpden*a%coefficient*a%coefficient*a%p*dtinv
                               
      !Momentum force term     
      elexm(1:ndime) = elexm(1:ndime) + gpden*a%coefficient*a%vel(1:ndime)*dtinv

      !Energy force term
      elexe(1) =  elexe(1) + gpden*a%accph*a%coefficient*a%t*dtinv

   end subroutine

end module

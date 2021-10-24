subroutine nsi_ManufacturedBoundaryConditions(a)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   integer(ip) :: npoin,ipoin,ibopo,ndime
   real(rp), pointer :: exnor(:,:)
   real(rp), pointer :: coord(:)
   ! work variables
   real(rp)    :: x,y,z,pi,r,teta,rad,h,sigma,v0,vx,vy,vz
   
   if (a%ManufacturedBoundaryCondition == 1) then
      a%kfl_fixbo = -1
      a%kfl_fixno = -1
      a%bvess = 0
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      
      if(ndime==2)then
         call runend('nsi_ManufacturedBoundaryConditions: exact solution not ready for 2d case')   
      end if

      pi= 4.0_rp*atan(1.0_rp)        
      rad  = 1.0_rp
      h  = 1.0_rp
      v0 = 1.0_rp
      sigma = 0.2_rp

      do ipoin = 1,npoin
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if (ibopo /= 0) then
            a%kfl_fixno(ndime,ipoin) = 1_ip
            call a%Mesh%GetPointCoord(ipoin,coord)
            x = coord(1) 
            y = coord(2) 
            z = coord(3) 
            r = sqrt(x*x+y*y)
            if (h-abs(z) <= 1.0D-6) then 
                vz = v0 * (rad*rad/(2*pi*sigma*sigma))*exp(-(r*r/(2*sigma*sigma))) 
            elseif (abs(z) <= 1.0D-6) then 
                vz = v0 * (rad*rad/(2*pi*sigma*sigma))*exp(-(r*r/(2*sigma*sigma))) 
            else
                vx = 0.0_rp
                vy = 0.0_rp
                vz = 0.0_rp
                a%kfl_fixno(1,ipoin) = 1_ip
                a%kfl_fixno(2,ipoin) = 1_ip
                a%bvess(1,ipoin,1) = vx
                a%bvess(2,ipoin,1) = vy
            endif
            a%bvess(ndime,ipoin,1) = vz
         endif
      enddo   

   endif


end subroutine

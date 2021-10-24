module Mod_SupOperations
   use typre

contains

      subroutine sup_GradStress(ndime,auxtens,auxGrad,grsig,gpgrad)
         implicit none   
         integer(ip), intent(in)   :: ndime,auxtens,auxGrad
         real(rp), intent(in)      :: grsig(auxtens,ndime)
         real(rp), intent(out)     :: gpgrad(auxGrad)
         integer(ip)               :: idime
            
         if(ndime==2)then

            gpgrad(1)=grsig(1,1)   !sxx
            gpgrad(2)=grsig(1,2)
            gpgrad(3)=grsig(2,1)   !syy
            gpgrad(4)=grsig(2,2)
            gpgrad(5)=grsig(3,1)   !sxy
            gpgrad(6)=grsig(3,2)       
               
         elseif(ndime==3)then
            
            gpgrad(1)=grsig(1,1)   !sxx
            gpgrad(2)=grsig(1,2)
            gpgrad(3)=grsig(1,3)
            gpgrad(4)=grsig(2,1)   !syy
            gpgrad(5)=grsig(2,2)
            gpgrad(6)=grsig(2,3)
            gpgrad(7)=grsig(3,1)   !szz
            gpgrad(8)=grsig(3,2)
            gpgrad(9)=grsig(3,3)
            gpgrad(10)=grsig(4,1)  !szy
            gpgrad(11)=grsig(4,2)
            gpgrad(12)=grsig(4,3)   
            gpgrad(13)=grsig(5,1)  !sxz
            gpgrad(14)=grsig(5,2)
            gpgrad(15)=grsig(5,3)
            gpgrad(16)=grsig(6,1)  !sxy
            gpgrad(17)=grsig(6,2)
            gpgrad(18)=grsig(6,3)
               
         end if

      end subroutine sup_GradStress
      
      subroutine sup_GradVel(ndime,auxSGrad,grvel,gpsgrad)
         implicit none   
         integer(ip), intent(in)   :: ndime,auxSGrad
         real(rp), intent(in)      :: grvel(ndime,ndime)
         real(rp), intent(out)     :: gpsgrad(auxSGrad)
         integer(ip)               :: idime
            
         if(ndime==2)then
            gpsgrad(1)=grvel(1,1)                          !dudx
            gpsgrad(2)=0.5_rp*(grvel(1,2) + grvel(2,1))    !0.5*(dudy+dvdx)
            gpsgrad(3)=gpsgrad(2)                          !0.5*(dudy+dvdx)
            gpsgrad(4)=grvel(2,2)                          !dvdy
               
         elseif(ndime==3)then

            gpsgrad(1)=grvel(1,1)                          !dudx
            gpsgrad(2)=0.5_rp*(grvel(1,2) + grvel(2,1))    !0.5*(dudy+dvdx)
            gpsgrad(3)=0.5_rp*(grvel(1,3) + grvel(3,1))    !0.5*(dudy+dvdx)
            gpsgrad(4)=gpsgrad(2)                          !0.5*(dudz+dwdx)
            gpsgrad(5)=grvel(2,2)                          !dvdy
            gpsgrad(6)=0.5_rp*(grvel(2,3) + grvel(3,2))    !0.5*(dvdz+dwdy) 
            gpsgrad(7)=gpsgrad(3)
            gpsgrad(8)=gpsgrad(6)
            gpsgrad(9)=grvel(3,3)
              
         end if

      end subroutine sup_GradVel
      
      subroutine sup_GradVelNS(ndime,auxSGrad,grvel,gpsgrad)
         use typre
         implicit none   
         integer(ip), intent(in)   :: ndime,auxSGrad
         real(rp), intent(in)      :: grvel(ndime,ndime)
         real(rp), intent(out)     :: gpsgrad(auxSGrad)
         integer(ip)               :: idime
            
         if(ndime==2)then
            gpsgrad(1)=grvel(1,1)                          !dudx
            gpsgrad(2)=grvel(1,2)                          !dudy
            gpsgrad(3)=grvel(2,1)                          !dvdx
            gpsgrad(4)=grvel(2,2)                          !dvdy
               
         elseif(ndime==3)then
            
            gpsgrad(1)=grvel(1,1)                          !dudx
            gpsgrad(2)=grvel(1,2)                          !dudy
            gpsgrad(3)=grvel(1,3)                          !dudz
            gpsgrad(4)=grvel(2,1)                          !dvdx
            gpsgrad(5)=grvel(2,2)                          !dvdy
            gpsgrad(6)=grvel(2,3)                          !dvdz 
            gpsgrad(7)=grvel(3,1)                          !dwdx
            gpsgrad(8)=grvel(3,2)                          !dwdy
            gpsgrad(9)=grvel(3,3)                          !dwdz
               
         end if

      end subroutine sup_GradVelNS      
      
      subroutine sup_grsigRP(ndime,auxtens,auxGrad,gpgrad,grsig)
         implicit none   
         integer(ip), intent(in)   :: ndime,auxtens,auxGrad
         real(rp), intent(out)     :: grsig(auxtens,ndime)
         real(rp), intent(in)      :: gpgrad(auxGrad)
         integer(ip)               :: idime
            
         if(ndime==2)then

            grsig(1,1)=gpgrad(1)   !sxx
            grsig(1,2)=gpgrad(2)
            grsig(2,1)=gpgrad(3)   !syy
            grsig(2,2)=gpgrad(4)
            grsig(3,1)=gpgrad(5)   !sxy
            grsig(3,2)=gpgrad(6)
                  
         elseif(ndime==3)then
            
            grsig(1,1)=gpgrad(1)   !sxx
            grsig(1,2)=gpgrad(2)
            grsig(1,3)=gpgrad(3)
            grsig(2,1)=gpgrad(4)   !syy
            grsig(2,2)=gpgrad(5)
            grsig(2,3)=gpgrad(6)
            grsig(3,1)=gpgrad(7)   !szz
            grsig(3,2)=gpgrad(8)
            grsig(3,3)=gpgrad(9)
            grsig(4,1)=gpgrad(10)  !szy
            grsig(4,2)=gpgrad(11)
            grsig(4,3)=gpgrad(12)   
            grsig(5,1)=gpgrad(13)  !sxz
            grsig(5,2)=gpgrad(14)
            grsig(5,3)=gpgrad(15)
            grsig(6,1)=gpgrad(16)  !sxy
            grsig(6,2)=gpgrad(17)
            grsig(6,3)=gpgrad(18)
               
         end if   
      end subroutine sup_grsigRP  
      
      subroutine sup_grvelRP(ndime,auxSGrad,gpsgrad,grvel)
      
         implicit none   
         integer(ip), intent(in)   :: ndime,auxSGrad
         real(rp), intent(out)     :: grvel(ndime,ndime)
         real(rp), intent(in)      :: gpsgrad(auxSGrad)
         integer(ip)               :: idime
            
         if(ndime==2)then      
               
            grvel(1,1)=gpsgrad(1)
            grvel(1,2)=gpsgrad(2)
            grvel(2,1)=gpsgrad(3)
            grvel(2,2)=gpsgrad(4)
               
         elseif(ndime==3)then
          
            grvel(1,1)=gpsgrad(1)
            grvel(1,2)=gpsgrad(2)
            grvel(1,3)=gpsgrad(3)
            grvel(2,1)=gpsgrad(4)
            grvel(2,2)=gpsgrad(5)
            grvel(2,3)=gpsgrad(6)
            grvel(3,1)=gpsgrad(7)
            grvel(3,2)=gpsgrad(8)
            grvel(3,3)=gpsgrad(9)       
               
         end if       

      end subroutine sup_grvelRP   
      
      subroutine sup_Getgpconv(ndime,acden,gpvel,grvel,gpconv)
            use typre
            implicit none   
            integer(ip), intent(in)   :: ndime
            real(rp), intent(in)      :: gpvel(ndime),grvel(ndime,ndime),acden
            real(rp), intent(out)     :: gpconv(ndime)
            integer(ip)               :: idime,jdime
            
            if(ndime==2)then
              gpconv(1) = acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2))
              gpconv(2) = acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2))
            elseif(ndime==3)then
              gpconv(1) = acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2) + gpvel(3)*grvel(1,3))
              gpconv(2) = acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2) + gpvel(3)*grvel(2,3))
              gpconv(3) = acden*(gpvel(1)*grvel(3,1) + gpvel(2)*grvel(3,2) + gpvel(3)*grvel(3,3))
            end if
      end subroutine sup_Getgpconv
      
      
      subroutine sup_divergenceU(ndime,grvel,gpdiv)
            use typre
            implicit none   
            integer(ip), intent(in)   :: ndime
            real(rp), intent(in)      :: grvel(ndime,ndime)
            real(rp), intent(out)     :: gpdiv(1)
            integer(ip)               :: idime,jdime
            
            if(ndime==2)then
              gpdiv(1) = grvel(1,1) + grvel(2,2)
            elseif(ndime==3)then
              gpdiv(1) = grvel(1,1) + grvel(2,2) + grvel(3,3)
            end if
      end subroutine sup_divergenceU      
      
      
      subroutine sup_GetDivS(ndime,auxtens,grsig,gpdivs)
            use typre
            implicit none   
            integer(ip), intent(in)   :: ndime,auxtens
            real(rp), intent(in)      :: grsig(auxtens,ndime)
            real(rp), intent(out)     :: gpdivs(ndime)
            integer(ip)               :: idime
            
            if(ndime==2)then

               gpdivs(1)=grsig(1,1) + grsig(3,2)                          
               gpdivs(2)=grsig(3,1) + grsig(2,2)                     

            elseif(ndime==3)then
            
               gpdivs(1)=grsig(1,1) + grsig(6,2) + grsig(5,3)                          
               gpdivs(2)=grsig(6,1) + grsig(2,2) + grsig(4,3)                         
               gpdivs(3)=grsig(5,1) + grsig(4,2) + grsig(3,3)     
               
            end if

      end subroutine sup_GetDivS      
      

end module
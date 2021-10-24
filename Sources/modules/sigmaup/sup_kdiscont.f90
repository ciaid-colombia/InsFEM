module Mod_sup_Kdiscont
contains
   subroutine sup_Kdiscont(e,cshock,auxtens,visac,lambda,chale,grsig,grsigRPO,grsigRP,grvel,vnorm,freq3,kdisc)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: grsig(auxtens,e%ndime),grsigRPO(auxtens,e%ndime),grvel(e%ndime,e%ndime),grsigRP(auxtens,e%ndime)
      real(rp),    intent(in)    :: visac,lambda,cshock(2)
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: kdisc
      
      real(rp)                   :: freq1,freq2,freq3,freqT,gradnorm,gradnormRPO,gradnormV,delta,gradnormRP,freq21,freq22,c1,c2
      integer(ip)                :: idime,jdime,expo
      
   !Initialization
   freq1=0.0_rp
   freq2=0.0_rp
   freq3=0.0_rp
   delta=1e-8
   expo=1
      
      
   !second invariant of velocity gradient
      gradnormV=0.0_rp
      gradnorm=0.0_rp
      gradnormRPO=0.0_rp
      gradnormRP=0.0_rp
      
      do idime=1,e%ndime
         do jdime=1,e%ndime
            gradnormV=gradnormV + grvel(idime,jdime)&
             *(grvel(idime,jdime))
         end do
      end do     
      gradnormV=0.5_rp*(gradnormV)**(0.5_rp)      
      
      do idime=1,auxtens
         do jdime=1,e%ndime
            gradnorm=gradnorm + grsig(idime,jdime)&
             *(grsig(idime,jdime))          
         end do
      end do 

      do idime=1,auxtens
         do jdime=1,e%ndime             
            gradnormRP=gradnormRP + grsigRP(idime,jdime)&
             *(grsigRP(idime,jdime))             
         end do
      end do        
      
      do idime=1,auxtens
         do jdime=1,e%ndime             
            gradnormRPO=gradnormRPO + grsigRPO(idime,jdime)&
             *(grsigRPO(idime,jdime))             
         end do
      end do       
      
      gradnorm=0.5_rp*(gradnorm)**(0.5_rp)
      gradnormRP=0.5_rp*(gradnormRP)**(0.5_rp)         
      gradnormRPO=0.5_rp*(gradnormRPO)**(0.5_rp)      

      c1=cshock(2)
      c2=cshock(1)
      
      !antes chale(2)= chale(2)/npol

      freq21=c1*vnorm*chale(2)
      freq22=c2*gradnormV*chale(2)*chale(2)
      ! Characteristic velocity
      freq1 = lambda/(2.0_rp*visac)
      freq2 = (freq21 + freq22)
      freq3 = gradnormRPO/(gradnorm + delta)
      freqT = freq1*freq2*freq3        
      !initialization
      kdisc=0.0_rp
 
      !Compute the stability parameter
      kdisc=freqT      

   end subroutine
end module
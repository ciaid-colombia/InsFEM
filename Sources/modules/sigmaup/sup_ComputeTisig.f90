module Mod_sup_ComputeTisig
contains
   subroutine sup_ComputeTiSig(e,lambda,visac,vnorm,grvel,staco,chale,auxG,auxPTT,gpsig,tisig)   
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxPTT
      real(rp),    intent(in)    :: staco(4),grvel(e%ndime,e%ndime),auxG
      real(rp),    intent(in)    :: visac,lambda,gpsig((e%ndime-1)*(e%ndime-1)+2)
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: tisig      
      real(rp)                   :: freq1,freq2,freq3,freq4,freqT,gradnorm,signorm
      real(rp)                   :: auto1,auto2,raiz,auxlamb,auxchale,traza,expo,h1
      integer(ip)                :: idime,jdime,auxtens
      
   !Initialization
   auto1=0.0_rp
   auto2=0.0_rp
   freq1=0.0_rp
   freq2=0.0_rp
   freq3=0.0_rp
   freq4=0.0_rp
   traza = 0.0_rp

   auxtens=(e%ndime-1)*(e%ndime-1)+2
      
   expo=1
   !second invariant of velocity gradient
      gradnorm=0.0_rp
      do idime=1,e%ndime
         
         traza = traza + gpsig(idime)
         
         do jdime=1,e%ndime
            gradnorm=gradnorm + grvel(idime,jdime)&
             *(grvel(idime,jdime))
            
         end do
      end do     
      gradnorm=0.5_rp*(gradnorm)**(0.5_rp)
      
      signorm=0.0_rp
      
      if(e%ndime==2)then
         signorm = gpsig(1)*gpsig(1) + gpsig(2)*gpsig(2) + 2.0_rp*gpsig(3)*gpsig(3)
      elseif(e%ndime==3)then   
         signorm = gpsig(1)*gpsig(1) + gpsig(2)*gpsig(2) + gpsig(3)*gpsig(3) &
               + 2.0_rp*gpsig(4)*gpsig(4) + 2.0_rp*gpsig(5)*gpsig(5) + 2.0_rp*gpsig(6)*gpsig(6)
      endif
      
      signorm = 0.5_rp*(signorm)**(0.5_rp)
      
      signorm = signorm*auxPTT + (1_ip - auxPTT)*abs(traza)
      
!      Using the max-eigenvalue 
!       
!      raiz=(grvel(1,1)+grvel(2,2))*(grvel(1,1)+grvel(2,2)) - 4.0_rp*((grvel(1,1)+grvel(2,2))-(grvel(2,1)*grvel(1,2)))
!      auxlamb= 0.5_rp*(grvel(1,1)+grvel(2,2))
!      if(raiz>=0.0_rp)then 
!       auto1= auxlamb + sqrt(raiz) 
!       auto2= auxlamb - sqrt(raiz)    
!       auto1=auto1*auto1
!       auto2=auto2*auto2
!      elseif(raiz<0.0_rp)then
!       auto1=auxlamb*auxlamb + raiz*raiz
!       auto2=auxlamb*auxlamb + raiz*raiz
!      endif      
!       gradnorm=max(auto1,auto2)
!       gradnorm=sqrt(gradnorm)

      auxchale=max(chale(1),chale(2))
      h1=chale(1)
      !h1=0.1_rp
      
      ! Characteristic velocity
      freq1 = (1.0_rp/(staco(3)*2.0_rp*visac)) !to control grad_sym
      freq2 = staco(4)*(lambda/(2.0_rp*visac))*(vnorm/((h1/e%npol)**(1.0_rp)))
      freq3 = staco(4)*(lambda/visac)*(gradnorm)
      freq4 = staco(2)*auxG*(lambda/(2.0_rp*visac))*signorm        
      
      
      freqT = (freq1) + (freq2) + (freq3) + (freq4)
      !initialization
      tisig=0.0_rp
 
      !Compute the stability parameter
      tisig=freqT**(-expo)

   end subroutine
   
 
   
   subroutine sup_LCR_ComputeTisig(e,lambda,lambda0,auxG,beta,visac,vnorm,grvel,ExpGpPsi_Matrix,staco,chale,auxPTT,gpsig,freq2,freq3,tisig)   
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxPTT
      real(rp),    intent(in)    :: staco(4),grvel(e%ndime,e%ndime),ExpGpPsi_Matrix(e%ndime,e%ndime)
      real(rp),    intent(in)    :: visac,lambda,gpsig((e%ndime-1)*(e%ndime-1)+2),lambda0, auxG
      real(rp),    intent(in)    :: vnorm,chale(2),beta
      real(rp),    intent(out)   :: tisig, freq2, freq3     
      real(rp)                   :: freq1,freq2b,freq4,freqT,gradnorm,signorm, Omeganorm, exponorm,exponegnorm
      real(rp)                   :: auto1,auto2,raiz,auxlamb,auxchale,traza,h,vcarac,h1
      real(rp)                   :: ExpMatrix(e%ndime,e%ndime) 
      integer(ip)                :: idime,jdime,expo,auxtens,i
      
      !Initialization
      auto1=0.0_rp
      auto2=0.0_rp
      freq1=0.0_rp
      freq2=0.0_rp
      freq2b=0.0_rp
      freq3=0.0_rp
      freq4=0.0_rp
      expo=1
      auxtens=(e%ndime-1)*(e%ndime-1)+2
            
      !second invariant of velocity gradient
      exponorm=0.0_rp
      gradnorm=0.0_rp
      expMatrix = 0.0_rp

      do i=1,e%ndime
         ExpMatrix(i,i)=ExpGpPsi_Matrix(i,i)-1.0_rp
      end do   
      
      do idime=1,e%ndime
         do jdime=1,e%ndime
            
            exponorm  = exponorm  + ExpMatrix(idime,jdime)*ExpMatrix(idime,jdime)
            gradnorm=gradnorm + grvel(idime,jdime)&
             *(grvel(idime,jdime))
   
         end do
      end do   
      
      exponorm=0.5_rp*(exponorm)**(0.5_rp)
      gradnorm=0.5_rp*(gradnorm)**(0.5_rp)
      
      h=chale(1)
      
      ! Characteristic velocity
      freq1 = 1.0_rp/(2.0_rp*staco(3)*visac*(1-beta))
      freq2 = staco(4)*(lambda/(2.0_rp*visac*(1.0_rp-beta)))*(vnorm/((h/e%npol)**(1.0_rp)))
      freq3 = staco(4)*(lambda/(visac*(1.0_rp-beta)))*gradnorm
      freq4 = staco(2)*auxG*(lambda0/(visac*(1.0_rp-beta)))*exponorm   !GIESEKUSLCR 
      
      freqT = freq1 + freq2 + freq3 + freq4
      
      !initialization
      tisig=0.0_rp

      !Compute the stability parameter
      tisig=freqT**(-expo)
       
      
   end subroutine
   
end module

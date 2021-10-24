module Mod_BlasiusProfile
   use typre
   implicit none
   
   !Using: http://www.calpoly.edu/~kshollen/ME347/Handouts/Blasius.pdf

   type :: BlasiusProfileGenerator
      real(rp) :: UFree = 1.0_rp, nu = 1.0_rp, X = 1.0_rp
      
      real(rp) :: sqrtUdivNuX   !sqrt(U/nu*x)
      
      real(rp) :: EtaTable(17) = (/   0.0,  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0 /)
      real(rp) :: FTable(17)   = (/ 0.0000, 0.1659, 0.3298, 0.4868, 0.6298, 0.7513, 0.8461, 0.9131, 0.9555, 0.9795, 0.9916, 0.9969, 0.9990, 0.9997, 0.9999, 1.0000, 1.0000 /) 
contains
      
      procedure :: SetUFree
      procedure :: SetNu
      procedure :: SetDistanceXFromPlateStartingPoint
      procedure :: GetU
   
   
   end type
   
   
 
contains
   
   subroutine SetUFree(a,UFree)
      class(BlasiusProfileGenerator) :: a
      real(rp) :: UFree
      
      a%UFree = UFree;
      call ComputeSqrtUdivNux(a)
   end subroutine
   
   subroutine SetNu(a,nu)
      class(BlasiusProfileGenerator) :: a
      real(rp) :: nu
   
      a%nu = nu
      call ComputeSqrtUdivNux(a)
   end subroutine
   
   subroutine SetDistanceXFromPlateStartingPoint(a,X)
      class(BlasiusProfileGenerator) :: a
      real(rp) :: X
      
      a%X = X;
      call ComputeSqrtUdivNux(a)
   end subroutine
   
   subroutine ComputeSqrtUdivNux(a)
      class(BlasiusProfileGenerator) :: a
      
      a%sqrtUdivNuX = sqrt(a%UFree/(a%nu*a%X))
   end subroutine
   
   
   subroutine GetU(a,y,U)
      class(BlasiusProfileGenerator) :: a
      real(rp) :: y,u
      
      real(rp) :: eta
      integer(ip) :: pos
      
      real(rp) :: eta0,eta1,etadist,coeff0,coeff1,Fval0,Fval1
      
      eta = y*a%sqrtUdivNuX
      
      if (eta > 8) then
         U = a%UFree;
      elseif (eta < 0.0_rp) then
         U = 0.0_rp
      else
         pos = floor(eta/0.5_rp)+1
         
         etadist = a%EtaTable(pos+1)- a%EtaTable(pos)
         eta0 = a%EtaTable(pos)
         eta1 = a%EtaTable(pos+1)
         
         coeff0 = 1.0_rp-(eta-eta0)/etadist
         coeff1 = 1.0_rp-coeff0
         
         Fval0 = a%FTable(pos)
         Fval1 = a%Ftable(pos+1)
         
         U = a%UFree*(coeff0*Fval0+coeff1*Fval1)
         
         
         
         
      endif
   end subroutine
   
   
   
end module




 

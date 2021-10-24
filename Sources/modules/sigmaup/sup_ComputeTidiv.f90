  
   subroutine sup_ComputeTidiv(e,denac,visac,vnorm,staco,chale,tidiv)
   
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco
      real(rp),    intent(in)    :: denac,visac
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: tidiv
      
      real(rp)                   :: freq1,freq2 !,freto,freq2 ! for the viscoelastic case

      ! Characteristic velocity
      freq1 = staco*(2.0_rp*visac)            !to control grad_sym
      !initialization
      tidiv=0.0_rp
 
      !Compute the stability parameter
      tidiv=freq1

      
   end subroutine

   

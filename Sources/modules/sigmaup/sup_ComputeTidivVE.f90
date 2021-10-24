Module Mod_sup_ComputeTidivVE
contains
   subroutine sup_ComputeTidivVE(e,timom,staco,chale,tidiv)
   
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco(4)
      real(rp),    intent(in)    :: timom,chale(2)
      real(rp),    intent(out)   :: tidiv
      
      real(rp)                   :: freq1 !,freto,freq2 ! for the viscoelastic case

      ! Characteristic velocity
      freq1 = (1.0_rp*chale(2)*chale(2)/staco(1))/e%npol4
      !initialization
      tidiv=0.0_rp
 
      !Compute the stability parameter
      tidiv=freq1/timom

      
   end subroutine

end module

Module Mod_MaterialProperties
!-----------------------------------------------------
!
! This module contains the physical properties of the materials used 
!
!-----------------------------------------------------
   use typre
   use def_parame     
   implicit none
   private
   public MatProperties
   
   type MatProperties
      real(rp) ::&
         densi,&              ! Density (rho)
         visco,&              ! Viscosity (mu)
         LawViParam(10)       ! Viscosity Parameters       
      integer(ip) :: &
         lawvi                ! Viscosity Law 0=constant 1=Power Law 2=Carreau-Yasuda 
  
contains   
   end type

end module

module Mod_plcd_RotationFactory
   use typre
   use Mod_plcd_Rotation
   use Mod_plcd_Rotators
   implicit none 
   private
   public AllocatePLCDRotator

contains

   subroutine AllocatePLCDRotator(string,Rotator)
      character(5) :: string
      class(Rotation), pointer :: Rotator

      select case (string)
         
         case ('GLOBA') 
            allocate(Rotator_Not::Rotator)
         case ('LOCAL')
            allocate(Rotator_Yes::Rotator)
         case default
            call runend('Wrong type of material')
      end select

   end subroutine
 
end module
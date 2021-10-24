module Mod_plcd_ConstitutiveTensorFactory
   use typre
   use Mod_plcd_ConstitutiveTensor
   use Mod_plcd_IsotropConstitutiveTensor
   use Mod_plcd_TransIsotropConstitutiveTensor
   use Mod_plcd_OrthotropConstitutiveTensor
   use Mod_plcd_Isotrop2dPlainStrain
   use Mod_plcd_Isotrop2dPlainStress
   implicit none
   
contains

   subroutine AllocateConstitutiveTensor(ndime,LargeStrainsflag,string,string2,CT)
      integer(ip) :: ndime
      character(5) :: string,string2
      class(ConstitutiveTensor), pointer :: CT
      integer(ip) :: LargeStrainsflag
      
      select case (string)
         case('ISOTR')
            if (ndime == 3) then
               allocate(IsotropCT::CT)
            elseif (ndime == 2) then
               if (string2 == 'STRES') then
                  allocate(Isotrop2dPlainStress::CT)
               elseif (string2 == 'STRAI') then
                  allocate(Isotrop2dPlainStrain::CT)
               else
                  call runend('Wrong type of constitutive tensor')
               endif
            endif
            
         case('TRANS')
            allocate(TransIsotropCT::CT)              
            
         case('ORTHO')
            allocate(OrthotropCT::CT)
         
         case('NOTDE')
            !Do nothing
         case default
            call runend('Wrong type of constitutive tensor')
      end select
      
      !Things which only depend on the dimension
      call CT%Initialize(ndime,LargeStrainsflag)
      
   end subroutine 
   
   
end module   
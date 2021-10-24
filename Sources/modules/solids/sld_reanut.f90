subroutine sld_reanut(a,itask)
   use typre
   use Mod_Solids
   implicit none
   
   class (SolidsProblem) :: a
   integer(ip)           :: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !  Initializations (defaults)
           
       a%delta_force   = 1.0_rp  ! Force delta

       a%sld_dispRelax = 1.0_rp  ! Relaxation coefficients

       a%beta          = 0.25_rp ! For newmark time integration
       a%omega         = 0.5_rp

   elseif (itask == 1) then     !Inside the Numerical Treatment Block   

       
       if(words(1)=='RELDI') then

           a%sld_dispRelax=param(1)

       else if(words(1)=='DELTA') then

           a%delta_force=param(1)

       else if(words(1)=='BETA ') then

           a%beta = param(1)

       else if(words(1)=='OMEGA') then

           a%omega= param(1)

       endif

   elseif (itask == 100) then   !Finalize reading operations  

       !There is no Picard in solids, so we make default RHS linearization
       if(a%kfl_linea == 1) a%kfl_linea = 0

   endif

   call a%SolidSpecificReanut(itask)

end subroutine sld_reanut


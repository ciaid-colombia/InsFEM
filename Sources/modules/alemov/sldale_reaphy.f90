subroutine sldale_reaphy(a,itask)
   !-----------------------------------------------------------------------
   !> This routine sets the mesh to ALE. The displacement and velocity pointers
   !! of the mesh are also initialized. 
   !-----------------------------------------------------------------------
   use typre  
   use def_parame
   use Mod_sldAlemov
   implicit none

   class(sldAlemovProblem)  :: a
   integer(ip):: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()

   call a%Listener%getarrs(words,param,nnpar,nnwor)


   if (itask == 2) then

       if(words(1).eq.'YOUNG') then
           a%solid%young    = a%Listener%getrea('YOUNG',0.0_rp,'#Young')
       endif
       if(words(1).eq.'POISS') then
           a%solid%poisson  = a%Listener%getrea('POISS',0.0_rp,'#Poisson')
       endif
       if(words(1).eq.'RELDI') then
           a%sldale_dispRelax = param(1)
           a%sldale_dispRelax_max = a%sldale_dispRelax
       endif
       if(words(1)=='AITKE') then
            if(a%Listener%exists('YES  ')) then 
                a%kfl_doAitken= .true.
            end if
       endif
   endif

end subroutine sldale_reaphy

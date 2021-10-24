subroutine sldup_reanut(a,itask)
   use typre
   use Mod_UPSolids
   implicit none
   
   class (UPSolidsProblem) :: a
   integer(ip)             :: itask
   real(rp),     pointer   :: param(:) => NULL()
   character(5), pointer   :: words(:) => NULL()
   integer(ip),  pointer   :: nnpar => NULL(),nnwor => NULL()
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !  Initializations (defaults)
           
       a%kfl_repro     = 0       ! Res. projection not used
       a%kfl_trasg     = 0       ! Tracking of subscales
       a%mtrit         = 1       ! Tracking iterations
       a%kfl_tacsg     = 0       ! Time accuracy of subscales
       a%kfl_nolsg     = 0       ! Non-linearity of subscales

       
       a%tau_u         = 1.0_rp  ! For SUP stabilization
       a%tau_p         = 1.0_rp

       a%kfl_PrbCharL = 1.0_rp

       a%kfl_printResiduals   = .false.   !Print residual duh

   elseif (itask == 1) then     !Inside the Numerical Treatment Block   

       
       if(words(1)=='TYPEO') then

           if(words(2)=='RESID') then
               a%kfl_repro = 1 
           else if(words(2)=='TOTAL') then
               a%kfl_repro = 0 
           endif

       else if(words(1)=='TRACK') then

           if (words(2) == 'ON   ') then
               a%kfl_trasg = 1
           else
               a%kfl_trasg = 0
           end if

      else if(words(1)=='ACCUR') then

         if(a%Listener%exists('DYNAM')) a%kfl_tacsg=1
         if(a%Listener%exists('QUASI')) a%kfl_tacsg=0
    
       else if(words(1)=='STABI') then

           a%tau_u = param(1)
           a%tau_p = param(2)

       else if(words(1)=='CHARL') then

           a%kfl_PrbCharL = param(1)

       endif

   elseif (itask == 100) then   !Finalize reading operations  

       if (a%kfl_tacsg > 0) a%kfl_trasg = 1
      
   endif

end subroutine sldup_reanut

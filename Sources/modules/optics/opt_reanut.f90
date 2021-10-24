subroutine opt_reanut(a,itask)
   !Read Numerical Treatment
   use typre
   use Mod_Listen
   use Mod_Mesh
   use Mod_Optics
   implicit none
   class (OpticsProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: istab,ipart,rpalo

   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)

   !Initializations (defaults)
   if (itask == 0) then         

      a%kfl_dissi = 0   !Default dissipation scheme is Smagorinsky
      a%a_obukhov = 2.4
   
   !Inside the Numerical Treatment Block   
   elseif (itask == 1) then     
      if(words(1)=='DISSI') then
         if (a%Listener%exists('SMAGO')) then
            a%kfl_dissi = 0
            a%turbu(1)=a%Listener%param(2)  
         elseif (a%Listener%exists('WALE ')) then
            a%kfl_dissi = 2
            a%turbu(1)=a%Listener%param(2)    
            
         elseif (a%Listener%exists('EXTER')) then   
            a%kfl_dissi = 1
         endif
      elseif(words(1)=='OBUKH') then
         a%a_Obukhov = param(1)
      endif
   
   !Finalize reading operations  
   elseif (itask == 100) then   

   endif
      


 

end subroutine opt_reanut

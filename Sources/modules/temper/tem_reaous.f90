subroutine tem_reaous(a,itask)
   use typre
   use Mod_Listen
   use Mod_Temperature
   implicit none

   class(TemperatureProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
 
   !Initializations
   if (itask == 0) then
      
      a%kfl_dispa = 0
      
   !Output section   
   elseif(itask == 1) then
     
      !Step or time for post-process
      if(words(1)=='POSTP') then
         if(words(2)=='TEMPE') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for T')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
            end if
         else if(words(2)=='HEATF') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for q')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
            end if
         else if(words(2)=='ERROR') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for e')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
            end if  
         else if(words(2)=='DISSI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for e')
               if (a%npp_stepi(4) /= 0) a%kfl_dispa = 1
               if(a%Listener%exists('ATTIM')) then 
                  a%pos_times(1:10,4)=param(4:13)
                  a%kfl_dispa = 1
               endif
               
            end if     
         elseif(words(2)=='AVGTE') then
            if(words(3)=='STEPS') then
               a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Tempe')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
            endif  
            
         elseif(words(2)=='AVGDI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
            endif  

         elseif(words(2)=='GRADI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Temperature Gradient')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            endif    
         elseif(words(2)=='SHOCK') then
            if(words(3)=='STEPS') then
               a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Shock Capturing Viscosity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,8)=param(4:13)
            endif     
            
        elseif(words(2)=='SIGMA') then
            if(words(3)=='STEPS') then
               a%npp_stepi(9) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Sigma Term')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,9)=param(4:13)
            endif 
          else if(words(2)=='TESGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(10) = a%Listener%getint('STEPS',1,'#Postprocess step interval for tesgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,10)=param(4:13)
             end if
            
        end if
     
     !Direction for 1d average
     elseif(words(1)=='AVG1D') then
         if (a%Listener%exists('X    ')) then
            a%Avg1DIdime = 1
         elseif (a%Listener%exists('Y    ')) then
            a%Avg1DIdime = 2
         elseif (a%Listener%exists('Z    ')) then
            a%Avg1DIdime = 3
         endif 
  
     !Force dissipation computation
     elseif(words(1)=='DISCO') then
        if (a%Listener%exists('ON   ')) a%kfl_dispa = 2
     
      ! Forces and Moments.
     elseif(words(1)=='FORCE') then
        if (a%Listener%exists('ON   ')) a%kfl_outfm = 1
     end if
   !Finalizations   
   elseif (itask == 100) then 

   endif

end subroutine

subroutine nsi_reaous(a,itask)
   use typre
   use Mod_iofile
   use Mod_NavierStokes
   implicit none

   class(NavierStokesProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   integer(ip) :: ndime,ipara
   
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !Initializations.
      a%kfl_dispa = 0
      a%kfl_postBtract = 0      
      a%kfl_outfm = 0 

   elseif (itask == 1) then  

      ! Step or time for post-process
      if(words(1)=='POSTP') then
         if(words(2)=='VELOC') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for velocities')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
            end if
         else if(words(2)=='PRESS') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for pressure')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
            end if
         else if(words(2)=='STREA') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for streamlines')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
            end if
         else if(words(2)=='DENSI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Density')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,4)=param(4:13)
            end if
         else if(words(2)=='VISCO') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Viscosity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
            end if
         else if(words(2)=='STRES') then
            if(words(3)=='STEPS') then
               a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Stress')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
            end if
         else if(words(2)=='TURBU') then
            if(words(3)=='STEPS') then
               a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Turbulent viscosity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            end if
         else if(words(2)=='PERME') then
            if(words(3)=='STEPS') then
               a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Permeability')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,8)=param(4:13)
            end if
         else if(words(2)=='TAUSM') then
            if(words(3)=='STEPS') then
               a%npp_stepi(9) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Tausmoothing')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,9)=param(4:13)
            end if   
         else if(words(2)=='VORTI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(10) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Vorticity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,10)=param(4:13)
            end if
         else if(words(2)=='DIVER') then
            if(words(3)=='STEPS') then
               a%npp_stepi(11) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Velocity divergence')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,11)=param(4:13)
            end if
         else if(words(2)=='KEFEM') then
            if(words(3)=='STEPS') then
               a%npp_stepi(12) = a%Listener%getint('STEPS',1,'#Postprocess step interval for fem Kinetic energy')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,12)=param(4:13)
            end if
         else if(words(2)=='KESGS') then
            if(words(3)=='STEPS') then
               a%npp_stepi(13) = a%Listener%getint('STEPS',1,'#Postprocess step interval for sgs Kinetic energy')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,13)=param(4:13)
            end if
         else if(words(2)=='VESGS') then
            if(words(3)=='STEPS') then
               a%npp_stepi(14) = a%Listener%getint('STEPS',1,'#Postprocess step interval for vesgs ')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,14)=param(4:13)
            end if
         else if(words(2)=='PRSGS') then
            if(words(3)=='STEPS') then
               a%npp_stepi(26) = a%Listener%getint('STEPS',1,'#Postprocess step interval for prsgs ')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,26)=param(4:13)
            end if
         else if(words(2)=='DISSI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(15) = a%Listener%getint('STEPS',1,'#Postprocess step interval for dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,15)=param(4:13)
               if (a%npp_stepi(15) /= 0) a%kfl_dispa=1  
            end if
         elseif(words(2)=='AVGVE') then
            if (words(3) == 'STEPS') then
               a%npp_stepi(16) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_velocity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,16)=param(4:13)
            endif
            
         elseif(words(2)=='AVGDI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(17) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,17)=param(4:13)
            endif  
         !Residual
         elseif(words(2)=='RESID') then
            if(words(3)=='STEPS') then
               a%npp_stepi(18) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,18)=param(4:13)
            endif  
         
         elseif(words(2)=='REPRO') then
            if(words(3)=='STEPS') then
               a%npp_stepi(19) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,19)=param(4:13)
            endif  

         elseif(words(2)=='GPDIS') then
            if(words(3)=='STEPS') then
               a%npp_stepi(20) = a%Listener%getint('STEPS',1,'#Postprocess step interval for avg1D_Dissipation')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,20)=param(4:13)
            endif      

         else if(words(2)=='QFACT') then
            if(words(3)=='STEPS') then
               a%npp_stepi(21) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Q-factor')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,21)=param(4:13)
            end if   

         elseif(words(2)=='BOUND') then
            if(words(3)=='STEPS') then
               a%npp_stepi(22) = a%Listener%getint('STEPS',1,'#Postprocess step interval for boundary traction')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,22)=param(4:13)
            end if   
            if(a%npp_stepi(22).gt.0 ) then
                a%kfl_postBtract=1
            end if
            
          elseif(words(2)=='RELAX') then
            if(words(3)=='STEPS') then
               a%npp_stepi(23) = a%Listener%getint('STEPS',1,'#Postprocess step interval for boundary traction')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,23)=param(4:13)
            end if  
            
          else if(words(2)=='STATI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(24) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,24)=param(4:13)
            end if  
            if(words(4)=='START') then
               a%StatisticsStartTime = a%Listener%getint('START',0,'#Postprocess start time for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,24)=param(4:13)
            end if 
          
          !External Forces
          else if(words(2)=='EXTER') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(25) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,25)=param(4:13)
            end if  
            if(words(4)=='START') then
               a%StatisticsStartTime = a%Listener%getint('START',0,'#Postprocess start time for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,25)=param(4:13)
            end if 
   
          !Exact solution
          elseif(words(2) == 'EXACT') then
            if(words(3) == 'STEPS') then
               a%npp_stepi(27) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,27)=param(4:13)
            endif
            if(words(4)=='START') then
               a%StatisticsStartTime = a%Listener%getint('START',0,'#Postprocess start time for Statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,27)=param(4:13)
            endif
                
          elseif(words(2)=='TAUST') then
            if(words(3)=='make STEPS') then
               a%npp_stepi(28) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stabilization taus')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,27)=param(4:13)
            end if  

          elseif(words(2)=='PSIVA') then
            if(words(3)=='STEPS') then
               a%npp_stepi(29) = a%Listener%getint('STEPS',1,'#Postprocess step interval for internal variable psi')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,28)=param(4:13)
            end if  
            
            elseif(words(2)=='PSIRE') then
            if(words(3)=='STEPS') then
               a%npp_stepi(30) = a%Listener%getint('STEPS',1,'#Postprocess step interval for real variable psi')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,29)=param(4:13)
            end if  

         end if
      
      elseif(words(1)=='AVG1D') then
         if (a%Listener%exists('X    ')) then
            a%Avg1DIdime = 1
         elseif (a%Listener%exists('Y    ')) then
            a%Avg1DIdime = 2
         elseif (a%Listener%exists('Z    ')) then
            a%Avg1DIdime = 3
         endif 
          
      ! Step or time for post-process
      !Force dissipation computation
      elseif(words(1)=='DISCO') then
         if (a%Listener%exists('ON   ')) a%kfl_dispa = 2
         
      ! Forces and Moments.
      elseif(words(1)=='FORCE') then
         a%kfl_outfm = 1      
         a%adimf(1) = a%Listener%getrea('FX   ',0.0_rp,'#x-component of CF')
         a%adimf(2) = a%Listener%getrea('FY   ',0.0_rp,'#y-component of CF')
         a%adimf(3) = a%Listener%getrea('FZ   ',0.0_rp,'#z-component of CF')
         a%adimm(1) = a%Listener%getrea('MX   ',0.0_rp,'#x-component of CM')
         a%adimm(2) = a%Listener%getrea('MY   ',0.0_rp,'#y-component of CM')
         a%adimm(3) = a%Listener%getrea('MZ   ',0.0_rp,'#z-component of CM')
         a%origm(1) = a%Listener%getrea('OX   ',0.0_rp,'#x-component of ori')
         a%origm(2) = a%Listener%getrea('OY   ',0.0_rp,'#y-component of ori')
         a%origm(3) = a%Listener%getrea('OZ   ',0.0_rp,'#z-component of ori')
      end if
      
   elseif (itask == 100) then   
      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            a%fil_trap = adjustl(trim(a%namda))//'.nsi.tp'
            call iofile(zero,a%lun_trap,a%fil_trap,'NSI TRACKING OF POINTS')
         endif
      endif 
   endif

   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine nsi_reaous


subroutine lmn_reaous(a,itask)
   use typre
   use Mod_iofile
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   integer(ip) :: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   integer(ip) :: ipara
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !Initializations.
      
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
         elseif(words(2)=='TEMPE') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Temperature')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
            end if
         else if(words(2)=='DENSI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Density')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,4)=param(4:13)
            end if
         else if(words(2)=='PTHER') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for thermodynamic pressure')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
            end if
         else if(words(2)=='HEATF') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Heat flux')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
            end if
         else if(words(2)=='STREA') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for streamlines')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            end if
          else if(words(2)=='VESGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for vesgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,8)=param(4:13)
             end if
          else if(words(2)=='TESGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(9) = a%Listener%getint('STEPS',1,'#Postprocess step interval for tesgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,9)=param(4:13)
             end if
          else if(words(2)=='PRSGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(10) = a%Listener%getint('STEPS',1,'#Postprocess step interval for prsgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,10)=param(4:13)
             end if

         !Residual
         elseif(words(2)=='RESID') then
            if(words(3)=='STEPS') then
               a%npp_stepi(11) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Residual')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,11)=param(4:13)
            endif  
         elseif(words(2)=='REPRO') then
            if(words(3)=='STEPS') then
               a%npp_stepi(12) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Residual Projection')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,12)=param(4:13)
            endif  
         
          end if

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
         a%adimh    = a%Listener%getrea('FQ   ',0.0_rp,'#heat flux coefficient')
      end if
     
   elseif (itask == 100) then    
      if (a%kfl_trasg < 1) then
         if (a%npp_stepi(8)> 0 .or. a%npp_stepi(9)> 0 .or. a%npp_stepi(10)> 0) then
            a%npp_stepi(8)  = 0
            a%npp_stepi(9)  = 0
            a%npp_stepi(10) = 0
            write(*,*) 'WARNING: Turn on tracking of SGS to print SGS and residual'
         end if
      end if         
      if (a%kfl_repro == 0 .and. a%npp_stepi(11)>0) then
         a%npp_stepi(11) = 0
         write(*,*) 'WARNING: Pick OSS residual to print residual projection'
      end if         
      
      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            a%fil_trap = adjustl(trim(a%namda))//'.lmn.tp'
            call iofile(zero,a%lun_trap,a%fil_trap,'LMN TRACKING OF POINTS')
         endif
      endif 

   endif

   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine lmn_reaous


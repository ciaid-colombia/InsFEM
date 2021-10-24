subroutine nsc_reaous(a,itask)
   use typre
   use def_parame
   use Mod_Listen
   use Mod_iofile
   use Mod_NSCompressible
   implicit none

   class(NSCompressibleProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL(),ndime => NULL()
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !Initializations.
      a%kfl_outfm = 0       !Do not postprocess forces and moments

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
         else if(words(2)=='TEMPE') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for temperature')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
            end if
         else if(words(2)=='DENSI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for density')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,4)=param(4:13)
            end if
         else if(words(2)=='MOMEN') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for momentum')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
            end if
         else if(words(2)=='ENERG') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for energy')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
            end if
         else if(words(2)=='STREA') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for streamlines')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            end if
         else if(words(2)=='HEATF') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for heat flux')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,8)=param(4:13)
            end if
         else if(words(2)=='ARTDI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(9) = a%Listener%getint('STEPS',1,'#Postprocess step interval for shock diffusivity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,9)=param(4:13)
            end if
         else if(words(2)=='SGSDI') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(10) = a%Listener%getint('STEPS',1,'#Postprocess step interval for subgrid diffusivity')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,10)=param(4:13)
            endif  
         elseif(words(2)=='SUBSC') then
            if(words(3)=='STEPS') then
               a%npp_stepi(12) = a%Listener%getint('STEPS',1,'#Postprocess step interval for subgrid scales')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,12)=param(4:13)
            endif  
         elseif(words(2)=='RESID') then
            if(words(3)=='STEPS') then
               a%npp_stepi(13) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Gauss point residual')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,13)=param(4:13)
            endif  
         elseif(words(2)=='REPRO') then
            if(words(3)=='STEPS') then
               a%npp_stepi(14) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Gauss point residual projection')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,14)=param(4:13)
            endif  
         elseif(words(2)=='STATI') then
            if(words(3)=='STEPS') then
               a%npp_stepi(15) = a%Listener%getint('STEPS',1,'#Postprocess step interval for statistics')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,15)=param(4:13)
            endif  
         endif  
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
         a%adimf(4) = a%Listener%getrea('FQ   ',0.0_rp,'#heat flux coefficient')
         if(words(3)=='FORCE') then  
               a%npp_stepi(11) = a%Listener%getint('STEPS',1,'#Postprocess step interval for forces and moments fields')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,11)=param(4:13)
         end if 
      end if

   elseif (itask == 100) then   
      !CAREFULL WITH THIS; MAY ONLY WORK FOR EXPLICIT 
      !These calculations are made in endste, could affect fields.
      if (a%npp_stepi(9) > 0) a%kfl_trasg = 1
      if (a%npp_stepi(10) > 0) a%kfl_trasg = 1
      if (a%npp_stepi(12) > 0) a%kfl_trasg = 1
      if (a%npp_stepi(13) > 0) a%kfl_trasg = 1
      if (a%npp_stepi(14) > 0) a%kfl_trasg = 1

      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            a%fil_trap = adjustl(trim(a%namda))//'.'//trim(adjustl(a%exmod))//'.tp'
            call iofile(zero,a%lun_trap,a%fil_trap,trim(adjustl(a%exmod))//' TRACKING OF POINTS')
         endif
      endif
      
   endif


   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine nsc_reaous


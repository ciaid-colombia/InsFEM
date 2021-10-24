subroutine plcd_reaous(a,itask)
   use typre
   use Mod_Listen
   use Mod_iofile
   use Mod_PLCD
   implicit none

   class(PLCDProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:)
   character(5), pointer :: words(:)
   integer(ip), pointer  :: nnpar,nnwor
   integer(ip) :: ndime,ipara
   
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then  

   elseif (itask == 1) then  

      ! Step or time for post-process
      if(words(1)=='POSTP') then
         if(words(2)=='DISPL') then  
            if(words(3)=='STEPS') then
                  a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for displacements')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
            end if
            if (words(3) == 'ITERA' .and. words(4) == 'ON   ') then
               a%kfl_PostprocessDisplacementAtEachIteration = 1
            endif 
         elseif(words(2)=='STRES') then  
            if(words(3)=='STEPS') then
                  a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for stress')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
            end if
         elseif(words(2)=='STRAI') then  
            if(words(3)=='STEPS') then
                  a%npp_stepi(3) = a%Listener%getint('STEPS',1,'#Postprocess step interval for strain')
                  a%PostprocessStrain = .true.
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,3)=param(4:13)
            end if   
         elseif(words(2)=='MATDA') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(4) = a%Listener%getint('STEPS',1,'#Postprocess step interval for material')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,4)=param(4:13)
            end if
            if (words(3) == 'ITERA' .and. words(4) == 'ON   ') then
               a%kfl_PostprocessMatDataAtEachIteration = 1
            endif   
         elseif(words(2)=='SUBSC') then  
            if(words(3)=='STEPS') then
                  a%npp_stepi(5) = a%Listener%getint('STEPS',1,'#Postprocess step interval for subscales')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,5)=param(4:13)
            end if   
            if (a%npp_stepi(5) /= 0) a%UPStoreSubscales = .true.   
         elseif(words(2)=='TOPOL') then  
            if(words(3)=='STEPS') then
                  a%npp_stepi(6) = a%Listener%getint('STEPS',1,'#Postprocess step interval for subscales')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,6)=param(4:13)
            end if   
         elseif(words(2)=='VELOC') then
            if(words(3)=='STEPS') then
                  a%npp_stepi(7) = a%Listener%getint('STEPS',1,'#Postprocess step interval for velocities')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            end if
         elseif(words(2)=='ACCEL') then
            if(words(3)=='STEPS') then
                  a%npp_stepi(8) = a%Listener%getint('STEPS',1,'#Postprocess step interval for accelerations')
                  if(a%Listener%exists('ATTIM')) a%pos_times(1:10,7)=param(4:13)
            end if
         end if
          
          
      ! Forces and Moments.
!       elseif(words(1)=='FORCE') then
!          a%kfl_outfm = 1      
!          do ipara = 1,3
!             a%adimf(ipara) = param(1)
!             a%adimm(ipara) = param(2)
!             a%origm(ipara) = param(3)
!          end do         
         
      end if

   elseif (itask == 100) then   

      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            a%fil_trap = adjustl(trim(a%namda))//'.plcd.tp'
            call iofile(zero,a%lun_trap,a%fil_trap,'SLD TRACKING OF POINTS')
         endif
      endif 

   endif

   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine plcd_reaous

